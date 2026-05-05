import os
import subprocess
import pathlib
import platform

import sys
import json as _json
import numpy as _np

"""
Applied by verywise before loading the main renderer.
1. Disables orjson for big-endian mesh compatibility (existing patch)
2. Auto-detects the best kaleido/Chromium rendering flags for the current environment
"""

# ── 1. orjson patch (kaleido) ────────────────────────────────────────────────
# kaleido v1.x calls orjson.dumps() internally, which rejects big-endian
# numpy arrays from FreeSurfer mesh files. This shim replaces orjson with
# a stdlib-json wrapper before kaleido imports it, making write_image()
# work regardless of the user's environment.

class _OrjsonShim:
    OPT_NON_STR_KEYS    = 1
    OPT_SERIALIZE_NUMPY = 2
    OPT_INDENT_2        = 4

    class JSONDecodeError(ValueError):
        pass

    @staticmethod
    def dumps(obj, *args, option=None, **kwargs):
        def _default(o):
            if isinstance(o, _np.ndarray):
                return o.astype(o.dtype.newbyteorder('='), copy=False).tolist()
            if isinstance(o, _np.generic):
                return o.item()
            raise TypeError(
                f'Object of type {type(o).__name__} is not JSON serializable'
            )
        return _json.dumps(obj, default=_default).encode()

    @staticmethod
    def loads(s):
        return _json.loads(s)

try:
    import orjson as _orjson  # noqa: F401
    sys.modules['orjson'] = _OrjsonShim()
    try:
        import importlib
        from kaleido._kaleido_tab import _tab
        importlib.reload(_tab)
    except (ImportError, ModuleNotFoundError):
        pass
except (ImportError, ModuleNotFoundError):
    pass

# ── 2. kaleido rendering patch ────────────────────────────────────────────────
def _patch_kaleido():
    # Respect explicit user override
    if "KALEIDO_CHROME_ARGS" in os.environ:
        return

    system = platform.system()

    if system == "Darwin":
        os.environ["KALEIDO_CHROME_ARGS"] = "--headless=new --no-sandbox"
        return

    if system != "Linux":
        return  # Windows: leave defaults

    # Ordered by likelihood of success on headless Linux / HPC
    candidates = [
        "--headless=new --no-sandbox --disable-dev-shm-usage "
        "--use-gl=angle --use-angle=swiftshader",

        "--headless=new --no-sandbox --disable-dev-shm-usage "
        "--use-gl=egl --disable-gpu-sandbox",

        "--headless=new --no-sandbox --disable-dev-shm-usage "
        "--use-gl=swiftshader --enable-unsafe-swiftshader",
    ]

    # Build a Mesa-preferring env to pass to each probe
    probe_env = dict(os.environ)
    probe_env.update({"LIBGL_ALWAYS_SOFTWARE": "1", "GALLIUM_DRIVER": "llvmpipe"})
    mesa_egl = "/lib64/libEGL_mesa.so.0"
    if os.path.exists(mesa_egl):
        probe_env["LD_PRELOAD"] = (
            (probe_env.get("LD_PRELOAD", "") + ":" + mesa_egl).lstrip(":")
        )

    chrome = _find_kaleido_chrome()

    if chrome is None:
        # Can't probe — apply Mesa env vars + best-guess flags
        _apply_mesa_env(probe_env)
        os.environ["KALEIDO_CHROME_ARGS"] = candidates[1]
        return

    for flags in candidates:
        try:
            r = subprocess.run(
                [chrome, *flags.split(), "--virtual-time-budget=500", "about:blank"],
                capture_output=True, timeout=10, env=probe_env,
            )
            if r.returncode == 0:
                _apply_mesa_env(probe_env)
                os.environ["KALEIDO_CHROME_ARGS"] = flags
                return
        except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError):
            continue

    # Nothing probed — set best-effort
    _apply_mesa_env(probe_env)
    os.environ["KALEIDO_CHROME_ARGS"] = candidates[0]


def _find_kaleido_chrome():
    try:
        import kaleido
        hits = list(pathlib.Path(kaleido.__file__).parent.rglob("chrome"))
        if hits:
            return str(hits[0])
    except Exception:
        pass
    import shutil
    return shutil.which("chromium-browser") or shutil.which("chromium")


def _apply_mesa_env(probe_env):
    """Copy Mesa-related vars into the live environment."""
    for key in ("LIBGL_ALWAYS_SOFTWARE", "GALLIUM_DRIVER", "LD_PRELOAD"):
        if key in probe_env:
            os.environ[key] = probe_env[key]


_patch_kaleido()
