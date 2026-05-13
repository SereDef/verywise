import importlib
import json as _json
import os
import pathlib
import platform
import shutil
import subprocess
import sys
import warnings

import numpy as _np


# ── 1. orjson shim ──────────────────────────────────────────────────────────
# kaleido v1.x calls orjson.dumps() internally, which rejects big-endian numpy
# arrays produced by FreeSurfer mesh files. This shim replaces orjson with a
# stdlib-json drop-in that handles them transparently.

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
                f'Object of type {type(o).__name__} is not JSON serializable')

        return _json.dumps(obj, default=_default).encode()

    @staticmethod
    def loads(s):
        return _json.loads(s)


try:
    import orjson as _orjson  # only imported to confirm it exists
    sys.modules['orjson'] = _OrjsonShim()
    try:
        from kaleido._kaleido_tab import _tab
        importlib.reload(_tab)
    except (ImportError, ModuleNotFoundError):
        pass
except (ImportError, ModuleNotFoundError):
    pass


# ── 2. Display helper ─────────────────────────────────────────────────────────
# Users on HPC clusters must run the following BEFORE starting R:
#
#   module load Xvfb
#   Xvfb :99 -screen 0 1280x1024x24 &
#   export DISPLAY=:99
#
# If DISPLAY is already in the shell environment it will be inherited by R and
# is already in os.environ — nothing to do. If it is not set, we default to
# :99 as prescribed above.

def _ensure_display():
    """Propagate DISPLAY=:99 and silence D-Bus noise into Python's env."""
    os.environ.setdefault("DISPLAY", ":99")
    os.environ.setdefault("DBUS_SESSION_BUS_ADDRESS", "disabled:")

def _xvfb_running():
    """Return True if Xvfb is running on display :99 (lock file check)."""
    return os.path.exists("/tmp/.X99-lock")

# ── 3. Chrome binary helper ───────────────────────────────────────────────────

def _find_kaleido_chrome():
    """Return the path to the Chrome binary bundled with choreographer/kaleido."""
    choreo = pathlib.Path.home() / ".local/share/choreographer/deps/chrome-linux64/chrome"
    if choreo.is_file():
        return str(choreo)
    try:
        import kaleido
        hits = list(pathlib.Path(kaleido.__file__).parent.rglob("chrome"))
        if hits:
            return str(hits[0])
    except Exception:
        pass
    return shutil.which("chromium-browser") or shutil.which("chromium")


# ── 4. Kaleido rendering patch ────────────────────────────────────────────────

def _patch_kaleido():
    system = platform.system()

    if system == "Darwin":
        os.environ["KALEIDO_CHROME_ARGS"] = "--headless=new --no-sandbox"
        return

    if system != "Linux":
        return  # on Windows: leave kaleido defaults

    ## Ensure DISPLAY and DBUS are visible to Python
    _ensure_display()

    # Respect an explicit user override — nothing more to do
    if "KALEIDO_CHROME_ARGS" in os.environ:
        return
    
    # Xvfb is running on :99 — use the proven-working bare flags
    if _xvfb_running():
        os.environ["KALEIDO_CHROME_ARGS"] = "--no-sandbox --disable-dev-shm-usage"
        return
    
    # No Xvfb on :99 — warn and probe headless fallbacks
    warnings.warn(
        "\n[verywise] Xvfb is not running on display :99.\n"
        "  Before starting R, run:\n"
        "    module load Xvfb\n"
        "    Xvfb :99 -screen 0 1280x1024x24 &\n"
        "    export DISPLAY=:99\n"
        "  Attempting headless fallback (may fail on this node).",
        UserWarning,
        stacklevel=3,
    )

    chrome = _find_kaleido_chrome()

    candidates = [
        # ANGLE/SwiftShader — works on most HPC nodes without GPU
        ("--headless=new --no-sandbox --disable-dev-shm-usage "
         "--use-gl=angle --use-angle=swiftshader --ignore-gpu-blocklist",
         {}),
        # EGL fallback
        ("--headless=new --no-sandbox --disable-dev-shm-usage "
         "--use-gl=egl --disable-gpu-sandbox",
         {}),
        # Pure software Mesa
        ("--headless=new --no-sandbox --disable-dev-shm-usage "
         "--use-gl=swiftshader --enable-unsafe-swiftshader",
         {"LIBGL_ALWAYS_SOFTWARE": "1", "GALLIUM_DRIVER": "llvmpipe"}),
    ]

    if chrome:
        for flags, extra_env in candidates:
            try:
                r = subprocess.run(
                    [chrome, *flags.split(), "--virtual-time-budget=500", "about:blank"],
                    capture_output=True, timeout=10, env={**os.environ, **extra_env})
                if r.returncode == 0:
                    for k, v in extra_env.items():
                        os.environ[k] = v
                    os.environ["KALEIDO_CHROME_ARGS"] = flags
                    return
            except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError):
                continue

    # Nothing worked — set best-effort flags so kaleido fails with a clear
    # error rather than a cryptic Python traceback.
    os.environ["KALEIDO_CHROME_ARGS"] = candidates[0][0]


_patch_kaleido()
