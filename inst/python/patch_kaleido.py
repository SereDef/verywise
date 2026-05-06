import os
import subprocess
import pathlib
import platform
import shutil
import time
import sys
import json as _json
import numpy as _np
import warnings


# ── 1. orjson patch ──────────────────────────────────────────────────────────

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


# ── 2. Xvfb helper ───────────────────────────────────────────────────────────

def _ensure_xvfb():
    # Fix restricted TMPDIR first
    tmpdir = os.environ.get("TMPDIR", "")
    if not tmpdir or not os.access(tmpdir, os.W_OK) or "jenkins" in tmpdir:
        safe_tmp = f"/tmp/{os.environ.get('USER', 'verywise')}_{os.getpid()}"
        os.makedirs(safe_tmp, exist_ok=True)
        os.environ["TMPDIR"] = safe_tmp

    if os.environ.get("DISPLAY"):
        return True

    if platform.system() != "Linux":
        return False

    xvfb = shutil.which("Xvfb")
    if xvfb is None:
        for candidate in ["/usr/bin/Xvfb", "/usr/local/bin/Xvfb"]:
            if os.path.isfile(candidate):
                xvfb = candidate
                break

    if xvfb is None:
        warnings.warn(
            "\n[verywise] No display found and Xvfb is not available.\n"
            "  On HPC clusters, run this before starting R:\n"
            "    module load Xvfb\n"
            "    Xvfb :99 -screen 0 1280x1024x24 &\n"
            "    export DISPLAY=:99\n"
            "  Or ask your sysadmin to make Xvfb available as a module.\n"
            "  Falling back to headless rendering (may fail on some nodes).",
            UserWarning,
            stacklevel=4,
        )
        return False

    # Find a free display number
    display = None
    for display_num in range(99, 120):
        if not os.path.exists(f"/tmp/.X{display_num}-lock"):
            display = f":{display_num}"
            break
    if display is None:
        return False

    try:
        proc = subprocess.Popen(
            [xvfb, display, "-screen", "0", "1280x1024x24"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        time.sleep(1.5)
        if proc.poll() is not None:
            warnings.warn(
                f"\n[verywise] Xvfb found at {xvfb} but failed to start "
                f"on display {display}.\n"
                "  Try starting it manually before launching R:\n"
                f"    Xvfb {display} -screen 0 1280x1024x24 &\n"
                f"    export DISPLAY={display}",
                UserWarning,
                stacklevel=4,
            )
            return False
        os.environ["DISPLAY"] = display
        os.environ["DBUS_SESSION_BUS_ADDRESS"] = "disabled:"
        import atexit
        atexit.register(proc.terminate)
        return True
    except Exception as e:
        warnings.warn(
            f"\n[verywise] Failed to start Xvfb: {e}",
            UserWarning,
            stacklevel=4,
        )
        return False

# ── 3. kaleido rendering patch ───────────────────────────────────────────────

def _find_kaleido_chrome():
    """Find the Chrome binary bundled with choreographer/kaleido."""
    # choreographer stores Chrome here
    choreo_chrome = pathlib.Path.home() / ".local/share/choreographer/deps/chrome-linux64/chrome"
    if choreo_chrome.is_file():
        return str(choreo_chrome)
    try:
        import kaleido
        hits = list(pathlib.Path(kaleido.__file__).parent.rglob("chrome"))
        if hits:
            return str(hits[0])
    except Exception:
        pass
    return shutil.which("chromium-browser") or shutil.which("chromium")


def _apply_env(extra: dict):
    for k, v in extra.items():
        if v is not None:
            os.environ[k] = v


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

    # Always fix TMPDIR and try Xvfb first — cleanest path on HPC
    has_display = _ensure_xvfb()

    # Always suppress D-Bus noise
    os.environ.setdefault("DBUS_SESSION_BUS_ADDRESS", "disabled:")

    if has_display:
        # Real X display available (Xvfb or desktop) — no headless flags needed
        os.environ["KALEIDO_CHROME_ARGS"] = (
            "--no-sandbox --disable-dev-shm-usage"
        )
        return

    # No display — try headless rendering strategies
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
            probe_env = {**os.environ, **extra_env}
            try:
                r = subprocess.run(
                    [chrome, *flags.split(),
                     "--virtual-time-budget=500", "about:blank"],
                    capture_output=True, timeout=10, env=probe_env,
                )
                if r.returncode == 0:
                    _apply_env(extra_env)
                    os.environ["KALEIDO_CHROME_ARGS"] = flags
                    return
            except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError):
                continue

    # Nothing worked — best-effort ANGLE flags
    warnings.warn(
        "\n[verywise] Could not find a working Chrome rendering configuration.\n"
        "  If you are on an HPC cluster, the most reliable fix is:\n"
        "    module load Xvfb\n"
        "    Xvfb :99 -screen 0 1280x1024x24 &\n"
        "    export DISPLAY=:99\n"
        "  Then restart your R session.",
        UserWarning,
        stacklevel=3,
    )
    os.environ["KALEIDO_CHROME_ARGS"] = candidates[0][0]


_patch_kaleido()