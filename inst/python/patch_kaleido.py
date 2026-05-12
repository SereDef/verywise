import os
import subprocess
import pathlib
import platform
import shutil
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
                f'Object of type {type(o).__name__} is not JSON serializable')

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
    """
    Propagate display environment variables into the Python process so that
    kaleido/choreographer can find the display when launching Chromium.

    This function does NOT start Xvfb. Users on HPC clusters must run the
    following BEFORE starting R (e.g. in their job script or shell profile):

        Xvfb :99 -screen 0 1280x1024x24 &
        export DISPLAY=:99
        sleep 1

    If DISPLAY is already set in the shell environment, it will have been
    inherited by R and is already in os.environ — nothing to do.
    If it is not set, we hard-code :99 as the expected display, since that
    is what the instructions above prescribe.
    """
    if not os.environ.get("DISPLAY"):
        os.environ["DISPLAY"] = ":99"
    
    # Suppress Chromium's verbose D-Bus probing on headless nodes
    os.environ.setdefault("DBUS_SESSION_BUS_ADDRESS", "disabled:")

    # Pass sandbox-disable flags to kaleido's Chromium launcher
    # os.environ.setdefault(
    #     "KALEIDO_CHROME_ARGS",
    #     "--no-sandbox --disable-dev-shm-usage",
    # )

    # tmpdir = os.environ.get("TMPDIR", "")
    # if not tmpdir or not os.access(tmpdir, os.W_OK) or "jenkins" in tmpdir:
    #     safe_tmp = f"/tmp/{os.environ.get('USER', 'verywise')}_{os.getpid()}"
    #     os.makedirs(safe_tmp, exist_ok=True)
    #     os.environ["TMPDIR"] = safe_tmp

    # if platform.system() != "Linux":
    #     return False

    # xvfb = shutil.which("Xvfb")
    # if xvfb is None:
    #     for candidate in ["/usr/bin/Xvfb", "/usr/local/bin/Xvfb"]:
    #         if os.path.isfile(candidate):
    #             xvfb = candidate
    #             break

    # if xvfb is None:
    #     warnings.warn(
    #         "\n[verywise] No display found and Xvfb is not available.\n"
    #         "  On HPC clusters, run this before starting R:\n"
    #         "    module load Xvfb\n"
    #         "    Xvfb :99 -screen 0 1280x1024x24 &\n"
    #         "    export DISPLAY=:99\n"
    #         "  Or ask your sysadmin to make Xvfb available as a module.\n"
    #         "  Falling back to headless rendering (may fail on some nodes).",
    #         UserWarning,
    #         stacklevel=4,
    #     )
    #     return False

    # # Find a free display number
    # display = None
    # for display_num in range(99, 120):
    #     if not os.path.exists(f"/tmp/.X{display_num}-lock"):
    #         display = f":{display_num}"
    #         break
    # if display is None:
    #     return False

    # try:
    #     proc = subprocess.Popen(
    #         [xvfb, display, "-screen", "0", "1280x1024x24"],
    #         stdout=subprocess.DEVNULL,
    #         stderr=subprocess.DEVNULL,
    #     )
    #     time.sleep(1.5)
    #     if proc.poll() is not None:
    #         warnings.warn(
    #             f"\n[verywise] Xvfb found at {xvfb} but failed to start "
    #             f"on display {display}.\n"
    #             "  Try starting it manually before launching R:\n"
    #             f"    Xvfb {display} -screen 0 1280x1024x24 &\n"
    #             f"    export DISPLAY={display}",
    #             UserWarning,
    #             stacklevel=4,
    #         )
    #         return False
    #     os.environ["DISPLAY"] = display
    #     os.environ["DBUS_SESSION_BUS_ADDRESS"] = "disabled:"
    #     import atexit
    #     atexit.register(proc.terminate)
    #     return True
    # except Exception as e:
    #     warnings.warn(
    #         f"\n[verywise] Failed to start Xvfb: {e}",
    #         UserWarning,
    #         stacklevel=4,
    #     )
    #     return False

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
    system = platform.system()

    if system == "Darwin":
        os.environ["KALEIDO_CHROME_ARGS"] = "--headless=new --no-sandbox"
        return

    if system != "Linux":
        return  # Windows: leave defaults

    # On Linux (including HPC): propagate DISPLAY + sandbox flags.
    # Step 1: ensure DISPLAY and DBUS are set in Python's env
    _ensure_xvfb()

    # Step 2: respect explicit user override
    if "KALEIDO_CHROME_ARGS" in os.environ:
        return
    
    # Step 3: probe whether the display is actually reachable

    # Check display reachability via X lock file
    # (xdpyinfo is not available on all HPC clusters, but the lock file is reliable)
    display = os.environ.get("DISPLAY", "")
    display_num = display.lstrip(":")
    display_reachable = bool(display_num) and os.path.exists(f"/tmp/.X{display_num}-lock")

    if display_reachable:
        os.environ["KALEIDO_CHROME_ARGS"] = "--no-sandbox --disable-dev-shm-usage"
        return
    
    # display_reachable = False
    # xdpyinfo = shutil.which("xdpyinfo")
    # if xdpyinfo:
    #     try:
    #         r = subprocess.run(
    #             [xdpyinfo, "-display", os.environ["DISPLAY"]],
    #             capture_output=True, timeout=5)
    #         display_reachable = r.returncode == 0
    #     except (subprocess.TimeoutExpired, FileNotFoundError):
    #         pass

    # if display_reachable:
    #     # Real X display (Xvfb or desktop) is up — this is the working path.
    #     # No GL override flags needed; Chromium uses the X display directly.
    #     os.environ["KALEIDO_CHROME_ARGS"] = (
    #         "--no-sandbox --disable-dev-shm-usage")
    #     return
    
    # Step 4: no reachable display — try headless fallbacks
    warnings.warn(
        "\n[verywise] DISPLAY is set to '{}' but no X server responded.\n"
        "  If you are on an HPC cluster, make sure you ran:\n"
        "    Xvfb :99 -screen 0 1280x1024x24 &\n"
        "    export DISPLAY=:99\n"
        "    sleep 1\n"
        "  BEFORE calling library(verywise). Attempting headless "
        "fallback.".format(os.environ["DISPLAY"]),
        UserWarning,
        stacklevel=3)

    # Headless fallback candidates (tried in order)
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
                    capture_output=True, timeout=10, env=probe_env)
                if r.returncode == 0:
                    _apply_env(extra_env)
                    os.environ["KALEIDO_CHROME_ARGS"] = flags
                    return
            except (subprocess.TimeoutExpired, FileNotFoundError, PermissionError):
                continue

    # Nothing worked — set best-effort flags and let kaleido fail with a
    # clear error rather than a cryptic Python traceback.
    os.environ["KALEIDO_CHROME_ARGS"] = candidates[0][0]


_patch_kaleido()
