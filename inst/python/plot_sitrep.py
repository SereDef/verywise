import importlib
import os
import pathlib
import platform
import shutil
import socket
import subprocess
import sys


# ── Formatting helpers ────────────────────────────────────────────────────────

def _tty():
    return hasattr(sys.stdout, "isatty") and sys.stdout.isatty()

def _c(code, text):
    return f"\033[{code}m{text}\033[0m" if _tty() else text

def ok(msg):   return _c("32;1", f"✔  {msg}")
def warn(msg): return _c("33;1", f"⚠  {msg}")
def err(msg):  return _c("31;1", f"✘  {msg}")
def info(msg): return _c("36",   f"ℹ  {msg}")
def head(msg): return _c("1",    f"\n{'─' * 56}\n  {msg}\n{'─' * 56}")

def row(label, status):
    print(f"  {label:<30} {status}")


# ── Section checks ────────────────────────────────────────────────────────────

def _section_platform():
    print(head("1 · Platform"))
    row("OS",         info(f"{platform.system()} {platform.release()}"))
    row("Node",       info(platform.node()))
    row("Machine",    info(platform.machine()))
    row("Python",     info(sys.version.split()[0]))
    row("Python exe", info(sys.executable))


def _section_display():
    print(head("2 · X Display (Xvfb)"))

    display = os.environ.get("DISPLAY", "")
    if display:
        row("DISPLAY", ok(display))
    else:
        row("DISPLAY", err("not set — run: export DISPLAY=:99"))

    lock = "/tmp/.X99-lock"
    if os.path.exists(lock):
        row("Xvfb :99 lock file", ok(f"{lock} exists"))
    else:
        row("Xvfb :99 lock file", err(
            f"{lock} missing — run:\n"
            "                                 "
            "  Xvfb :99 -screen 0 1280x1024x24 &"))

    # Socket reachability
    display_num = display.lstrip(":")
    reached = False
    for path in (f"/tmp/.X11-unix/X{display_num}",
                 f"\0/tmp/.X11-unix/X{display_num}"):
        try:
            s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            s.settimeout(1)
            s.connect(path)
            s.close()
            reached = True
            break
        except OSError:
            pass
    row("X11 socket reachable",
        ok("yes") if reached else warn("no (lock file is the primary check)"))

    dbus = os.environ.get("DBUS_SESSION_BUS_ADDRESS", "")
    row("DBUS_SESSION_BUS_ADDRESS",
        ok(dbus[:50]) if dbus else warn("not set (Chromium may print D-Bus noise)"))


def _section_chrome():
    print(head("3 · Chrome binary"))

    choreo = pathlib.Path.home() / ".local/share/choreographer/deps/chrome-linux64/chrome"
    system_bins = ["chromium-browser", "chromium", "google-chrome"]

    found = None

    if choreo.is_file():
        row("choreographer bundled", ok(str(choreo)))
        found = str(choreo)
    else:
        row("choreographer bundled", warn("not found"))

    for name in system_bins:
        p = shutil.which(name)
        if p:
            row(name, ok(p))
            if found is None:
                found = p
        else:
            row(name, warn("not found"))

    # kaleido bundled
    try:
        import kaleido
        hits = list(pathlib.Path(kaleido.__file__).parent.rglob("chrome"))
        if hits:
            row("kaleido bundled", ok(str(hits[0])))
            if found is None:
                found = str(hits[0])
        else:
            row("kaleido bundled", warn("not found"))
    except ImportError:
        row("kaleido bundled", warn("kaleido not importable"))

    if not found:
        print(f"\n  {err('No Chrome binary found — static export will fail')}")
        return

    # Live launch test
    flags = os.environ.get("KALEIDO_CHROME_ARGS", "--no-sandbox --disable-dev-shm-usage")
    print(f"\n  Launch test: {info(flags)}")
    try:
        r = subprocess.run(
            [found, *flags.split(), "--virtual-time-budget=500", "about:blank"],
            capture_output=True, timeout=10, env={**os.environ})
        snip = (r.stderr or b"")[:200].decode(errors="replace")
        canvas_fail = any(k in snip.lower() for k in ("canvas", "webgl", "error 525"))
        if r.returncode == 0 and not canvas_fail:
            row("Chrome launch (about:blank)", ok(f"exit 0"))
        else:
            row("Chrome launch (about:blank)",
                err(f"exit {r.returncode} | {snip[:120]}"))
    except subprocess.TimeoutExpired:
        row("Chrome launch (about:blank)", warn("timed out (10 s)"))
    except Exception as exc:
        row("Chrome launch (about:blank)", err(str(exc)))


def _section_packages():
    print(head("4 · Python packages"))

    packages = [
        "kaleido", "plotly", "orjson",
        "numpy", "nibabel", "nilearn",
        "matplotlib", "pandas", "scipy",
        "choreographer", "scikit-learn",
    ]
    for pkg in packages:
        try:
            mod = importlib.import_module(pkg)
            ver = getattr(mod, "__version__", "installed")
            row(pkg, ok(ver))
        except ImportError:
            row(pkg, warn("not installed"))

    # orjson shim status
    if "orjson" in sys.modules:
        orj = sys.modules["orjson"]
        is_shim = not hasattr(orj, "__spec__") or hasattr(orj, "OPT_NON_STR_KEYS")
        row("orjson shim active",
            ok("yes — big-endian FreeSurfer arrays safe") if is_shim
            else warn("no — orjson will handle serialisation (may fail on HPC)"))


def _section_kaleido():
    print(head("5 · Kaleido / Plotly config"))

    chrome_args = os.environ.get("KALEIDO_CHROME_ARGS", "")
    if chrome_args:
        row("KALEIDO_CHROME_ARGS", ok(chrome_args))
    else:
        row("KALEIDO_CHROME_ARGS", warn("not set — kaleido uses its own defaults"))

    # Flag conflict: --disable-gpu with a real display causes Error 525
    display = os.environ.get("DISPLAY", "")
    if display and chrome_args:
        bad = [f for f in ("--disable-gpu", "--use-gl=swiftshader") if f in chrome_args]
        if bad:
            row("Flag conflict", err(
                f"{', '.join(bad)} + DISPLAY set → likely cause of Error 525"))

    try:
        import kaleido
        row("kaleido path",
            info(str(pathlib.Path(kaleido.__file__).parent)))
        row("kaleido version",
            info(getattr(kaleido, "__version__", "unknown")))
    except ImportError:
        row("kaleido", err("not importable"))

    try:
        import plotly.io as pio
        row("plotly renderer", info(str(pio.renderers.default) or "(empty)"))
    except Exception as exc:
        row("plotly.io", warn(str(exc)))


def _section_export():
    print(head("6 · End-to-end Plotly → PNG export"))
    import tempfile
    try:
        import plotly.express as px
        fig = px.scatter(x=[1, 2, 3], y=[1, 4, 9], title="verywise sitrep test")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmp = f.name
        fig.write_image(tmp, format="png")
        size = pathlib.Path(tmp).stat().st_size
        pathlib.Path(tmp).unlink(missing_ok=True)
        if size > 1000:
            row("write_image (PNG)", ok(f"{size:,} bytes — success"))
        else:
            row("write_image (PNG)", warn(f"only {size} bytes — may be corrupt"))
    except Exception as exc:
        row("write_image (PNG)", err(str(exc)[:250]))


# ── Public entry point ────────────────────────────────────────────────────────

def plot_sitrep():
    """
    Print a full sit-rep of the plotting environment for verywise.

    Covers: platform, Xvfb/display, Chrome binary (+ live launch test),
    Python package versions, Kaleido config, and an end-to-end PNG export test.

    Call from R:
        reticulate::py_run_string(
          "from plot_sitrep import plot_sitrep; plot_sitrep()")
    """
    print(_c("1", "\n╔════════════════════════════════════════════════╗"))
    print(_c("1",   "║   verywise · Plotting Sit-Rep                  ║"))
    print(_c("1",   "╚════════════════════════════════════════════════╝"))

    _section_platform()
    _section_display()
    _section_chrome()
    _section_packages()
    _section_kaleido()
    _section_export()

    print(f"\n{'─' * 56}")
    print(_c("1", "  Done."))
    print("  Tip: if Error 525 persists after all checks pass,")
    print("  restart R — kaleido caches its Chrome process at startup.")
    print(f"{'─' * 56}\n")


if __name__ == "__main__":
    plot_sitrep()
