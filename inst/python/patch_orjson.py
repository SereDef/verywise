import sys
import json as _json
import numpy as _np

# ── Patch kaleido to use stdlib json instead of orjson ──────────────────
# kaleido v1.x calls orjson.dumps() internally, which rejects big-endian
# numpy arrays from FreeSurfer mesh files. This shim replaces orjson with
# a stdlib-json wrapper before kaleido imports it, making write_image()
# work regardless of the user's environment.

class _OrjsonShim:
    OPT_NON_STR_KEYS   = 1
    OPT_SERIALIZE_NUMPY = 2
    OPT_INDENT_2        = 4

    class JSONDecodeError(ValueError):
        pass

    @staticmethod
    def dumps(obj, *args, **kwargs):
        def _default(o):
            if isinstance(o, _np.ndarray):
                return o.astype(o.dtype.newbyteorder('='), copy=False).tolist()
            if isinstance(o, _np.generic):
                return o.item()
            raise TypeError(f'Object of type {type(o).__name__} is not JSON serializable')
        return _json.dumps(obj, default=_default).encode()

    @staticmethod
    def loads(s):
        return _json.loads(s)

# Only install the shim when orjson is actually present (it's the source of
# the problem); if absent kaleido already falls back to stdlib json itself.
try:
    import orjson as _orjson  # noqa: F401
    sys.modules['orjson'] = _OrjsonShim()
    try:
        # force kaleido's _tab module to pick up the replacement
        import importlib
        from kaleido._kaleido_tab import _tab
        importlib.reload(_tab)
    except (ImportError, ModuleNotFoundError):
        pass  # kaleido structure changed — shim is still in sys.modules, may still work
except (ImportError, ModuleNotFoundError):
    pass  # orjson not installed: nothing to do

