from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent
_MODEL_DIR = _PACKAGE_DIR / "Models"

apkp_dir = str(_MODEL_DIR / "AP_KP") + "/"
_msis_dir = _MODEL_DIR / "nrlmsise00"
_msis_candidates = [
    _msis_dir / "libnrlmsise-00.so",
    _msis_dir / "libnrlmsise-00.dll",
]
_msis_lib = next((p for p in _msis_candidates if p.exists()), _msis_candidates[0])
msislib_path = str(_msis_lib)
model_dir = str(_MODEL_DIR) + "/"
plot_dir = str(_PACKAGE_DIR.parent / "Plots") + "/"

if not _msis_lib.exists():
    raise FileNotFoundError(
        f"MSIS library not found. Checked: {', '.join(str(p) for p in _msis_candidates)}. "
        f"Please go into the {_msis_dir} directory and run 'make' to build the library."
        )

__all__ = ["apkp_dir", "msislib_path", "model_dir", "plot_dir"]
