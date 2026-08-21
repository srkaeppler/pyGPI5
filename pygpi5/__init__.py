from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent
_MODEL_DIR = _PACKAGE_DIR / "Models"

apkp_dir = str(_MODEL_DIR / "AP_KP") + "/"
msislib_path = str(_MODEL_DIR / "nrlmsise00" / "libnrlmsise-00.so")
model_dir = str(_MODEL_DIR) + "/"
plot_dir = str(_PACKAGE_DIR.parent / "Plots") + "/"

if not Path(msislib_path).exists():
    raise FileNotFoundError(
        f"MSIS library not found at {msislib_path}. Please go into the "
        f"{Path(msislib_path).parent} directory and run 'make' to build the library."
        )

__all__ = ["apkp_dir", "msislib_path", "model_dir", "plot_dir"]
