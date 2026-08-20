from pathlib import Path

# Resolve paths from the installed package location.
_PACKAGE_DIR = Path(__file__).resolve().parent
_MODEL_DIR = _PACKAGE_DIR / "Models"

# Keep names aligned with config.cfg keys for easy reuse.
apkp_dir = str(_MODEL_DIR / "AP_KP") + "/"
msislib_path = str(_MODEL_DIR / "nrlmsise00" / "libnrlmsise-00.so")
model_dir = str(_MODEL_DIR) + "/"
plot_dir = str(_PACKAGE_DIR.parent / "Plots") + "/"

__all__ = ["apkp_dir", "msislib_path", "model_dir", "plot_dir"]
