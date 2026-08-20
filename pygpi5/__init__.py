from pathlib import Path

# Resolve paths from the installed package location.
_PACKAGE_DIR = Path(__file__).resolve().parent
_MODEL_DIR = _PACKAGE_DIR / "Models"

# Keep names aligned with config.cfg keys for easy reuse.
APKP_Path = str(_MODEL_DIR / "AP_KP") + "/"
MSISLib_Path = str(_MODEL_DIR / "nrlmsise00" / "libnrlmsise-00.so")
Model_Path = str(_MODEL_DIR) + "/"

__all__ = ["APKP_Path", "MSISLib_Path", "Model_Path"]
