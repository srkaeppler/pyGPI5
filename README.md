# pyGPI5
A simplified Python time dependent D-region and E-region ionospheric chemistry model based on the Glukov, Pasko, and Inan (GPI) implementation (1991) that was updated by Lehtinen in 2007.

Note: this model is valid below altitudes of 200 km where molecular chemistry dominates in the ionosphere.  Atomic Oxygen is considered in the model, but diffusion is not calculated expliclity.  If you use this code for F-region physics, you use it at your own risk.  This is your disclaimer.  F-region chemistry and diffusion may eventually be included in future releases.

# Version 
Version: 2.1
Updated: 10 August 2026

Please use Version 2.1 of the code

# Installation
1. Create a Python environment, for example:
```bash
python -m venv .venv
source .venv/bin/activate
```

2. Clone the repository:
```bash
git clone https://github.com/srkaeppler/pyGPI5
cd pyGPI5
```

3. Install the pygpi5 package in your virtual environment:
```bash
pip install -e .
```

or if you need to install with exact dependencies run:
```bash
pip install -e ".[dev]"
```

If you use uv, the equivalent workflow is:
```bash
uv sync --no-dev
```

or if you need to install with exact dependencies run:
```bash
uv sync --group dev
```

Now you need to get the MSIS wrapper working.

4. Make sure you have `gcc` installed and navigate to `pygpi5/Models/nrlmsise00`.

5. Run:
```bash
make
./nrlmsise-test
```
You should see either a `.so` or `.dll` generated depending on the operating system.

6. Update AP/KP as needed. In `pygpi5/Models/AP_KP`, download/update files from:
`https://amisr.com/geophys_params/`

7. Path handling is now package-relative in `pygpi5/__init__.py` (`apkp_dir`, `msislib_path`, `model_dir`, `plot_dir`), so manual absolute-path edits in `config.cfg` are no longer required for these paths.

8. You can now run `tutorials/RunExamples.ipynb`, which should reproduce Figures 1-4 from in the [pyGPI5 Paper](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.1028042/full) and save them to the `plots/` directory (and overwrite existing plots).

## Notes on Packages and other items
Dependencies are managed in `pyproject.toml` (not `requirements.txt`) and include `numpy>=1.26.3`, `scipy>=1.15.3`, `iricore>=1.9.0`, and `pymsis>=0.12.0`.

`jupyter` is still not included by default and is strongly recommended for tutorials.

`iricore` is from `https://github.com/MIST-Experiment/iricore` and pyMSIS corresponds to `https://github.com/SWxTREC/pymsis`.

There are both IRI and MSIS wrappers in the ```Models``` directory.  However in 2026, I discovered that newer versions of ```numpy``` required the use of ```meson``` over what I had done with ```f2py``` which were the original wrappers at least for IRI.  Instead of actually learning ```meson``` I decided to look for a different IRI wrapper and found ```iricore``` which is a drop in replacement.  If you have interest in using the original wrappers for IRI, I recommend you contact me.

For MSIS, we are still using the wrapper I got/update/wrote many years ago, although in the ```Models``` directory the file ```testNewMSISWrapper.py``` uses ```pymsis```.  Further validation is really needed and more carefully implementing, which is a TODO item.   Also, I am still using MSIS00 which is not the most current version of MSIS.

# Paper and Citation
If you do use the code, please cite the [pyGPI5 Paper](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.1028042/full) 





