# pyGPI5
Simple python D-region Ionospheric Chemistry model

# Version 
Version: 2.0
Updated: 19 June 2026

Please use Version 2.0 of the code

# Installation
1. make a python enviornment: ```python -m venv Name_Of_Enviornment``` 

2. Clone the repository: ```git clone https://github.com/srkaeppler/pyGPI5```

3. Activate into your python enviornment and install using pip: ```pip install -r requirements.txt```

Now you need to get the MSIS wrapper working.

4. Make sure you have ```gcc``` installed and navigate to the ```nrlmsise00``` subdirectory in the ```Models``` directory.

5. Run ```make ``` you should see the code compile.  You should also notice either a ```.so``` or ```.dll``` file be produced depending on the operating system (I mainly use linux and mac).  You can also ```./nrlmsise-test```. If you produce text outputs and not errors, that mean it compiled correctly.

6. Update AP/KP as needed.  In the ```Model``` directory, you should download the AP/KP files which you can find at: ```https://amisr.com/geophys_params/``` This is particularly important if you are doing something very current, i.e., in 2026. I periodically update these files. 

7. Update the config file.  In the main directory, you will now need to update ```config.cfg``` file.  This is a personal choice I am imposing onto you, but I give the full path to the location of the ```libnrlmsise-00.so``` file, the directory containing the AP/KP files, and the ```Model``` directory.  There are probably better ways to do this and feel free to contact me with those better ways, but this works for me.  

8. You can now run the ```RunExample.ipynb``` which should reproduce Figures 1-4 in the directory.  BE AWARE THAT YOU WILL OVERWRITE THE FIGURES!  ```RunExample_test.ipynb``` does the same thing and is a little bit safer in that you will not overwrite the figures. This reproduces the Figures in the[pyGPI5 Paper](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.1028042/full)

## Notes on Packages and other items
Note the requirements are a minimal set of requirements and does not include ```jupyter``` which is strongly recommended.  Also note that the ```numpy==1.26.3```, ```scipy==1.15.3```, and ```iricore==1.9.0``` which do not conflict with each other.  Additionally ```pymsis==0.12``` is also included but I have not gotten pyMSIS properly wrapped yet(as of June 2026).  ```iricore``` is from 
```https://github.com/MIST-Experiment/iricore``` and pyMSIS corresponds to ```https://github.com/SWxTREC/pymsis```

There are both IRI and MSIS wrappers in the ```Models``` directory.  However in 2026, I discovered that newer versions of ```numpy``` required the use of ```meson``` over what I had done with ```f2py``` which were the original wrappers at least for IRI.  Instead of actually learning ```meson``` I decided to look for a different IRI wrapper and found ```iricore``` which is a drop in replacement.  If you have interest in using the original wrappers for IRI, I recommend you contact me.

For MSIS, we are still using the wrapper I got/update/wrote many years ago, although in the ```Models``` directory the file ```testNewMSISWrapper.py``` uses ```pymsis```.  Further validation is really needed and more carefully implementing, which is a TODO item.   Also, I am still using MSIS00 which is not the most current version of MSIS.

# Paper and Citation
If you do use the code, please cite the [pyGPI5 Paper](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.1028042/full) 





