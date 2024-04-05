# iEMmod

This repository contains the code that comes with the publication

> Hunziker, J., J. Thorbecke, J. Brackenhoff, and E. Slob, 2016,  
> Inversion of controlled-source electromagnetic reflection responses:  
> Geophysics, 81(5), F49-F57;  
> DOI: [10.1190/geo2015-0320.1](https://doi.org/10.1190/geo2015-0320.1);  

The original code was published with the article in a static archive by the SEG, see [wiki.seg.org/wiki/Software:iemmod](https://wiki.seg.org/wiki/Software:iemmod). The version in the SEG archive is (c) 2016 by the Society of Exploration Geophysicists, for more info consult the file [wiki.seg.org/wiki/Software:Disclaimer](https://wiki.seg.org/wiki/Software:Disclaimer).

We[^1] release our source code here under the CPL-1.0 license, see the file `Common_Public_License.txt`. Please cite the above article if you use the code in your research.



## LICENSES

THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THIS COMMON PUBLIC LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT.

A copy of this license can be found in the file `Common_Public_License.txt` in the directory where you have found this README.

http://www.opensource.org/licenses/cpl1.0.php

The following routines are from Seismic Unix:

```
atopkge.c
getpars.c
docpkge.c
par.h
```

For those the `SU LEGAL_STATEMENT`, which is included in the source code, applies and not the common public license. 

The following routines are from slatec: http://www.netlib.org/slatec/

```
fdump.f
i1mach.f
j4save.f
xercnt.f
xerhlt.f
xermsg.f
xerprn.f
xersve.f
xgetua.f
zqk61n.f
```

For those, the slatec guidelines apply (http://www.netlib.org/slatec/guide) and not the common public license. Routines from slatec include one of the following lines in the corresponding source file: 
```
LIBRARY   SLATEC
```

Note: Some of these routines are slightly modified. Therefore, using the original library from slatec will not work with iEMmod.


## COMPILATION-CYGWIN

COMPILE THE CODE USING WINDOWS

1. Go to www.cygwin.org and download the newest version of cygwin

2. Run the setup-file, which will install cygwin. During the installation, you can select packages to be installed. Next to the basic install, also select the following packages from the devel (development) section:

   - make: The GNU version of the 'make' utility
   - gcc-core: GNU compiler collection (C, OpenMP)
   - gcc-fortran: GNU compiler collection (Fortran)

   If you have cygwin installed already, but not these packages, just run the setup-file again and select the missing packages. 

3. In the directory, in which you have installed cygwin, you will find a home-directory. This home directory contains another directory with your user name. Copy all files related to the iEMmod code into a folder in the directory with your user name. 

4. Open the cygwin terminal and follow the steps described below on how to "compile the code using a linux or a linux-like system".


## COMPILATION-LINUX

COMPILE THE CODE USING A LINUX OR A LINUX-LIKE SYSTEM

1. To compile and link the code you first have to set the IEMMODROOT variable in the `Make_include` file which can be found in the directory where you have found this README. To use the code and demo-scripts you must also define the environment variables:

   ```
   export IEMMODROOT='your-root-iemmod-directory' (for bash/sh/ksh)
   export PATH=$IEMMODROOT/bin:$PATH (for bash/sh/ksh)
   
   setenv IEMMODROOT 'your-root-iemmod-directory' (for csh/tcsh)
   setenv PATH $IEMMODROOT/bin:$PATH (for csh/tcsh)
   ```

2. Check the compiler, `CFLAGS` and `FFLAGS` options in the file `Make_include` and adapt to the system you are using. The default options are set for the GNU C and Fortran compiler on a Linux system. The Makefile has been tested with GNU make.

3. If the compiler options are set in the `Make_include` file you can type 

```
make clean
make
```

Depending on the compiler an error about missing object files for the Bessel functions J0, J1 can be encountered during linking on the iemmod executable. These missing symbols are defined in `libm.a`; the mathematical C library. Some C-compilers add this library automatically during linking, others don't.  When you encounter this error you can solve it by adding `-lm` to line 27 in the `Make_include` file: 

```
  LIBS = -lgfortran -lm
```

## RUNNING EXAMPLES

If the compilation has finished without errors and produced an executable called iemmod you can run one of the following demo programs found in the subdirectory iemmod by running 

```
./rTE_only.scr
./refTM_refTE_newmod_iso_01.scr
```

These two scripts and variations thereof were used to compute the results of the conjugate-gradient algorithm presented in the accompanying GEOPHYSICS paper "Inversion of CSEM reflection responses" by Jürg Hunziker, Jan Thorbecke, Joeri Brackenhoff and Evert Slob. The first script `rTE_only.scr` reproduces the data shown in Figures 5.c) and 6.b). The second script is not actually used in the paper, but was added to illustrate how to invert several datasets jointly. Note that running these scripts can take several hours up to days. If you wish to only run a quick test, change the line `maxit=500` to for example `maxit=3`. This will reduce the maximum amount of iterations from 500 to 3. 

Information about the parameters are found in the file `manual.pdf` in the subdirectory doc. 

In the subfolder matlab, a Matlab script can be found to plot the convergence history of the inversion process and the model found. 

- If you get the error `iemmod: command not found`, you forgot to set the `PATH` variable as mentioned above (under 1). 

- If you encounter the error: `./rTE_only.scr: Permission denied` you have to set the execution permission on the script files with the command:

```
chmod +x *.scr.
```

[^1]: 2024-04-05: The code was uploaded to
  [github.com/emsig/iEMmod](https://github.com/emsig/iEMmod) by
  [@prisae](https://github.com/prisae),
  upon the request of the main author Jürg Hunziker.
