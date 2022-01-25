# BEACON
BEACON is a Python tool that finds chemo-kinematic patterns in resolved stellar systems both observed and simulated.

## License and Referencing
This code is released under a BSD 2-clause license.

If you find this code useful for your research, please cosider citing [del Pino et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.3708D/abstract).

## Usage
At the moment, BEACON does not have an installator. Please download the code and execute it locally in your machine. BEACON runs automatically with the parameters found in the options file. An example of BEACON's execution is:

$ ./python launch_beacon.py Options_file_example.py

## Input files
BEACON requires from two files to be executed. The first one is a small python file containing all the required parameter definitions and options. An example of this file with a short description for all the options can be found in "Options_file_example.py".

The second file is the table containing the data to be analyzed. This data has to be cleaned from contaminants and poorly measured stars prior to the execution of BEACON. The file should be in ascii format with its columns arranged as:  

#x   y   vx   evx   vy   evy   vz   evz   Fe_H   eFe_H ...[Extra columns]

where x and y are the sky-projected coordinates of the stars and vx, vy the velocities. Columns preceded by "e" are the corresponding uncertainties. If those coordinates are expresed in the WCS, one can choose to use

#RA   Dec   vRA   evRA  vDec   evDec   vLOS    evLOS    Fe_H     eFe_H ...[Extra columns]

Spatial coordinates and velocities are required by BEACON. The rest of the features (columns) are optional. BEACON can use as many features as the user provides:

#RA   Dec   vRA   evRA  vDec   evDec   vLOS    evLOS    Coo1   eCoo1   Coo2   eCoo2 ...[Extra columns]

The table has to contain the columns for the sky-projected velocities, even if these are unknown:

#RA  Dec  vRA   evRA  vDec   evDec   vLOS   evLOS   Fe_H   eFe_H

39.73225  -34.51991  0.00  1.00  0.00  1.00  49.4  3.6  -1.55  0.10
...

If more than one quemical abundance is going to be used (for example: Fe/H, O/Fe, Mg/Fe...), the user should set the option "Extra_Coo_is_FeH" to "False" in the options file.

## Requirements

This code makes use of Numpy, Scipy, networkx, and Matplotlib among other scientific libraries.
