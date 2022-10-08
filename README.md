# Energy-direct Tree Search

ConfMaker.f is the fortran program to generate a full conformational space.

EDTS_autorun_part1-3.pl are perl scripts to implement the Energy-direct tree search (EDTS) algorithm to explore the comformational space of molecule and idetify the global minimum conformer.

__DISCLAMIER__

The following scripts are based on, updated and bug fixed versions of two seperate scripts; ConfMaker.f and Confsearch.pl.  

Original Author(s): C.Y. Lin, E.I. Izgorodina, D. Brittain, and K. Zhang

Original script can be found at: https://rsc.anu.edu.au/~mcoote/scripts.php
   
When publishing findings using this codes please using the following citation:
* E. I. Izgorodina , C.Y. Lin, M. L. Coote. Phys. Chem. Chem. Phys. 9, 2507 (2007)

## Compatibility

Program | Version 
--------- | ----------
Slurm Workload Manager | == 17.02.1-2
Gaussian | == 16
GCC | == 7.1.0-2.28
perl | == 5.24.0

* Support for GAMESS, NWChem and orca is currently unavailable.

## Installation 

Once the repository has been cloned run the included config file.

``` shell
cd ./edts
bash config
cd ..
rm -rf ./edts
```
## Usage 

### ConfMaker.f

``` shell
# to execute Confmaker
mkdir g16
cd ./g16
module load GCC/7.1.0-2.28
ConfMaker $mol
```
* This requires the optimised geometry in a zmat format and saved as $mol.zmat and a $mol.input file containing all the desired rotated diheral angles.
* The .input format can be found in ConfMaker.f

### EDTS_autorun_part1-3.pl

``` shell
# to execute EDTS
# NOTE from the directory above g16
module load Perl/5.24.0-GCC-5.4.0-2.26 
EDTS $mol 
```
* all data files should be saved in the subdirectory ./g16

### Restarting EDTS after a gaussain error.
* Once you have fixed the gaussian error and have a compplete .out file
``` shell
# to execute restart EDTS
# NOTE from the directory containing the restart file or give the full path to the directory contain the restart file
module load Perl/5.24.0-GCC-5.4.0-2.26 
RESTART EDTS $<optional path> 
```

