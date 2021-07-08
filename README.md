# RICH_Matrices

**RICH performance matrices computation framework**:
 - **UserEvents** to fill trees of **K0, Lambda, phi and rho**.
 - **Fitting framework = `RooRarFit`**
 - **Software to:**
   - **histogram K0,Lambda,phi and rho, w/ and w/o RICH-based pi,K,p PID,**
   - **simultaneously fit them, w/ line shapes common to PID flavours,**
   - **build RICH matrices**

## Contents
### `./userevents`:
 - `UserEvent13010.cc`: K0, Lambda (*to be processed w/* `normal/fit_table`).
 - `UserEvent13060.cc`: Inclusive phi (*to be processed w/* `normal/fit_table`).
 - `UserEvent7102.cc`:  Inclusive rho.
 - `UserEvent103`: K0, Lambda, (Incl.|Excl.)phi, ... (*w/* `CSEvent/fit_table`).
### `./normal`:
 - Original fitting software for K0, Lambda, Incl. phi.
 - Execution driven by options in `option_fit.dat`.
 - Executable: `fit_table`. Get help w/ `fit_table -h`.
### `./CSEvent`:
 - Alternative sotware.
 - Same operating mode as *supra*.
### `root.macros`:
 - Macros to process ROOT file outputs.

## Setup & Build
### UserEvents
```
cp <contents_of_userevents> $PHAST/user # $PHAST is PHAST installation base dir.
cd $PHAST
make
```
### Fitting software
```
cd RooRarFit
make [ROOTBUILD=debug]
cd ../normal # As an example
make [DEBUG=1]
```

## Execution

**Add `RooRarFit` to `LD_LIBRARY_PATH`. E.g. for csh:**
```
setenv RICH_MATRICES $HOME/analyse/RICH/rich_matrices
setenv LD_LIBRARY_PATH .:$RICH_MATRICES/RooRarFit/tmp:$LD_LIBRARY_PATH
setenv ROOT_INCLUDE_PATH $RICH_MATRICES/RooRarFit/tmp/RooRarFit
```
**See README in sub-directories**
