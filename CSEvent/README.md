## git pull

## git add -A // -A depuis la racine du projet

## git commit -m "..." // remplacer ... par la description de la modification

## git push
# CSEvent

**Alternative selection/fit software**

## Introduction:
 - Input data from `../userevents/UserEvent103`.
 - Operating mode is the same as that of `../normal`. It follows the following
  steps:
  1. Selection: `fit_table` executed w/ arg. `plots`.
  2. Fit: &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;`fit_table` executed w/ arg. `fit`.
 - The execution is driven by an options file, <i>e.g.</i>: `./options_fit.dat`.
 - Execute: `fit_table -h` for more help.

### Step 1.: `plots`.
 - Produces ROOT files of invariant mass distributions for **hadrons K0,
  Lambda and phi**.  
   Two distinct **selections** for the latter: inclusive and exclusive.
 - ROOT files are per hadron selection; 4 in all hence. 
 - Mass distributions come in **two subsets: positive and negative**, depending on
 which **polarity of the decay particle** they are meant to **examine**.  
   Indeed, when it comes to evaluate the RICH response to some particle type (<it>e.g.</it> $\pi$+) one is allowed to ID its counterpart (<it>e.g.</it> $\pi$&minus;).  
   This, neglecting the possible correlation between $\pi$+ID and $\pi$&minus;ID.
 - For each polarity, distributions come in **five species**, depending upon the ID of the
 particle examined, this time.  
   <i>Viz.</i>: `a,pi,K,p,u` for all (no ID), $\pi$ID, $K$ID, $p$ID and unID'd (RICH response is ambiguous).
 - Each species is further subdivided in **bins of momentum and angle** of the
 examined particle.  
  (<i>The angle being that of the incidence on the RICH.</i>)
 - The histogram naming refects this whole subdivision. <i>E.g.</i>:  
  `h_K0_pim_K_10_3` &nbsp; for &nbsp; `K0- K+ 25<P<27 0.12<#theta<0.30`, &nbsp; <i>i.e.</i>  
 $K0$ w/ $\pi$&minus;ID for examining $\pi$+ mis-indentification as $K$+ in the
 bin [25,27] GeV &times; [.12,.3] mrd.

### Step 2.: `fit`:
 - Simultaneously fits all 5 histos, to extract **efficiency** and
  **mis-identification** probabilities.
 - Output in the shape of ROOT (<i>recommended</i>) or PDF files:  
 `rich.root`: `TGraphErrors` of efficiency and mis-identification.  
 `test_\*.root`: `TCanvas` of `RooPlots`.

## Contents
### `libCSEvent.so`:
 - One can directly look into the input files <i>via</i>:  
 ```
 root -l
 root [0] gSystem->Load("libCSEvent.so");
 root [1] TFile *_file0 = TFile::Open(<path-to-input-data>);
 ```

### Executables:
 - Interactive: `fit_table`.
 - Batch submission: `fit_table.csh`.

### Options files:
 - Default file: `options_fit.dat`.
 - One can create dedicated options files, for processing subsets of data or else.  
  And feed them to `fit_table` from the command line.

## Setup & Build
### `UserEvent103`:
```
cp <contents_of_userevents> $PHAST/user # $PHAST is PHAST installation base dir.
cd $PHAST
make
```
### `fit_table`
```
make [DEBUG=1]
```

## Execution
 - Add `RooRarFit` to `LD_LIBRARY_PATH`. <i>E.g.</i> under csh:
```
setenv RICH_MATRICES $HOME/analyse/RICH/rich_matrices
setenv LD_LIBRARY_PATH .:$RICH_MATRICES/RooRarFit/tmp:$LD_LIBRARY_PATH
setenv ROOT_INCLUDE_PATH $RICH_MATRICES/RooRarFit/tmp/RooRarFit
```
 - Can be done <i>via</i> the script `../installRichMatrices.[c]sh`:
```
 source ../installRichMatrices.csh
```

## ROOT macros
 - Nota Bene: Local `.rootrc` overrides any global ROOT setting.
