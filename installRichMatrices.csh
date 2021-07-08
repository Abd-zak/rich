#!/bin/csh

########## INSTALL ENVIRONMENT FOR rich_matrices

##### FATAL ERROR IF SCRIPT IS NOT SOURCED
# $0 contains the location of the script it if it was run.
if ( `echo $0 | grep -c installRichMatrices.csh` != 0 ) then 
    echo "$0 : This script has to be sourced !"
    exit 1
endif

##### RETRIEVE LOCATION OF SCRIPT
# In tcsh, $_ at the beginning of the script will contain the location if the file was sourced and
# But this is not true if the script is sourced from another script.
# =>=> Instead look into the currrent thread "$$".
set script_path = `ls -l /proc/$$/fd | sed -e 's/^[^/]*//' | grep "/installRichMatrices.csh" | sed -e 's\/installRichMatrices.*\\'`

##### $RICH_MATRICES = DIRECTORY HOUSING SCRIPT 
set rich_matrices = 0; set isMacOS = 0
if ( $?TERM_PROGRAM != 0 ) then
    if ( $TERM_PROGRAM == "Apple_Terminal" ) then  # If Mac OS or not
	set isMacOS = 1
	set rich_matrices = `realpath -m $script_path`
    endif
endif
if ( $rich_matrices == 0 ) then
    set rich_matrices = `readlink -m $script_path`
endif
setenv RICH_MATRICES $rich_matrices
echo ============ $RICH_MATRICES

##### UPDATE ENVIRONMENT
setenv LD_LIBRARY_PATH $rich_matrices/CSEvent:$rich_matrices/RooRarFit/tmp:$LD_LIBRARY_PATH
setenv ROOT_INCLUDE_PATH $rich_matrices/RooRarFit/tmp/RooRarFit
