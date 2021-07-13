#! /bin/bash

########## INSTALL ENVIRONMENT FOR rich_matrices

##### FATAL ERROR IF SCRIPT IS NOT SOURCED
if [ $0 == "$BASH_SOURCE" ]; then

        echo "$_ : This script has to be sourced !"
        exit 1
fi

##### RETRIEVE LOCATION OF SCRIPT
[ "$TERM_PROGRAM" == "Apple_Terminal" ] && export READLINK="realpath -m" || export READLINK="readlink -m" # If Mac OS or not

##### $RICH_MATRICES = DIRECTORY HOUSING SCRIPT 
export RICH_MATRICES=$(dirname $($READLINK $BASH_SOURCE))
#echo ============ $RICH_MATRICES

##### UPDATE ENVIRONMENT
export LD_LIBRARY_PATH="$RICH_MATRICES/CSEvent:$RICH_MATRICES/RooRarFit/tmp:$LD_LIBRARY_PATH"
export ROOT_INCLUDE_PATH="$RICH_MATRICES/RooRarFit/tmp/RooRarFit"
echo "$CORAL"
return;
