#!/bin/csh

# Execute "fit_table" on batch (at Lyon)
# To be submitted w/:
# qsub -P P_compass -V -q long -l os=cl7 -cwd -l sps=1 -o phi.log -j yes fit_table.csh [-f <optFile>] [-v] (plots|fit) 

set debug = 0
if ( $debug != 0 ) then
    echo Args. $#
    if ( $# > 0 ) then
	echo 1 $1
	if ( $# > 1 ) then
	    echo 2 $2
	endif
    endif
    echo =================
    echo $LD_LIBRARY_PATH
    echo =================
endif

##### PARSE COMMAND LINE FOR OPTIONS
set verbeux = " "; set optFile = 0; set period = 0
while ( $# > 0 && `printf "%s" $1 | egrep -c "^-"` != 0 )
    set arg = `printf "%s" $1 | sed -e {s/-//}`   # Strip arg of its leading "-"
    if      ( $arg == h ) then
 printf ' * fit_table.csh: Execute "./fit_table" w/ proper environment.\n'
 printf 'Usage:  fit_table [-f <optFile>] [ -p <Pij>] [-v] <mode>\n'
 printf '  -p <Pij>: Outputs copied to "hist-(K0|Lambda|iphi).<Pij>.root".\n'
 printf '            Default <optFile> = "options_fit.Pij.dat".\n'
 printf 'Batch submission @ Lyon: qsub -P P_compass -V -q long -l os=cl7 -cwd -l sps=1 -o phi.log -j yes fit_table.csh [-f <optFile>] [-v] (plots|fit)\n'
	exit 1
    else if ( $arg == v ) then
	set verbeux = "-v"
    else if ( $arg == f || $arg == p ) then
	if ( $# < 2 ) then
	    printf '** fit_table:\a Wrong arg list!\n\n'; exit 1
	endif
	shift
	if      ( $arg == f ) then
	    set optFile = $1
	else if ( $arg == p ) then
	    set period = $1
	endif
    else
	printf '** fit_table:\a Wrong option "-%s" in arg list\!\n\n' $arg
	goto usage
    endif
    shift
end

##### PAR COMMAND LINE FOR ARG.
if ( $# != 1 ) then
    printf '** fit_table:\a Wrong arg list!\n\n'; exit 1
endif
set mode = $1

##### DEFAULT optFile FOR GIVEN period
# If optFile is independently specified, it takes precedence
if ( $period != 0 && $optFile == 0 ) then
    set optFile = options_fit.$period.dat
endif

##### SET ENVIRONMENT FOR "./fit_table"
set rich_matrices = $HOME/analyse/RICH/rich_matrices

setenv LD_LIBRARY_PATH .:$rich_matrices/RooRarFit/tmp:$LD_LIBRARY_PATH
setenv ROOT_INCLUDE_PATH $rich_matrices/RooRarFit/tmp/RooRarFit

##### EXECUTE "./fit_table" IN PROPER DIRECTORY
cd $rich_matrices/CSEvent

if ( $optFile == 0 ) then
    ./fit_table $verbeux $mode
    set exeStatus = $status
else
    ./fit_table $verbeux -f $optFile $mode
    set exeStatus = $status
endif
if ( $exeStatus != 0 ) then
    printf '** fit_table: "./fit_table" fails => Exit...\n'
    exit 1
endif

if ( $period != 0 ) then
    foreach K0Lphi ( K0 Lambda iphi ephi )
	set outFile = hist_$K0Lphi.root
	if ( -e $outFile ) then
	    set outFil2 = hist_$K0Lphi.$period.root
	    printf ' * fit_table.csh: Moving "%s" -> "%s"\n' $outFile $outFil2
	    \mv $outFile $outFil2
	endif
    end
endif
exit 0
