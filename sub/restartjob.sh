#!/bin/bash

###### current Perl job re-submission. Run after correcting any previous Gaussian calcs fails

# Usage:     restart <mol directory>  or if missing uses current working directory to locate restart file

#echo $1
if [ $# -eq 0 ] 
then
    cmdlne="restart"
else
    cmdlne="$1/restart"
fi

#echo $cmdlne

if [ -f $cmdlne ]
then
    echo "Restarting jobstep " 
    #cat $cmdlne
    /bin/bash $cmdlne
else
    echo "No restart point found";
fi

