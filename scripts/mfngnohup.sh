#!/bin/sh

# Runs mfng with nohup. It does not halt the program when we log out.

if test $# -ge 1
then FILE=$1
else FILE=nohupout
fi

echo "This script runs the mfng with nohup (no hang up)."
echo "If we run it on a server and log out the process will no hang up."
echo "Results will be in project*/runs.* files, details in ${FILE}_details.txt and errors in ${FILE}_err.txt."
nohup ./mrun.py >> ${FILE}_details.txt 2>> ${FILE}_err.txt < /dev/null &
