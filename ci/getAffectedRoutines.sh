#! /usr/bin/env bash

# TODO: which directory is this run from?
# TODO: also test if a YODA ref file is updated, since autobinning could have been broken

OUTFILE=$1

function getAffectedFiles {
    fn="$@";
    if grep -Fxq "$fn" alreadyChecked.txt; then
        echo "[INFO] $fn already checked, skip it"
    else
        echo "[INFO] Checking which analyses depend on $fn";
        echo $fn > affectedFiles.txt
        # echo $fn >> alreadyChecked.txt

        headerName=`basename ${p%.*}`.hh
        grep -iRl "$headerName" src/ include/ test/ >> affectedFiles.txt

        echo "[INFO] $fn affects these analyses:"
        cat affectedFiles.txt | grep -E "^analyses/.*\.(hh$|cc)$"
        cat affectedFiles.txt | grep -E "^analyses/.*\.(hh$|cc)$" >> $OUTFILE

        > filesToScan.txt
        cat affectedFiles.txt | grep -E "^src/" >> filesToScan.txt
        cat affectedFiles.txt | grep -E "\.hh$" >> filesToScan.txt
        # cat filesToScan.txt
        echo "[INFO] Look recursively for '$fn' dependencies"
        cat filesToScan.txt | while read q; do
            echo "$q"
            getAffectedFiles $q
        done
    fi
}



> $OUTFILE
> alreadyChecked.txt

git diff-tree --no-commit-id --name-only -r $CI_COMMIT_SHA | while read p; do
    getAffectedFiles $p
done
sort --unique -o $OUTFILE $OUTFILE

echo
echo "Analyses to compile:"
cat $OUTFILE

rm -f alreadyChecked.txt filesToScan.txt
