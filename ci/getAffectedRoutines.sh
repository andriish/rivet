#! /usr/bin/env bash

# TODO: which directory is this run from?
# TODO: also test if a YODA ref file is updated, since autobinning could have been broken

## Define output as stdout or write to file
OUTFILE=$1
if [[ -z "$OUTFILE" || "$OUTFILE" = "-" ]]; then
    OUTFILE=/dev/stdout
fi

AFFECTEDFILES=/tmp/$$_affected.txt
CHECKEDFILES=/tmp/$$_checked.txt
SCANFILES=/tmp/$$_scanned.txt

function getAffectedFiles {
    fn="$@";
    if grep -Fxq "$fn" $CHECKEDFILES; then
        : #pass
        #echo "[INFO] $fn already checked, skip it"
    else
        echo "[INFO] Checking which analyses depend on $fn";
        echo $fn > $AFFECTEDFILES
        echo $fn >> $CHECKEDFILES

        headerName=`basename ${p%.*}`
        grep -iRl "$headerName" $CI_PROJECT_DIR/src/ $CI_PROJECT_DIR/include/ $CI_PROJECT_DIR/test/ >> $AFFECTEDFILES

        echo "[INFO] $fn affects these analyses:"
        cat $AFFECTEDFILES | grep -E "^analyses/.*\.(hh$|cc)$"
        cat $AFFECTEDFILES | grep -E "^analyses/.*\.(hh$|cc)$" | while read p; do
          echo $CI_PROJECT_DIR/$p $OUT
        done

        > $SCANFILES
        cat $AFFECTEDFILES | grep -E "^src/" >> $SCANFILES
        cat $AFFECTEDFILES | grep -E "\.hh$" >> $SCANFILES
        # cat $SCANFILES
        echo "[INFO] Look recursively for '$fn' dependencies"
        cat $SCANFILES | while read q; do
            #echo "$q"
            getAffectedFiles $q
        done
    fi
}



> $OUTFILE
> $CHECKEDFILES

if [[ -n "$CI_COMMIT_SHA" ]]; then
    git diff-tree --no-commit-id --name-only -r $CI_COMMIT_SHA | while read p; do
        getAffectedFiles $p
    done
fi
sort --unique -o $OUTFILE $OUTFILE

echo
echo "Analyses to compile:"
cat $OUTFILE

rm -f $CHECKEDFILES $SCANFILES
