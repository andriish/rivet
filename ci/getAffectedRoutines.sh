function reset { touch "$@" ; rm "$@" ; touch "$@" ;}

reset $1
reset alreadyChecked.txt


function getAffectedFiles {
 fn="$@"; 
 if grep -Fxq "$fn" alreadyChecked.txt;
   then echo "[INFO] $fn already checked, skip it"
   else
   echo $fn >> alreadyChecked.txt
   headerName=`basename ${p%.*}`.hh
   echo "[INFO] Checking which analyses depend on $p";
   grep -iRl $headerName | tee affectedFiles.txt
   cat affectedFiles.txt| grep analyses  >> $1
   echo "[INFO] $fn affects these analyses:"
   cat affectedFiles.txt| grep analyses
   reset filesToScan.txt
   cat affectedFiles.txt| grep src >> filesToScan.txt
   cat affectedFiles.txt| grep hh >> filesToScan.txt
   cat filesToScan.txt
   echo "[INFO] Look recursively '$fn' dependencies"
   while read q ; do 
     getAffectedFiles $q
   done  < filesToScan.txt;
 fi;
   }


git diff-tree --no-commit-id --name-only -r $CI_COMMIT_SHA  >> modifiedFiles.txt

while read p ; do
  getAffectedFiles $p
done < modifiedFiles.txt

echo "========"
echo Analyses to compile 
echo "========"
cat $1 | sort --unique > $1
rm alreadyChecked.txt
rm modifiedFiles.txt
rm filesToScan.txt
cat $1
