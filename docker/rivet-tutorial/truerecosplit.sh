OUTTRUE=${1/.yoda/-true.yoda}
yoda2yoda -m ".*_true" $1 - | sed -e s/_true// > $OUTTRUE

OUTRECO=${1/.yoda/-reco.yoda}
yoda2yoda -m ".*_reco" $1 - | sed -e s/_reco// > $OUTRECO
