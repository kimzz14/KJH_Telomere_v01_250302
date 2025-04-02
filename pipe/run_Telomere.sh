############################################################################################
threadN=$1
prefix=$2
############################################################################################
python script/01.countGapN.py     -p ${prefix} -t ${threadN}
python script/02.find_telomere.py -p ${prefix} -t ${threadN}
python script/03.draw.py          -p ${prefix}
python script/04.teloCounter.py   -p ${prefix}
