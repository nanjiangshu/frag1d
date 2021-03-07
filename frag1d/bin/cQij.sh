#!/bin/bash
# Created 2009-11-08 
pssmPath=../noLowComFilteredVariedEvalue/
seqmapPath=../new_chain6105Seqmap 
outQijPath=data/qijmatrix


pssmfilelistFile=pssmfilelist.$$.list
find $pssmPath -name "*.pssm" > $pssmfilelistFile
pssm2Qij --idtype 1  --bkfile data/background_freq_aa.txt  --seqmap $seqmapPath --sad -a AVLIPFMKRHGSTCYNEWDQ -d $outQijPath --typeprofile 1 --wb -list $pssmfilelistFile 

rm -f $pssmfilelistFile
