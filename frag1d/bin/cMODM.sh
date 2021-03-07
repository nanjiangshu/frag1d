#!/bin/bash
# Created 2009-11-08 
pssmPath=../noLowComFilteredVariedEvalue/
seqmapPath=../new_chain6105Seqmap 
outPath=data/modmatrix


pssmfilelistFile=pssmfilelist.$$.list
find $pssmPath -name "*.pssm" > $pssmfilelistFile
pssm2modm -t per --idtype 1  --seqmap $seqmapPath --sad -a AVLIPFMKRHGSTCYNEWDQ -d $outPath --typeprofile 1 --wb -list $pssmfilelistFile 

rm -f $pssmfilelistFile
