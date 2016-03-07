#!/bin/bash 
# copyright (c) Nanjiang Shu, 
# Department of Biochemistry and Biophysics, 
# Stockholm Bioinformatics Center
# Sweden.
# Email: nanjiang.shu@mmk.su.se
# 
# Description:
#   Predicting one dimensional structure of proteins, i.e. three-state
#   secondary structure, eight-state Shape Strings and three-state Shape
#   Strings. 
# 
# Version: 1.3
#
# Reference:
#  Tuping Zhou*, Nanjiang Shu* and Sven Hovmoller. A Novel Method for
#  Accurate One-dimensional Protein Structure Prediction based on Fragment
#  Matching, Bioinformatics, 2010;26(4):470-477. (*Co-first author)
#   
# Usage: ./predzinc.sh [options] sequence-file-in-fasta-format
#
# Updated 2011-10-14
AddAbsolutePath(){ #$path#{{{
    local var=$1
    if [ "${var:0:1}" != "/" ];then
        var=$PWD/$var # add the absolut path
    fi
    echo $var
    return 0
}
#}}}
IsProgExist(){ #{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting $0" >&2; exit 1; }
}
#}}}
IsPathExist(){ #{{{
# supply the effective path of the program 
    if ! test -d $1; then
        echo "Directory $1 does not exist. Aborting $0" >&2
        exit
    fi
}
#}}}
IsFileExist(){ #{{{
    if ! test -f $1; then 
        echo "File $1 does not exist. Aborting $0"   >&2 
        exit 
    fi
}
#}}}
exec_cmd() { #{{{
    if [ $isPrintVerboseInfo -eq 1 ]; then 
        echo "$*"
    fi
    eval "$*"
}
#}}}

FRAG1D=`dirname $0`
FRAG1D=`AddAbsolutePath $FRAG1D`
export FRAG1D

LD_LIBRARY_PATH=$FRAG1D/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH
export BINPATH=$FRAG1D/bin
export DATADIR=$FRAG1D/data

export BLASTBIN=$BINPATH
export BLASTMAT=$DATADIR
blastbin=$BLASTBIN
blastdbname=nr
dbname=nr30
dbtype=1  #reading dumped database

usage="
`cat $FRAG1D/version`

Usage: frag1d.sh  [-cpu INT] [-outpath DIR] [-blastdb STR] [-db STR]
                  [-pssm FILE] [-pssmfilelist LISTFILE]
                  [-showexample] [-verbose] [-version] [-not-clean] [-h]
                  [-outname STR]
                   FILE

Note: the supplied sequence should be in FASTA format

Options:
 -cpu           INT  Set the number of cores to be used to run blastpgp
                     (default: 8)
 -outpath       DIR  Output the result to the specified path, (default:./)
 -outname       STR  Output name, (default: query)
 -blastdb      FILE  Database for psi-blast, (default: nr)
 -db            STR  Database for Frag1D, (default: nr30)
 -pssm         FILE  Supply pssm file in PSI-BLAST -Q flag output format,
                     If supplied, PSI-BLAST will not run
 -pssmfilelist FILE  Supply a file with a list of pssm files for batch mode
                     prediction
 -not-clean          If supplied, the intermediate files will be kept
 -verbose            Print verbose information 
 -version            Print version 
 -showexample        Print How to use with examples
 -h, --help          Print this help message and exit

Created 2009-05-07, updated 2011-10-19, Nanjiang Shu

"
examples="
Examples:
# In the subfolder \"test\" 
# Carry out the prediction by supplying a single sequence file in FASTA format
    ../frag1d.sh test1.aa

# Carry out the prediction by supplying a single pssm file
    ../frag1d.sh --pssm test2.pssm

# Carry out the prediction by a list file with a number of sequence files
    ../frag1d.sh --seqfilelist test3.seqfilelist
"
debugoptioninfo="
 --tmpdir  DIR     Set the directory to store temporary files, 
                   by default it will e created under /tmp by mktemp
 -dbtype   INT     Set the database type, (default: 1)
                   0: individual 1: dumped
"
VERSION="FRAG1D version 1.2, updated 2011-10-14"
if [ -s $FRAG1D/version ]; then 
    VERSION=`cat $FRAG1D/version`
fi

MIN_LENGTH_SEQ=10
MAX_LENGTH_SEQ=10000
SEQ_MODE_AASEQ=0
SEQ_MODE_PSSM=1

CheckSequence(){     #$file $mode#{{{
    local file=$1
    local mode=$2
    local lengthseq=0
    local cntCHDE=0
    if [ $mode -eq $SEQ_MODE_PSSM ]; then
        lengthseq=`awk 'function isnum(x){return(x==x+0)}BEGIN{cnt=0}{if (isnum(substr($1,1,1))) {cnt+=1}} END{print cnt}'  $file ` 
    elif [ $mode -eq $SEQ_MODE_AASEQ ] ; then 
        lengthseq=`awk 'BEGIN{cnt=0}/^[^>]/{str=$0;sub(/ \r\t\n/,"",str); cnt+=length(str);}END{print cnt}'  $file ` 
    else
        echo "Unrecognized sequence mode for $file. Ignore" >&2
        echo "1"
        return
    fi

    if [ $lengthseq -lt  $MIN_LENGTH_SEQ -o $lengthseq -gt $MAX_LENGTH_SEQ ];  then
        echo "Length of the sequence for $file ($lengthseq) is out of range ($MIN_LENGTH_SEQ - $MAX_LENGTH_SEQ). Ignore."  >&2
        echo "1"
        return
    else
        echo "0"
        return 
    fi
}
#}}}
PrintHelp(){ #{{{
    echo "$usage"
}
#}}}
CheckErrMsg(){ #$errFile#{{{
    local errFile=$1
    if [ -s "$errFile" ]; then 
        cat $errFile >&2 
        exit
    fi  
}
#}}}
PrintVersion(){ #{{{
    echo "$VERSION"
}
#}}}
BuildProfileFromAASeq(){ #seqfile outpath pssmfilelist #{{{
    local seqfile="$1"
    local outpath="$2"
    local pssmfilelist=$3

    local basename=`basename $seqfile`
    local rootname=${basename%.*}
    local pssmfile=$outpath/${rootname}.pssm
    local blastfile=$outpath/${rootname}.blast
    local chkfile=$outpath/${rootname}.chk
    echo "Running PSI-BLAST for $seqfile..."
    exec_cmd "$blastbin/blastpgp -a $numCPU -d "$blastdbname"  -j 3 -h 0.001 -i $seqfile -o $blastfile -m 9 -Q $pssmfile -C $chkfile"
    if [ -s $pssmfile ]; then 
        echo $pssmfile >> $pssmfilelist
    fi
}
#}}}
RunFrag1D(){ #$outpath $pssmFileListFile $frag1dParaFile #{{{
    local outpath=$1
    local pssmFileListFile=$2
    local frag1dParaFile=$3

    local aaBakFreqFile=$DATADIR/$dbname/background_freq_aa.txt
    #step 1, get pssmfilelist
    testIDListFile=$TMPDIR/test.idlist
    cat $pssmFileListFile | $BINPATH/rootname > $testIDListFile

    local testQijPath=$TMPDIR/Qij
    local testModmPath=$TMPDIR/modm
    local testFragAccPath=$TMPDIR/fragacc
    mkdir -p $testQijPath $testModmPath $testFragAccPath

    # ========  step 2: create modmatrix and qijmatrix for the given pssmFiles
    exec_cmd "$BINPATH/pssm2Qij --wb --idtype 1  --bkfile $aaBakFreqFile --newseq --sad -a AVLIPFMKRHGSTCYNEWDQ -d $testQijPath -list $pssmFileListFile 1> /dev/null 2> $errFile"
    CheckErrMsg $errFile

    exec_cmd "$BINPATH/pssm2modm -t per --wb --idtype 1 --newseq --sad -a AVLIPFMKRHGSTCYNEWDQ -d $testModmPath -list $pssmFileListFile  1> /dev/null 2> $errFile"
    CheckErrMsg $errFile

    for file in $(find $testQijPath -name "*.Qijbin");do
        local bname1=`basename $file`
        local rtname1=${bname1%.*}
        exec_cmd "ln -sf '$file' '$TMPDIR/fragacc/${rtname1}.fragaccbin'"
    done

    local trainQijPath=$DATADIR/$dbname/train_Qij
    local trainModmPath=$DATADIR/$dbname/train_modm
    local trainFragAccPath=$DATADIR/$dbname/train_fragacc
    local shapeTripleFile=$DATADIR/$dbname/train.triple.freq
    local trainIDListFile=$DATADIR/$dbname/train.idlist

    local blosumFile=$DATADIR/blosum62_tu.txt

    local qijformat=1
    local modmformat=1
    local fragaccformat=1
    local NPER=0

    IsFileExist $trainIDListFile
    IsFileExist $shapeTripleFile
    IsFileExist $aaBakFreqFile
    IsFileExist $blosumFile

    ## ========= Step 3 ==========================
    # Step 3 :  segment matching for finding candidate fragments
    local fragPath=$TMPDIR/modmFragRatio3
    echo "Searching for candidate fragments..."
    exec_cmd "$BINPATH/search_new -dbtype $dbtype --train $trainIDListFile --test $testIDListFile --train-qij $trainQijPath --test-qij $testQijPath --qijformat $qijformat --modm $trainModmPath --test-modm $testModmPath  --modmformat $modmformat  --fragacc $trainFragAccPath --test-fragacc $testFragAccPath --fragaccformat $fragaccformat  --bkfile $aaBakFreqFile --mergeside 0 --ratioscheme 3 --topN 100 --speed 0 -w 9  --result $fragPath --rb --NPer $NPER --scoretype 0 --wb --fragformat 0  1> /dev/null 2> $errFile"
    CheckErrMsg $errFile

    ## ========= Step 4 ==========================
    # Step 4: doing 1D structure prediction
    local dsspMapMethod=1
    local databaseTrain=$trainIDListFile.withnumber
    local NPER=40
    local NPte=50

    # 4.1 checkfirst.sh
    local chkFstRound1Path=$TMPDIR/chkFstRound1
    echo "Predicting 1D structure, round 1..."
    exec_cmd "$BINPATH/checkfirst  -dbtype $dbtype --dsspmap $dsspMapMethod --specialhead 4 --specialtail 4 --tshift 0 --rb --ratioscheme 0 -ltst $testIDListFile -qtst $testQijPath -mtst $testModmPath  -sloc  $fragPath  -ltrn $databaseTrain -qtrn $trainQijPath -mtrn $trainModmPath -atrn $trainFragAccPath -nptr $NPER -phsr 0 --outpath $chkFstRound1Path  1> /dev/null 2> $errFile"
    CheckErrMsg $errFile


    # 4.2 checkresult
    local chkRstRound1Path=$TMPDIR/chkRstRound1
    exec_cmd "$BINPATH/checkresult  -dbtype $dbtype --dsspmap $dsspMapMethod --specialhead 4 --specialtail 4 --tshift 0  --rb --ratioscheme 0 -ltst $testIDListFile -qtst $testQijPath -mtst $testModmPath  -sloc  $fragPath  -ltrn $databaseTrain -qtrn $trainQijPath -mtrn $trainModmPath -atrn $trainFragAccPath -nptr $NPER --round1-pred $chkFstRound1Path  -blosum $blosumFile  --outpath $chkRstRound1Path 1> /dev/null 2> $errFile"
    CheckErrMsg $errFile

    # 4.3 buildHSRFrag
    echo "Building structural profile for the predicted 1D structure..."
    HSRProfilePath=$TMPDIR/HSRprofile_chkRstRound1
    exec_cmd "$BINPATH/build_hsrfrag -dbtype $dbtype -ltst $testIDListFile -qtst $testQijPath -rloc $chkRstRound1Path -bloc $blosumFile  -ltrn $databaseTrain -qtrn $trainQijPath  --outpath $HSRProfilePath 1> /dev/null 2> $errFile"
    CheckErrMsg $errFile

    # 4.4 checkfirstRound2.sh
    local chkFstRound2Path=$TMPDIR/chkFstRound2
    echo "Predicting 1D structure, round 2..."
    exec_cmd "$BINPATH/checkfirst   -dbtype $dbtype --dsspmap $dsspMapMethod --minconf 0 --maxconf 200 --specialhead 4 --specialtail 4  --tshift 0 --rb --ratioscheme 0 -ltst $testIDListFile -qtst $testQijPath -mtst $testModmPath  -npte $NPte -sloc  $fragPath  -ltrn $databaseTrain -qtrn $trainQijPath -mtrn $trainModmPath -atrn $trainFragAccPath -nptr $NPER -phsr 1 -ploc $HSRProfilePath --outpath $chkFstRound2Path 1>  /dev/null 2> $errFile"
    CheckErrMsg $errFile

    # 4.5 checkresultRound2.sh
    local chkRstRound2Path=$TMPDIR/chkRstRound2
    exec_cmd "$BINPATH/checkresult   -dbtype $dbtype  --dsspmap $dsspMapMethod --specialhead 4 --specialtail 4 --tshift 0   --rb --ratioscheme 0 -ltst $testIDListFile -qtst $testQijPath -mtst $testModmPath  -sloc  $fragPath  -ltrn $databaseTrain -qtrn $trainQijPath -mtrn $trainModmPath -atrn $trainFragAccPath -nptr $NPER --round1-pred $chkFstRound2Path  -blosum $blosumFile --outpath $chkRstRound2Path 1>  /dev/null 2> $errFile "
    CheckErrMsg $errFile

    # ==== step 5 ================
    # report the result
    exec_cmd "$BINPATH/reportFrag1D --datapath $TMPDIR --outpath $outpath --outname $outname -l $testIDListFile 1> /dev/null 2> $errFile"
    CheckErrMsg $errFile
    echo "Frag1D succeeded. Results has been output to"
    echo "  $outpath/$outname.predfrag1d"
    if [ -d $TMPDIR/psiblast ]; then 
        /bin/cp -f $TMPDIR/psiblast/*.pssm $outpath
        /bin/cp -f $TMPDIR/psiblast/*.chk $outpath
    fi
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit 1
fi

TMPDIR=
isClean=1
outpath=./
outname=query
seqFile=
pssmFile=
pssmFileListFile=
numCPU=8

isPrintVerboseInfo=0

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        idList="$idList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit 1 ;;
            -showexample|--showexample) echo "$examples"; exit 0;;
            --version|-version) PrintVersion; exit 1;;
            --pssm|-pssm) pssmFile=$2;shift;;
            --pssmfilelist|-pssmfilelist) pssmFileListFile=$2;shift;;
            --blastdb|-blastdb) blastdbname=$2;shift;;
            -db|--db|-dbname|--dbname) dbname=$2;shift;;
            --outpath|-outpath) outpath=$2;shift;;
            -outname|--outname) outname=$2;shift;;
            -cpu|--cpu) numCPU=$2;shift;;
            -verbose|--verbose) isPrintVerboseInfo=1;;
            -nc|--nc|-not-clean|--not-clean) isClean=0;;
            -tmpdir|--tmpdir) TMPDIR=$2;shift;;
            -dbtype|--dbtype) dbtype=$2;shift;;
            -q|--q|--quiet|-quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1"; exit 2;;
        esac
    else
        seqFile=$1
    fi
    shift
done

#check for necessary programs and files
#-----------------------------------------------------------
IsProgExist $blastbin/blastpgp
IsPathExist $DATADIR
#-----------------------------------------------------------
if [ $isPrintVerboseInfo -eq 1 ]; then
    echo "BLASTMAT=$BLASTMAT"
    echo "BLASTBIN=$BLASTBIN"
    echo "BLASTDB=$BLASTDB"
fi
if [ ! -d "$outpath" ];then
    mkdir -p $outpath
fi

# create the temporary dirs
if [ "$TMPDIR" == "" ]; then
    TMPDIR=$(mktemp -d /tmp/tmpdir.frag1d.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }
else
    mkdir -p $TMPDIR
fi

TMPDIR=`AddAbsolutePath "$TMPDIR"`
errFile=$TMPDIR/frag1d.err

# create frag1dParaFile, parameters
blast_version=`$blastbin/blastpgp - | awk '/^blastpgp/{print $2}'`
trainingSet="2727 PDB chains"
blastdb_para="10,083,419 sequences"
frag1dParaFile=$TMPDIR/frag1dparafile.txt
echo "Training set: $trainingSet" >$frag1dParaFile
echo "blastpgp version: $blast_version" >>$frag1dParaFile
echo "blast nr database: $blastdb_para" >>$frag1dParaFile

OS=`uname -o`
case $OS in 
    Cygwin*)
        blastdbname=`cygpath -w "$blastdbname" `
        export BLASTMAT=`cygpath -w "$BLASTMAT"`
        ;;
esac

tmpSeqFileListFile0=$TMPDIR/query.0.seqfilelist
tmpPSSMFileListFile0=$TMPDIR/query.0.pssmfilelist

tmpPSSMFileListFile1=$TMPDIR/query.1.pssmfilelist

if [  "$seqFile" != "" ] ; then #the input is aa sequence
    if [ -s "$seqFile" ]; then
        numSeq=`grep "^>" $seqFile | wc -l`
        if [ $numSeq -gt 0 ]; then 
            echo "$numSeq sequences recognized in the input file $seqFile"
            $BINPATH/splitfasta.py $seqFile -outpath $TMPDIR/seq -q
            cntValidSeq=0
            for file in $(find $TMPDIR/seq -name "*.aa"); do
                v=`CheckSequence $file $SEQ_MODE_AASEQ`
                if [ $v -eq 0 ]; then 
                    ((cntValidSeq++))
                    echo "$file" >> $tmpSeqFileListFile0 
                fi
            done

            if [ $cntValidSeq -gt 0 ]; then

                if [ -f "$blastdbname.phr"  -o -f "$blastdbname.00.phr" ]; then 
                    echo "Using blastdb under current directory."
                    export BLASTDB=$PWD
                elif [ -f "$DATADIR/blastdb/$blastdbname.phr"  -o  -f "$DATADIR/blastdb/$blastdbname.00.phr" ]; then 
                    echo "Using blastdb at $DATADIR/data/blastdb."
                    export BLASTDB=$DATADIR/data/blastdb
                elif [ "$BLASTDB" != "" -a  \( -f "$BLASTDB/$blastdbname.phr" -o  -f "$BLASTDB/$blastdbname.00.phr" \) ]; then 
                    echo "Env BLASTDB is set and using blastdb at $BLASTDB/$blastdbname."
                else
                    echo "Fatal! Could not find blastdb $blastdbname" >&2 
                    exit
                fi

                echo "$cntValidSeq sequences satisfy the requirement and will be processed further..."
                echo "Running PSI-BLAST for $cntValidSeq sequences."
                echo "Be patient with PSI-BLAST for building sequence profiles..."
                echo
                mkdir -p $TMPDIR/psiblast
                for file in $(cat $tmpSeqFileListFile0); do
                    BuildProfileFromAASeq $file $TMPDIR/psiblast $tmpPSSMFileListFile1
                done
            else
                echo "No sequences satisfy the requirement. Ignore..."
            fi
        else
            echo "Sequence file $seqFile format error (not Fasta). Ignore!" >&2
        fi
    else
        echo "Sequence file $seqFile does not exists or is empty. Ignore!" >&2
    fi
fi

if [ "$pssmFile" != ""  ]; then
    if [ ! -s "$pssmFile"  ]; then
        echo "PSSM file $pssmFile does not exist or empty. Ignore."  >&2 
    else 
        echo "1 PSSM file \"$pssmFile\" recognized."
        echo $pssmFile >> $tmpPSSMFileListFile0
    fi
fi

if [ "$pssmFileListFile" != "" ] ; then
    if [ ! -s "$pssmFileListFile"  ]; then
        echo "PSSM list file $pssmFileListFile does not exist or empty. Ignore."  >&2 
    else 
        cntfile=0
        for file in $(cat $pssmFileListFile); do 
            echo $file >> $tmpPSSMFileListFile0
            ((cntfile++))
        done 
        echo "$cntfile PSSM files recognized from the list file \"$pssmFileListFile\"."
    fi
fi

if [ -s $tmpPSSMFileListFile0 ]; then
    for file in $(cat $tmpPSSMFileListFile0); do
        v=`CheckSequence $file $SEQ_MODE_PSSM` 
        cntValidSeq=0
        mkdir -p $TMPDIR/pssm
        if [ $v == 0 ]; then 
            ((cntValidSeq++))
            echo $file  >> $tmpPSSMFileListFile1 
        fi
    done
fi

NFile=`cat $tmpPSSMFileListFile1 | wc -l`
echo "Run Frag1D for $NFile sequence(s)"
echo
RunFrag1D $outpath $tmpPSSMFileListFile1 $frag1dParaFile

# clean temporary files
if [ "$isCleanUp" == "true" ]; then
    rm -rf $TMPDIR
fi
