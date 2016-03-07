#!/bin/bash
# 2009-11-08
# build accfrag profile

usage="
Usage:  run-database-build.sh argument-list
  test argument passing
Options:
  -l          file: set the idListFile
  -fragacc  path  : outpath for the fragacc matrices
  -frag     path  : outpath 
  -triple file    : supply shape triple frequency file
  -qij            : qij path
  -workdir        : specify the working dir
  -h|--help       : print this help message and exit
Created 2009-07-27, updated 2009-11-08, Nanjiang
"
function PrintHelp()
{
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

qijformat=1  #format $id.Qij
qijpath=data/qijmatrix
tripleFile=data/chain5860.triple.freq
result=data/dbFrag
resultacc=data/accfrag
idListFile=

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        idList="$idList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -qij|--qij) qijpath=$2;shift;;
            -fragacc|--fragacc) resultacc=$2;shift;;
            -frag|--frag) result=$2;shift;;
            -workdir|--workdir) WORKINGDIR=$2;shift;;
            -triple|--triple) tripleFile=$2;shift;;
            -l) idListFile=$2;shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        echo "Error! Wrong argument: $1";
        exit
    fi
    shift
done


if [ ! -d "$qijpath" ]; then 
    echo "Error! qijpath = \"$qijpath\" does not exist. Exit..."
    exit
fi
if [ ! -f "$tripleFile" ]; then 
    echo "Error! tripleFile = \"$tripleFile\" does not exist. Exit..."
    exit
fi
if [ ! -f "$idListFile" ]; then 
    echo "Error! idListFile = \"$idListFile\" does not exist. Exit..."
    exit
fi

if [  "$result"  == "" ]; then 
    echo "Error! resultpath = \"$result\" does not set. Exit..."
    exit
fi
if [ "$resultacc" == "" ]; then 
    echo "Error! resultaccpath = \"$resultacc\" does not set. Exit..."
    exit
fi

echo "database_build --rb --wb --triple  $tripleFile --qijpath $qijpath --qijformat $qijformat --result $result --resultacc $resultacc $idListFile "
database_build --rb --wb --triple  $tripleFile --qijpath $qijpath --qijformat $qijformat --result $result --resultacc $resultacc $idListFile 
