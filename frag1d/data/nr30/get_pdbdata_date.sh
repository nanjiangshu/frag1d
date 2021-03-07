#!/bin/bash

# get the newest date of the pdb entry of the database
usage="
usage: $0 pdblistfile
"

if [ "$1" == "" ];then
    echo "$usage"
    exit 1
fi

pdblistfile=$1

basename=`basename $pdblistfile`
rootname=${basename%.*}

tmpdir=$(mktemp -d /tmp/tmpdir.get_pdbdata_date.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }


trap 'rm -rf $tmpdir' INT EXIT TERM

getpdbfilepath -l $pdblistfile | sort -u > $tmpdir/pdbfilelist.txt
(for file in $(cat $tmpdir/pdbfilelist.txt); do
    head -n 1 $file
done) > $tmpdir/pdb_head1.txt

awk '{print substr($0,51,9)}' $tmpdir/pdb_head1.txt | convertdate-ddmmmyy-to-yyyymmdd.awk | sort -rg |head -n 1  > $rootname.latest.date.txt

cp $tmpdir/pdb_head1.txt  $rootname.pdbhead1.txt

latest_date=$(cat $rootname.latest.date.txt)
echo "Latest date of PDB entry for $pdblistfile: $latest_date"

rm -rf $tmpdir
