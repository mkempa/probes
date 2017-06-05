#!/bin/bash

#Script simply selects those sequences and their names from reference which are 
#present in names parameter
names="$1"
reference="$2"
outfile="$3"

if [ ! -f $names ];
then
   	echo "File $names does not exist."
	exit
fi

if [ ! -f $reference ];
then
   	echo "File $reference does not exist."
	exit
fi

for l in $( cat $names ); do
	echo $( grep -A1 $l $reference ) >> $outfile
done
