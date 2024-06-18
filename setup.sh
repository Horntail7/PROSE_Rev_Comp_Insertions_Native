#!/bin/bash

htname=2024-06-07_HT23
temp=Bp
bars=Bp
nsample=100
score_cut_off=0.95
ncores=40


dir=`pwd`
master='/data/PROSE_QC_Multi_fast_Master_Rev_Native_template'

if [ ! -d $htname ]
then
    mkdir $htname
else
    echo $htname "already exits"
    exit
fi
cd $htname

cp -r $master/Scripts .

if [ ! -d Inputs ]
then
    mkdir Inputs
fi


cp $master/Inputs/$temp.fasta ./Inputs/Templates.fasta
cp $master/Inputs/$bars.txt ./Inputs/barcodes.txt

sed "s|HTNAME|$htname|g" $master/snakefile.temp | sed "s|CUTOFF|$score_cut_off|g" | sed "s|NSAMPLE|$nsample|g" > snakefile
sed "s|NCORES|$ncores|g" $master/run.temp > run.sh

echo "Done"



