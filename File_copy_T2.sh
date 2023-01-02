#!/bin/bash
FILE=$1

export src_dir="/dpm/indiacms.res.in/home/cms/store/user/abala/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/221224_084643/0001/"

export t2_dir="root://se01.indiacms.res.in//store/user/abala/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/221224_084643/0001/"

numFile=`eval rfdir $src_dir | awk '{print $9}' | grep "rootuple_" -wc`
#numFile=`eval rfdir $src_dir | awk '{print $9}' | grep "output_" -wc`
echo $numFile Files

for file in `rfdir $src_dir | awk '{print $9}' | grep \\rootuple_`; do
	echo $file
	echo $t2_dir$file >> ${FILE}
#	rfrm $src_dir$file
done

#ls $dest_dir | grep "root" -wc	 
