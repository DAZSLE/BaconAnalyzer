#!/bin/bash

#argv[1] => input bacon file name                                                                                                                                                          
#argv[2] => dataset type: "mc", "mczbb", "mczcc", "mcwcs", "mcvlf", "data"                                                                                                              
#argv[3] => JSON file for run-lumi filtering of data, specify "none" for MC or no filtering                                                                                             
#argv[4] => cross section (pb), ignored for data                                                                                                                                       
#argv[5] => total weight, ignored for data       

ismc=$1
dir=$2
label=$3
xs=$4
xdir=$5
exec=$6
option1=$7
option2=$8

options1="Output.root --passSumEntries 5:Events  -a 2:"${ismc}" -a 3:none  -n 7000 -q 2nw4cores"

if [ $ismc = "data" ]; then
     options1="Output.root -a 5:1  -a 2:"${ismc}" -a 3:$PWD/../data/Cert_271036-282037_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt -n 7000 -q 1nh"
fi

scandir=$dir$label
submitted=0
total=0

function scan { 
    local scandir=$1
    echo $scandir
    submitted=0
    for x in `ls $scandir |  awk '{print $1}'`; do
        hasfailed=`echo $scandir | grep -c "failed"`
	if [ "$hasfailed" -eq 1 ]; then
	    break
	fi
        hasdot=`echo $x | grep -c "\."`
        if [ "$hasdot" -eq 0 ]; then
            echo "Scanning: "$scandir/$x
            scan  $scandir/$x
        fi
        hasroot=`echo $x | grep -c .root`
        if [ "$hasroot" -ne 0 ]; then
            if [ "$submitted" -eq 0 ]; then
		python baconBatch.py $exec  $options1 -a "4:"$xs -d "1:"$scandir  -o $PWD/../${xdir}/${label}_${ismc}   $option1 $option2
		submitted=$((submitted+1))
		total=$((total+1))
            fi
        fi
    done
}

for x in `eos ls $scandir |  awk '{print $scandir}'`; do
    hasrootdir=`echo $x | grep -c .root`
    if [ "$hasrootdir" -eq 1 ]; then
	break
    fi
done
if [ "$hasrootdir" -ne 0 ]; then
    submitted=0	
    if [ "$submitted" -eq 0 ]; then
	python baconBatch.py $exec  $options1 -a "4:"$xs -d "1:"$scandir  -o $PWD/../${xdir}/${label}_${ismc}   $option1 $option2
	submitted=$((submitted+1))
        total=$((total+1))	
    fi
else
    scan $scandir
fi
exit
