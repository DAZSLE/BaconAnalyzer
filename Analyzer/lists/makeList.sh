#!/bin/bash
# Usage: source makeList.sh [ListOfDirectories]. Eg: source makeList.sh 11.txt

eos='eos root://cmseos.fnal.gov'
#eos=''

function print_to_file {
    local printdir=$1
    local writefile=$2
    xrdfs root://cmseos.fnal.gov ls -u $printdir | grep '\.root' >> $writefile
    #find ../eos/cms/$printdir -maxdepth 1 -type f -name "*.root" >> $writefile
}

function scan {
    local scandir=$1
    local outfile=$2
    written=0
    nfiles=`$eos ls $scandir/*.root | wc -l`
    echo Scanning: $scandir

    if [ "$nfiles" -eq 0 ]
    then
        for x in `$eos ls $scandir |  awk '{print $1}'`
        do
            hasfailed=`echo $x | grep -c "failed"`
            hasdot=`echo $x | grep -c "\."`
            # skip the failed directory
            if [ "$hasfailed" -eq 1 ]
            then
	            continue
            else
                # recursive call to subdir if there is no root file in the current directory
                if [ "$hasdot" -eq 0 ]
                then
                    scan  $scandir/$x $outfile
                fi
            fi
        done
    else
        # print to the list if there are root files in the directory  
        print_to_file $scandir $outfile
        echo Print $scandir to $outfile
        written=1
    fi
}

# Get the version number from input file. Eg: 11.txt -> version = 11
version=`echo $1 | sed 's/\.txt//'`
echo $version
if [ ! -d "production${version}" ]
then
    echo "Making production$version"
    mkdir production$version
fi

while read -r line
do
    if [ -f "production${version}/${line}.txt" ]
    then
        rm production${version}/${line}.txt
    fi

    #scan  /store/cmst3/group/monojet/production/${version}/$line production${version}/${line}.txt
    scan  /store/group/lpcbacon/${version}/$line production${version}/${line}.txt    
    if [ "$written" == 1 ]
    then
        #sed -i 's/\.\.\/eos\/cms/root:\/\/eoscms.cern.ch/g' production${version}/${line}.txt 
        sed -i 's/\.\.\/eos\/cms/root:\/\/cmseos.fnal.gov/g' production${version}/${line}.txt
    fi
done < $1
