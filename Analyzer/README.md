# Analyzer

 * Depends on BaconAna and BaconProd packages
 * Place BaconAnalyzer in `$CMSSW_BASE/src` area

Setup
-------------
 * Setup $CMSSW_BASE area
 * Setup extra packages
 * Setup BaconProd, BaconAna, and add the missing python module
 * Setup BaconAnalyzer, Development Packages
 * Compile
 * For tag:14 of BaconProd/BaconAna, CMSSW_VERSION = CMSSW_9_4_7
```
cmsrel CMSSW_9_4_7
cd CMSSW_9_4_7/src
cmsenv
git clone https://github.com/ksung25/BaconProd
git clone https://github.com/ksung25/BaconAna
cd BaconAna
git remote add jmgd https://github.com/jmduarte/BaconAna
git fetch jmgd
git merge jmgd/add_python
cd ..
git clone https://github.com/DAZSLE/BaconAnalyzer
git fetch
git checkout -b 90x origin/90x
scram b clean
scram b -j 10
cd BaconAnalyzer
```

Define Analysis
----------
Define $YOURANALYSIS in bin/run$YOURANALYSIS.cpp - don't forget to include it in bin/BuildFile.xml

Modify it according to your analysis preselection.

Triggers are listed in 	    
```
/src/BaconAna/DataFormats/data/HLTFile_25ns
```

For other types of Jets  e.g:
```
    VJetLoader      *fVJet15CHS     = 0;
    fVJet15CHS    = new VJetLoader    (lTree,"CA15CHS","AddCA15CHS");
    fVJet15CHS   ->setupTree      (lOut,"bst15_CHSjet");
```

After modifications compile before running.

Baconbits production
-----------
Define samples on submitZprime.py

Run makeList.sh on $PROD.txt to list bacon files.

Run submitZprime.py for a given SAMPLE and TAG as following:

```
python submitZprime.py -s SAMPLE -t TAG --monitor sub # (to submit)
python submitZprime.py -s SAMPLE -t TAG --monitor check # (to check status - also bjobs)
python submitZprime.py -s SAMPLE -t TAG --monitor resub # (to resubmit)
```

Some temporary instructions for cmslpc
-----------
EVERYTIME after compiling: re-tar CMSSW environment and data/ dir

For v14 (BaconProd,BaconAna)
CMSMSW_VERSION = CMSSW_9_4_7

```
cd ../../
tar --exclude-caches-all --exclude-vcs --exclude-caches-all --exclude-vcs -cvzf CMSSW_9_4_7.tgz CMSSW_9_4_7 --exclude=src --exclude=tmp
cd CMSSW_9_4_7/src/BaconAnalyzer/Analyzer/
tar -zcvf data.tgz data
```

Those are transferred to condor. Follow same instructions to produce submission files.

Also, initialize your proxy before submitting jobs.
```
voms-proxy-init --voms cms --valid 168:00
```

To update lists
-----------
In cmslpc:

```
(Before doing cmsenv)
cd lists/
eosls /store/group/lpcbacon/${PROD} > ${PROD}.txt
bash makeList.sh ${PROD}.txt
```
