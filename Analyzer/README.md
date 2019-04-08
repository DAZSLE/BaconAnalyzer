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
 * For tag:14 of BaconProd/BaconAna, CMSSW_VERSION = CMSSW_10_2_6
```
cmsrel CMSSW_10_2_6
cd CMSSW_10_2_6/src
cmsenv
git clone https://github.com/ksung25/BaconProd
git clone https://github.com/ksung25/BaconAna
git clone https://github.com/DAZSLE/BaconAnalyzer
git fetch
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

In general, you can define a new object and include it in bin/run$YOURANALYSIS.cpp, e.g. for VJets:
```
VJetLoader      *fVJet15CHS     = 0;
fVJet15CHS    = new VJetLoader    (lTree,"CA15CHS","AddCA15CHS");
fVJet15CHS   ->setupTree      (lOut,"CHSjet");
```

IMPORTANT
```
After modifications compile before running.
Also, if you are submitting to condor, make sure you re-tar your environment after compiling.
```

Baconbits production
-----------

For e.g. `Zprime` samples submission.

Define samples in `production/submitZprime.py`

List all samples in `$PROD.txt` e.g. `ls $EOSDIR > $PROD.txt`.

Run `makeList.sh` on `$PROD` to list bacon files: `bash makeList.sh $PROD`. This will produce a `production${PROD}` folder inside list.

Run `submitZprime.py` for a given SAMPLE and TAG as following:

```
python submitZprime.py -s SAMPLE -t TAG (to create submission files)
python submitZprime.py -s SAMPLE -t TAG --monitor sub # (to submit)
```

LPC instructions
-----------
Since 90X BaconAnalyzer is configured to run @ CMSLPC cluster. It should also work on lxplus but you need to modify eosdir input and output...

EVERYTIME after compiling: re-tar `CMSSW environment` and `data/` dir

For v15 (BaconProd,BaconAna): `CMSMSW_VERSION = CMSSW_10_2_6`

```
cd ../../
tar --exclude-caches-all --exclude-vcs --exclude-caches-all --exclude-vcs -cvzf CMSSW_10_2_6.tgz CMSSW_10_2_6 --exclude=src --exclude=tmp --exclude="*.scram" --exclude="*.SCRAM"
cd CMSSW_10_2_6/src/BaconAnalyzer/Analyzer/
tar -zcvf data.tgz data
```

Those are transferred to condor. Follow same instructions to produce submission files.

Also, DO NOT FORGET TO initialize your proxy before submitting jobs.
```
voms-proxy-init --voms cms --valid 168:00
```

To update lists (i.e. if BaconProd changes or new samples are added)
-----------
In cmslpc:

```
(Before doing cmsenv)
cd lists/
eosls /store/group/lpcbacon/${PROD} > ${PROD}.txt
bash makeList.sh ${PROD}
```

Also, add the new samples to `submitZprime.py`.

To produce skim
-----------
Skim is a version of baconbits with usually a simple pT cut e.g. 350 GeV. In cmslpc do:

```
e.g. EOSDIR = /eos/uscms//store/user/lpcbacon/dazsle/zprimebits-v15.01/
python skim.py -s SAMPLE -i EOSDIR --files-to-hadd NFILES
```

Skimmed files should appear on sub-dir of baconbits.