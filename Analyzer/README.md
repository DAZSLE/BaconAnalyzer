# Analyzer

 * Depends on BaconAna and BaconProd packages
 * Place BaconAnalyzer in `$CMSSW_BASE/src` area

Setup
-------------
 * Setup $CMSSW_BASE area
 * Setup extra packages
 * Pull new b-tagger 
 * Setup BaconProd, BaconAna, and add the missing python module
 * Setup BaconAnalyzer, Development Packages
 * Compile
```
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv
git cms-addpkg DataFormats/BTauReco
git cms-addpkg PhysicsTools/PatAlgos
git cms-addpkg PhysicsTools/PatUtils
git cms-addpkg RecoBTag/SecondaryVertex
git cms-addpkg RecoBTau/JetTagComputer
git cms-addpkg RecoMET/METPUSubtraction
git pull official-cmssw pull/16224/head
git clone https://github.com/jmduarte/ShowerDeconstruction 
git clone https://github.com/ksung25/BaconProd
git clone https://github.com/ksung25/BaconAna
cd BaconAna
git remote add jmgd https://github.com/jmduarte/BaconAna
git fetch jmgd
git merge jmgd/add_python
cd ..
git clone https://github.com/DAZSLE/BaconAnalyzer
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
1) Define list of samples on production/submit$YOURANALYSIS.sh along with xsec

2) Make output directory

   `mkdir $YOURANALYSISbits` e.g. `mkdir zprimebits`

3) Make configuration files for the $OPTION samples (or All) specified in `production/submit$YOURANALYSIS.sh`
```
cd production/
./submit$YOURANALYSIS.sh $OPTION
```

4) After compiling, submit jobs to Batch as 
```
./submit$YOURANALYSIS.sh $OPTION --monitor sub
```

5) When production is done combine files using
```
./combine$YOURANALYSIS.sh $OPTION
```

Baconbits are stored in $YOURANALYSISbits/*.root
