# Analyzer

 * Depends on BaconAna and BaconProd packages
 * Place BaconAnalyzer in `$CMSSW_BASE/src` area

Setup
-------------
 * Setup $CMSSW_BASE area
 * Setup BaconProd, BaconAna, and add the missing python module
	
    `git clone https://github.com/ksung25/BaconProd.git`

    `git clone https://github.com/ksung25/BaconAna.git`

    `cp -r /afs/cern.ch/work/m/mcremone/public/CMSSW_7_6_2/src/BaconAna/Utils/python BaconAna/Utils/`

 * Setup BaconAnalyzer, Development Packages

    `git clone https://github.com/DAZSLE/BaconAnalyzer.git`

    `cd BaconAnalyzer`

   
 * Compile

   scram b -j 10

Define Analysis
----------
Define $YOURANALYSIS in bin/run$YOURANALYSIS.cpp - don't forget to include it in bin/BuildFile.xml

Modify it according to your analysis preselection.

Triggers are listed in 	    
	 `/src/BaconAna/DataFormats/data/HLTFile_25ns`	

For other types of Jets  e.g:
    `VJetLoader      *fVJet15CHS     = 0;`
    `fVJet15CHS    = new VJetLoader    (lTree,"CA15CHS","AddCA15CHS");`
    `fVJet15CHS   ->setupTree      (lOut,"bst15_CHSjet");`

After modifications compile before running.

Baconbits production
-----------
1) Define list of samples on production/submit$YOURANALYSIS.sh along with xsec

2) Make output directory

   `mkdir $YOURANALYSISbits`

3) Make configuration files for the $OPTION samples (or All)

   `cd production/`
   
   `./submit$YOURANALYSIS.sh $OPTION`

4) After compiling, submit jobs to Batch as 
   `./submit$YOURANALYSIS.sh $OPTION --monitor sub`

5) When production is done combine files using
   `./combine$YOURANALYSIS.sh $OPTION`

Baconbits are stored in $YOURANALYSISbits/*.root

