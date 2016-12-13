#ifndef BTagWeightLoader_H
#define BTagWeightLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "MonoXUtils.hh"

// B-tag scale factors
#include "BTagUnc.hh"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "BTagCalibrationStandalone.h"

using namespace baconhep;

class BTagWeightLoader { 
public:
  BTagWeightLoader(TTree *iTree,std::string btagScaleFactorFilename = "/afs/cern.ch/work/c/cmantill/public/Bacon/CMSSW_8_0_20/src/BaconAnalyzer/Analyzer/data/CSVv2_ichep.csv");
  ~BTagWeightLoader();
  void reset();
  void setupTree(TTree *iTree, std::string iJetLabel);
  void addBTag(std::string iHeader,TTree *iTree,std::string iLabel,std::vector<std::string> iLabels,int iN,std::vector<float> &iVals);
  void fillBTag(std::vector<const TJet*> iObjects);
  
  const double CSVL = 0.460; // CSVv2 WP
  const double CSVM = 0.800;
  const double CSVT = 0.935;
  
protected:
  std::vector<std::string> fLabels = {"CENT", "MISTAGUP","MISTAGDO","BTAGUP","BTAGDO"};   // nominal, mistagup, mistagdown, btagup and btagdown
  std::vector<std::string> measurementTypes = {"comb", "incl"};                           // measurements type
  std::vector<std::string> variationTypes = {"central", "up", "down"};                    // systematics type
  std::vector<std::string> flavorTypes = {"Ms", "Bs"};                                    // nominal, mistag and btag
  std::vector<std::string> wpTypes = {"L","M","T"};                                       // working points                           
  std::vector<float>      fVars;
  BTagCalibration          *fJetCalib;
  std::vector<BTagCalibrationReader*> freadersL;
  std::vector<BTagCalibrationReader*> freadersM;
  std::vector<BTagCalibrationReader*> freadersT;
  std::vector<std::vector<BTagCalibrationReader*>> freaders;
  
};
#endif
