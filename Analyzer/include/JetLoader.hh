#ifndef JetLoader_H
#define JetLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "Utils.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JECLoader.hh"
#include "TRandom3.h"

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree,  bool iData=false,  std::string iLabel="2017");
  ~JetLoader();
  void reset();
  void setupTree(TTree *iTree, std::string iJetLabel);
  void load(int iEvent);
  void selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets, double iRho, unsigned int runNum);
  std::vector<TJet*> fTightJets;
  std::vector<const TJet*> fGoodJets, selectedJets15, selectedJets8;
  std::vector<TLorentzVector> selectedJets;

  //Fillers
  void fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum);
  void addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals);
  void fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, std::vector<TLorentzVector> iVJets, double iRho, unsigned int runNum);
  
  TRandom3* r;

  // default is 2017
  double DEEPCSVL;
  double DEEPCSVM;
  double DEEPCSVT;

  int           fNJetsPt30;
  int           fNJetsPt30jesUp;
  int           fNJetsPt30jesDown;
  int           fNJetsPt30jerUp;
  int           fNJetsPt30jerDown;
  double        MetXCorrjesUp;
  double        MetYCorrjesUp;
  double        MetXCorrjesDown;
  double        MetYCorrjesDown;
  double        MetXCorrjerUp;
  double        MetYCorrjerUp;
  double        MetXCorrjerDown;
  double        MetYCorrjerDown;

  int           fNVars;
  int           fNOtherVars; 
protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;
  TClonesArray *fJetsCHS;
  TBranch      *fJetBrCHS;

  TTree        *fTree;

  JECLoader    *fJEC;

  int           fNFwdPt30;
  int           fNBTagsLPt30;
  int           fNBTagsMPt30;
  int           fNBTagsTPt30;

  std::vector<int>           fNJetsPt30dR08;
  std::vector<int>           fNJetsPt30dR08jesUp;
  std::vector<int>           fNJetsPt30dR08jesDown;
  std::vector<int>           fNJetsPt30dR08jerDown;
  std::vector<int>           fNJetsPt30dR08jerUp;
  std::vector<int>           fNBTagsLPt50dR08;
  std::vector<int>           fNBTagsMPt50dR08;
  std::vector<int>           fNBTagsTPt50dR08;
  std::vector<int>           fNBTagsLPt100dR08;
  std::vector<int>           fNBTagsMPt100dR08;
  std::vector<int>           fNBTagsTPt100dR08;
  std::vector<int>           fNBTagsLPt150dR08;
  std::vector<int>           fNBTagsMPt150dR08;
  std::vector<int>           fNBTagsTPt150dR08;

  int  fN;
  int  fNV;
  bool isData;
  std::vector<double>      fVars;

  // Gaussian random numbers (one for each jet)
  std::vector<double> x1List;

  // b-tag WP
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
  double fdeepCSVL2016 = 0.2217;
  double fdeepCSVM2016 = 0.6321;
  double fdeepCSVT2016 = 0.8953;

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
  double fdeepCSVL2017 = 0.1522;
  double fdeepCSVM2017 = 0.4941;
  double fdeepCSVT2017 = 0.8001;

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
  double fdeepCSVL2018 = 0.1241;
  double fdeepCSVM2018 = 0.4184;
  double fdeepCSVT2018 = 0.7527;

  std::string label2016 = "2016";
  std::string label2017 = "2017";
  std::string label2018 = "2018";

};
#endif
