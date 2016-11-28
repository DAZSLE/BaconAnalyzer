#ifndef VJetLoader_H
#define VJetLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "MonoXUtils.hh"

// B-tag scale factors
#include "SJBTagUnc.hh"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "BTagCalibrationStandalone.h"

using namespace baconhep;

class VJetLoader { 
public:
  VJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS="AK8CHS",std::string iAddJetCHS="AddAK8CHS",int iN=1,std::string subjetbtagScaleFactorFilename ="/afs/cern.ch/work/c/cmantill/public/Bacon/CMSSW_8_0_10/src/BaconAnalyzer/Analyzer/data/subjet_CSVv2_ichep.csv");
  ~VJetLoader();
  double correction(TJet &iJet,double iRho);
  void reset();
  void resetZprime();
  void resetCHS();
  void resetDoubleB();
  void resetMonoX();
  void resetSubJetBTag();
  void setupTree(TTree *iTree,std::string iJetLabel);
  void setupTreeZprime(TTree *iTree,std::string iJetLabel);
  void setupTreeCHS(TTree *iTree,std::string iJetLabel);
  void setupTreeMonoX(TTree *iTree,std::string iJetLabel);
  void setupTreeSubJetBTag(TTree *iTree, std::string iJetLabel);
  void load(int iEvent);
  void selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho);
  void selectVJetsByDoubleBCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho);  
  void selectVJetsCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho);
  void fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho = 0);
  int getMatchedCHSJetIndex(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR);
  void matchJet(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex);
  void matchJet15(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex);
  void fillVJetCHS(TJet *iJet, int jIndex);
  TAddJet *getAddJet(TJet *iJet);
  TAddJet *getAddJetCHS(TJet *iJet);
  void addSubJetBTag(std::string iHeader,TTree *iTree,std::string iLabel,std::vector<std::string> &iLabels,int iN,std::vector<float> &iVals);
  void fillSubJetBTag(const TClonesArray* iGens, std::vector<TLorentzVector> iObjects);

  
  double dPhi(TLorentzVector v1, TLorentzVector v2, TLorentzVector v3);
  
  
  double ftopSize, ftopMatching, fvSize, fvMatching,fRatioPt;
  int fisHadronicTop, fisHadronicV;
  int fisTightVJet, fisTightVJetCHS;
  int fpartonFlavor, fhadronFlavor, fnbHadrons, fncHadrons, fnCharged, fnNeutrals, fnParticles;
  float fVMT,fdR_sj0dR,fdPhi_sj0dPhi, fdPhiJRF_sj0dPhiJRF;

  std::vector<TJet*> fLooseVJets, fLooseVJetsCHS;
  std::vector<TLorentzVector> fGoodVSubJets, selectedVJets, fGoodVSubJetsCHS, selectedVJetsCHS;
  std::vector<double> fdoublecsvCHS, fdoublesubCHS, fptCHS, fetaCHS, fphiCHS;
  std::vector<TJet*> fLooseVJetsByDoubleB, fLooseVJetsCHSByDoubleB;
  std::vector<TLorentzVector> selectedVJetsByDoubleB, selectedVJetsCHSByDoubleB;

  const double CSVL = 0.460; // CSVv2SubJet WP 
  const double CSVM = 0.800;

protected: 
  TClonesArray *fVJets;
  TBranch      *fVJetBr;
  TClonesArray *fVAddJets;
  TBranch      *fVAddJetBr;
  TClonesArray *fVJetsCHS;
  TBranch      *fVJetBrCHS;
  TClonesArray *fVAddJetsCHS;
  TBranch      *fVAddJetBrCHS;

  TTree        *fTree;

  int           fNLooseVJets, fNLooseVJetsCHS;
  int           fNTightVJets, fNTightVJetsCHS;

  std::vector<std::string> fBtagLabels = {"CENT", "MISTAGUP","MISTAGDO","BTAGUP","BTAGDO"};  // nominal, mistagup, mistagdown, btagup and btagdown
  std::vector<std::string> measurementTypes = {"lt", "incl"};                                // measurements type
  std::vector<std::string> variationTypes = {"central", "up", "down"};                       // systematics type 
  std::vector<std::string> flavorTypes = {"Ms", "Bs"};                                       // nominal, mistag and btag
  std::vector<std::string> wpTypes = {"L","M"};                                              // working points 
  std::vector<double> fVars, fVarsZprime;
  std::vector<float>  fSubJetBTagVars;
  std::vector<std::string> fLabels, fLabelsZprime;
  std::vector<std::string> fTrigString;

  BTagCalibration          *fSubJetCalib;
  std::vector<BTagCalibrationReader*> fSubJetreadersL;
  std::vector<BTagCalibrationReader*> fSubJetreadersM;
  std::vector<std::vector<BTagCalibrationReader*>> fSubJetreaders;

  int           fN;
};
#endif
