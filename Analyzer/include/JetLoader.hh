#ifndef JetLoader_H
#define JetLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "MonoXUtils.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "SimpleJetResolution.h"
#include "TRandom3.h"

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree);
  ~JetLoader();
  void reset();
  void setupTree(TTree *iTree, std::string iJetLabel);
  void load(int iEvent);
  void selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets, double iRho);
  std::vector<TJet*> fLooseJets;
  std::vector<const TJet*> fGoodJets, selectedJets15, selectedJets8;
  std::vector<TLorentzVector> selectedJets;

  //Fillers
  void fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho = 0);
  void addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals);
  void fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, std::vector<TLorentzVector> iVJets, double iRho);

  
  // JEC tools
  std::vector<FactorizedJetCorrector*> getJetCorrector() { return JetCorrector; }
  std::vector<std::pair<int,int> > getJetCorrectionsIOV() { return JetCorrectionsIOV; }
  std::vector<SimpleJetResolution*> getJetResolutionCalculator() { return JetResolutionCalculator; }
  double getJecUnc( float pt, float eta, int run );
  double getJerSF( float eta, int nsigma);
  TRandom3* r;

  const double CSVL = 0.460; // CSVv2 WP
  const double CSVM = 0.800;
  const double CSVT = 0.935;

  int           fNJetsPt30;

  int           fNVars;
  int           fNOtherVars;
protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;
  TClonesArray *fJetsCHS;
  TBranch      *fJetBrCHS;

  TTree        *fTree;
  int           fNFwdPt30;
  int           fNBTagsLPt30;
  int           fNBTagsMPt30;
  int           fNBTagsTPt30;

  std::vector<int>           fNJetsPt30dR08;
  std::vector<int>           fNBTagsLPt50dR08;
  std::vector<int>           fNBTagsMPt50dR08;
  std::vector<int>           fNBTagsTPt50dR08;
  std::vector<int>           fNBTagsLPt100dR08;
  std::vector<int>           fNBTagsMPt100dR08;
  std::vector<int>           fNBTagsTPt100dR08;
  std::vector<int>           fNBTagsLPt150dR08;
  std::vector<int>           fNBTagsMPt150dR08;
  std::vector<int>           fNBTagsTPt150dR08;

  int           fN;
  int           fNV;
  
  std::vector<double>      fVars;
  FactorizedJetCorrector   *fJetCorr;

  std::string cmsswPath;
  
  // for jet energy corrections
  void loadCMSSWPath();
  void loadJECs(bool isData);
  void loadJECs_Rereco(bool isData);
  std::vector<std::vector<JetCorrectorParameters> > correctionParameters;
  std::vector<JetCorrectorParameters*> JetResolutionParameters;
  std::vector<FactorizedJetCorrector*> JetCorrector;
  std::vector<JetCorrectionUncertainty*> jecUnc;
  std::vector<SimpleJetResolution*> JetResolutionCalculator;
  std::vector<std::pair<int,int> > JetCorrectionsIOV;

};
#endif
