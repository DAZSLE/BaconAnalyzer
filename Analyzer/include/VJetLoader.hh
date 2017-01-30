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
#include "Utils.hh"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// B-tag scale factors
#include "SJBTagUnc.hh"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "BTagCalibrationStandalone.h"

#include "TRandom3.h"

using namespace baconhep;

class VJetLoader { 
public:
  VJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS="AK8CHS",std::string iAddJetCHS="AddAK8CHS",int iN=1, bool iData=false);
  ~VJetLoader();
  double correction(TJet &iJet,double iRho);
  void reset();
  void resetZprime();
  void resetCHS();
  void resetDoubleB();
  void setupTree(TTree *iTree,std::string iJetLabel);
  void setupTreeZprime(TTree *iTree,std::string iJetLabel);
  void setupTreeCHS(TTree *iTree,std::string iJetLabel);
  void load(int iEvent);
  void selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum);
  void selectVJetsByDoubleBCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum);  
  void selectVJetsCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum);
  void fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum);
  int getMatchedCHSJetIndex(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR);
  void matchJet(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex);
  void matchJet15(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR);
  void fillVJetCHS(TJet *iJet, int jIndex);
  void fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum);
  TAddJet *getAddJet(TJet *iJet);
  TAddJet *getAddJetCHS(TJet *iJet);

  double fvSize, fvMatching,fRatioPt;
  int fisHadronicV;
  std::vector<int> fisTightVJet, fisTightVJetCHS;
  int fpartonFlavor, fhadronFlavor, fnbHadrons, fncHadrons, fnCharged, fnNeutrals, fnParticles;

  std::vector<TJet*> fLooseVJets, fLooseVJetsCHS;
  std::vector<TLorentzVector> selectedVJets, selectedVJetsCHS;
  std::vector<double> fdoublecsvCHS, fdoublesubCHS, fptCHS, fetaCHS, fphiCHS;
  std::vector<TJet*> fLooseVJetsByDoubleB, fLooseVJetsCHSByDoubleB;
  std::vector<TLorentzVector> selectedVJetsByDoubleB, selectedVJetsCHSByDoubleB;


  const double CSVL = 0.460; // CSVv2SubJet WP 
  const double CSVM = 0.800;
  // JEC tools
  std::vector<FactorizedJetCorrector*> getJetCorrector() { return JetCorrector; }
  std::vector<std::pair<int,int> > getJetCorrectionsIOV() { return JetCorrectionsIOV; }
  double getJecUnc( float pt, float eta, int run );
  double JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
				    double rho, double jetArea,
				    int run,
				    std::vector<std::pair<int,int> > JetCorrectionsIOV,
				    std::vector<FactorizedJetCorrector*> jetcorrector,  
				    int jetCorrectionLevel = -1,
				    bool printDebug = false);
  double JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
				    double rho, double jetArea,					  
				    FactorizedJetCorrector* jetcorrector,  
				    int jetCorrectionLevel = -1,
				    bool printDebug = false);
  TRandom3* r;

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

  std::vector<double> fVars, fVarsZprime;
  std::vector<std::string> fLabels, fLabelsZprime;
  std::vector<std::string> fTrigString;

  int           fN;
  
  std::string cmsswPath;
  
  // for jet energy corrections
  bool isData;
  void loadCMSSWPath();
  void loadJECs(bool isData);
  void loadJECs_Rereco(bool isData);
  std::vector<std::vector<JetCorrectorParameters> > correctionParameters;
  std::vector<FactorizedJetCorrector*> JetCorrector;
  std::vector<JetCorrectionUncertainty*> jecUnc;
  std::vector<std::pair<int,int> > JetCorrectionsIOV;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor resolution_sf;
  // Gaussian random numbers (one for each jet)
  std::vector<double> x1List;
  std::vector<double> x2List;
  std::vector<double> x3List;  
};
#endif
