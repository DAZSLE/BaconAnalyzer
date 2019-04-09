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
#include "JECLoader.hh"
#include "TRandom3.h"

using namespace baconhep;

class VJetLoader { 
public:
  VJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,int iN=1, bool iData=false, std::string iLabel="2017");
  ~VJetLoader();
  double correction(TJet &iJet,double iRho);
  void reset();
  void resetZprime();
  void resetDoubleB();
  void setupTree(TTree *iTree,std::string iJetLabel,bool iHWW=false);
  void setupTreeZprime(TTree *iTree,std::string iJetLabel);
  void load(int iEvent);
  void selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum, bool iHWW=false);
  void selectVJetsByDoubleB(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum);  
  void fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum, bool iHWW=false);
  void fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum);
  TAddJet *getAddJet(TJet *iJet);

  std::vector<double> fvSize, fvMatching,fRatioPt;
  std::vector<int> fisHadronicV, fisTightVJet, fisMatchedVJet;
  std::vector<int> fpartonFlavor, fhadronFlavor, fnbHadrons, fncHadrons, fnCharged, fnNeutrals, fnParticles;

  std::vector<TJet*> fLooseVJets;
  std::vector<TLorentzVector> selectedVJets;
  std::vector<TJet*> fLooseVJetsByDoubleB; 
  std::vector<TLorentzVector> selectedVJetsByDoubleB; 

  TRandom3* r;

protected: 
  TClonesArray *fVJets;
  TBranch      *fVJetBr;
  TClonesArray *fVAddJets;
  TBranch      *fVAddJetBr;

  TTree        *fTree;

  JECLoader    *fJEC;

  int           fNLooseVJets;
  int           fNTightVJets; 

  std::vector<double> fVars, fVarsZprime;
  std::vector<std::string> fLabels, fLabelsZprime;
  std::vector<std::string> fTrigString;
  std::string fYear;

  int  fN;
  bool isData;
  // Gaussian random numbers (one for each jet)
  std::vector<double> x1List;
};
#endif
