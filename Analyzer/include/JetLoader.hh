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

// Razor utils
#include "RazorUtils.hh"

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree, std::string iJet = "AK4Puppi");
  ~JetLoader();
  void reset();
  void setupTree(TTree *iTree, std::string iJetLabel);
  void setupTreeDiJet(TTree *iTree, std::string iJetLabel);
  void setupTreeRazor(TTree *iTree, std::string iJetLabel);
  void load(int iEvent);
  void selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets,float iMet, float iMetPhi,float iFMet,float iFMetPhi);
  void fillGoodJets(std::vector<TLorentzVector> iVJets,double dR, std::vector<const TJet*> &iGoodJets);
  std::vector<TJet*> fLooseJets;
  std::vector<const TJet*> fGoodJets, selectedJets15, selectedJets8;
  std::vector<TLorentzVector> selectedJets;

  //Fillers
  void addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals);
  void fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals);
  void addDijet(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals);
  void fillDiJet();
  void fillRazor(std::vector<TLorentzVector> iJets,float iMet, float iFMet);
  double correction(TJet &iJet,float iRho);

  const double CSVL = 0.460; // CSVv2 WP
  const double CSVM = 0.800;
  const double CSVT = 0.935;

  int           fNJets;
  int           fNJetsAbove80GeV;
  float         fMT;
  double        fMinDPhi;
  double        fMinDFPhi;
  float         fMR, fdeltaPhi;
  float         fHT, fMHT, fMT2;
  float         falphaT, fdPhiMin;
  float         fRsq; 

  enum RazorBox { //boxes for razor inclusive analysis
    MuEle = 0,  
    MuMu = 1,
    EleEle = 2,
    MuSixJet = 3,
    MuFourJet = 4,
    MuJet = 5,
    EleSixJet = 6,
    EleFourJet = 7,
    EleJet = 8,
    LooseLeptonSixJet = 9,
    LooseLeptonFourJet = 10, 
    SixJet = 11, 
    FourJet = 12, 
    LooseLeptonDiJet = 13, 
    DiJet = 14,         
    TwoBJet = 15, 
    OneBJet = 16, 
    ZeroBJet = 17, 
    MuMultiJet = 18, 
    EleMultiJet = 19, 
    LooseLeptonMultiJet = 20, 
    MultiJet = 21, 
    NONE = 999 
  };
  
protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;
  TTree        *fTree;
  int           fNBTagsL;
  int           fNBTagsM;
  int           fNBTagsT;

  int           fNJetsdR2;
  int           fNBTagsLdR2;
  int           fNBTagsMdR2;
  int           fNBTagsTdR2;

  int           fNJetsdR08;
  int           fNBTagsLdR08;
  int           fNBTagsMdR08;
  int           fNBTagsTdR08;

  int           fNJetsdR15;
  int           fNBTagsLdR15;
  int           fNBTagsMdR15;
  int           fNBTagsTdR15;

  int           fN;
  std::vector<double>      fVars;
  FactorizedJetCorrector   *fJetCorr;
};
#endif
