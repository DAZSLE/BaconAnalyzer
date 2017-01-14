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

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree);
  ~JetLoader();
  void reset();
  void setupTree(TTree *iTree, std::string iJetLabel);
  void load(int iEvent);
  void selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets);
  std::vector<TJet*> fLooseJets;
  std::vector<const TJet*> fGoodJets, selectedJets15, selectedJets8;
  std::vector<TLorentzVector> selectedJets;

  //Fillers
  void addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals);
  void fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, std::vector<TLorentzVector> iVJets);

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
};
#endif
