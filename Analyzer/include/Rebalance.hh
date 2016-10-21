#ifndef Rebalance_H
#define Rebalance_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "MonoXUtils.hh"
#include "RBZero.hh"

using namespace baconhep;

class Rebalance { 
public:
  Rebalance(TTree *iTree,bool isData);
  ~Rebalance();
  void reset();
  void setupTree(TTree *iTree);
  void load(int iEvent);
  bool selectJets(float &iMet,float &iMetPhi,std::vector<TLorentzVector> &iVetoes);
  void fillJets(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals,std::vector<TLorentzVector> &iVec,TLorentzVector &iMet);
  float dPhi(TLorentzVector &iMet,std::vector<TLorentzVector> iJets,double &iId,int iN=0);
  TLorentzVector met(std::vector<TLorentzVector> iJets,int iN=0);
  std::vector<TJet*> fSelJets;

protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;
  TTree        *fTree;

  double fHT       ;
  double fSHT      ;
  double fSHTD     ;
  double fSHTU     ;
  double fNJets    ;
  double fNFail    ;
  double fNBTag    ;
  double fChi2     ; 
  double fProb     ; 
  double fMet      ;
  double fZMet     ;
  double fSMet     ;
  double fSMetU    ;
  double fSMetD    ;
  double fMetPhi   ;
  double fZMetPhi  ;
  double fSMetPhi  ;
  double fSMetUPhi ;
  double fSMetDPhi ;
  double fMetDPhi  ;
  double fZMetDPhi ;
  double fSMetMDPhi ;
  double fSMetUDPhi ;
  double fSMetDDPhi ;

  double fIMetDPhi  ;
  double fIZMetDPhi ;
  double fISMetMDPhi ;
  double fISMetUDPhi ;
  double fISMetDDPhi ;
  RBZero *fRBZero;
  std::vector<double> fVars;
  std::vector<std::string> fLabels;
  int           fN;
};
#endif
