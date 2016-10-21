#ifndef gFatJetLoader_H
#define gFatJetLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TgFatJet.hh"
#include "MonoXUtils.hh"

using namespace baconhep;

class gFatJetLoader { 
public:
  gFatJetLoader(TTree *iTree);
  ~gFatJetLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load(int iEvent);
  bool selectgFatJets(std::vector<TLorentzVector> &iVetoes);
  std::vector<TgFatJet*> fSelgFatJets;

protected: 
  TClonesArray *fgFatJets;
  TBranch      *fgFatJetBr;
  TTree        *fTree;
  int           fNgFatJets;
  std::vector<double> fVars;
  int           fN;
};
#endif
