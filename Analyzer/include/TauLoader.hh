#ifndef TauLoader_H
#define TauLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "Utils.hh"

using namespace baconhep;

class TauLoader { 
public:
  TauLoader(TTree *iTree,std::string iLabel="2017");
  ~TauLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load(int iEvent);
  void selectTaus(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons);
  std::vector<TTau*> fSelTaus;
  int           fNTaus;
  int           fNTausTight;

protected: 
  TClonesArray *fTaus;
  TBranch      *fTauBr;
  TTree        *fTree;
  int           fN;
  std::string   fYear;
};
#endif
