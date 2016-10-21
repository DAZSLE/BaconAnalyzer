#ifndef ElectronLoader_H
#define ElectronLoader_H
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "MonoXUtils.hh"

using namespace baconhep;

class ElectronLoader { 
public:
  ElectronLoader(TTree *iTree);
  ~ElectronLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load(int iEvent);
  void selectElectrons(double iRho,std::vector<TLorentzVector> &iElectrons);
  std::vector<TElectron*> fLooseElectrons, fTightElectrons;
  std::vector<TLorentzVector> looseElectrons, tightElectrons;
  int           fNElectronsLoose;
  int           fNElectronsTight;
  int           fisele0Tight;
  int           fisele1Tight;

protected: 
  TClonesArray *fElectrons;
  TBranch      *fElectronBr;
  TTree        *fTree;
  std::vector<double>     fVars;
  int           fN;

};
#endif
