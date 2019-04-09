#ifndef ElectronLoader_H
#define ElectronLoader_H
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "Utils.hh"

using namespace baconhep;

class ElectronLoader { 
public:
  ElectronLoader(TTree *iTree,std::string iLabel="2017");
  ~ElectronLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load(int iEvent);
  void selectElectrons(double iRho,double iMet, std::vector<TLorentzVector> &iElectrons);
  std::vector<TElectron*> fLooseElectrons, fTightElectrons, fHEEPElectrons;
  std::vector<TLorentzVector> looseElectrons, tightElectrons, HEEPElectrons;
  int           fNElectronsLoose;
  int           fNElectronsTight;
  int           fNElectronsHEEP;
  int           fisele0Tight;
  int           fisele1Tight;

protected: 
  TClonesArray *fElectrons;
  TBranch      *fElectronBr;
  TTree        *fTree;
  std::vector<double>     fVars;
  int           fN;
  std::string   fYear;
};
#endif
