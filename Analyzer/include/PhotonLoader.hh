#ifndef PhotonLoader_H
#define PhotonLoader_H
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "Utils.hh"

using namespace baconhep;

class PhotonLoader { 
public:
  PhotonLoader(TTree *iTree,std::string iLabel="2017");
  ~PhotonLoader();
  void reset();
  void setupTree(TTree *iTree);
  void setupTreeQbert(TTree *iTree);
  void load(int iEvent);
  void selectPhotons(double iRho,std::vector<TLorentzVector> &iVetoes,std::vector<TLorentzVector> &iPhotons);

  std::vector<TPhoton*> fTightPhotons;
  int           fNPhotonsLoose;
  int           fNPhotonsTight;
  int           fispho0Tight;

protected: 
  TClonesArray *fPhotons;
  TBranch      *fPhotonBr;
  TTree        *fTree;
  
  std::vector<double> fVars;
  int           fN;
  double        fIso;
  std::string   fYear;
};
#endif
