#include "../include/PhotonLoader.hh"
#include <cmath>
#include <iostream> 

using namespace baconhep;

PhotonLoader::PhotonLoader(TTree *iTree,std::string iLabel) { 
  fPhotons  = new TClonesArray("baconhep::TPhoton");
  iTree->SetBranchAddress("Photon",       &fPhotons);
  fPhotonBr  = iTree->GetBranch("Photon");
  fN = 1;
  fYear = iLabel;
  for(int i0 = 0; i0 < fN*3.; i0++) {double pVar = 0; fVars.push_back(pVar);}
}
PhotonLoader::~PhotonLoader() { 
  delete fPhotons;
  delete fPhotonBr;
}
void PhotonLoader::reset() { 
  fNPhotonsLoose  = 0; 
  fNPhotonsTight  = 0;
  fispho0Tight    = 0;
  fTightPhotons.clear();
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
}
void PhotonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("nphoLoose",   &fNPhotonsLoose, "fNPhotonsLoose/I");  // loose photon multiplicity 
  fTree->Branch("nphoTight",   &fNPhotonsTight, "fNPhotonsTight/I");  // tight photon multiplicity
  fTree->Branch("ispho0Tight", &fispho0Tight,   "fispho0Tight/I");
  setupNtuple("vpho",iTree,fN,fVars);                                 // pho0_pt,_eta,_phi (1*3=3)
}
void PhotonLoader::setupTreeQbert(TTree *iTree) {
  reset();
  fTree = iTree;
  setupNtuple("vpho",iTree,fN,fVars);        
}
void PhotonLoader::load(int iEvent) { 
  fPhotons  ->Clear();
  fPhotonBr ->GetEntry(iEvent);
}
void PhotonLoader::selectPhotons(double iRho,std::vector<TLorentzVector> &iVetoes,std::vector<TLorentzVector> &iPhotons) {
  reset();
  int lCount(0), lTCount(0);
  for  (int i0 = 0; i0 < fPhotons->GetEntriesFast(); i0++) {
    TPhoton *pPhoton = (TPhoton*)((*fPhotons)[i0]);

    if(pPhoton->pt        <=  15)                       continue;
    if(fabs(pPhoton->eta) >=  2.5)                      continue;
    if(passVeto(pPhoton->eta,pPhoton->phi,0.4,iVetoes)) continue;
    if(!passPhoLooseSel(pPhoton,iRho,fYear))            continue;
    lCount++;

    if(pPhoton->pt        <= 175)                       continue;
    if(fabs(pPhoton->eta) >= 1.4442)                    continue;
    if(!passPhoTightSel(pPhoton, iRho,fYear))           continue;
    lTCount++;

    if(lCount==1) fispho0Tight = 1;
    addPhoton(pPhoton,fTightPhotons); // sort photons by pt
  }
  addVPhoton(fTightPhotons,iPhotons,0.); // add cleaning (tight) photons

  if(fVars.size() > 0) fillPhoton(fN,fTightPhotons,fVars); // save leading leptons

  fNPhotonsLoose = lCount;
  fNPhotonsTight = lTCount;
}
