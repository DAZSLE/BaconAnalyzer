#include "../include/TauLoader.hh"
#include <cmath>
#include <iostream> 

using namespace baconhep;

TauLoader::TauLoader(TTree *iTree) { 
  fTaus  = new TClonesArray("baconhep::TTau");
  iTree->SetBranchAddress("Tau",       &fTaus);
  fTauBr  = iTree->GetBranch("Tau");
  fN = 1;
}
TauLoader::~TauLoader() { 
  delete fTaus;
  delete fTauBr;
}
void TauLoader::reset() { 
  fNTaus = 0; 
  fSelTaus.clear();
}
void TauLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("ntau",&fNTaus,"fNTaus/I");                 // tau multiplicity
}
void TauLoader::load(int iEvent) { 
  fTaus   ->Clear();
  fTauBr ->GetEntry(iEvent);
}
void TauLoader::selectTaus(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons) {
  reset(); 
  int lCount = 0; 
  for  (int i0 = 0; i0 < fTaus->GetEntriesFast(); i0++) { 
    TTau *pTau = (TTau*)((*fTaus)[i0]);
    if(passVeto(pTau->eta,pTau->phi,0.4,iElectrons)) continue; 
    if(passVeto(pTau->eta,pTau->phi,0.4,iMuons))     continue;
    if(pTau->pt        <=  18)                       continue;
    if(fabs(pTau->eta) >=  2.3)                      continue;
    if(!passTauSel(pTau))                            continue;
    addTau(pTau,fSelTaus);
    lCount++;
  }
  fNTaus = lCount;
}
