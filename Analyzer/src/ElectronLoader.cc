#include "../include/ElectronLoader.hh"

#include "TFile.h"
#include <cmath>
#include <iostream> 
#include <sstream> 

using namespace baconhep;

ElectronLoader::ElectronLoader(TTree *iTree){
  fElectrons  = new TClonesArray("baconhep::TElectron");
  iTree->SetBranchAddress("Electron",       &fElectrons);
  fElectronBr  = iTree->GetBranch("Electron");
  fN = 2;

  for(int i0 = 0; i0 < fN*3.; i0++) {double pVar = 0; fVars.push_back(pVar);}
  for(int i0 = 0; i0 <     4; i0++) {double pVar = 0; fVars.push_back(pVar);}
}
ElectronLoader::~ElectronLoader() { 
  delete fElectrons;
  delete fElectronBr;
}
void ElectronLoader::reset() { 
  fNElectronsLoose = 0; 
  fNElectronsTight = 0;
  fisele0Tight     = 0;
  fisele1Tight     = 0;
  fLooseElectrons.clear();
  fTightElectrons.clear();
  fHEEPElectrons.clear();
  looseElectrons.clear();
  tightElectrons.clear();
  HEEPElectrons.clear();
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
}
void ElectronLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("neleLoose"   ,&fNElectronsLoose ,"fNElectronsLoose/I");
  fTree->Branch("neleTight"   ,&fNElectronsTight ,"fNElectronsTight/I");
  fTree->Branch("neleHEEP"    ,&fNElectronsHEEP  ,"fNElectronsHEEP/I");
  setupNtuple  ("veleLoose"   ,iTree,fN,fVars);        // 2 electrons ele*_pt,ele*_eta,ele*_phi (2*4=8)
}
void ElectronLoader::load(int iEvent) { 
  fElectrons  ->Clear();
  fElectronBr ->GetEntry(iEvent);
}
void ElectronLoader::selectElectrons(double iRho, double iMet, std::vector<TLorentzVector> &iElectrons) {
  reset();
  int lCount = 0,lTCount=0, lHEEPCount =0; 

  // Electrons selection
  for  (int i0 = 0; i0 < fElectrons->GetEntriesFast(); i0++) { 
    TElectron *pElectron = (TElectron*)((*fElectrons)[i0]);
    if(pElectron->pt        <=  10)                                            continue;
    if(fabs(pElectron->eta) >=  2.5)                                           continue;
    if(fabs(pElectron->eta) > 1.4442 && fabs(pElectron->eta) < 1.566)          continue;
    if(!passEleLooseSel(pElectron, iRho))                                      continue;
    lCount++;
    addElectron(pElectron,fLooseElectrons);

    if(passEleTightSel(pElectron,iRho) && pElectron->pt>40 && fabs(pElectron->eta)<2.5){
      if(lCount==1) fisele0Tight = 1;
      if(lCount==2) fisele1Tight = 1;
      lTCount++;
      addElectron(pElectron,fTightElectrons);
    }

    if(passEleHEEPSel(pElectron, iRho, iMet) && pElectron->pt>40 && fabs(pElectron->eta)<2.5){
      lHEEPCount++;
      addElectron(pElectron,fHEEPElectrons);
    }

  }
  addVElectron(fLooseElectrons,looseElectrons,0.000511);
  addVElectron(fTightElectrons,tightElectrons,0.000511);
  addVElectron(fHEEPElectrons,HEEPElectrons,0.000511);

  fNElectronsLoose = lCount;
  fNElectronsTight = lTCount;
  fNElectronsHEEP = lHEEPCount;

  if(fVars.size() > 0) fillElectron(fN,fLooseElectrons,fVars);

  // Cleaning iElectrons
  if(fTightElectrons.size()>1){
    iElectrons.push_back(tightElectrons[0]); // save first tight electron
    iElectrons.push_back(looseElectrons[1]); // if leading lepton is tight, save the subleading loose one
  }
}
