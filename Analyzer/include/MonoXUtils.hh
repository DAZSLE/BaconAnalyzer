//================================================================================================
//
// Various functions for Mono-X analysis 
//
//________________________________________________________________________________________________
#ifndef MONOXUTILS_HH
#define MONOXUTILS_HH

// bacon object headers
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
//Root
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
//STL
#include <vector>
#include <cassert>
#include <iostream>

//=== FUNCTION DECLARATIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
bool   passJet04Sel           (const baconhep::TJet *jet);
bool   passJetLooseSel        (const baconhep::TJet *jet);
bool   passJetTightLepVetoSel (const baconhep::TJet *jet);
double eleEffArea             (const double eta);
double phoEffArea             (const double eta, const int type);
bool   passEleSel             (const baconhep::TElectron *electron, const double rho);
bool   passEleTightSel        (const baconhep::TElectron *electron, const double rho);
bool   passMuonLooseSel       (const baconhep::TMuon *muon);
bool   passMuonTightSel       (const baconhep::TMuon *muon);
bool   passTauSel             (const baconhep::TTau *tau);
bool   passPhoSel             (const baconhep::TPhoton *photon, const double rho);
bool   passPhoLooseSel        (const baconhep::TPhoton *photon, const double rho);
bool   passPhoMediumSel       (const baconhep::TPhoton *photon, const double rho);
double eleIso                 (const baconhep::TElectron *electron, const double rho);
double phoEffAreaHighPt       (const double eta, const int type);
double getVal                 (TH1D*h,double val);
double getVal2D               (TH2D*h,double val1, double val2);
bool   passVeto               (double iEta,double iPhi,double idR,std::vector<TLorentzVector> &iVetoes);
void   setupNtuple            (std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals);
void   setupNtuple            (std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals,int iHead,std::vector<std::string> &iLabels);
void   setupNtuple            (std::string iHeader,TTree *iTree,int iN,std::vector<float> &iVals,std::vector<std::string> &iLabels);
void   setupNtupleVector      (std::string iHeader,TTree *iTree,std::vector<double> &pt, std::vector<double> &eta, std::vector<double> &phi);
void   setupNtupleVector      (std::string iHeader,TTree *iTree,std::vector< std::vector<double> > &iValVectors,std::vector<std::string> &iLabels);
void   setupNtupleVector      (std::string iHeader,TTree *iTree,std::vector< std::vector<float> > &iValVectors,std::vector<std::string> &iLabels);

template<class T> void addObject(T *iObject,std::vector<T*> &iObjects) {
  bool lFill = false;
  for(typename std::vector<T*>::iterator pIter = iObjects.begin(); pIter != iObjects.end(); pIter++) {
    if((*pIter)->pt > iObject->pt) continue;
    iObjects.insert(pIter,iObject);
    lFill = true;
    break;
  }
  if(!lFill)  iObjects.push_back(iObject);
}
template<class T> void fillObject(int iN,std::vector<T*> &iObjects,std::vector<double> &iVals) { 
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    iVals[i0*3+0] = iObjects[i0]->pt;
    iVals[i0*3+1] = iObjects[i0]->eta;
    iVals[i0*3+2] = iObjects[i0]->phi;
  }
}
template<class T> void addVeto(std::vector<T*> &iObjects,std::vector<TLorentzVector> &iVetoes,double iMass) { 
  for(typename std::vector<T*>::iterator pIter = iObjects.begin(); pIter != iObjects.end(); pIter++) {
    TLorentzVector lVec; lVec.SetPtEtaPhiM((*pIter)->pt,(*pIter)->eta,(*pIter)->phi,iMass);
    iVetoes.push_back(lVec);
  }
}
template<class T> void addVetoV(std::vector<T*> &iObjects,std::vector<TLorentzVector> &iVetoes) {
  for(typename std::vector<T*>::iterator pIter = iObjects.begin(); pIter != iObjects.end(); pIter++) {
    TLorentzVector lVec; lVec.SetPtEtaPhiM((*pIter)->pt,(*pIter)->eta,(*pIter)->phi,(*pIter)->mass);
    iVetoes.push_back(lVec);
  }
}



#define addElectron  addObject<baconhep::TElectron>
#define addMuon      addObject<baconhep::TMuon>
#define addTau       addObject<baconhep::TTau>
#define addJet       addObject<baconhep::TJet>
#define addPhoton    addObject<baconhep::TPhoton>

#define fillElectron  fillObject<baconhep::TElectron>
#define fillMuon      fillObject<baconhep::TMuon>
#define fillTau       fillObject<baconhep::TTau>
#define fillJet       fillObject<baconhep::TJet>
#define fillPhoton    fillObject<baconhep::TPhoton>

#define addVElectron   addVeto<baconhep::TElectron>
#define addVMuon       addVeto<baconhep::TMuon>
#define addVTau        addVeto<baconhep::TTau>
#define addVJet        addVetoV<baconhep::TJet>
#define addVPhoton     addVeto<baconhep::TPhoton>

#endif
