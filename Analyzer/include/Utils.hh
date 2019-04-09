#ifndef UTILS_HH
#define UTILS_HH

#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>

#include <vector>
#include <cassert>
#include <iostream>

bool   passJetTightSel        (const baconhep::TJet *jet,std::string iLabel="2017");
bool   passJetTightLepVetoSel (const baconhep::TJet *jet,std::string iLabel="2017");

bool   passEleVetoSel         (const baconhep::TElectron *electron, const double rho,std::string iLabel="2017");
bool   passEleLooseSel        (const baconhep::TElectron *electron, const double rho,std::string iLabel="2017");
bool   passEleTightSel        (const baconhep::TElectron *electron, const double rho,std::string iLabel="2017");
bool   passEleHEEPSel         (const baconhep::TElectron *electron, const double rho, const double met,std::string iLabel="2017");

bool   passPhoLooseSel        (const baconhep::TPhoton *photon, const double rho,std::string iLabel="2017");
bool   passPhoMediumSel       (const baconhep::TPhoton *photon, const double rho,std::string iLabel="2017");
bool   passPhoTightSel        (const baconhep::TPhoton *photon, const double rho,std::string iLabel="2017");

double eleEffArea2016         (const double eta);
double eleEffArea2017         (const double eta);
double phoEffArea2016         (const double eta, const int type);
double phoEffArea2017         (const double eta, const int type);

bool   passMuonLooseSel       (const baconhep::TMuon *muon,std::string iLabel="2017");
bool   passMuonMediumSel      (const baconhep::TMuon *muon,std::string iLabel="2017");
bool   passMuonTightSel       (const baconhep::TMuon *muon,std::string iLabel="2017");
bool   passMuonSoftSel        (const baconhep::TMuon *muon,std::string iLabel="2017");
bool   passMuonHighPtSel      (const baconhep::TMuon *muon,std::string iLabel="2017");

bool   passTauSel             (const baconhep::TTau *tau,std::string iLabel="2017");
bool   passTauTightSel        (const baconhep::TTau *tau,std::string iLabel="2017");

double getVal                 (TH1D*h,double val);
double getVal2D               (TH2D*h,double val1, double val2);
bool   passVeto               (double iEta,double iPhi,double idR,std::vector<TLorentzVector> &iVetoes);
void   setupNtuple            (std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals);
void   setupNtuple            (std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals,int iHead,std::vector<std::string> &iLabels);
void   setupNtuple            (std::string iHeader,TTree *iTree,int iN,std::vector<float> &iVals,std::vector<std::string> &iLabels);
void   setupNtupleVector      (std::string iHeader,TTree *iTree,std::vector<double> &pt, std::vector<double> &eta, std::vector<double> &phi);
void   setupNtupleVector      (std::string iHeader,TTree *iTree,std::vector< std::vector<double> > &iValVectors,std::vector<std::string> &iLabels);
void   setupNtupleVector      (std::string iHeader,TTree *iTree,std::vector< std::vector<float> > &iValVectors,std::vector<std::string> &iLabels);

extern std::string label2016;
extern std::string label2017;
extern std::string label2018;

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
