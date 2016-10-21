#include "../include/MonoXUtils.hh"
//================================================================================================
//
// Various functions for Mono-X analysis 
//
//________________________________________________________________________________________________

// bacon object headers
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include "../include/LeptonEffUtils.hh"

#include <string>
#include <sstream>
#include <vector>

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
bool passJet04Sel(const baconhep::TJet *jet)
{
  // Loose PFJet ID
  if(jet->neuHadFrac >= 0.99) return false;
  if(jet->neuEmFrac  >= 0.99) return false;
  if(jet->nParticles <= 1)    return false;
  if(jet->muonFrac   >= 0.8)  return false;
  if(fabs(jet->eta)<2.4) {
    if(jet->chHadFrac == 0)    return false;
    if(jet->nCharged  == 0)    return false;
    if(jet->chEmFrac  >= 0.99) return false;
  }

  // PU Jet ID
  if     (0    <= fabs(jet->eta) && fabs(jet->eta) < 2.5  && jet->mva < -0.63) return false;
  else if(2.5  <= fabs(jet->eta) && fabs(jet->eta) < 2.75 && jet->mva < -0.60) return false;
  else if(2.75 <= fabs(jet->eta) && fabs(jet->eta) < 3    && jet->mva < -0.55) return false;
  else if(3    <= fabs(jet->eta) && fabs(jet->eta) < 5    && jet->mva < -0.45) return false;

  return true;
}
//--------------------------------------------------------------------------------------------------
bool passJetLooseSel(const baconhep::TJet *jet)
{
  // Loose PFJet ID
  if(jet->neuHadFrac >= 0.99) return false;
  if(jet->neuEmFrac  >= 0.99) return false;
  if(jet->nParticles <= 1)    return false;
  if(fabs(jet->eta)<2.4) {
    if(jet->chHadFrac == 0)    return false;
    if(jet->nCharged  == 0)    return false;
    if(jet->chEmFrac  >= 0.99) return false;
  }
  return true;
}
//--------------------------------------------------------------------------------------------------
bool passJetTightLepVetoSel(const baconhep::TJet *jet)
{
  // Tight PFJet ID
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
  if(jet->neuHadFrac >= 0.90) return false;
  if(jet->neuEmFrac  >= 0.90) return false;
  if(jet->nParticles <= 1)    return false;
  if(jet->muonFrac   >= 0.8)  return false;
  if(fabs(jet->eta)<2.4) {
    if(jet->chHadFrac == 0)    return false;
    if(jet->nCharged  == 0)    return false;
    if(jet->chEmFrac  >= 0.90) return false;
  }
  return true;
}
//--------------------------------------------------------------------------------------------------
extern double eleEffArea(const double eta)
{
  // effective area for PU correction
  // (see slide 4 of https://indico.cern.ch/event/370494/contribution/2/material/slides/0.pdf)

  if     (fabs(eta) >= 0.0 && fabs(eta) < 0.8) { return 0.1752; }
  else if(fabs(eta) >= 0.8 && fabs(eta) < 1.3) { return 0.1862; }
  else if(fabs(eta) >= 1.3 && fabs(eta) < 2.0) { return 0.1411; }
  else if(fabs(eta) >= 2.0 && fabs(eta) < 2.2) { return 0.1534; }
  else if(fabs(eta) >= 2.2 && fabs(eta) < 2.3) { return 0.1903; }
  else if(fabs(eta) >= 2.3 && fabs(eta) < 2.4) { return 0.2243; }
  else                                         { return 0.2687; }
}

//--------------------------------------------------------------------------------------------------
extern double phoEffArea(const double eta, const int type)
{
  // effective area for PU correction
  // (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Pointers_for_PHYS14_selection_im)

  const int kCH_HAD  = 0;
  const int kNEU_HAD = 1;
  const int kPHOTON  = 2;

  if(type==kCH_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0157; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0143; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0115; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0094; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0095; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0068; }
    else                                             { return 0.0053; }

  } else if(type==kNEU_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0143; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0210; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0147; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0082; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0124; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0186; }
    else                                             { return 0.0320; }

  } else if(type==kPHOTON) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0725; }//{ return 0.0780; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0604; }//{ return 0.0629; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0320; }//{ return 0.0264; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0512; }//{ return 0.0462; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0766; }//{ return 0.0740; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0949; }//{ return 0.0924; }
    else                                             { return 0.1160; }//{ return 0.1484; }

  } else { assert(0); }
}
extern double phoEffAreaHighPt(const double eta, const int type)
{
  // effective area for PU correction
  // (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Pointers_for_PHYS14_selection_im)

  const int kPHOTON  = 2;

  if(type==kPHOTON) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.17; }//{ return 0.0780; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.14; }//{ return 0.0629; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.11; }//{ return 0.0264; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.14; }//{ return 0.0462; }
    else                                             { return 0.22; }//{ return 0.1484; }
  } else { assert(0); }
}
bool passEleSel(const baconhep::TElectron *electron, const double rho)
{
  // Phys14 PU20 bx25 cut-based veto ID
  // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_PHYS14_sample

  if(electron->isConv) return false;

  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );

  if(fabs(electron->scEta)<1.479) {
    if(iso >= 0.126*(electron->pt)) return false;

    if(electron->sieie              >= 0.01140)                        return false;
    if(fabs(electron->dEtaIn)       >= 0.01520)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.21600)                        return false;
    if(electron->hovere             >= 0.18100)                        return false;
    if(fabs(1.0 - electron->eoverp) >= 0.20700*(electron->ecalEnergy)) return false;
    if(fabs(electron->d0)           >= 0.05640)                        return false;
    if(fabs(electron->dz)           >= 0.47200)                        return false;
    if(electron->nMissingHits       >  2)                               return false;

  } else {
    if(iso >= 0.144*(electron->pt)) return false;

    if(electron->sieie              >= 0.03520)                        return false;
    if(fabs(electron->dEtaIn)       >= 0.01130)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.23700)                        return false;
    if(electron->hovere             >= 0.11600)                        return false;
    if(fabs(1.0 - electron->eoverp) >= 0.17400*(electron->ecalEnergy)) return false;
    if(fabs(electron->d0)           >= 0.22200)                        return false;
    if(fabs(electron->dz)           >= 0.92100)                        return false;
    if(electron->nMissingHits       >  3)                              return false;
  }

  return true;
}
double eleIso(const baconhep::TElectron *electron, const double rho) {
  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );
  return iso;
}
bool passEleTightSel(const baconhep::TElectron *electron, const double rho)
{
  // Phys14 PU20 bx25 cut-based veto ID
  // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_PHYS14_sample
  // if(electron->pt < 40 || fabs(electron->eta) > 2.5) return false;

  if(electron->isConv) return false;

  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );

  if(fabs(electron->scEta)<1.479) {
    if(iso >= 0.0354*(electron->pt)) return false;

    if(electron->sieie              >= 0.01010)                        return false;
    if(fabs(electron->dEtaIn)       >= 0.00926)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.03360)                        return false;
    if(electron->hovere             >= 0.05970)                        return false;
    if(fabs(1.0 - electron->eoverp) >= 0.01200*(electron->ecalEnergy)) return false;
    if(fabs(electron->d0)           >= 0.01110)                        return false;
    if(fabs(electron->dz)           >= 0.04660)                        return false;
    if(electron->nMissingHits       >  2)                              return false;
  } else {
    if(iso >= 0.0646*(electron->pt)) return false;

    if(electron->sieie              >= 0.02790)                        return false;
    if(fabs(electron->dEtaIn)       >= 0.00724)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.09180)                        return false;
    if(electron->hovere             >= 0.06150)                        return false;
    if(fabs(1.0 - electron->eoverp) >= 0.00999*(electron->ecalEnergy)) return false;
    if(fabs(electron->d0)           >= 0.03510)                        return false;
    if(fabs(electron->dz)           >= 0.41700)                        return false;
    if(electron->nMissingHits       >  1)                              return false;
  }

  return true;
}

//-------------------------------------------------------------------------------------------------
bool passMuonLooseSel(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGLooseMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.2*(muon->pt)) return false;

  return true;
}

//-------------------------------------------------------------------------------------------------
bool passMuonTightSel(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGTightMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.12*(muon->pt)) return false;
  //if(muon->pt < 20 || fabs(muon->eta) > 2.4) return false;
  return true;
}

//--------------------------------------------------------------------------------------------------
bool passTauSel(const baconhep::TTau *tau)
{
  if(!(tau->hpsDisc & baconhep::kByDecayModeFinding)) return false;
  if(tau->rawIso3Hits > 5.)                           return false;

  return true;
}

//--------------------------------------------------------------------------------------------------
bool passPhoLooseSel(const baconhep::TPhoton *photon, const double rho)
{
  // Loose photon ID (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#PHYS14_selections_PU20_bunch_cro)
  //if(!(photon->passElectronVeto)) return false;  // conversion safe electron veto

  double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea(photon->scEta, 0), (double)0.);
  double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea(photon->scEta, 1), (double)0.);
  double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea(photon->scEta, 2), (double)0.);

  if(fabs(photon->scEta) <= 1.479) {
    if(photon->sthovere > 0.05)                                        return false;
    if(photon->sieie    > 0.0103)                                      return false;
    if(chHadIso         > 2.44)                                        return false;
    if(neuHadIso        > 2.57 + TMath::Exp(0.0044*photon->pt+0.5809)) return false;
    if(phoIso           > 1.92 + 0.0043*photon->pt)                    return false;

  } else {
    if(photon->sthovere > 0.05)                                        return false;
    if(photon->sieie    > 0.0277)                                      return false;
    if(chHadIso         > 1.84)                                        return false;
    if(neuHadIso        > 4.00 + TMath::Exp(0.0040*photon->pt+0.9402)) return false;
    if(phoIso           > 2.15 + 0.0041*photon->pt)                    return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
bool passPhoMediumSel(const baconhep::TPhoton *photon, const double rho)
{
  // Medium photon ID (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#PHYS14_selections_PU20_bunch_cro)
  //if(!(photon->passElectronVeto)) return false;  // conversion safe electron veto

  double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea(photon->scEta, 0), (double)0.);
  double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea(photon->scEta, 1), (double)0.);
  double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea(photon->scEta, 2), (double)0.);

  if(fabs(photon->scEta) <= 1.479) {
    if(photon->sthovere > 0.05)                                        return false;
    if(photon->sieie    > 0.0100)                                      return false;
    if(chHadIso         > 1.31)                                        return false;
    if(neuHadIso        > 0.60 + TMath::Exp(0.0044*photon->pt+0.5809)) return false;
    if(phoIso           > 1.33 + 0.0043*photon->pt)                    return false;
  } else {
    if(photon->sthovere > 0.05)                                        return false;
    if(photon->sieie    > 0.0267)                                      return false;
    if(chHadIso         > 1.25)                                        return false;
    if(neuHadIso        > 1.65 + TMath::Exp(0.0040*photon->pt+0.9402)) return false;
    if(phoIso           > 1.02 + 0.0041*photon->pt)                    return false;
  }

  return true;
}
///tools
//--------------------------------------------------------------------------------------------------
bool passVeto(double iEta,double iPhi,double idR, std::vector<TLorentzVector> &iVetoes) { 
  bool pMatch = false;
  for(unsigned int i1 = 0; i1 < iVetoes.size(); i1++) { 
    double pDEta = iEta - iVetoes[i1].Eta();
    double pDPhi = iPhi - iVetoes[i1].Phi();
    if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
    if(sqrt(pDPhi*pDPhi+pDEta*pDEta) > idR) continue;
    if(iVetoes[i1].Pt() < 0) continue;
    pMatch = true;
  }
  return pMatch;
}
//--------------------------------------------------------------------------------------------------
void setupNtuple(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals) { 
  for(int i0 = 0; i0 < iN; i0++) { 
    int iBase = i0*3;
    std::stringstream pSPt,pSEta,pSPhi;  
    pSPt  << iHeader << i0 << "_pt";
    pSEta << iHeader << i0 << "_eta";
    pSPhi << iHeader << i0 << "_phi";
    iTree->Branch(pSPt .str().c_str(),&iVals[iBase+0],(pSPt .str()+"/D").c_str()); 
    iTree->Branch(pSEta.str().c_str(),&iVals[iBase+1],(pSEta.str()+"/D").c_str());
    iTree->Branch(pSPhi.str().c_str(),&iVals[iBase+2],(pSPhi.str()+"/D").c_str());
  }
}
//--------------------------------------------------------------------------------------------------
void setupNtuple(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals,int iHead,std::vector<std::string> &iLabels) { 
  int lBase  = iHead;
  int lCount = 0;
  for(int i0 = 0; i0 < iN*(int(iLabels.size())); i0++) { 
    std::stringstream pVal;
    pVal  << iHeader << lCount  << "_" << iLabels[i0 % iLabels.size()];
    iTree->Branch(pVal .str().c_str(),&iVals[lBase],(pVal .str()+"/D").c_str());
    if(i0 % int(iLabels.size()) == 0 && i0 > 0) lCount++; 
    lBase++;
  }
}
//--------------------------------------------------------------------------------------------------
void setupNtuple(std::string iHeader,TTree *iTree,int iN,std::vector<float> &iVals,std::vector<std::string> &iLabels) {
  int lBase  = 0;
  int lCount = 0;
  for(int i0 = 0; i0 < iN*(int(iLabels.size())); i0++) {
    std::stringstream pVal;
    pVal  << iHeader << lCount  << "_" << iLabels[i0 % iLabels.size()];
    iTree->Branch(pVal .str().c_str(),&iVals[lBase],(pVal .str()+"/F").c_str());
    if(i0 % int(iLabels.size()) == 0 && i0 > 0) lCount++;
    lBase++;
  }
}
//--------------------------------------------------------------------------------------------------
double getVal(TH1D*h,double val) {
  return h->GetBinContent(h->FindBin(val));
}
//--------------------------------------------------------------------------------------------------
double getVal2D(TH2D*h,double val1, double val2) {
  val2 = (val2 > h->GetYaxis()->GetXmax()) ? h->GetYaxis()->GetBinCenter(h->GetNbinsY()) : val2;
  val2 = (val2 < h->GetYaxis()->GetXmin()) ? h->GetYaxis()->GetBinCenter(1) : val2;
  val1 = (val1 > h->GetXaxis()->GetXmax()) ? h->GetXaxis()->GetBinCenter(h->GetNbinsX()) : val1;
  val1 = (val1 < h->GetXaxis()->GetXmin()) ? h->GetXaxis()->GetBinCenter(1) : val1;
  return h->GetBinContent(h->FindBin(val1,val2));
}
