#include "../include/Utils.hh"
//================================================================================================
//
// Various functions using 80X recommendations
//
//________________________________________________________________________________________________

#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include <string>
#include <sstream>
#include <vector>

//=== FUNCTION DEFINITIONS ======================================================================================

// JETS
// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
//--------------------------------------------------------------------------------------------------
/*bool passJetLooseSel(const baconhep::TJet *jet)
{
  if(fabs(jet->eta)<= 2.7){
    if(jet->neuHadFrac >= 0.99) return false;
    if(jet->neuEmFrac  >= 0.99) return false;
    if(jet->nParticles <= 1)    return false;
  }
  if(fabs(jet->eta)<= 2.4) {
    if(jet->chHadFrac == 0)     return false;
    if(jet->nCharged  == 0)     return false;
    if(jet->chEmFrac  >= 0.99)  return false;
  }
  if(fabs(jet->eta) > 2.7 && fabs(jet->eta) <= 3.0) {
    if(jet->neuEmFrac <= 0.01)  return false;
    if(jet->neuHadFrac >= 0.98) return false;
    if(jet->nNeutrals <= 2)     return false;
  }
  if(fabs(jet->eta) > 3.0) {
    if(jet->neuEmFrac >= 0.90)  return false;
    if(jet->nNeutrals <= 10)    return false;
  }
  return true;
}*/
//--------------------------------------------------------------------------------------------------
bool passJetTightSel(const baconhep::TJet *jet)
{
  if(fabs(jet->eta)<= 2.7){
    if(jet->neuHadFrac >= 0.90) return false;
    if(jet->neuEmFrac  >= 0.90) return false;
    if(jet->nParticles <= 1)    return false;
    if(fabs(jet->eta)<= 2.4) {
      if(jet->chHadFrac == 0)     return false;
      if(jet->nCharged  == 0)     return false;
    }
  }
  if(fabs(jet->eta) > 2.7 && fabs(jet->eta) <= 3.0) {
    if(jet->neuEmFrac <= 0.02)  return false;
    if(jet->neuHadFrac >= 0.99) return false;
    if(jet->nNeutrals <= 2)     return false;
  }
  if(fabs(jet->eta) > 3.0) {
    if(jet->neuEmFrac >= 0.90)  return false;
    if(jet->neuHadFrac <= 0.02) return false;
    if(jet->nNeutrals <= 10)    return false;
  }
  return true;
}
//--------------------------------------------------------------------------------------------------
bool passJetTightLepVetoSel(const baconhep::TJet *jet)
{
  if(fabs(jet->eta)<= 2.7){
    if(jet->neuHadFrac >= 0.90) return false;
    if(jet->neuEmFrac  >= 0.90) return false;
    if(jet->nParticles <= 1)    return false;
    if(jet->muonFrac   >= 0.8)  return false;
    if(fabs(jet->eta)<= 2.4) {
      if(jet->chHadFrac == 0)     return false;
      if(jet->nCharged  == 0)     return false;
      if(jet->chEmFrac  >= 0.80)  return false;
    }
  }
  return true;
}
// ELECTRONS
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
//--------------------------------------------------------------------------------------------------
bool passEleVetoSel(const baconhep::TElectron *electron, const double rho)
{
  if(electron->isConv) return false;

  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );

  if(fabs(electron->scEta)<=1.479) {
    if(iso >= 0.168*(electron->pt)) return false;

    if(electron->sieie              >= 0.01280)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00523)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.15900)                        return false;
    if(electron->hovere             >= 0.05 + 1.12/(electron->ecalEnergy) + 0.036*rho/(electron->ecalEnergy))                                          return false;
    if(fabs(1.0 - electron->eoverp) >= 0.19300*(electron->ecalEnergy)) return false;
    if(electron->nMissingHits       >  2)                              return false;
  } else {
    if(iso >= 0.185*(electron->pt)) return false;

    if(electron->sieie              >= 0.04450)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00984)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.15700)                        return false;
    if(electron->hovere             >= 0.05 + 0.5/(electron->ecalEnergy) + 0.201*rho/(electron->ecalEnergy))                                           return false;
    if(fabs(1.0 - electron->eoverp) >= 0.09620*(electron->ecalEnergy)) return false;
    if(electron->nMissingHits       >  3)                              return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
bool passEleLooseSel(const baconhep::TElectron *electron, const double rho)
{
  if(electron->isConv) return false;

  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );

  if(fabs(electron->scEta)<=1.479) {
    if(iso >= 0.133*(electron->pt)) return false;

    if(electron->sieie              >= 0.0105)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00387)                       return false;
    if(fabs(electron->dPhiIn)       >= 0.0761)                        return false;
    if(electron->hovere             >= 0.05 + 1.12/(electron->ecalEnergy) + 0.036*rho/(electron->ecalEnergy))                                         return false;
    if(fabs(1.0 - electron->eoverp) >= 0.129*(electron->ecalEnergy))  return false;
    if(electron->nMissingHits       >  1)                             return false;
  } else {
    if(iso >= 0.146*(electron->pt)) return false;

    if(electron->sieie              >= 0.03560)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00720)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.14700)                        return false;
    if(electron->hovere             >= 0.0414 + 0.5/(electron->ecalEnergy) + 0.201*rho/(electron->ecalEnergy))                                         return false;
    if(fabs(1.0 - electron->eoverp) >= 0.08750*(electron->ecalEnergy)) return false;
    if(electron->nMissingHits       >  1)                              return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
bool passEleMediumSel(const baconhep::TElectron *electron, const double rho)
{
  if(electron->isConv) return false;

  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );

  if(fabs(electron->scEta)<=1.479) {
    if(iso >= 0.0718*(electron->pt)) return false;

    if(electron->sieie              >= 0.01050)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00365)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.05880)                        return false;
    if(electron->hovere             >= 0.026 + 1.12/(electron->ecalEnergy) + 0.036*rho/(electron->ecalEnergy))                                         return false;
    if(fabs(1.0 - electron->eoverp) >= 0.03270*(electron->ecalEnergy)) return false;
    if(electron->nMissingHits       >  1)                              return false;
  } else {
    if(iso >= 0.143*(electron->pt)) return false;

    if(electron->sieie              >= 0.03090)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00625)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.03550)                        return false;
    if(electron->hovere             >= 0.026 + 0.5/(electron->ecalEnergy) + 0.201*rho/(electron->ecalEnergy))                                          return false;
    if(fabs(1.0 - electron->eoverp) >= 0.03350*(electron->ecalEnergy)) return false;
    if(electron->nMissingHits       >  1)                              return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
bool passEleTightSel(const baconhep::TElectron *electron, const double rho)
{
  if(electron->isConv) return false;

  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );

  if(fabs(electron->scEta)<=1.479) {
    if(iso >= 0.0361*(electron->pt)) return false;

    if(electron->sieie              >= 0.01040)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00353)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.04990)                        return false;
    if(electron->hovere             >= 0.026 + 1.12/(electron->ecalEnergy) + 0.036*rho/(electron->ecalEnergy))                                         return false;
    if(fabs(1.0 - electron->eoverp) >= 0.02780*(electron->ecalEnergy)) return false;
    if(electron->nMissingHits       >  1)                              return false;
  } else {
    if(iso >= 0.094*(electron->pt)) return false;

    if(electron->sieie              >= 0.03050)                        return false;
    if(fabs(electron->dEtaInSeed)   >= 0.00567)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.01650)                        return false;
    if(electron->hovere             >= 0.026 + 0.5/(electron->ecalEnergy) + 0.201*rho/(electron->ecalEnergy))                                          return false;
    if(fabs(1.0 - electron->eoverp) >= 0.01580*(electron->ecalEnergy)) return false;
    if(electron->nMissingHits       >  1)                              return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------                                                            
// HEEP Electrons
// https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2#Selection_Cuts_HEEP_V7_0_recomme
bool passEleHEEPSel(const baconhep::TElectron *electron, const double rho, const double met)
{
  if(fabs(electron->scEta)<1.442) {
    if(met                             <= 35)                                              return false;
    if(!(electron->typeBits & baconhep::kEcalDriven))                                      return false;
    if(fabs(electron->dEtaInSeed)      >= 0.00400)                                         return false;
    if(fabs(electron->dPhiIn)          >= 0.06000)                                         return false;
    if(electron->hovere                >= (1/(electron->ecalEnergy))+0.05)                 return false;
    if((electron->e2x5/electron->e5x5) <= 0.94 || (electron->e1x5/electron->e1x5) <= 0.83) return false;
    if(electron->hcalDepth1Iso         >= 2.0+0.03*met+0.28*rho)                           return false;
    if(electron->trkIso                >= 5)                                               return false;
    if(electron->nMissingHits          >  1)                                               return false;
    if(electron->d0                    >= 0.02)                                            return false;
    return true;
  } 
  else if(fabs(electron->scEta) > 1.566 && fabs(electron->scEta) < 2.5){
    if(met                             <= 35)                                              return false;
    if(!(electron->typeBits & baconhep::kEcalDriven))                                      return false;
    if(fabs(electron->dEtaInSeed)      >= 0.00600)                                         return false;
    if(fabs(electron->dPhiIn)          >= 0.06000)                                         return false;
    if(electron->hovere                >= (5/(electron->ecalEnergy))+0.05)                 return false;
    if(electron->sieie                 >= 0.03)                                            return false;
    if(met < 50){
      if(electron->hcalDepth1Iso       >= 2.5+0.28*rho)                                    return false;
    }
    else{
      if(electron->hcalDepth1Iso       >= 2.5+0.03*(met-50)+0.28*rho)                      return false;
    }
    if(electron->trkIso                >= 5)                                               return false;
    if(electron->nMissingHits          >  1)                                               return false;
    if(electron->d0                    >= 0.05)                                            return false;
    return true;
  }
  else{
    return false;
  }
  return true;
}
// PHOTONS
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Recommended_Working_points_for_2
//--------------------------------------------------------------------------------------------------
bool passPhoLooseSel(const baconhep::TPhoton *photon, const double rho)
{
  double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea(photon->scEta, 0), (double)0.);
  double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea(photon->scEta, 1), (double)0.);
  double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea(photon->scEta, 2), (double)0.);

  double maxneuHadIsoB = 9.188+0.0126*photon->pt+0.000026*photon->pt*photon->pt;
  double maxneuHadIsoE = 10.471+0.0119*photon->pt+0.000025*photon->pt*photon->pt;

  if(fabs(photon->scEta) <= 1.479) {
    if(photon->hovere   > 0.1050)                                      return false;
    if(photon->sieie    > 0.01030)                                     return false;
    if(chHadIso         > 2.839)                                       return false;
    if(neuHadIso        > maxneuHadIsoB)                               return false;                                  
    if(phoIso           > 2.956 + 0.0035*photon->pt)                   return false;

  } else {
    if(photon->hovere   > 0.0290)                                      return false;
    if(photon->sieie    > 0.02760)                                     return false;
    if(chHadIso         > 2.150)                                       return false;
    if(neuHadIso        > maxneuHadIsoE)                               return false;
    if(phoIso           > 4.895 + 0.0040*photon->pt)                   return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
bool passPhoMediumSel(const baconhep::TPhoton *photon, const double rho)
{
  double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea(photon->scEta, 0), (double)0.);
  double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea(photon->scEta, 1), (double)0.);
  double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea(photon->scEta, 2), (double)0.);

  double maxneuHadIsoB = 2.491+0.0126*photon->pt+0.000026*photon->pt*photon->pt;
  double maxneuHadIsoE = 9.131+0.0119*photon->pt+0.000025*photon->pt*photon->pt;

  if(fabs(photon->scEta) <= 1.479) {
    if(photon->hovere   > 0.0350)                                      return false;
    if(photon->sieie    > 0.01030)                                     return false;
    if(chHadIso         > 1.416)                                       return false;
    if(neuHadIso        > maxneuHadIsoB)                               return false;
    if(phoIso           > 2.952 + 0.0040*photon->pt)                   return false;

  } else {
    if(photon->hovere   > 0.0270)                                      return false;
    if(photon->sieie    > 0.02710)                                     return false;
    if(chHadIso         > 1.012)                                       return false;
    if(neuHadIso        > maxneuHadIsoE)                               return false;
    if(phoIso           > 4.095 + 0.0040*photon->pt)                   return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
bool passPhoTightSel(const baconhep::TPhoton *photon, const double rho)
{
  double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea(photon->scEta, 0), (double)0.);
  double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea(photon->scEta, 1), (double)0.);
  double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea(photon->scEta, 2), (double)0.);

  double maxneuHadIsoB = 1.267+0.0126*photon->pt+0.000026*photon->pt*photon->pt;
  double maxneuHadIsoE = 8.916+0.0119*photon->pt+0.000025*photon->pt*photon->pt;

  if(fabs(photon->scEta) <= 1.479) {
    if(photon->hovere   > 0.0200)                                      return false;
    if(photon->sieie    > 0.01030)                                     return false;
    if(chHadIso         > 1.158)                                       return false;
    if(neuHadIso        > maxneuHadIsoB)                               return false;
    if(phoIso           > 2.065 + 0.0035*photon->pt)                   return false;

  } else {
    if(photon->hovere   > 0.0250)                                      return false;
    if(photon->sieie    > 0.02710)                                     return false;
    if(chHadIso         > 0.575)                                       return false;
    if(neuHadIso        > maxneuHadIsoE)                               return false;
    if(phoIso           > 3.272 + 0.0040*photon->pt)                   return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
// Effective Areas
// https://github.com/lsoffi/cmssw/blob/CMSSW_9_2_X_TnP/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt
extern double eleEffArea(const double eta)
{
  if     (fabs(eta) >= 0.0000 && fabs(eta) < 1.0000) { return 0.1566; }
  else if(fabs(eta) >= 1.0000 && fabs(eta) < 1.4790) { return 0.1626; }
  else if(fabs(eta) >= 1.4790 && fabs(eta) < 2.0000) { return 0.1073; }
  else if(fabs(eta) >= 2.0000 && fabs(eta) < 2.2000) { return 0.0854; }
  else if(fabs(eta) >= 2.2000 && fabs(eta) < 2.3000) { return 0.1051; }
  else if(fabs(eta) >= 2.3000 && fabs(eta) < 2.4000) { return 0.1204; }
  else                                               { return 0.1524; }
}
//--------------------------------------------------------------------------------------------------
//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Offline_selection_criteria
extern double phoEffArea(const double eta, const int type)
{
  const int kCH_HAD  = 0;
  const int kNEU_HAD = 1;
  const int kPHOTON  = 2;

  if(type==kCH_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0385; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0468; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0435; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0378; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0338; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0314; }
    else                                             { return 0.0269; }

  } else if(type==kNEU_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0636; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.1103; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0759; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0236; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0151; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.00007; }
    else                                             { return 0.0132; }

  } else if(type==kPHOTON) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.1240; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.1093; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0631; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0779; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0999; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.1155; }
    else                                             { return 0.1373; }

  } else { assert(0); }
}
// MUONS
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2 
//-------------------------------------------------------------------------------------------------
bool passMuonLooseSel(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGLooseMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;

  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonMediumSel(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGMediumMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;

  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonTightSel(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGTightMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.15*(muon->pt)) return false;
  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonSoftSel(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGSoftMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;
  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonHighPtSel(const baconhep::TMuon *muon)
{
  if(!(muon->pogIDBits & baconhep::kPOGHighPtMuon)) return false;

  // PF-isolation with Delta-beta correction                                                                                  
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;
  return true;
}
// TAUS https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#2017v2_discriminators
//--------------------------------------------------------------------------------------------------
bool passTauSel(const baconhep::TTau *tau)
{
  if(!(tau->hpsDisc & baconhep::kByDecayModeFinding)) return false;
  if(tau->rawIso3Hits > 4.5)                           return false;

  return true;
}
// old decay mode finding? + tight MVA isolation + loose anti-electron MVA + tight anti-muon cut
bool passTauTightSel(const baconhep::TTau *tau)
{
  if(!(tau->hpsDisc & baconhep::kByDecayModeFinding)) return false;
  if(!(tau->hpsDisc & baconhep::kByVTightIsolationMVA3oldDMwLT)) return false;
  if(!(tau->hpsDisc & baconhep::kByMVA6LooseElectronRejection)) return false;
  if(!(tau->hpsDisc & baconhep::kByTightMuonRejection3)) return false;
  return true;
}
// Tools
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
    if(i0 % int(iLabels.size()) == int(iLabels.size())-1 && i0 > 0) lCount++; 
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
void setupNtupleVector(std::string iHeader,TTree *iTree,std::vector<double> &pt, std::vector<double> &eta, std::vector<double> &phi) {
  std::stringstream pSPt,pSEta,pSPhi;  
  pSPt  << iHeader << "_pt";
  pSEta << iHeader << "_eta";
  pSPhi << iHeader  << "_phi";
  iTree->Branch(pSPt .str().c_str(), &pt);
  iTree->Branch(pSEta.str().c_str(), &eta);
  iTree->Branch(pSPhi.str().c_str(), &phi);
}
//--------------------------------------------------------------------------------------------------
void setupNtupleVector(std::string iHeader,TTree *iTree,std::vector< std::vector<double> > &iValVectors,std::vector<std::string> &iLabels) { 
  for(int i0 = 0; i0 < int(iLabels.size()); i0++) { 
    std::stringstream pVal;
    pVal  << iHeader << "_" << iLabels[i0];
    iTree->Branch(pVal .str().c_str(),&iValVectors[i0]);
  }
}
//--------------------------------------------------------------------------------------------------
void setupNtupleVector(std::string iHeader,TTree *iTree,std::vector< std::vector<float> > &iValVectors,std::vector<std::string> &iLabels) { 
  for(int i0 = 0; i0 < int(iLabels.size()); i0++) { 
    std::stringstream pVal;
    pVal  << iHeader << "_" << iLabels[i0];
    iTree->Branch(pVal .str().c_str(),&iValVectors[i0]);
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
