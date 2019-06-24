#include "../include/Utils.hh"
//================================================================================================
//
// Various functions using 8X/9X/10X recommendations
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

std::string label2016="2016";
std::string label2017="2017";
std::string label2018="2018";

// JETS
bool passJetTightSel(const baconhep::TJet *jet,std::string iLabel)
{
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
  if(iLabel==label2016){
    if(fabs(jet->eta)<= 2.7){
      if(jet->neuHadFrac >= 0.90) return false;
      if(jet->neuEmFrac  >= 0.90) return false;
      if(jet->nParticles <= 1)    return false;
      if(fabs(jet->eta)<= 2.4) {
	if(jet->chHadFrac == 0)     return false;
	if(jet->nCharged  == 0)     return false;
	if(jet->chEmFrac  >= 0.90)  return false;
      }
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
  }
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
  if(iLabel==label2017){
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
  // this is very prelim and no puppi? https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
  if(iLabel==label2018){
    if(fabs(jet->eta)<= 2.6){
      if(jet->neuHadFrac >= 0.90) return false;
      if(jet->neuEmFrac  >= 0.90) return false;
      if(jet->nParticles <= 1)    return false;
      if(jet->chHadFrac == 0)     return false;
      if(jet->nCharged  == 0)     return false;
    }
    if(fabs(jet->eta) > 2.6 && fabs(jet->eta) <= 2.7) {
      if(jet->neuHadFrac >= 0.90) return false;
      if(jet->neuEmFrac  >= 0.99) return false;
      if(jet->nCharged  == 0)     return false;
    }
    if(fabs(jet->eta) > 2.7 && fabs(jet->eta) <= 3.0) {
      if(jet->neuEmFrac <= 0.02)  return false;
      if(jet->neuHadFrac >= 0.99) return false;
      if(jet->nNeutrals <= 2)     return false;
    }
    if(fabs(jet->eta) > 3.0) {
      if(jet->neuEmFrac >= 0.90)  return false;
      if(jet->neuHadFrac <= 0.2) return false;
      if(jet->nNeutrals <= 10)    return false;
    }
    return true;
  }
  return false;
}
//--------------------------------------------------------------------------------------------------
bool passJetTightLepVetoSel(const baconhep::TJet *jet,std::string iLabel)
{
  if(iLabel==label2016){
    if(fabs(jet->eta)<= 2.7){
      if(jet->neuHadFrac >= 0.90) return false;
      if(jet->neuEmFrac  >= 0.90) return false;
      if(jet->nParticles <= 1)    return false;
      if(jet->muonFrac   >= 0.8)  return false;
    }
    if(fabs(jet->eta)<= 2.4) {
      if(jet->chHadFrac == 0)     return false;
      if(jet->nCharged  == 0)     return false;
      if(jet->chEmFrac  >= 0.90)  return false;
    }
    return true;
  }
  if(iLabel==label2017){
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
  if(iLabel==label2018){
    if(fabs(jet->eta)<= 2.6){
      if(jet->neuHadFrac >= 0.90) return false;
      if(jet->neuEmFrac  >= 0.90) return false;
      if(jet->nParticles <= 1)    return false;
      if(jet->muonFrac   >= 0.8)  return false;
      if(jet->chHadFrac == 0)     return false;
      if(jet->nCharged  == 0)     return false;
      if(jet->chEmFrac  >= 0.80)  return false;
    }
    if(fabs(jet->eta) > 2.6 && fabs(jet->eta) <= 2.7) {
      if(jet->neuHadFrac >= 0.90) return false;
      if(jet->neuEmFrac  >= 0.99) return false;
      if(jet->muonFrac   >= 0.8)  return false;
    if(jet->nCharged  == 0)     return false;
    if(jet->chEmFrac  >= 0.80)  return false;
    }
    return true;
  }
  return false;
}
// ELECTRONS
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
// 2018=2017
//--------------------------------------------------------------------------------------------------
bool passEleVetoSel(const baconhep::TElectron *electron, const double rho,std::string iLabel)
{  
  if(iLabel==label2016){
    if(electron->isConv) return false;
    double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea2016(electron->eta)) );
    if(fabs(electron->scEta)<1.479) {
      if(iso >= 0.175*(electron->pt)) return false;
      if(electron->sieie              >= 0.01150)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00749)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.22800)                        return false;
      if(electron->hovere             >= 0.35600)                        return false;
      if(fabs(1.0 - electron->eoverp) >= 0.29900*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  2)                              return false;
    } else {
      if(iso >= 0.159*(electron->pt)) return false;
      if(electron->sieie              >= 0.03700)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00895)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.21300)                        return false;
      if(electron->hovere             >= 0.21100)                        return false;
      if(fabs(1.0 - electron->eoverp) >= 0.15000*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  3)                              return false;
    }
    return true;
  }
  if(iLabel==label2017||iLabel==label2018){
    if(electron->isConv) return false;
    double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea2017(electron->eta)) );
    if(fabs(electron->scEta)<=1.479) {
      if(iso >= 0.198*(electron->pt)) return false;
      if(electron->sieie              >= 0.01260)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00463)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.14800)                        return false;
      if(electron->hovere             >= 0.05 + 1.16/(electron->ecalEnergy) + 0.0324*rho/(electron->ecalEnergy))                                          return false;
      if(fabs(1.0 - electron->eoverp) >= 0.209*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  2)                              return false;
    } else {
      if(iso >= 0.203*(electron->pt)) return false;
      if(electron->sieie              >= 0.04570)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00814)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.19000)                        return false;
      if(electron->hovere             >= 0.05 + 2.54/(electron->ecalEnergy) + 0.183*rho/(electron->ecalEnergy))                                           return false;
      if(fabs(1.0 - electron->eoverp) >= 0.132*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  3)                              return false;
    }
    return true;
  }
  return false;
}
//--------------------------------------------------------------------------------------------------
bool passEleLooseSel(const baconhep::TElectron *electron, const double rho,std::string iLabel)
{
  if(iLabel==label2016){
    if(electron->isConv) return false;
    double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea2016(electron->eta)) );
    if(fabs(electron->scEta)<1.479) {
      if(iso >= 0.0994*(electron->pt)) return false;
      if(electron->sieie              >= 0.01100)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00477)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.22200)                        return false;
      if(electron->hovere             >= 0.29800)                        return false;
      if(fabs(1.0 - electron->eoverp) >= 0.24100*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  1)                              return false;
    } else {
      if(iso >= 0.107*(electron->pt)) return false;
      if(electron->sieie              >= 0.03140)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00868)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.21300)                        return false;
      if(electron->hovere             >= 0.10100)                        return false;
      if(fabs(1.0 - electron->eoverp) >= 0.14000*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  1)                              return false;
    }
    return true;
  }
  if(iLabel==label2017||iLabel==label2018){
    if(electron->isConv) return false;
    double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea2017(electron->eta)) );
    if(fabs(electron->scEta)<=1.479) {
      if(iso >= 0.133*(electron->pt)) return false;
      if(electron->sieie              >= 0.0112)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00377)                       return false;
      if(fabs(electron->dPhiIn)       >= 0.0884)                        return false;
      if(electron->hovere             >= 0.05 + 1.16/(electron->ecalEnergy) + 0.0324*rho/(electron->ecalEnergy))                                         return false;
      if(fabs(1.0 - electron->eoverp) >= 0.112*(electron->ecalEnergy))  return false;
      if(electron->nMissingHits       >  1)                             return false;
    } else {
      if(iso >= 0.146*(electron->pt)) return false;
      if(electron->sieie              >= 0.0425)                         return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00674)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.169)                          return false;
      if(electron->hovere             >= 0.0441 + 2.54/(electron->ecalEnergy) + 0.183*rho/(electron->ecalEnergy))                                         return false;
      if(fabs(1.0 - electron->eoverp) >= 0.108*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  1)                              return false;
    }
    return true;
  }
  return false;
}
bool passEleTightSel(const baconhep::TElectron *electron, const double rho,std::string iLabel)
{
  if(iLabel==label2016){
    if(electron->isConv) return false;
    double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea2016(electron->eta)) );
    if(fabs(electron->scEta)<1.479) {
      if(iso >= 0.0588*(electron->pt)) return false;
      if(electron->sieie              >= 0.00998)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00308)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.08160)                        return false;
      if(electron->hovere             >= 0.04140)                        return false;
      if(fabs(1.0 - electron->eoverp) >= 0.12900*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  1)                              return false;
    } else {
      if(iso >= 0.0571*(electron->pt)) return false;
      if(electron->sieie              >= 0.02920)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00605)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.03940)                        return false;
      if(electron->hovere             >= 0.06410)                        return false;
      if(fabs(1.0 - electron->eoverp) >= 0.12900*(electron->ecalEnergy)) return false;
      if(electron->nMissingHits       >  1)                              return false;
    }
    return true;
  }
  if(iLabel==label2017||iLabel==label2018){
    if(electron->isConv) return false;
    double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea2017(electron->eta)) );
    if(fabs(electron->scEta)<1.479) {
      if(iso >= 0.0287*(electron->pt)) return false;
      if(electron->sieie              >= 0.01040)                        return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00255)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.02200)                        return false;
      if(electron->hovere             >= 0.026 + 1.15/(electron->ecalEnergy) + 0.0324*rho/(electron->ecalEnergy))                                         return false;
      if(fabs(1.0 - electron->eoverp) >= 0.159*(electron->ecalEnergy))   return false;
      if(electron->nMissingHits       >  1)                              return false;
    } else {
      if(iso >= 0.0445*(electron->pt)) return false;
      if(electron->sieie              >= 0.0353)                         return false;
      if(fabs(electron->dEtaInSeed)   >= 0.00501)                        return false;
      if(fabs(electron->dPhiIn)       >= 0.0236)                         return false;
      if(electron->hovere             >= 0.0188 + 2.06/(electron->ecalEnergy) + 0.183*rho/(electron->ecalEnergy))                                         return false;
      if(fabs(1.0 - electron->eoverp) >= 0.0197*(electron->ecalEnergy))  return false;
      if(electron->nMissingHits       >  1)                              return false;
    }
    return true;
  }
  return false;
}
//--------------------------------------------------------------------------------------------------
// HEEP Electrons
// https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2#Selection_Cuts_HEEP_V7_0_recomme
// 2016=2017=2018
bool passEleHEEPSel(const baconhep::TElectron *electron, const double rho, const double met,std::string iLabel)
{
  if(fabs(electron->scEta)<1.442) {
    if(met                             <= 35)                                              return false;
    if(!(electron->typeBits & baconhep::kEcalDriven))                                      return false;
    if(fabs(electron->dEtaInSeed)      >= 0.00400)                                         return false;
    if(fabs(electron->dPhiIn)          >= 0.06000)                                         return false;
    if(electron->hovere                >= (1/(electron->ecalEnergy))+0.05)                 return false;
    if((electron->e2x5/electron->e5x5) <= 0.94 || (electron->e1x5/electron->e5x5) <= 0.83) return false;
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
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
// 2018=2017 for now
//--------------------------------------------------------------------------------------------------
bool passPhoLooseSel(const baconhep::TPhoton *photon, const double rho,std::string iLabel)
{
  if(iLabel==label2016){
    double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea2016(photon->scEta, 0), (double)0.);
    double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea2016(photon->scEta, 1), (double)0.);
    double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea2016(photon->scEta, 2), (double)0.);
    
    double maxneuHadIsoB = 10.910+0.0148*photon->pt+0.000017*photon->pt*photon->pt;
    double maxneuHadIsoE = 5.931+0.0163*photon->pt+0.000014*photon->pt*photon->pt;

    if(fabs(photon->scEta) <= 1.479) {
      if(photon->hovere   > 0.0597)                                      return false;
      if(photon->sieie    > 0.01031)                                     return false;
      if(chHadIso         > 1.295)                                       return false;
      if(neuHadIso        > maxneuHadIsoB)                               return false;                                  
      if(phoIso           > 3.630 + 0.0047*photon->pt)                   return false;
      
    } else {
      if(photon->hovere   > 0.0481)                                      return false;
      if(photon->sieie    > 0.03013)                                     return false;
      if(chHadIso         > 1.011)                                       return false;
      if(neuHadIso        > maxneuHadIsoE)                               return false;
      if(phoIso           > 6.641 + 0.0034*photon->pt)                   return false;
    }
    return true;
  }
  if(iLabel==label2017||iLabel==label2018){
    double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea2017(photon->scEta, 0), (double)0.);
    double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea2017(photon->scEta, 1), (double)0.);
    double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea2017(photon->scEta, 2), (double)0.);
    
    double maxneuHadIsoB = 24.032+0.01512*photon->pt+0.00002259*photon->pt*photon->pt;
    double maxneuHadIsoE = 19.722+0.0117*photon->pt+0.000023*photon->pt*photon->pt;

    if(fabs(photon->scEta) <= 1.479) {
      if(photon->hovere   > 0.04596)                                     return false;
      if(photon->sieie    > 0.0106)                                      return false;
      if(chHadIso         > 1.694)                                       return false;
      if(neuHadIso        > maxneuHadIsoB)                               return false;                                  
      if(phoIso           > 2.876 + 0.004017*photon->pt)                 return false;
    } else {
      if(photon->hovere   > 0.0590)                                      return false;
      if(photon->sieie    > 0.0272)                                      return false;
      if(chHadIso         > 2.089)                                       return false;
      if(neuHadIso        > maxneuHadIsoE)                               return false;
      if(phoIso           > 4.162 + 0.0037*photon->pt)                   return false;
    }
    return true;
  }
  return false;
}
//--------------------------------------------------------------------------------------------------
bool passPhoMediumSel(const baconhep::TPhoton *photon, const double rho,std::string iLabel)
{
  if(iLabel==label2016){
    double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea2016(photon->scEta, 0), (double)0.);
    double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea2016(photon->scEta, 1), (double)0.);
    double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea2016(photon->scEta, 2), (double)0.);
    
    double maxneuHadIsoB = 2.725+0.0148*photon->pt+0.000017*photon->pt*photon->pt;
    double maxneuHadIsoE = 1.715+0.0163*photon->pt+0.000014*photon->pt*photon->pt;
    
    if(fabs(photon->scEta) <= 1.479) {
      if(photon->hovere   > 0.0396)                                      return false;
      if(photon->sieie    > 0.01022)                                     return false;
      if(chHadIso         > 0.441)                                       return false;
      if(neuHadIso        > maxneuHadIsoB)                               return false;
      if(phoIso           > 2.571 + 0.0047*photon->pt)                   return false;
    } else {
      if(photon->hovere   > 0.0219)                                      return false;
      if(photon->sieie    > 0.03001)                                     return false;
      if(chHadIso         > 0.442)                                       return false;
      if(neuHadIso        > maxneuHadIsoE)                               return false;
      if(phoIso           > 3.863 + 0.0034*photon->pt)                   return false;
    }
    return true;
  }
  if(iLabel==label2017||iLabel==label2018){
    double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea2017(photon->scEta, 0), (double)0.);
    double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea2017(photon->scEta, 1), (double)0.);
    double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea2017(photon->scEta, 2), (double)0.);
    
    double maxneuHadIsoB = 1.189+0.01512*photon->pt+0.00002259*photon->pt*photon->pt;
    double maxneuHadIsoE = 2.718+0.0117*photon->pt+0.000023*photon->pt*photon->pt;
    
    if(fabs(photon->scEta) <= 1.479) {
      if(photon->hovere   > 0.02197)                                     return false;
      if(photon->sieie    > 0.01015)                                     return false;
      if(chHadIso         > 1.141)                                       return false;
      if(neuHadIso        > maxneuHadIsoB)                               return false;
      if(phoIso           > 2.08 + 0.004017*photon->pt)                  return false;
    } else {
      if(photon->hovere   > 0.0326)                                      return false;
      if(photon->sieie    > 0.0272)                                      return false;
      if(chHadIso         > 1.051)                                       return false;
      if(neuHadIso        > maxneuHadIsoE)                               return false;
      if(phoIso           > 3.867 + 0.0037*photon->pt)                   return false;
    }
    return true;
  }
  return false;
}
//--------------------------------------------------------------------------------------------------
bool passPhoTightSel(const baconhep::TPhoton *photon, const double rho, std::string iLabel)
{
  if(iLabel==label2016){
    double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea2016(photon->scEta, 0), (double)0.);
    double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea2016(photon->scEta, 1), (double)0.);
    double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea2016(photon->scEta, 2), (double)0.);
    
    double maxneuHadIsoB = 0.264+0.0148*photon->pt+0.000017*photon->pt*photon->pt;
    double maxneuHadIsoE = 0.586+0.0163*photon->pt+0.000014*photon->pt*photon->pt;
    
    if(fabs(photon->scEta) <= 1.479) {
      if(photon->hovere   > 0.0269)                                      return false;
      if(photon->sieie    > 0.00994)                                     return false;
      if(chHadIso         > 0.202)                                       return false;
      if(neuHadIso        > maxneuHadIsoB)                               return false;
      if(phoIso           > 2.362 + 0.0047*photon->pt)                   return false;
    } else {
      if(photon->hovere   > 0.0213)                                      return false;
      if(photon->sieie    > 0.03000)                                     return false;
      if(chHadIso         > 0.034)                                       return false;
      if(neuHadIso        > maxneuHadIsoE)                               return false;
      if(phoIso           > 2.617 + 0.0034*photon->pt)                   return false;
    }
    return true;
  }
  if(iLabel==label2017||iLabel==label2018){
    double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea2017(photon->scEta, 0), (double)0.);
    double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea2017(photon->scEta, 1), (double)0.);
    double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea2017(photon->scEta, 2), (double)0.);
    
    double maxneuHadIsoB = 0.317+0.01512*photon->pt+0.0000229*photon->pt*photon->pt;
    double maxneuHadIsoE = 2.716+0.0117*photon->pt+0.000023*photon->pt*photon->pt;
    
    if(fabs(photon->scEta) <= 1.479) {
      if(photon->hovere   > 0.02148)                                     return false;
      if(photon->sieie    > 0.00996)                                     return false;
      if(chHadIso         > 0.65)                                        return false;
      if(neuHadIso        > maxneuHadIsoB)                               return false;
      if(phoIso           > 2.044 + 0.004017*photon->pt)                 return false;
    } else {
      if(photon->hovere   > 0.0321)                                      return false;
      if(photon->sieie    > 0.0271)                                      return false;
      if(chHadIso         > 0.517)                                       return false;
      if(neuHadIso        > maxneuHadIsoE)                               return false;
      if(phoIso           > 3.032 + 0.0037*photon->pt)                   return false;
    }
    return true;
  }
  return false;
}
//--------------------------------------------------------------------------------------------------
// Effective Areas
//https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
extern double eleEffArea2016(const double eta)
{
  if     (fabs(eta) >= 0.0000 && fabs(eta) < 1.0000) { return 0.1703; }
  else if(fabs(eta) >= 1.0000 && fabs(eta) < 1.4790) { return 0.1715; }
  else if(fabs(eta) >= 1.4790 && fabs(eta) < 2.0000) { return 0.1213; }
  else if(fabs(eta) >= 2.0000 && fabs(eta) < 2.2000) { return 0.1230; }
  else if(fabs(eta) >= 2.2000 && fabs(eta) < 2.3000) { return 0.1635; }
  else if(fabs(eta) >= 2.3000 && fabs(eta) < 2.4000) { return 0.1937; }
  else                                               { return 0.2393; }
}
// https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
extern double eleEffArea2017(const double eta)
{
  if     (fabs(eta) >= 0.0000 && fabs(eta) < 1.0000) { return 0.1440; }
  else if(fabs(eta) >= 1.0000 && fabs(eta) < 1.4790) { return 0.1562; }
  else if(fabs(eta) >= 1.4790 && fabs(eta) < 2.0000) { return 0.1032; }
  else if(fabs(eta) >= 2.0000 && fabs(eta) < 2.2000) { return 0.0859; }
  else if(fabs(eta) >= 2.2000 && fabs(eta) < 2.3000) { return 0.1116; }
  else if(fabs(eta) >= 2.3000 && fabs(eta) < 2.4000) { return 0.1321; }
  else                                               { return 0.1654; }
}
//--------------------------------------------------------------------------------------------------
//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
extern double phoEffArea2016(const double eta, const int type)
{
  const int kCH_HAD  = 0;
  const int kNEU_HAD = 1;
  const int kPHOTON  = 2;

  if(type==kCH_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0360; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0377; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0306; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0283; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0254; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0217; }
    else                                             { return 0.0167; }

  } else if(type==kNEU_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0597; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0807; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0629; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0197; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0184; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0284; }
    else                                             { return 0.0591; }

  } else if(type==kPHOTON) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.1210; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.1107; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0699; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.1056; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.1457; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.1719; }
    else                                             { return 0.1998; }

  } else { assert(0); }
}
extern double phoEffArea2017(const double eta, const int type)
{
  const int kCH_HAD  = 0;
  const int kNEU_HAD = 1;
  const int kPHOTON  = 2;

  if(type==kCH_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0112; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0108; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0106; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.01002; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0098; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0089; }
    else                                             { return 0.0087; }

  } else if(type==kNEU_HAD) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.0668; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.1054; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0786; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0233; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.0078; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.0028; }
    else                                             { return 0.0137; }

  } else if(type==kPHOTON) {
    if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.1113; }
    else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.0953; }
    else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.0619; }
    else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.0837; }
    else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.1070; }
    else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.1212; }
    else                                             { return 0.1466; }

  } else { assert(0); }
}
// MUONS
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2 
//-------------------------------------------------------------------------------------------------
bool passMuonLooseSel(const baconhep::TMuon *muon,std::string iLabel)
{
  if(!(muon->pogIDBits & baconhep::kPOGLooseMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;

  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonMediumSel(const baconhep::TMuon *muon,std::string iLabel)
{
  if(!(muon->pogIDBits & baconhep::kPOGMediumMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;

  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonTightSel(const baconhep::TMuon *muon,std::string iLabel)
{
  if(!(muon->pogIDBits & baconhep::kPOGTightMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.15*(muon->pt)) return false;
  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonSoftSel(const baconhep::TMuon *muon,std::string iLabel)
{
  if(!(muon->pogIDBits & baconhep::kPOGSoftMuon)) return false;

  // PF-isolation with Delta-beta correction
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;
  return true;
}
//-------------------------------------------------------------------------------------------------
bool passMuonHighPtSel(const baconhep::TMuon *muon,std::string iLabel)
{
  if(!(muon->pogIDBits & baconhep::kPOGHighPtMuon)) return false;

  // PF-isolation with Delta-beta correction                                                                                  
  double iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso), double(0));
  if(iso >= 0.25*(muon->pt)) return false;
  return true;
}
// TAUS https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#2017v2_discriminators
//--------------------------------------------------------------------------------------------------
bool passTauSel(const baconhep::TTau *tau, std::string iLabel)
{
  if(!(tau->hpsDisc & baconhep::kByDecayModeFinding)) return false;
  if(tau->rawIso3Hits > 4.5)                           return false;

  return true;
}
// old decay mode finding? + tight MVA isolation + loose anti-electron MVA + tight anti-muon cut
bool passTauTightSel(const baconhep::TTau *tau, std::string iLabel)
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
//--------------------------------------------------------------------------------------------------
double deltaR2(double iEta, double iPhi, double jEta, double jPhi) {
  double pDEta = iEta - jEta;
  double pDPhi = iPhi - jPhi;
  if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
  return pDPhi*pDPhi + pDEta*pDEta;
}
