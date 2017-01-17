#include "../include/VJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>

using namespace baconhep;

VJetLoader::VJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS,std::string iAddJetCHS,int iN) { 
  fVJets         = new TClonesArray("baconhep::TJet");
  fVAddJets      = new TClonesArray("baconhep::TAddJet");
  fVJetsCHS      = new TClonesArray("baconhep::TJet");
  fVAddJetsCHS   = new TClonesArray("baconhep::TAddJet");

  iTree->SetBranchAddress(iJet.c_str(),       &fVJets);
  iTree->SetBranchAddress(iAddJet.c_str(),    &fVAddJets);
  iTree->SetBranchAddress(iJetCHS.c_str(),    &fVJetsCHS);
  iTree->SetBranchAddress(iAddJetCHS.c_str(), &fVAddJetsCHS);

  fVJetBr        = iTree->GetBranch(iJet.c_str());
  fVAddJetBr     = iTree->GetBranch(iAddJet.c_str());
  fVJetBrCHS     = iTree->GetBranch(iJetCHS.c_str());
  fVAddJetBrCHS  = iTree->GetBranch(iAddJetCHS.c_str());

  fN = iN;

  r = new TRandom3(1988);

  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
}
VJetLoader::~VJetLoader() { 
  delete fVJets;
  delete fVJetBr;
  delete fVAddJets;
  delete fVAddJetBr;
  delete fVJetsCHS;
  delete fVJetBrCHS;
  delete fVAddJetsCHS;
  delete fVAddJetBrCHS;
}
void VJetLoader::reset() { 
  fNLooseVJets        = 0;
  fNTightVJets        = 0;  
  for(int i0 = 0; i0 < int(fisTightVJet.size()); i0++) fisTightVJet[i0] = -999;
  selectedVJets.clear();
  fLooseVJets.clear();
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
  resetZprime();
}
void VJetLoader::resetCHS() {
  fNLooseVJetsCHS     = 0;
  fNTightVJetsCHS     = 0;
  for(int i0 = 0; i0 < int(fisTightVJetCHS.size()); i0++) fisTightVJetCHS[i0] = -999;
  selectedVJetsCHS.clear();
  fLooseVJetsCHS.clear();
  fdoublecsvCHS.clear();
  fdoublesubCHS.clear();
  fptCHS.clear();
  fetaCHS.clear();
  fphiCHS.clear();
}
void VJetLoader::resetDoubleB() {
  fLooseVJetsByDoubleB.clear();
  selectedVJetsByDoubleB.clear();
  fLooseVJetsCHSByDoubleB.clear();
  selectedVJetsCHSByDoubleB.clear();
}
void VJetLoader::resetZprime() {
  fvSize              = 999;
  fvMatching          = 999;
  fisHadronicV        = 0;
  fRatioPt =0;
}
void VJetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();
  fLabels.clear();
  fLabels.push_back("mass");
  fLabels.push_back("csv");
  fLabels.push_back("CHF");
  fLabels.push_back("NHF");
  fLabels.push_back("NEMF");
  fLabels.push_back("tau21");
  fLabels.push_back("tau32");
  fLabels.push_back("msd");
  fLabels.push_back("rho");
  fLabels.push_back("minsubcsv");
  fLabels.push_back("maxsubcsv");
  fLabels.push_back("doublecsv");
  fLabels.push_back("doublesub");
  fLabels.push_back("ptraw");
  fLabels.push_back("genpt");
  fLabels.push_back("e2_b1"); // Correlation function inputs beta=1
  fLabels.push_back("e3_b1");
  fLabels.push_back("e3_v1_b1");
  fLabels.push_back("e3_v2_b1");
  fLabels.push_back("e4_v1_b1");
  fLabels.push_back("e4_v2_b1");
  fLabels.push_back("e2_b2"); // Correlation function inputs beta=2
  fLabels.push_back("e3_b2");
  fLabels.push_back("e3_v1_b2");
  fLabels.push_back("e3_v2_b2");
  fLabels.push_back("e4_v1_b2");
  fLabels.push_back("e4_v2_b2");
  fLabels.push_back("e2_sdb1"); // Correlation function inputs beta=1 soft-dropped 
  fLabels.push_back("e3_sdb1");
  fLabels.push_back("e3_v1_sdb1");
  fLabels.push_back("e3_v2_sdb1");
  fLabels.push_back("e4_v1_sdb1");
  fLabels.push_back("e4_v2_sdb1");
  fLabels.push_back("e2_sdb2"); // Correlation function inputs beta=2 soft-dropped 
  fLabels.push_back("e3_sdb2");
  fLabels.push_back("e3_v1_sdb2");
  fLabels.push_back("e3_v2_sdb2");
  fLabels.push_back("e4_v1_sdb2");
  fLabels.push_back("e4_v2_sdb2");
  fLabels.push_back("N2sdb1"); // 2-prong ECFs observables
  fLabels.push_back("N2sdb2");
  fLabels.push_back("M2sdb1");
  fLabels.push_back("M2sdb2");
  fLabels.push_back("D2sdb1");
  fLabels.push_back("D2sdb2");
  fLabels.push_back("N2b1");
  fLabels.push_back("N2b2");
  fLabels.push_back("M2b1");
  fLabels.push_back("M2b2");
  fLabels.push_back("D2b1");
  fLabels.push_back("D2b2");
  fLabels.push_back("pt_cen");
  fLabels.push_back("pt_JESUp");
  fLabels.push_back("pt_JESDown");
  fLabels.push_back("pt_JERUp");
  fLabels.push_back("pt_JERDown");

  std::stringstream pSNJ;   pSNJ << "n" << iJetLabel << "s";
  fTree = iTree;
  for(int i0 = 0; i0 < fN*4.;                    i0++) {double pVar = 0; fVars.push_back(pVar);} // declare array of vars
  for(int i0 = 0; i0 < fN*(int(fLabels.size())); i0++) {double pVar = 0; fVars.push_back(pVar);} 
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars);                                                 // from Utils.cc => fN =1 *_pt,*_eta,*_phi for vjet0 (3*1=3)
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars,fN*3,fLabels);
  fTree->Branch(pSNJ.str().c_str() ,&fNLooseVJets         ,(pSNJ.str()+"/I").c_str());  
  for(int i0 = 0; i0 < fN;                    i0++) fisTightVJet.push_back(-999);
  for(int i0 = 0; i0 < fN;                    i0++) {
    std::stringstream pSTJ;   pSTJ << iJetLabel << i0 << "_isTightVJet";
    fTree->Branch(pSTJ.str().c_str() ,&fisTightVJet[i0]         ,(pSTJ.str()+"/I").c_str());
  }
}
void VJetLoader::setupTreeZprime(TTree *iTree, std::string iJetLabel) {
  resetZprime();
  std::stringstream pSiV;   pSiV << iJetLabel << "0_isHadronicV";
  std::stringstream pSVM;   pSVM << iJetLabel << "0_vMatching";
  std::stringstream pSVS;   pSVS << iJetLabel << "0_vSize";
  std::stringstream pSpF;   pSpF << iJetLabel << "0_partonFlavor";
  std::stringstream pShF;   pShF << iJetLabel << "0_hadronFlavor";
  std::stringstream pSnC;   pSnC << iJetLabel << "0_nCharged";
  std::stringstream pSnN;   pSnN << iJetLabel << "0_nNeutrals";
  std::stringstream pSnP;   pSnP << iJetLabel << "0_nParticles";
  std::stringstream pSratio; pSratio << iJetLabel << "0_ratioCA15_04";

  fTree = iTree;
  fTree->Branch(pSiV.str().c_str() ,&fisHadronicV         ,(pSiV.str()+"/I").c_str());
  fTree->Branch(pSVM.str().c_str() ,&fvMatching           ,(pSVM.str()+"/D").c_str());
  fTree->Branch(pSVS.str().c_str() ,&fvSize               ,(pSVS.str()+"/D").c_str());
  fTree->Branch(pSpF.str().c_str() ,&fpartonFlavor        ,(pSpF.str()+"/I").c_str());
  fTree->Branch(pShF.str().c_str() ,&fhadronFlavor        ,(pShF.str()+"/I").c_str());
  fTree->Branch(pSnC.str().c_str() ,&fnCharged            ,(pSnC.str()+"/I").c_str());
  fTree->Branch(pSnN.str().c_str() ,&fnNeutrals           ,(pSnN.str()+"/I").c_str());
  fTree->Branch(pSnP.str().c_str() ,&fnParticles          ,(pSnP.str()+"/I").c_str());
  fTree->Branch(pSratio.str().c_str() ,&fRatioPt          ,(pSratio.str()+"/D").c_str());
}
void VJetLoader::setupTreeCHS(TTree *iTree, std::string iJetLabel) {
  resetCHS();  
  fTree = iTree;
  for(int i0 = 0; i0 < fN; i0++) {
    fdoublecsvCHS.push_back(-999);
    fdoublesubCHS.push_back(-999);
    fptCHS.push_back(-999);
    fetaCHS.push_back(-999);
    fphiCHS.push_back(-999);
    fisTightVJetCHS.push_back(-999);
  }
  for(int i0 = 0; i0 < fN; i0++) {
    std::stringstream pSdc;   pSdc << iJetLabel << i0 << "_doublecsv";
    std::stringstream pSds;   pSds << iJetLabel << i0 << "_doublesub";    
    std::stringstream pSpt;   pSpt << iJetLabel << i0 << "_pt";    
    std::stringstream pSeta;   pSeta << iJetLabel << i0 << "_eta";  
    std::stringstream pSphi;   pSphi << iJetLabel << i0 << "_phi";
    std::stringstream pSTJ;   pSTJ << iJetLabel << i0 << "_isTightVJet";    
    fTree->Branch(pSTJ.str().c_str() ,&fisTightVJetCHS[i0]         ,(pSTJ.str()+"/I").c_str());
    fTree->Branch(pSdc.str().c_str() ,&fdoublecsvCHS.at(i0)        ,(pSdc.str()+"/D").c_str());
    fTree->Branch(pSds.str().c_str() ,&fdoublesubCHS.at(i0)       ,(pSds.str()+"/D").c_str());
    fTree->Branch(pSpt.str().c_str() ,&fptCHS.at(i0)       ,(pSpt.str()+"/D").c_str());
    fTree->Branch(pSeta.str().c_str() ,&fetaCHS.at(i0)       ,(pSeta.str()+"/D").c_str());
    fTree->Branch(pSphi.str().c_str() ,&fphiCHS.at(i0)       ,(pSphi.str()+"/D").c_str());
  }
}
void VJetLoader::load(int iEvent) { 
  fVJets       ->Clear();
  fVJetBr      ->GetEntry(iEvent);
  fVAddJets    ->Clear();
  fVAddJetBr   ->GetEntry(iEvent);
  fVJetsCHS    ->Clear();
  fVJetBrCHS   ->GetEntry(iEvent);
  fVAddJetsCHS ->Clear();
  fVAddJetBrCHS->GetEntry(iEvent);
}
void VJetLoader::selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho){
  reset(); 
  int lCount(0), lCountT(0);
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    TJet *pVJet = (TJet*)((*fVJets)[i0]);
    if(pVJet->pt        <=  350)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJets);
    lCount++;

    if(!passJetTightLepVetoSel(pVJet))                                     continue;
    lCountT++;
  }
  addVJet(fLooseVJets,selectedVJets);

  for  (int i0 = 0; i0 < int(selectedVJets.size()); i0++) { 
    if(passJetTightLepVetoSel(fLooseVJets[i0])) fisTightVJet[i0] = 1;
  }
  fNLooseVJets = lCount;
  fNTightVJets = lCountT;

  fillJet( fN,fLooseVJets,fVars);
  fillVJet(fN,fLooseVJets,fVars,iRho);
}

void VJetLoader::selectVJetsByDoubleBCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho){
  // first do Puppi jets (pT > 500 GeV)
  resetDoubleB(); 
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    TJet *pVJet = (TJet*)((*fVJets)[i0]);
    if(pVJet->pt        <=  500)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJetsByDoubleB);
    if(!passJetTightLepVetoSel(pVJet))                                     continue;
  }
  addVJet(fLooseVJetsByDoubleB,selectedVJetsByDoubleB);

  // now do CHS jets (pT > 400 GeV)
  for  (int i0 = 0; i0 < fVJetsCHS->GetEntriesFast(); i0++) {
    TJet *pVJet = (TJet*)((*fVJetsCHS)[i0]);
    if(pVJet->pt        <=  400)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJetsCHSByDoubleB);
    if(!passJetTightLepVetoSel(pVJet))                                     continue;
  }
  addVJet(fLooseVJetsCHSByDoubleB,selectedVJetsCHSByDoubleB);
  
  std::vector<int> indexCHS;
  std::vector<double> doubleBCHS;
  for (int i0 = 0; i0 < int(selectedVJetsByDoubleB.size()); i0++) {    
    int iCHSJet = -999;
    double dbCHSJet = -999;
    iCHSJet = getMatchedCHSJetIndex(selectedVJetsCHSByDoubleB,selectedVJetsByDoubleB[i0],0.8);
    if (iCHSJet > -999) {
      //selectedVJetsCHS[iCHSJet];      
      TAddJet *pAddJetCHS = getAddJetCHS(fLooseVJetsCHSByDoubleB[iCHSJet]);  
      dbCHSJet = pAddJetCHS->doublecsv;      
    }
    std::cout << "index CHS    = " << iCHSJet << std::endl;
    std::cout << "double-b CHS = " << dbCHSJet << std::endl;
    indexCHS.push_back(iCHSJet);
    doubleBCHS.push_back(dbCHSJet);
  }

  
}
void VJetLoader::selectVJetsCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho){
  resetCHS();
  int lCount(0), lCountT(0);
  for  (int i0 = 0; i0 < fVJetsCHS->GetEntriesFast(); i0++) {
    TJet *pVJet = (TJet*)((*fVJetsCHS)[i0]);
    if(pVJet->pt        <=  300)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJetsCHS);
    lCount++;

    if(!passJetTightLepVetoSel(pVJet))                                     continue;
    lCountT++;
  }
  addVJet(fLooseVJetsCHS,selectedVJetsCHS);

  
  for  (int i0 = 0; i0 < int(selectedVJetsCHS.size()); i0++) { 
    if(passJetTightLepVetoSel(fLooseVJetsCHS[i0])) fisTightVJetCHS[i0] = 1;
  }
  fNLooseVJetsCHS = lCount;
  fNTightVJetsCHS = lCountT;
}
void VJetLoader::fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho){ 
  int lBase = 3.*fN;
  int lMin = iObjects.size();
  int lNLabel = int(fLabels.size());
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    TAddJet *pAddJet = getAddJet(iObjects[i0]);
    iVals[lBase+i0*lNLabel+0]  = iObjects[i0]->mass;
    iVals[lBase+i0*lNLabel+1]  = iObjects[i0]->csv;
    iVals[lBase+i0*lNLabel+2]  = iObjects[i0]->chHadFrac;
    iVals[lBase+i0*lNLabel+3]  = iObjects[i0]->neuHadFrac;
    iVals[lBase+i0*lNLabel+4]  = iObjects[i0]->neuEmFrac;
    iVals[lBase+i0*lNLabel+5]  = (pAddJet->tau2/pAddJet->tau1);
    iVals[lBase+i0*lNLabel+6]  = (pAddJet->tau3/pAddJet->tau2);
    iVals[lBase+i0*lNLabel+7]  = pAddJet->mass_sd0;
    iVals[lBase+i0*lNLabel+8]  = log((pAddJet->mass_sd0*pAddJet->mass_sd0)/iObjects[i0]->pt);
    iVals[lBase+i0*lNLabel+9]  = TMath::Min(pAddJet->sj1_csv,pAddJet->sj2_csv);
    iVals[lBase+i0*lNLabel+10] = TMath::Max(TMath::Max(pAddJet->sj1_csv,pAddJet->sj2_csv),TMath::Max(pAddJet->sj3_csv,pAddJet->sj4_csv));
    iVals[lBase+i0*lNLabel+11] = pAddJet->doublecsv;
    iVals[lBase+i0*lNLabel+12] = pAddJet->Double_sub;
    iVals[lBase+i0*lNLabel+13] = iObjects[i0]->ptRaw;
    iVals[lBase+i0*lNLabel+14] = iObjects[i0]->genpt;
    iVals[lBase+i0*lNLabel+15] = pAddJet->e2_b1;
    iVals[lBase+i0*lNLabel+16] = pAddJet->e3_b1;
    iVals[lBase+i0*lNLabel+17] = pAddJet->e3_v1_b1;
    iVals[lBase+i0*lNLabel+18] = pAddJet->e3_v2_b1;
    iVals[lBase+i0*lNLabel+19] = pAddJet->e4_v1_b1;
    iVals[lBase+i0*lNLabel+20] = pAddJet->e4_v2_b1;
    iVals[lBase+i0*lNLabel+21] = pAddJet->e2_b2;
    iVals[lBase+i0*lNLabel+22] = pAddJet->e3_b2;
    iVals[lBase+i0*lNLabel+23] = pAddJet->e3_v1_b2;
    iVals[lBase+i0*lNLabel+24] = pAddJet->e3_v2_b2;
    iVals[lBase+i0*lNLabel+25] = pAddJet->e4_v1_b2;
    iVals[lBase+i0*lNLabel+26] = pAddJet->e4_v2_b2;
    iVals[lBase+i0*lNLabel+27] = pAddJet->e2_sdb1;
    iVals[lBase+i0*lNLabel+28] = pAddJet->e3_sdb1;
    iVals[lBase+i0*lNLabel+29] = pAddJet->e3_v1_sdb1;
    iVals[lBase+i0*lNLabel+30] = pAddJet->e3_v2_sdb1;
    iVals[lBase+i0*lNLabel+31] = pAddJet->e4_v1_sdb1;
    iVals[lBase+i0*lNLabel+32] = pAddJet->e4_v2_sdb1;
    iVals[lBase+i0*lNLabel+33] = pAddJet->e2_sdb2;
    iVals[lBase+i0*lNLabel+34] = pAddJet->e3_sdb2;
    iVals[lBase+i0*lNLabel+35] = pAddJet->e3_v1_sdb2;
    iVals[lBase+i0*lNLabel+36] = pAddJet->e3_v2_sdb2;
    iVals[lBase+i0*lNLabel+37] = pAddJet->e4_v1_sdb2;
    iVals[lBase+i0*lNLabel+38] = pAddJet->e4_v2_sdb2;
    iVals[lBase+i0*lNLabel+39] = pAddJet->e3_v2_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+40] = pAddJet->e3_v2_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    iVals[lBase+i0*lNLabel+41] = pAddJet->e3_v1_sdb1/(pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+42] = pAddJet->e3_v1_sdb2/(pAddJet->e2_sdb2);
    iVals[lBase+i0*lNLabel+43] = pAddJet->e3_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+44] = pAddJet->e3_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    iVals[lBase+i0*lNLabel+45] = pAddJet->e3_v2_b1/(pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+46] = pAddJet->e3_v2_b2/(pAddJet->e2_b2*pAddJet->e2_b2);
    iVals[lBase+i0*lNLabel+47] = pAddJet->e3_v1_b1/(pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+48] = pAddJet->e3_v1_b2/(pAddJet->e2_b2);
    iVals[lBase+i0*lNLabel+49] = pAddJet->e3_b1/(pAddJet->e2_b1*pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+50] = pAddJet->e3_b2/(pAddJet->e2_b2*pAddJet->e2_b2*pAddJet->e2_b2);
    //JEC    
    double x1 = r->Gaus();
    double x2 = r->Gaus();
    double x3 = r->Gaus();
    double sf = getJerSF(iObjects[i0]->eta,0);
    double sfUp = getJerSF(iObjects[i0]->eta,1);
    double sfDown = getJerSF(iObjects[i0]->eta,-1);
    double sigma_MC = 0.05; // assume 5% for now
    double jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    double jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
    double jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;

    double unc = iObjects[i0]->unc;
    double jetCorrPt = iObjects[i0]->pt;
    double jetCorrPtSmear = (iObjects[i0]->pt)*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;
    /*
      std::cout << "VJet" << std::endl;
    std::cout << "i0 =" << i0 << std::endl;
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;
    std::cout << "x3 = " << x2 << std::endl;
    std::cout << "unc = " << unc << std::endl;
    std::cout << "sf = " << sf << std::endl;
    std::cout << "pt = " << jetCorrPt << std::endl;
    std::cout << "ptcen = " << jetCorrPtSmear << std::endl;
    std::cout << "jesup = " << jetPtJESUp << std::endl;
    std::cout << "jesdown = " << jetPtJESDown << std::endl;
    std::cout << "jerup = " << jetPtJERUp << std::endl;
    std::cout << "jerdown = " << jetPtJERDown << std::endl;
    */
    iVals[lBase+i0*lNLabel+51] = jetCorrPtSmear;
    iVals[lBase+i0*lNLabel+52] = jetPtJESUp;
    iVals[lBase+i0*lNLabel+53] = jetPtJESDown;
    iVals[lBase+i0*lNLabel+54] = jetPtJERUp;
    iVals[lBase+i0*lNLabel+55] = jetPtJERDown;

    fpartonFlavor   = iObjects[0]->partonFlavor;
    fhadronFlavor   = iObjects[0]->hadronFlavor;
    fnCharged       = iObjects[0]->nCharged;
    fnNeutrals      = iObjects[0]->nNeutrals;
    fnParticles     = iObjects[0]->nParticles;

  }
}
int VJetLoader::getMatchedCHSJetIndex(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR) {  
  TLorentzVector iJet1;
  int iJet1id(0), nmatched(0);
  float mindR = dR;
  for(int i0 = 0; i0 < int(iJets1.size()); i0++) {
    if ((iJets1[i0].DeltaR(iJet2) < mindR) && (fabs(iJets1[i0].Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))) {
      nmatched++;
      iJet1 = iJets1[i0];
      iJet1id = i0;
      mindR= iJets1[i0].DeltaR(iJet2);	
    }
  }
  if (nmatched >0 && (iJet1.DeltaR(iJet2) < dR) && (fabs(iJet1.Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))){    
    return iJet1id;
  }
  return -999;
}
void VJetLoader::matchJet(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex){
  TLorentzVector iJet1;
  int iJet1id(0), nmatched(0);
  float mindR = dR;
  for(int i0 = 0; i0 < int(iJets1.size()); i0++) {
    if ((iJets1[i0].DeltaR(iJet2) < mindR) && (fabs(iJets1[i0].Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))) {
      nmatched++;
      iJet1 = iJets1[i0];
      iJet1id = i0;
      mindR= iJets1[i0].DeltaR(iJet2);	
    }
  }
  if (nmatched >0 && (iJet1.DeltaR(iJet2) < dR) && (fabs(iJet1.Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))){
    fillVJetCHS(fLooseVJetsCHS[iJet1id], jIndex);
  }
}

void VJetLoader::matchJet15(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR){
  TLorentzVector iJet1;
  int nmatched(0);
  float mindR = dR;  
  for(int i0 = 0; i0 < int(iJets1.size()); i0++) {
    if ((iJets1[i0].DeltaR(iJet2) < mindR)) {
      nmatched++;
      iJet1 = iJets1[i0];
      mindR= iJets1[i0].DeltaR(iJet2);
    }
  }
  if (nmatched >0 && (iJet1.DeltaR(iJet2) < dR)){
    fRatioPt = iJet2.Pt()/iJet1.Pt();
  }
}
void VJetLoader::fillVJetCHS(TJet *iJet, int jIndex){
  TAddJet *pAddJet = getAddJetCHS(iJet);
  fdoublecsvCHS[jIndex] = double(pAddJet->doublecsv);
  fdoublesubCHS[jIndex] = double(pAddJet->Double_sub);
  fptCHS[jIndex] = double(iJet->pt);
  fetaCHS[jIndex] = double(iJet->eta);
  fphiCHS[jIndex] = double(iJet->phi);
}
TAddJet *VJetLoader::getAddJet(TJet *iJet) { 
  int lIndex = -1;
  TAddJet *lJet = 0; 
  for(int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    if((*fVJets)[i0] == iJet) { lIndex = i0; break;}
  }
  if(lIndex == -1) return 0;
  for  (int i0 = 0; i0 < fVAddJets->GetEntriesFast(); i0++) { 
    TAddJet *pJet = (TAddJet*)((*fVAddJets)[i0]);
    if(pJet->index == fabs(lIndex)) { lJet = pJet; break;}
  }
  return lJet;
}
TAddJet *VJetLoader::getAddJetCHS(TJet *iJet) {
  int lIndex = -1;
  TAddJet *lJet = 0;
  for(int i0 = 0; i0 < fVJetsCHS->GetEntriesFast(); i0++) {
    if((*fVJetsCHS)[i0] == iJet) { lIndex = i0; break;}
  }
  if(lIndex == -1) return 0;
  for  (int i0 = 0; i0 < fVAddJetsCHS->GetEntriesFast(); i0++) {
    TAddJet *pJet = (TAddJet*)((*fVAddJetsCHS)[i0]);
    if(pJet->index == fabs(lIndex)) { lJet = pJet; break;}
  }
  return lJet;
}


// Retrieve jet energy resolution scale factor as a function of eta
double VJetLoader::getJerSF( float eta, int nsigma) {
  double sf[3] = {1.0, 1.0, 1.0};
  if (fabs(eta) < 0.5) {
    sf[0] = 1.109; 
    sf[1] = 1.101; 
    sf[2] = 1.117;
  }
  else if (fabs(eta) >= 0.5 && fabs(eta) < 0.8) {
    sf[0] = 1.138;
    sf[1] = 1.125;
    sf[2] = 1.151;
  }
  else if (fabs(eta) >= 0.8 && fabs(eta) < 1.1) {
    sf[0] = 1.114;
    sf[1] = 1.101; 
    sf[2] = 1.127;
  }
  else if (fabs(eta) >= 1.1 && fabs(eta) < 1.3) {
    sf[0] = 1.123;
    sf[1] = 1.099;
    sf[2] = 1.147;
  }
  else if (fabs(eta) >= 1.3 && fabs(eta) < 1.7) {
    sf[0] = 1.084;
    sf[1] = 1.073;
    sf[2] = 1.095;
  }
  else if (fabs(eta) >= 1.7 && fabs(eta) < 1.9) {
    sf[0] = 1.082;
    sf[1] = 1.047;
    sf[2] = 1.117;
  }
  else if (fabs(eta) >= 1.9 && fabs(eta) < 2.1) {
    sf[0] = 1.140;
    sf[1] = 1.093;
    sf[2] = 1.187;
  }
  else if (fabs(eta) >= 2.1 && fabs(eta) < 2.3) {
    sf[0] = 1.067;
    sf[1] = 1.014;
    sf[2] = 1.120;
  }
  else if (fabs(eta) >= 2.3 && fabs(eta) < 2.5) {
    sf[0] = 1.177;
    sf[1] = 1.136;
    sf[2] = 1.218;
  }
  else if (fabs(eta) >= 2.5 && fabs(eta) < 2.8) {
    sf[0] = 1.364;
    sf[1] = 1.325;
    sf[2] = 1.403;
  }
  else if (fabs(eta) >= 2.8 && fabs(eta) < 3.0) {
    sf[0] = 1.857;
    sf[1] = 1.786;
    sf[2] = 1.928;
  }
  else if (fabs(eta) >= 3.0 && fabs(eta) < 3.2) {
    sf[0] = 1.328;
    sf[1] = 1.306;
    sf[2] = 1.350;
  }
  else if (fabs(eta) >= 3.2 && fabs(eta) < 4.7) {
    sf[0] = 1.160;
    sf[1] = 1.131;
    sf[2] = 1.189;
  }

  if (nsigma == 0) {
    return sf[0];
  }
  else if (nsigma == -1) {
    return sf[1];
  }
  else if (nsigma == 1) {
    return sf[2];
  }
  // if nsigma ! = 0, -1, 1
  std::cout << "Warning: nsigma = " << nsigma << " is not -1, 0, 1" << std::endl;
  return 1;
}
