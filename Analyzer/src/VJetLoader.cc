#include "../include/VJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>

using namespace baconhep;

VJetLoader::VJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS,std::string iAddJetCHS,int iN, bool iData) { 
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

  isData = iData;  
  loadJECs(isData);  

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
void VJetLoader::selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum){
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

  fillJetCorr( fN,fLooseVJets,fVars,iRho,runNum);
  fillVJet(fN,fLooseVJets,fVars,iRho,runNum);
}
void VJetLoader::fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum){ 
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  
  TLorentzVector vPJet;
  
  
  for(int i0 = 0; i0 < lMin; i0++) {
    
    double JEC_old = (iObjects[i0]->pt)/(iObjects[i0]->ptRaw);
    vPJet.SetPtEtaPhiM(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, (iObjects[i0]->mass)/JEC_old);
    double jetE = vPJet.E();
    double JEC = JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
    					   iRho, iObjects[i0]->area, 
    					   runNum,
   					   JetCorrectionsIOV,JetCorrector);
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    double unc = getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default  
    iVals[i0*3+0] = iObjects[i0]->pt;
    iVals[i0*3+1] = iObjects[i0]->eta;
    iVals[i0*3+2] = iObjects[i0]->phi;
  }
}
void VJetLoader::selectVJetsByDoubleBCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum){
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
void VJetLoader::selectVJetsCHS(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum){
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
// Retrieve jet energy uncertainty as a function of pt and eta
double VJetLoader::getJecUnc( float pt, float eta , int run) {

  int foundIndex = -1;
  for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }

  jecUnc[foundIndex]->setJetPt(pt);
  jecUnc[foundIndex]->setJetEta(eta);
  return jecUnc[foundIndex]->getUncertainty(true);
}
void VJetLoader::fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum){ 
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
    
    std::cout << "AK8 pt = " << iObjects[i0]->pt << ", eta = " << iObjects[i0]->eta << ", rho = " << iRho << std::endl; 
    JME::JetParameters parameters = {{JME::Binning::JetPt, iObjects[i0]->pt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, iRho}};
    float sigma_MC = resolution.getResolution(parameters);
    float sf = resolution_sf.getScaleFactor(parameters);
    float sfUp = resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = resolution_sf.getScaleFactor(parameters, Variation::DOWN);

    std::cout << "AK8 sf = " << sf << ", " <<  sfUp << ", " << sfDown << std::endl;
    std::cout << "AK8 sigma_MC = " << sigma_MC << std::endl;
    
    double jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    double jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
    double jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;
    
    
    double JEC_old = (iObjects[i0]->pt)/(iObjects[i0]->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, (iObjects[i0]->mass)/JEC_old);
    double jetE = vPJet.E();
    
    double JEC = JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
    					   iRho, iObjects[i0]->area, 
    					   runNum,
   					   JetCorrectionsIOV,JetCorrector);
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    
    double unc_old = iObjects[i0]->unc;
    double unc = getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default
    
    std::cout << "runNum = " << runNum << std::endl;
    std::cout << "AK8 unc_old = " << JEC_old << std::endl;
    std::cout << "AK8 unc = " << JEC << std::endl;
    std::cout << "AK8 JEC_old = " << JEC_old << std::endl;
    std::cout << "AK8 JEC = " << JEC << std::endl;
    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
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

//2016 Prompt Reco
void VJetLoader::loadJECs(bool isData) {
    std::cout << "VJetLoader: loading jet energy correction constants" << std::endl;
    // initialize
    loadCMSSWPath();
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    resolution = JME::JetResolution(Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_PtResolution_AK8PFPuppi.txt",jecPathname.c_str()));
    resolution_sf = JME::JetResolutionScaleFactor(Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_SF_AK8PFPuppi.txt",jecPathname.c_str()));
    
    if (isData) {      
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));

      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_MC/Spring16_25nsV6_MC_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}
void VJetLoader::loadJECs_Rereco(bool isData) {
    std::cout << "VJetLoader: loading Rereco jet energy correction constants" << std::endl;
    // initialize
    loadCMSSWPath();
    std::string jecPathname = cmsswPath + "/src/BaconAnalyzer/Analyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    resolution = JME::JetResolution(Form("%s/Spring16_25nsV10_MC/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt",jecPathname.c_str()));
    resolution_sf = JME::JetResolutionScaleFactor(Form("%s/Spring16_25nsV10_MC/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",jecPathname.c_str()));
 
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016BCDV3_DATA/Spring16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016BCDV3_DATA/Spring16_23Sep2016BCDV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016BCDV3_DATA/Spring16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016BCDV3_DATA/Spring16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Spring16_23Sep2016BCDV3_DATA/Spring16_23Sep2016BCDV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016EFV3_DATA/Spring16_23Sep2016EFV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016EFV3_DATA/Spring16_23Sep2016EFV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016EFV3_DATA/Spring16_23Sep2016EFV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016EFV3_DATA/Spring16_23Sep2016EFV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Spring16_23Sep2016EFV3_DATA/Spring16_23Sep2016EFV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);

      correctionParameters.push_back(correctionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));

      //IOV: 2016G
      std::vector<JetCorrectorParameters> correctionParametersG = std::vector<JetCorrectorParameters> ();
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016GV3_DATA/Spring16_23Sep2016GV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016GV3_DATA/Spring16_23Sep2016GV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016GV3_DATA/Spring16_23Sep2016GV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016GV3_DATA/Spring16_23Sep2016GV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorG = new FactorizedJetCorrector(correctionParametersG);
      std::string jecUncPathG = jecPathname+"/Spring16_23Sep2016GV3_DATA/Spring16_23Sep2016GV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncG = new JetCorrectionUncertainty(jecUncPathG);

      correctionParameters.push_back(correctionParametersG);
      JetCorrector.push_back( JetCorrectorG );
      jecUnc.push_back(jecUncG);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 280385 ));

      //IOV: 2016H
      std::vector<JetCorrectorParameters> correctionParametersH = std::vector<JetCorrectorParameters> ();
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016HV3_DATA/Spring16_23Sep2016HV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016HV3_DATA/Spring16_23Sep2016HV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016HV3_DATA/Spring16_23Sep2016HV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016HV3_DATA/Spring16_23Sep2016HV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorH = new FactorizedJetCorrector(correctionParametersH);
      std::string jecUncPathH = jecPathname+"/Spring16_23Sep2016HV3_DATA/Spring16_23Sep2016HV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncH = new JetCorrectionUncertainty(jecUncPathH);

      correctionParameters.push_back(correctionParametersH);
      JetCorrector.push_back( JetCorrectorH );
      jecUnc.push_back(jecUncH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 280919, 99999999 ));

    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V3_MC/Spring16_23Sep2016V3_MC_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V3_MC/Spring16_23Sep2016V3_MC_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V3_MC/Spring16_23Sep2016V3_MC_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Spring16_23Sep2016V3_MC/Spring16_23Sep2016V3_MC_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
  
}

void VJetLoader::loadCMSSWPath() {
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in VJetLoader::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        cmsswPath = "";
    }
    cmsswPath = std::string(cmsswPathChar);
}


//Jet Energy Corrections
double VJetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 int run,
						 std::vector<std::pair<int,int> > JetCorrectionsIOV,
						 std::vector<FactorizedJetCorrector*> jetcorrector,
						 int jetCorrectionLevel,
						 bool printDebug) {

  int foundIndex = -1;
  for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }

  if (!jetcorrector[foundIndex]) {
    std::cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }

  jetcorrector[foundIndex]->setJetEta(jetEta);
  jetcorrector[foundIndex]->setJetPt(jetRawPt);
  jetcorrector[foundIndex]->setJetPhi(jetPhi);
  jetcorrector[foundIndex]->setJetE(jetE);
  jetcorrector[foundIndex]->setRho(rho);
  jetcorrector[foundIndex]->setJetA(jetArea);

  std::vector<float> corrections;
  corrections = jetcorrector[foundIndex]->getSubCorrections();

  if (printDebug) std::cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) std::cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) std::cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}

//Jet Energy Corrections
double VJetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 FactorizedJetCorrector *jetcorrector,
						 int jetCorrectionLevel,
						 bool printDebug) {
  if (!jetcorrector) {
    std::cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }

  jetcorrector->setJetEta(jetEta);
  jetcorrector->setJetPt(jetRawPt);
  jetcorrector->setJetPhi(jetPhi);
  jetcorrector->setJetE(jetE);
  jetcorrector->setRho(rho);
  jetcorrector->setJetA(jetArea);

  std::vector<float> corrections;
  corrections = jetcorrector->getSubCorrections();

  if (printDebug) std::cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) std::cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) std::cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}
