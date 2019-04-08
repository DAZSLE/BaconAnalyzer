#include "../include/VJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>

using namespace baconhep;

VJetLoader::VJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,int iN, bool iData, std::string iLabel) { 
  fVJets         = new TClonesArray("baconhep::TJet");
  fVAddJets      = new TClonesArray("baconhep::TAddJet");

  iTree->SetBranchAddress(iJet.c_str(),       &fVJets);
  iTree->SetBranchAddress(iAddJet.c_str(),    &fVAddJets);

  fVJetBr        = iTree->GetBranch(iJet.c_str());
  fVAddJetBr     = iTree->GetBranch(iAddJet.c_str());

  fN = iN;

  fJEC = new JECLoader(iData,iLabel,"AK8PFPuppi");
  r = new TRandom3(1988);
  isData=iData;
}
VJetLoader::~VJetLoader() { 
  delete fVJets;
  delete fVJetBr;
  delete fVAddJets;
  delete fVAddJetBr;
  delete fJEC;
}
void VJetLoader::reset() { 
  fNLooseVJets        = 0;
  fNTightVJets        = 0;  
  for(int i0 = 0; i0 < int(fisTightVJet.size()); i0++) fisTightVJet[i0] = -999;
  for(int i0 = 0; i0 < int(fisMatchedVJet.size()); i0++) fisMatchedVJet[i0] = -999;
  selectedVJets.clear();
  fLooseVJets.clear();
  x1List.clear();
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = -999;
  resetZprime();  
}
void VJetLoader::resetDoubleB() {
  fLooseVJetsByDoubleB.clear();
  selectedVJetsByDoubleB.clear();
}
void VJetLoader::resetZprime() {
  fisHadronicV.clear();
  fvMatching.clear();
  fvSize.clear();
  fpartonFlavor.clear();
  fhadronFlavor.clear();
  fnCharged.clear();
  fnNeutrals.clear();
  fnParticles.clear();
}
void VJetLoader::setupTree(TTree *iTree, std::string iJetLabel, bool iHWW) { 
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
  fLabels.push_back("e2_sdb1"); // Correlation function inputs beta=1 soft-dropped 
  fLabels.push_back("e3_sdb1");
  fLabels.push_back("e3_v1_sdb1");
  fLabels.push_back("e3_v2_sdb1");
  fLabels.push_back("e4_v1_sdb1");
  fLabels.push_back("e4_v2_sdb1");
  fLabels.push_back("N2sdb1"); // 2-prong ECFs observables
  fLabels.push_back("M2sdb1");
  fLabels.push_back("D2sdb1");
  fLabels.push_back("N2b1");
  fLabels.push_back("M2b1");
  fLabels.push_back("D2b1");
  fLabels.push_back("pt_old");
  fLabels.push_back("pt_JESUp");
  fLabels.push_back("pt_JESDown");
  fLabels.push_back("pt_JERUp");
  fLabels.push_back("pt_JERDown");
  fLabels.push_back("deepdoubleb");
  fLabels.push_back("deepdoublec");
  fLabels.push_back("deepdoublecvb");
  fLabels.push_back("deepdoubleb_nomasssculptpen");
  fLabels.push_back("deepdoublec_nomasssculptpen");
  fLabels.push_back("deepdoublecvb_nomasssculptpen");

  if(iHWW) {
    fLabels.push_back("lepCPt");
    fLabels.push_back("lepCEta");
    fLabels.push_back("lepCPhi");
    fLabels.push_back("lepCId");
    fLabels.push_back("lsfCInc");
    fLabels.push_back("lsfC_2");
    fLabels.push_back("lsfC_3");
    fLabels.push_back("lsfC_4");
    fLabels.push_back("lmdCInc");
    fLabels.push_back("lmdC_2");
    fLabels.push_back("lmdC_3");
    fLabels.push_back("lmdC_4");
    fLabels.push_back("lsfC_3_sj1_pt");
    fLabels.push_back("lsfC_3_sj1_eta");
    fLabels.push_back("lsfC_3_sj1_phi");
    fLabels.push_back("lsfC_3_sj1_m");
    fLabels.push_back("lsfC_3_sj2_pt");
    fLabels.push_back("lsfC_3_sj2_eta");
    fLabels.push_back("lsfC_3_sj2_phi");
    fLabels.push_back("lsfC_3_sj2_m");
    fLabels.push_back("lsfC_3_sj3_pt");
    fLabels.push_back("lsfC_3_sj3_eta");
    fLabels.push_back("lsfC_3_sj3_phi");
    fLabels.push_back("lsfC_3_sj3_m");
    fLabels.push_back("sj1_pt");
    fLabels.push_back("sj1_eta");
    fLabels.push_back("sj1_phi");
    fLabels.push_back("sj1_m");
    fLabels.push_back("sj2_pt");
    fLabels.push_back("sj2_eta");
    fLabels.push_back("sj2_phi");
    fLabels.push_back("sj2_m");
    fLabels.push_back("sj3_pt");
    fLabels.push_back("sj3_eta");
    fLabels.push_back("sj3_phi");
    fLabels.push_back("sj3_m");
    fLabels.push_back("sj4_pt");
    fLabels.push_back("sj4_eta");
    fLabels.push_back("sj4_phi");
    fLabels.push_back("sj4_m");
    fLabels.push_back("tau31");
    fLabels.push_back("tau42");
  }

  /*
  if(iGen) {
    fLabels.push_back("genmsd");
    fLabels.push_back("geneta");
    fLabels.push_back("genphi");
    }*/

  fTree = iTree;
  for(int i0 = 0; i0 < fN*3.;                    i0++) {double pVar = 0; fVars.push_back(pVar);} // declare array of vars
  for(int i0 = 0; i0 < fN*(int(fLabels.size())); i0++) {double pVar = 0; fVars.push_back(pVar);} 

  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars);                                                 // from Utils.cc => fN =1 *_pt,*_eta,*_phi for vjet0 (3*1=3)
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars,fN*3,fLabels);

  std::stringstream pSNJ;   pSNJ << "n" << iJetLabel << "s";
  fTree->Branch(pSNJ.str().c_str() ,&fNLooseVJets         ,(pSNJ.str()+"/I").c_str());  
  for(int i0 = 0; i0 < fN;                       i0++) {
    fisTightVJet.push_back(-999);
    fisMatchedVJet.push_back(-999);
  }
  for(int i0 = 0; i0 < fN;                    i0++) {
    std::stringstream pSTJ;   pSTJ << iJetLabel << i0 << "_isTightVJet";
    std::stringstream pSMJ;   pSMJ << iJetLabel << i0 << "_isMatchedVJet";
    fTree->Branch(pSTJ.str().c_str() ,&fisTightVJet[i0]         ,(pSTJ.str()+"/I").c_str());
    fTree->Branch(pSMJ.str().c_str() ,&fisMatchedVJet[i0]       ,(pSMJ.str()+"/I").c_str());
  }
}
void VJetLoader::setupTreeZprime(TTree *iTree, std::string iJetLabel) {
  resetZprime();
  for(int i0 = 0; i0 < fN; i0++) {
    double pVar = -1;
    fisHadronicV.push_back(pVar);
    fvMatching.push_back(pVar);
    fvSize.push_back(pVar);
    fpartonFlavor.push_back(pVar);
    fhadronFlavor.push_back(pVar);
    fnCharged.push_back(pVar);
    fnNeutrals.push_back(pVar);
    fnParticles.push_back(pVar);
  }
  for(int i0 = 0; i0 < fN; i0++) {
    std::stringstream pSiV;   pSiV << iJetLabel << i0 << "_isHadronicV";
    std::stringstream pSVM;   pSVM << iJetLabel << i0 << "_vMatching";
    std::stringstream pSVS;   pSVS << iJetLabel << i0 << "_vSize";
    std::stringstream pSpF;   pSpF << iJetLabel << i0 << "_partonFlavor";
    std::stringstream pShF;   pShF << iJetLabel << i0 << "_hadronFlavor";
    std::stringstream pSnC;   pSnC << iJetLabel << i0 << "_nCharged";
    std::stringstream pSnN;   pSnN << iJetLabel << i0 << "_nNeutrals";
    std::stringstream pSnP;   pSnP << iJetLabel << i0 << "_nParticles";

    fTree->Branch(pSiV.str().c_str() ,&fisHadronicV[i0]     ,(pSiV.str()+"/I").c_str());
    fTree->Branch(pSVM.str().c_str() ,&fvMatching[i0]       ,(pSVM.str()+"/D").c_str());
    fTree->Branch(pSVS.str().c_str() ,&fvSize[i0]           ,(pSVS.str()+"/D").c_str());
    fTree->Branch(pSpF.str().c_str() ,&fpartonFlavor[i0]    ,(pSpF.str()+"/I").c_str());
    fTree->Branch(pShF.str().c_str() ,&fhadronFlavor[i0]    ,(pShF.str()+"/I").c_str());
    fTree->Branch(pSnC.str().c_str() ,&fnCharged[i0]        ,(pSnC.str()+"/I").c_str());
    fTree->Branch(pSnN.str().c_str() ,&fnNeutrals[i0]       ,(pSnN.str()+"/I").c_str());
    fTree->Branch(pSnP.str().c_str() ,&fnParticles[i0]      ,(pSnP.str()+"/I").c_str());
  }
}
void VJetLoader::load(int iEvent) { 
  fVJets       ->Clear();
  fVJetBr      ->GetEntry(iEvent);
  fVAddJets    ->Clear();
  fVAddJetBr   ->GetEntry(iEvent);
}
void VJetLoader::selectVJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum,bool iHWW){
  reset();  
  int lCount(0), lCountT(0);
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    TJet *pVJet = (TJet*)((*fVJets)[i0]);    
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    
    double JEC_old = (pVJet->pt)/(pVJet->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(pVJet->ptRaw, pVJet->eta, pVJet->phi, (pVJet->mass)/JEC_old);
    double jetE = vPJet.E();
    double JEC = fJEC->JetEnergyCorrectionFactor(pVJet->ptRaw, pVJet->eta, pVJet->phi, jetE, 
						 iRho, pVJet->area, 
						 runNum);
    double jetCorrPt = JEC*(pVJet->ptRaw);
    double jetCorrE = JEC*(vPJet.E());
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, pVJet->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; 
    float sigma_MC = fJEC->resolution.getResolution(parameters);
    float sf = fJEC->resolution_sf.getScaleFactor(parameters);
    float sfUp = fJEC->resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = fJEC->resolution_sf.getScaleFactor(parameters, Variation::DOWN);
    double x1 = r->Gaus();
    double jetEnergySmearFactor = 1.0; 
    double jetEnergySmearFactorUp = 1.0; 
    double jetEnergySmearFactorDown = 1.0;    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
      jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x1;
      jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x1;
    }    
    double unc = fJEC->getJecUnc( jetCorrPt, pVJet->eta, runNum ); //use run=999 as default
    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;
    
    double jetCorrESmear = jetCorrE*jetEnergySmearFactor*(1+unc);    
    double jetEJESUp = jetCorrE*jetEnergySmearFactor*(1+unc);
    double jetEJESDown = jetCorrE*jetEnergySmearFactor/(1+unc);
    double jetEJERUp = jetCorrE*jetEnergySmearFactorUp;
    double jetEJERDown = jetCorrE*jetEnergySmearFactorDown;
    
    TLorentzVector thisJet;  thisJet.SetPtEtaPhiE(jetCorrPtSmear, pVJet->eta, pVJet->phi, jetCorrESmear);
    TLorentzVector thisJetJESUp;  thisJetJESUp.SetPtEtaPhiE(jetPtJESUp, pVJet->eta, pVJet->phi, jetEJESUp);
    TLorentzVector thisJetJESDown; thisJetJESDown.SetPtEtaPhiE(jetPtJESDown, pVJet->eta, pVJet->phi, jetEJESDown);
    TLorentzVector thisJetJERUp;  thisJetJERUp.SetPtEtaPhiE(jetPtJERUp,  pVJet->eta, pVJet->phi, jetEJERUp);
    TLorentzVector thisJetJERDown; thisJetJERDown.SetPtEtaPhiE(jetPtJERDown,  pVJet->eta, pVJet->phi, jetEJERDown);
    
    if(jetCorrPtSmear   <=  200)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(!passJetTightSel(pVJet))                                            continue;

    addJet(pVJet,fLooseVJets);
    lCount++;
    x1List.push_back(x1);
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
  fillVJet(fN,fLooseVJets,fVars,iRho,runNum,iHWW);
}
void VJetLoader::fillJetCorr(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum){ 
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  
  for(int i0 = 0; i0 < lMin; i0++) {
    
    double JEC_old = (iObjects[i0]->pt)/(iObjects[i0]->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, (iObjects[i0]->mass)/JEC_old);
    double jetE = vPJet.E();
    double JEC = fJEC->JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
						 iRho, iObjects[i0]->area, 
						 runNum);
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    double x1 = x1List[i0];
    double jetEnergySmearFactor = 1.0;    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; 
    float sigma_MC = fJEC->resolution.getResolution(parameters);
    float sf = fJEC->resolution_sf.getScaleFactor(parameters);    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    }    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    iVals[i0*3+0] = jetCorrPtSmear;
    iVals[i0*3+1] = iObjects[i0]->eta;
    iVals[i0*3+2] = iObjects[i0]->phi;
  }
}
void VJetLoader::selectVJetsByDoubleB(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, double dR, double iRho, unsigned int runNum){
  // first do Puppi jets (pT > 500 GeV)
  resetDoubleB(); 
  for  (int i0 = 0; i0 < fVJets->GetEntriesFast(); i0++) { 
    TJet *pVJet = (TJet*)((*fVJets)[i0]);
    if(pVJet->pt        <=  500)                                           continue;
    if(fabs(pVJet->eta) >=  2.5)                                           continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iElectrons))                      continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iMuons))                          continue;
    if(passVeto(pVJet->eta,pVJet->phi,dR,iPhotons))                        continue;
    if(!passJetTightSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJetsByDoubleB);
    if(!passJetTightLepVetoSel(pVJet))                                     continue;
  }
  addVJet(fLooseVJetsByDoubleB,selectedVJetsByDoubleB);
}
void VJetLoader::fillVJet(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, double iRho, unsigned int runNum, bool iHWW){ 
  int lBase = 3.*fN;
  int lMin = iObjects.size();
  int lNLabel = int(fLabels.size());
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    TAddJet *pAddJet = getAddJet(iObjects[i0]);
    //JEC    
    double x1 = x1List[i0];
    double JEC_old = (iObjects[i0]->pt)/(iObjects[i0]->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, (iObjects[i0]->mass)/JEC_old);
    double jetE = vPJet.E();
    
    double JEC = fJEC->JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
    					   iRho, iObjects[i0]->area, 
    					   runNum);
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    
    double unc = fJEC->getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default
    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; 
    float sigma_MC = fJEC->resolution.getResolution(parameters);
    float sf = fJEC->resolution_sf.getScaleFactor(parameters);
    float sfUp = fJEC->resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = fJEC->resolution_sf.getScaleFactor(parameters, Variation::DOWN);

    double jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    double jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x1;
    double jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x1;
    //std::cout << "JER SF " << sf << " up " << sfUp <<  " dn " << sfDown << " sigma_MC " << sigma_MC << " x1 " << x1 << std::endl;
    //std::cout << "smearup " << jetEnergySmearFactorUp << " dn "<< jetEnergySmearFactorDown << std::endl;

    //double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;

    /*if (true){
      std::cout << "VJet" << std::endl;
      std::cout << "i0 =" << i0 << std::endl;
      std::cout << "ptraw = " << iObjects[i0]->ptRaw << std::endl;
      std::cout << "eta = " << iObjects[i0]->eta << std::endl;
      std::cout << "rho = " << iRho << std::endl;
      std::cout << "runNum = " << runNum << std::endl;
      std::cout << "sf = " << sf << ", " <<  sfUp << ", " << sfDown << std::endl;
      std::cout << "sigma_MC = " << sigma_MC << std::endl;    
      std::cout << "runNum = " << runNum << std::endl;
      std::cout << "unc = " << unc << std::endl;
      std::cout << "JEC_old = " << JEC_old << std::endl;
      std::cout << "JEC = " << JEC << std::endl;
      std::cout << "x1 = " << x1 << std::endl;
      std::cout << "ptcorr = " << jetCorrPt << std::endl;
      std::cout << "jesup = " << jetPtJESUp << std::endl;
      std::cout << "jesdown = " << jetPtJESDown << std::endl;
      std::cout << "jerup = " << jetPtJERUp << std::endl;
      std::cout << "jerdown = " << jetPtJERDown << std::endl;
      }*/

    iVals[lBase+i0*lNLabel+0]  = JEC*jetEnergySmearFactor*(iObjects[i0]->mass);
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
    iVals[lBase+i0*lNLabel+21] = pAddJet->e2_sdb1;
    iVals[lBase+i0*lNLabel+22] = pAddJet->e3_sdb1;
    iVals[lBase+i0*lNLabel+23] = pAddJet->e3_v1_sdb1;
    iVals[lBase+i0*lNLabel+24] = pAddJet->e3_v2_sdb1;
    iVals[lBase+i0*lNLabel+25] = pAddJet->e4_v1_sdb1;
    iVals[lBase+i0*lNLabel+26] = pAddJet->e4_v2_sdb1;
    iVals[lBase+i0*lNLabel+27] = pAddJet->e3_v2_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+28] = pAddJet->e3_v1_sdb1/(pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+29] = pAddJet->e3_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    iVals[lBase+i0*lNLabel+30] = pAddJet->e3_v2_b1/(pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+31] = pAddJet->e3_v1_b1/(pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+32] = pAddJet->e3_b1/(pAddJet->e2_b1*pAddJet->e2_b1*pAddJet->e2_b1);
    iVals[lBase+i0*lNLabel+33] = iObjects[i0]->pt;
    iVals[lBase+i0*lNLabel+34] = jetPtJESUp;
    iVals[lBase+i0*lNLabel+35] = jetPtJESDown;
    iVals[lBase+i0*lNLabel+36] = jetPtJERUp;
    iVals[lBase+i0*lNLabel+37] = jetPtJERDown;
    iVals[lBase+i0*lNLabel+38] = pAddJet->deepdoubleb;
    iVals[lBase+i0*lNLabel+39] = pAddJet->deepdoublec;
    iVals[lBase+i0*lNLabel+40] = pAddJet->deepdoublecvb;
    iVals[lBase+i0*lNLabel+41] = pAddJet->deepdoubleb_nomasssculptpen;
    iVals[lBase+i0*lNLabel+42] = pAddJet->deepdoublec_nomasssculptpen;
    iVals[lBase+i0*lNLabel+43] = pAddJet->deepdoublecvb_nomasssculptpen;

    int iN = 43;

    if (iHWW) {
      iVals[lBase+i0*lNLabel+iN+1] = pAddJet->lepCPt;
      iVals[lBase+i0*lNLabel+iN+2] = pAddJet->lepCEta;
      iVals[lBase+i0*lNLabel+iN+3] = pAddJet->lepCPhi;
      iVals[lBase+i0*lNLabel+iN+4] = pAddJet->lepCId;
      iVals[lBase+i0*lNLabel+iN+5] = pAddJet->lsfCInc;
      iVals[lBase+i0*lNLabel+iN+6] = pAddJet->lsfC_2;
      iVals[lBase+i0*lNLabel+iN+7] = pAddJet->lsfC_3;
      iVals[lBase+i0*lNLabel+iN+8] = pAddJet->lsfC_4;
      iVals[lBase+i0*lNLabel+iN+9] = pAddJet->lmdCInc;
      iVals[lBase+i0*lNLabel+iN+10] = pAddJet->lmdC_2;
      iVals[lBase+i0*lNLabel+iN+11] = pAddJet->lmdC_3;
      iVals[lBase+i0*lNLabel+iN+12] = pAddJet->lmdC_4;
      iVals[lBase+i0*lNLabel+iN+13] = pAddJet->lsfC_3_sj1_pt;
      iVals[lBase+i0*lNLabel+iN+14] = pAddJet->lsfC_3_sj1_eta;
      iVals[lBase+i0*lNLabel+iN+15] = pAddJet->lsfC_3_sj1_phi;
      iVals[lBase+i0*lNLabel+iN+16] = pAddJet->lsfC_3_sj1_m;
      iVals[lBase+i0*lNLabel+iN+17] = pAddJet->lsfC_3_sj2_pt;
      iVals[lBase+i0*lNLabel+iN+18] = pAddJet->lsfC_3_sj2_eta;
      iVals[lBase+i0*lNLabel+iN+19] = pAddJet->lsfC_3_sj2_phi;
      iVals[lBase+i0*lNLabel+iN+20] = pAddJet->lsfC_3_sj2_m;
      iVals[lBase+i0*lNLabel+iN+21] = pAddJet->lsfC_3_sj3_pt;
      iVals[lBase+i0*lNLabel+iN+22] = pAddJet->lsfC_3_sj3_eta;
      iVals[lBase+i0*lNLabel+iN+23] = pAddJet->lsfC_3_sj3_phi;
      iVals[lBase+i0*lNLabel+iN+24] = pAddJet->lsfC_3_sj3_m;
      iVals[lBase+i0*lNLabel+iN+25] = pAddJet->sj1_pt;
      iVals[lBase+i0*lNLabel+iN+26] = pAddJet->sj1_eta;
      iVals[lBase+i0*lNLabel+iN+27] = pAddJet->sj1_phi;
      iVals[lBase+i0*lNLabel+iN+28] = pAddJet->sj1_m;
      iVals[lBase+i0*lNLabel+iN+29] = pAddJet->sj2_pt;
      iVals[lBase+i0*lNLabel+iN+30] = pAddJet->sj2_eta;
      iVals[lBase+i0*lNLabel+iN+31] = pAddJet->sj2_phi;
      iVals[lBase+i0*lNLabel+iN+32] = pAddJet->sj2_m;
      iVals[lBase+i0*lNLabel+iN+33] = pAddJet->sj3_pt;
      iVals[lBase+i0*lNLabel+iN+34] = pAddJet->sj3_eta;
      iVals[lBase+i0*lNLabel+iN+35] = pAddJet->sj3_phi;
      iVals[lBase+i0*lNLabel+iN+36] = pAddJet->sj3_m;
      iVals[lBase+i0*lNLabel+iN+37] = pAddJet->sj4_pt;
      iVals[lBase+i0*lNLabel+iN+38] = pAddJet->sj4_eta;
      iVals[lBase+i0*lNLabel+iN+39] = pAddJet->sj4_phi;
      iVals[lBase+i0*lNLabel+iN+40] = pAddJet->sj4_m;
      iVals[lBase+i0*lNLabel+iN+41]  = (pAddJet->tau3/pAddJet->tau1);
      iVals[lBase+i0*lNLabel+iN+42]  = (pAddJet->tau4/pAddJet->tau2);
      iN = iN+42;
    }

    /*
    bool iGen = true;
    if (iGen) {
      iVals[lBase+i0*lNLabel+iN+1] = iObjects[i0]->genmsd;
      iVals[lBase+i0*lNLabel+iN+2] = iObjects[i0]->geneta;
      iVals[lBase+i0*lNLabel+iN+3] = iObjects[i0]->genphi;
      iN = iN+3;
      }*/

    fpartonFlavor[i0]   = iObjects[i0]->partonFlavor;
    fhadronFlavor[i0]   = iObjects[i0]->hadronFlavor;
    fnCharged[i0]       = iObjects[i0]->nCharged;
    fnNeutrals[i0]      = iObjects[i0]->nNeutrals;
    fnParticles[i0]     = iObjects[i0]->nParticles;

  }
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
