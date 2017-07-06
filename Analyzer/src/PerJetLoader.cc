#include "../include/PerJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>
#include <unordered_set>

#define NCPF 50
#define NNPF 50
#define NSV 5

using namespace baconhep;

PerJetLoader::PerJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS,std::string iAddJetCHS,int iN, bool iData) { 
  fVJets         = new TClonesArray("baconhep::TJet");
  fVAddJets      = new TClonesArray("baconhep::TAddJet");
  fGens         = new TClonesArray("baconhep::TGenParticle");
  fPFs = new TClonesArray("baconhep::TPFPart");
  fSVs = new TClonesArray("baconhep::TSVtx");

  iTree->SetBranchAddress(iJet.c_str(),       &fVJets);
  iTree->SetBranchAddress(iAddJet.c_str(),    &fVAddJets);
  iTree->SetBranchAddress("GenParticle", &fGens);
  iTree->SetBranchAddress("PFPart", &fPFs);
  iTree->SetBranchAddress("SV", &fSVs);

  fVJetBr        = iTree->GetBranch(iJet.c_str());
  fVAddJetBr     = iTree->GetBranch(iAddJet.c_str());
  fGenBr         = iTree->GetBranch("GenParticle");
  fPFBr = iTree->GetBranch("PFPart");
  fSVBr = iTree->GetBranch("SV");

  fN = iN;

  isData = iData;  
  loadJECs_Rereco(isData);

  r = new TRandom3(1993);

  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
}

PerJetLoader::~PerJetLoader() { 
  delete fVJets;
  delete fVJetBr;
  delete fVAddJets;
  delete fVAddJetBr;
  delete fGens;
  delete fGenBr;
  delete fPFs;
  delete fPFBr;
  delete fSVs;
  delete fSVBr;
  for (auto &iter : fCPFArrs) 
    delete iter.second;
  for (auto &iter : fNPFArrs) 
    delete iter.second;
  for (auto &iter : fSVArrs) 
    delete iter.second;
}

void PerJetLoader::reset() { 
  fNLooseVJets        = 0;
  fNTightVJets        = 0;  
  for(int i0 = 0; i0 < int(fisTightVJet.size()); i0++) fisTightVJet[i0] = -999;
  selectedVJets.clear();
  fLooseVJets.clear();
  x1List.clear();
  x2List.clear();
  x3List.clear();    
  for (auto &iter : fSingletons) {
    fSingletons[iter.first] = 0;
  }
  for (auto &iter : fCPFArrs) {
    for (unsigned i = 0; i != NCPF; ++i) 
      fCPFArrs[iter.first][i] = 0;
  }
  for (auto &iter : fNPFArrs) {
    for (unsigned i = 0; i != NNPF; ++i) 
      fNPFArrs[iter.first][i] = 0;
  }
  for (auto &iter : fSVArrs) {
    for (unsigned i = 0; i != NSV; ++i) 
      fSVArrs[iter.first][i] = 0;
  }
  resetZprime();  
}

void PerJetLoader::resetZprime() {
  fvSize              = -999;
  fvMatching          = -999;
  fisHadronicV        = 0;
}

void PerJetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();

  fSingletons.clear();
  fCPFArrs.clear();
  fNPFArrs.clear();
  fSVArrs.clear();

  fSingletons["pt"] = 0;
  fSingletons["eta"] = 0;
  fSingletons["phi"] = 0;
  fSingletons["mass"] = 0;
  fSingletons["csv"] = 0;
  fSingletons["CHF"] = 0;
  fSingletons["NHF"] = 0;
  fSingletons["NEMF"] = 0;
  fSingletons["tau21"] = 0;
  fSingletons["tau32"] = 0;
  fSingletons["msd"] = 0;
  fSingletons["rho"] = 0;
  fSingletons["minsubcsv"] = 0;
  fSingletons["maxsubcsv"] = 0;
  fSingletons["doublecsv"] = 0;
  fSingletons["doublesub"] = 0;
  fSingletons["ptraw"] = 0;
  fSingletons["genpt"] = 0;
  fSingletons["e2_b1"] = 0; // Correlation function inputs beta=1
  fSingletons["e3_b1"] = 0;
  fSingletons["e3_v1_b1"] = 0;
  fSingletons["e3_v2_b1"] = 0;
  fSingletons["e4_v1_b1"] = 0;
  fSingletons["e4_v2_b1"] = 0;
  fSingletons["e2_b2"] = 0; // Correlation function inputs beta=2
  fSingletons["e3_b2"] = 0;
  fSingletons["e3_v1_b2"] = 0;
  fSingletons["e3_v2_b2"] = 0;
  fSingletons["e4_v1_b2"] = 0;
  fSingletons["e4_v2_b2"] = 0;
  fSingletons["e2_sdb1"] = 0; // Correlation function inputs beta=1 soft-dropped 
  fSingletons["e3_sdb1"] = 0;
  fSingletons["e3_v1_sdb1"] = 0;
  fSingletons["e3_v2_sdb1"] = 0;
  fSingletons["e4_v1_sdb1"] = 0;
  fSingletons["e4_v2_sdb1"] = 0;
  fSingletons["e2_sdb2"] = 0; // Correlation function inputs beta=2 soft-dropped 
  fSingletons["e3_sdb2"] = 0;
  fSingletons["e3_v1_sdb2"] = 0;
  fSingletons["e3_v2_sdb2"] = 0;
  fSingletons["e4_v1_sdb2"] = 0;
  fSingletons["e4_v2_sdb2"] = 0;
  fSingletons["N2sdb1"] = 0; // 2-prong ECFs observables
  fSingletons["N2sdb2"] = 0;
  fSingletons["M2sdb1"] = 0;
  fSingletons["M2sdb2"] = 0;
  fSingletons["D2sdb1"] = 0;
  fSingletons["D2sdb2"] = 0;
  fSingletons["N2b1"] = 0;
  fSingletons["N2b2"] = 0;
  fSingletons["M2b1"] = 0;
  fSingletons["M2b2"] = 0;
  fSingletons["D2b1"] = 0;
  fSingletons["D2b2"] = 0;
  fSingletons["pt_old"] = 0;
  fSingletons["pt_JESUp"] = 0;
  fSingletons["pt_JESDown"] = 0;
  fSingletons["pt_JERUp"] = 0;
  fSingletons["pt_JERDown"] = 0;
  fSingletons["e2_sdb05"] = 0; // Correlation function inputs beta=0.5 soft-dropped 
  fSingletons["e3_sdb05"] = 0;
  fSingletons["e3_v1_sdb05"] = 0;
  fSingletons["e3_v2_sdb05"] = 0;
  fSingletons["e4_v1_sdb05"] = 0;
  fSingletons["e4_v2_sdb05"] = 0;
  fSingletons["e2_sdb4"] = 0; // Correlation function inputs beta=4 soft-dropped 
  fSingletons["e3_sdb4"] = 0;
  fSingletons["e3_v1_sdb4"] = 0;
  fSingletons["e3_v2_sdb4"] = 0;
  fSingletons["e4_v1_sdb4"] = 0;
  fSingletons["e4_v2_sdb4"] = 0;
  fSingletons["flavour"] = 0;
  fSingletons["nbHadrons"] = 0;
  fSingletons["nSV"] = 0;
  fSingletons["jetNTracks"] = 0;
  fSingletons["tau_flightDistance2dSig_1"] = 0;
  fSingletons["SubJet_csv"] = 0;
  fSingletons["z_ratio"] = 0;
  fSingletons["trackSipdSig_3"] = 0;
  fSingletons["trackSipdSig_2"] = 0;
  fSingletons["trackSipdSig_1"] = 0;
  fSingletons["trackSipdSig_0"] = 0;
  fSingletons["trackSipdSig_1_0"] = 0;
  fSingletons["trackSipdSig_0_0"] = 0;
  fSingletons["trackSipdSig_1_1"] = 0;
  fSingletons["trackSipdSig_0_1"] = 0;
  fSingletons["trackSip2dSigAboveCharm_0"] = 0;
  fSingletons["trackSip2dSigAboveBottom_0"] = 0;
  fSingletons["trackSip2dSigAboveBottom_1"] = 0;
  fSingletons["tau1_trackEtaRel_0"] = 0;
  fSingletons["tau1_trackEtaRel_1"] = 0;
  fSingletons["tau1_trackEtaRel_2"] = 0;
  fSingletons["tau0_trackEtaRel_0"] = 0;
  fSingletons["tau0_trackEtaRel_1"] = 0;
  fSingletons["tau0_trackEtaRel_2"] = 0;
  fSingletons["tau_vertexMass_0"] = 0;
  fSingletons["tau_vertexEnergyRatio_0"] = 0;
  fSingletons["tau_vertexDeltaR_0"] = 0;
  fSingletons["tau_flightDistance2dSig_0"] = 0;
  fSingletons["tau_vertexMass_1"] = 0;
  fSingletons["tau_vertexEnergyRatio_1"] = 0;
  fSingletons["nProngs"] = 0;

  fCPFArrs["cpfPt"] = new float[NCPF]; // example - add more

  fSVArrs["svPt"] = new float[NSV];

  fTree = iTree;

  for (auto &iter : fSingletons) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    fTree->Branch(bname.str().c_str(), &(iter.second), (bname.str()+"/f").c_str());
  }
  for (auto &iter : fCPFArrs) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname.str() << "[" << NCPF << "]/f";
    fTree->Branch(bname.str().c_str(), &(iter.second), bname2.str().c_str());
  }
  for (auto &iter : fNPFArrs) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname.str() << "[" << NNPF << "]/f";
    fTree->Branch(bname.str().c_str(), &(iter.second), bname2.str().c_str());
  }
  for (auto &iter : fSVArrs) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname.str() << "[" << NSV << "]/f";
    fTree->Branch(bname.str().c_str(), &(iter.second), bname2.str().c_str());
  }

}

void PerJetLoader::setupTreeZprime(TTree *iTree, std::string iJetLabel) {
  resetZprime();
  std::stringstream pSiV;   pSiV << iJetLabel << "0_isHadronicV";
  std::stringstream pSVM;   pSVM << iJetLabel << "0_vMatching";
  std::stringstream pSVS;   pSVS << iJetLabel << "0_vSize";
  std::stringstream pSpF;   pSpF << iJetLabel << "0_partonFlavor";
  std::stringstream pShF;   pShF << iJetLabel << "0_hadronFlavor";
  std::stringstream pSnC;   pSnC << iJetLabel << "0_nCharged";
  std::stringstream pSnN;   pSnN << iJetLabel << "0_nNeutrals";
  std::stringstream pSnP;   pSnP << iJetLabel << "0_nParticles";
  std::stringstream pSvF;   pSvF << iJetLabel << "0_vertexFlavor"; 
  std::stringstream pSvFI;   pSvFI << iJetLabel << "0_vertexFlavorInfo";

  fTree = iTree;
  fTree->Branch(pSiV.str().c_str() ,&fisHadronicV         ,(pSiV.str()+"/I").c_str());
  fTree->Branch(pSVM.str().c_str() ,&fvMatching           ,(pSVM.str()+"/D").c_str());
  fTree->Branch(pSVS.str().c_str() ,&fvSize               ,(pSVS.str()+"/D").c_str());
  fTree->Branch(pSpF.str().c_str() ,&fpartonFlavor        ,(pSpF.str()+"/I").c_str());
  fTree->Branch(pShF.str().c_str() ,&fhadronFlavor        ,(pShF.str()+"/I").c_str());
  fTree->Branch(pSnC.str().c_str() ,&fnCharged            ,(pSnC.str()+"/I").c_str());
  fTree->Branch(pSnN.str().c_str() ,&fnNeutrals           ,(pSnN.str()+"/I").c_str());
  fTree->Branch(pSnP.str().c_str() ,&fnParticles          ,(pSnP.str()+"/I").c_str());
  fTree->Branch(pSvF.str().c_str() ,&fnVtxFlavor          ,(pSvF.str()+"/I").c_str());
  fTree->Branch(pSvFI.str().c_str() ,&fnVtxFlavInfo, (pSvFI.str()+"/I").c_str());
}

void PerJetLoader::load(int iEvent) { 
  fVJets       ->Clear();
  fVJetBr      ->GetEntry(iEvent);
  fVAddJets    ->Clear();
  fVAddJetBr   ->GetEntry(iEvent);
  fGens        ->Clear();
  fGenBr       ->GetEntry(iEvent);
  fPFs        ->Clear();
  fPFBr       ->GetEntry(iEvent);
  fSVs        ->Clear();
  fSVBr       ->GetEntry(iEvent);
}

void PerJetLoader::selectVJets(std::vector<TLorentzVector> &iElectrons, 
                               std::vector<TLorentzVector> &iMuons, 
                               std::vector<TLorentzVector> &iPhotons, 
                               double dR, 
                               double iRho, 
                               unsigned int runNum)
{
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
    double JEC = JetEnergyCorrectionFactor(pVJet->ptRaw, pVJet->eta, pVJet->phi, jetE, 
                 iRho, pVJet->area, 
                 runNum,
                JetCorrectionsIOV,JetCorrector);
    double jetCorrPt = JEC*(pVJet->ptRaw);
    double jetCorrE = JEC*(vPJet.E());
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, pVJet->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; // max 44.30 for Spring16_25nsV6_MC JER (CHANGE ONCE UPDATED)
    float sigma_MC = resolution.getResolution(parameters);
    float sf = resolution_sf.getScaleFactor(parameters);
    float sfUp = resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = resolution_sf.getScaleFactor(parameters, Variation::DOWN);
    double x1 = r->Gaus();
    double x2 = r->Gaus();
    double x3 = r->Gaus();
    double jetEnergySmearFactor = 1.0; 
    double jetEnergySmearFactorUp = 1.0; 
    double jetEnergySmearFactorDown = 1.0;    
    if (!isData) {      
      jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
      jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
      jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;
    }    
    double unc = getJecUnc( jetCorrPt, pVJet->eta, runNum ); //use run=999 as default
    
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
    if(!passJetLooseSel(pVJet))                                            continue;
    addJet(pVJet,fLooseVJets);
    lCount++;
    x1List.push_back(x1);
    x2List.push_back(x2);
    x3List.push_back(x3);
    if(!passJetTightLepVetoSel(pVJet))                                     continue;
    lCountT++;
  }
  addVJet(fLooseVJets,selectedVJets);

  for  (int i0 = 0; i0 < int(selectedVJets.size()); i0++) { 
    if(passJetTightLepVetoSel(fLooseVJets[i0])) fisTightVJet[i0] = 1;
  }
  fNLooseVJets = lCount;
  fNTightVJets = lCountT;

  fillVJet(fN,fLooseVJets,dR,iRho,runNum);
}


void PerJetLoader::fillVJet(int iN,
                            std::vector<TJet*> &iObjects,
                            double dR, 
                            double iRho, 
                            unsigned int runNum)
{ 
  int lBase = 3.*fN;
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    
    //JEC    
    double x1 = x1List[i0];
    double x2 = x2List[i0];
    double x3 = x3List[i0];
    
    double JEC_old = (iObjects[i0]->pt)/(iObjects[i0]->ptRaw);
    TLorentzVector vPJet;
    vPJet.SetPtEtaPhiM(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, (iObjects[i0]->mass)/JEC_old);
    double jetE = vPJet.E();
    double JEC = JetEnergyCorrectionFactor(iObjects[i0]->ptRaw, iObjects[i0]->eta, iObjects[i0]->phi, jetE, 
                 iRho, iObjects[i0]->area, 
                 runNum,
                JetCorrectionsIOV,JetCorrector);    
    double jetCorrPt = JEC*(iObjects[i0]->ptRaw);
    double unc = getJecUnc( jetCorrPt, iObjects[i0]->eta, runNum ); //use run=999 as default    
    JME::JetParameters parameters = {{JME::Binning::JetPt, jetCorrPt}, {JME::Binning::JetEta, iObjects[i0]->eta}, {JME::Binning::Rho, TMath::Min(iRho,44.30)}}; // max 44.30 for Spring16_25nsV6_MC JER (CHANGE ONCE UPDATED)
    float sigma_MC = resolution.getResolution(parameters);
    float sf = resolution_sf.getScaleFactor(parameters);    
    float sfUp = resolution_sf.getScaleFactor(parameters, Variation::UP);
    float sfDown = resolution_sf.getScaleFactor(parameters, Variation::DOWN);

    double jetEnergySmearFactor = 1.0 + sqrt(sf*sf - 1.0)*sigma_MC*x1;
    double jetEnergySmearFactorUp = 1.0 + sqrt(sfUp*sfUp - 1.0)*sigma_MC*x2;
    double jetEnergySmearFactorDown = 1.0 + sqrt(sfDown*sfDown - 1.0)*sigma_MC*x3;
    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;

    fSingletons["pt"] = jetCorrPtSmear;
    fSingletons["eta"] = iObjects[i0]->eta;
    fSingletons["phi"] = iObjects[i0]->phi;
    
    TAddJet *pAddJet = getAddJet(iObjects[i0]);

    

    fSingletons["mass"]  = JEC*jetEnergySmearFactor*(iObjects[i0]->mass);
    fSingletons["csv"]  = iObjects[i0]->csv;
    fSingletons["CHF"]  = iObjects[i0]->chHadFrac;
    fSingletons["NHF"]  = iObjects[i0]->neuHadFrac;
    fSingletons["NEMF"]  = iObjects[i0]->neuEmFrac;
    fSingletons["tau21"]  = (pAddJet->tau2/pAddJet->tau1);
    fSingletons["tau32"]  = (pAddJet->tau3/pAddJet->tau2);
    fSingletons["msd"]  = pAddJet->mass_sd0;
    fSingletons["rho"]  = log((pAddJet->mass_sd0*pAddJet->mass_sd0)/iObjects[i0]->pt);
    fSingletons["minsubscv"]  = TMath::Min(pAddJet->sj1_csv,pAddJet->sj2_csv);
    fSingletons["maxsubscv"] = TMath::Max(TMath::Max(pAddJet->sj1_csv,pAddJet->sj2_csv),TMath::Max(pAddJet->sj3_csv,pAddJet->sj4_csv));
    fSingletons["doublecsv"] = pAddJet->doublecsv;
    fSingletons["doublesub"] = pAddJet->Double_sub;
    fSingletons["ptraw"] = iObjects[i0]->ptRaw;
    fSingletons["genpt"] = iObjects[i0]->genpt;
    fSingletons["e2_b1"] = pAddJet->e2_b1;
    fSingletons["e3_b1"] = pAddJet->e3_b1;
    fSingletons["e3_v1_b1"] = pAddJet->e3_v1_b1;
    fSingletons["e3_v2_b1"] = pAddJet->e3_v2_b1;
    fSingletons["e4_v1_b1"] = pAddJet->e4_v1_b1;
    fSingletons["e4_v2_b1"] = pAddJet->e4_v2_b1;
    fSingletons["e2_b2"] = pAddJet->e2_b2;
    fSingletons["e3_b2"] = pAddJet->e3_b2;
    fSingletons["e3_v1_b2"] = pAddJet->e3_v1_b2;
    fSingletons["e3_v2_b2"] = pAddJet->e3_v2_b2;
    fSingletons["e4_v1_b2"] = pAddJet->e4_v1_b2;
    fSingletons["e4_v2_b2"] = pAddJet->e4_v2_b2;
    fSingletons["e2_sdb1"] = pAddJet->e2_sdb1;
    fSingletons["e3_sdb1"] = pAddJet->e3_sdb1;
    fSingletons["e3_v1_sdb1"] = pAddJet->e3_v1_sdb1;
    fSingletons["e3_v2_sdb1"] = pAddJet->e3_v2_sdb1;
    fSingletons["e4_v1_sdb1"] = pAddJet->e4_v1_sdb1;
    fSingletons["e4_v2_sdb1"] = pAddJet->e4_v2_sdb1;
    fSingletons["e2_sdb2"] = pAddJet->e2_sdb2;
    fSingletons["e3_sdb2"] = pAddJet->e3_sdb2;
    fSingletons["e3_v1_sdb2"] = pAddJet->e3_v1_sdb2;
    fSingletons["e3_v2_sdb2"] = pAddJet->e3_v2_sdb2;
    fSingletons["e4_v1_sdb2"] = pAddJet->e4_v1_sdb2;
    fSingletons["e4_v2_sdb2"] = pAddJet->e4_v2_sdb2;
    fSingletons["N2sdb1"] = pAddJet->e3_v2_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    fSingletons["N2sdb2"] = pAddJet->e3_v2_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    fSingletons["M2sdb1"] = pAddJet->e3_v1_sdb1/(pAddJet->e2_sdb1);
    fSingletons["M2sdb2"] = pAddJet->e3_v1_sdb2/(pAddJet->e2_sdb2);
    fSingletons["D2sdb1"] = pAddJet->e3_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    fSingletons["D2sdb2"] = pAddJet->e3_sdb2/(pAddJet->e2_sdb2*pAddJet->e2_sdb2*pAddJet->e2_sdb2);
    fSingletons["N2b1"] = pAddJet->e3_v2_b1/(pAddJet->e2_b1*pAddJet->e2_b1);
    fSingletons["N2b2"] = pAddJet->e3_v2_b2/(pAddJet->e2_b2*pAddJet->e2_b2);
    fSingletons["M2b1"] = pAddJet->e3_v1_b1/(pAddJet->e2_b1);
    fSingletons["M2b2"] = pAddJet->e3_v1_b2/(pAddJet->e2_b2);
    fSingletons["D2b1"] = pAddJet->e3_b1/(pAddJet->e2_b1*pAddJet->e2_b1*pAddJet->e2_b1);
    fSingletons["D2b2"] = pAddJet->e3_b2/(pAddJet->e2_b2*pAddJet->e2_b2*pAddJet->e2_b2);
    fSingletons["pt_old"] = iObjects[i0]->pt;
    fSingletons["jetPtJESUp"] = jetPtJESUp;
    fSingletons["jetPtJESDown"] = jetPtJESDown;
    fSingletons["jetPtJERUp"] = jetPtJERUp;
    fSingletons["jetPtJERDown"] = jetPtJERDown;
    fSingletons["e2_sdb05"] = pAddJet->e2_sdb05;
    fSingletons["e3_sdb05"] = pAddJet->e3_sdb05;
    fSingletons["e3_v1_sdb05"] = pAddJet->e3_v1_sdb05;
    fSingletons["e3_v2_sdb05"] = pAddJet->e3_v2_sdb05;
    fSingletons["e4_v1_sdb05"] = pAddJet->e4_v1_sdb05;
    fSingletons["e4_v2_sdb05"] = pAddJet->e4_v2_sdb05;
    fSingletons["e2_sdb4"] = pAddJet->e2_sdb4;
    fSingletons["e3_sdb4"] = pAddJet->e3_sdb4;
    fSingletons["e3_v1_sdb4"] = pAddJet->e3_v1_sdb4;
    fSingletons["e3_v2_sdb4"] = pAddJet->e3_v2_sdb4;
    fSingletons["e4_v1_sdb4"] = pAddJet->e4_v1_sdb4;
    fSingletons["e4_v2_sdb4"] = pAddJet->e4_v2_sdb4;
    fSingletons["flavour"] = pAddJet->flavour;
    fSingletons["nbHadrons"] = pAddJet->nbHadrons;
    fSingletons["nSV"] = pAddJet->nSV;
    fSingletons["jetNTracks"] = pAddJet->jetNTracks;
    fSingletons["tau_flightDistance2dSig_1"] = pAddJet->tau_flightDistance2dSig_1;
    fSingletons["SubJet_csv"] = pAddJet->SubJet_csv;
    fSingletons["z_ratio"] = pAddJet->z_ratio;
    fSingletons["trackSipdSig_3"] = pAddJet->trackSipdSig_3;
    fSingletons["trackSipdSig_2"] = pAddJet->trackSipdSig_2;
    fSingletons["trackSipdSig_1"] = pAddJet->trackSipdSig_1;
    fSingletons["trackSipdSig_0"] = pAddJet->trackSipdSig_0;
    fSingletons["trackSipdSig_1_0"] = pAddJet->trackSipdSig_1_0;
    fSingletons["trackSipdSig_0_0"] = pAddJet->trackSipdSig_0_0;
    fSingletons["trackSipdSig_1_1"] = pAddJet->trackSipdSig_1_1;
    fSingletons["trackSipdSig_0_1"] = pAddJet->trackSipdSig_0_1;
    fSingletons["trackSip2dSigAboveCharm_0"] = pAddJet->trackSip2dSigAboveCharm_0;
    fSingletons["trackSip2dSigAboveBottom_0"] = pAddJet->trackSip2dSigAboveBottom_0;
    fSingletons["trackSip2dSigAboveBottom_1"] = pAddJet->trackSip2dSigAboveBottom_1;
    fSingletons["tau1_trackEtaRel_0"] = pAddJet->tau1_trackEtaRel_0;
    fSingletons["tau1_trackEtaRel_1"] = pAddJet->tau1_trackEtaRel_1;
    fSingletons["tau1_trackEtaRel_2"] = pAddJet->tau1_trackEtaRel_2;
    fSingletons["tau0_trackEtaRel_0"] = pAddJet->tau0_trackEtaRel_0;
    fSingletons["tau0_trackEtaRel_1"] = pAddJet->tau0_trackEtaRel_1;
    fSingletons["tau0_trackEtaRel_2"] = pAddJet->tau0_trackEtaRel_2;
    fSingletons["tau_vertexMass_0"] = pAddJet->tau_vertexMass_0;
    fSingletons["tau_vertexEnergyRatio_0"] = pAddJet->tau_vertexEnergyRatio_0;
    fSingletons["tau_vertexDeltaR_0"] = pAddJet->tau_vertexDeltaR_0;
    fSingletons["tau_flightDistance2dSig_0"] = pAddJet->tau_flightDistance2dSig_0;
    fSingletons["tau_vertexMass_1"] = pAddJet->tau_vertexMass_1;
    fSingletons["tau_vertexEnergyRatio_1"] = pAddJet->tau_vertexEnergyRatio_1;

    unsigned nG = fGens->GetEntriesFast();
    unsigned nP = 0;
    double dR2 = dR * dR;
    double threshold = 0.2 * iObjects[i0]->pt;
    std::unordered_set<TGenParticle*> partons; // avoid double-counting
    for (unsigned iG = 0; iG != nG; ++iG) {
      TGenParticle *part = (TGenParticle*)((*fGens)[iG]);
      unsigned apdgid = abs(part->pdgId);
      if (apdgid > 5 &&
          apdgid != 21 &&
          apdgid != 15 &&
          apdgid != 11 &&
          apdgid != 13) 
        continue;

      TGenParticle *parent = part;
      bool foundParent = false;
      while (parent->parent > 0) {
        parent = (TGenParticle*)((*fGens)[parent->parent]);
        if (partons.find(parent) != partons.end()) {
          foundParent = true;
          break;
        }
      }
      if (foundParent)
        continue;

      if (part->pt < threshold)
        continue;
      if (deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, part->eta, part->phi) > dR2)
        continue;
      partons.insert(part);
      ++nP;
    }  
    fSingletons["nProngs"] = nP;

    // fill neutral and charged PF candidates
    std::vector<TPFPart*> jetPFs;
    for (auto idx : iObjects[i0]->pfCands) {
      jetPFs.push_back( (TPFPart*)(fPFs->At(idx)) );
    }
    std::sort(jetPFs.begin(),
              jetPFs.end(),
              [](TPFPart *x, TPFPart *y) {return x->pt > y->pt;});
    unsigned iCPF=0, iNPF=0;
    for (auto *pf : jetPFs) {
      if (pf->q && iCPF < NCPF) { // charged PF
        fCPFArrs["cpfPt"][iCPF] = pf->pt;
        iCPF++;
      } else if (pf->q == 0 && iNPF < NNPF) {
        iNPF++;
      }
    }

    // fill PF 
    std::vector<TSVtx*> jetSVs;
    for (auto idx : pAddJet->svtx) {
      jetSVs.push_back( (TSVtx*)(fSVs->At(idx)) );
    }
    std::sort(jetSVs.begin(),
              jetSVs.end(),
              [](TSVtx *x, TSVtx *y) {return x->pt > y->pt;});
    unsigned iSV=0;
    for (auto *sv : jetSVs) {
      if (iSV == NSV)
        break;
      fSVArrs["svPt"][iSV] = sv->pt;
      iSV++;
    }

    fpartonFlavor   = iObjects[0]->partonFlavor;
    fhadronFlavor   = iObjects[0]->hadronFlavor;
    fnCharged       = iObjects[0]->nCharged;
    fnNeutrals      = iObjects[0]->nNeutrals;
    fnParticles     = iObjects[0]->nParticles;
    fnVtxFlavor     = iObjects[0]->vtxFlavor;
    fnVtxFlavInfo   = iObjects[0]->vtxFlavInfo;

    fTree->Fill(); 

  }
}

void PerJetLoader::matchJet(std::vector<TLorentzVector> iJets1, TLorentzVector iJet2, double dR, int jIndex){
  TLorentzVector iJet1;
  int nmatched(0);
  float mindR = dR;
  for(int i0 = 0; i0 < int(iJets1.size()); i0++) {
    if ((iJets1[i0].DeltaR(iJet2) < mindR) && (fabs(iJets1[i0].Pt()-iJet2.Pt())<0.35*fabs(iJet2.Pt()))) {
      nmatched++;
      iJet1 = iJets1[i0];
      mindR= iJets1[i0].DeltaR(iJet2);  
    }
  }
}

TAddJet *PerJetLoader::getAddJet(TJet *iJet) { 
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

//2016 Prompt Reco
void PerJetLoader::loadJECs(bool isData) {
    std::cout << "PerJetLoader: loading jet energy correction constants" << std::endl;
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
void PerJetLoader::loadJECs_Rereco(bool isData) {
    std::cout << "PerJetLoader: loading Rereco jet energy correction constants" << std::endl;
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
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Summer16_23Sep2016BCDV3_DATA/Summer16_23Sep2016BCDV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Summer16_23Sep2016EFV3_DATA/Summer16_23Sep2016EFV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);

      correctionParameters.push_back(correctionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));

      //IOV: 2016G
      std::vector<JetCorrectorParameters> correctionParametersG = std::vector<JetCorrectorParameters> ();
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorG = new FactorizedJetCorrector(correctionParametersG);
      std::string jecUncPathG = jecPathname+"/Summer16_23Sep2016GV3_DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncG = new JetCorrectionUncertainty(jecUncPathG);

      correctionParameters.push_back(correctionParametersG);
      JetCorrector.push_back( JetCorrectorG );
      jecUnc.push_back(jecUncG);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 280385 ));

      //IOV: 2016H
      std::vector<JetCorrectorParameters> correctionParametersH = std::vector<JetCorrectorParameters> ();
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorH = new FactorizedJetCorrector(correctionParametersH);
      std::string jecUncPathH = jecPathname+"/Summer16_23Sep2016HV3_DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncH = new JetCorrectionUncertainty(jecUncPathH);

      correctionParameters.push_back(correctionParametersH);
      JetCorrector.push_back( JetCorrectorH );
      jecUnc.push_back(jecUncH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 280919, 99999999 ));

    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK8PFPuppi.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK8PFPuppi.txt", jecPathname.c_str())));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK8PFPuppi.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);

      correctionParameters.push_back(correctionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
  
}

void PerJetLoader::loadCMSSWPath() {
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in PerJetLoader::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        cmsswPath = "";
    }
    cmsswPath = std::string(cmsswPathChar);
}



// Retrieve jet energy uncertainty as a function of pt and eta
double PerJetLoader::getJecUnc( float pt, float eta , int run) {

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


//Jet Energy Corrections
double PerJetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
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
double PerJetLoader::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
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

