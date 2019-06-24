#include "../include/PerJetLoader.hh"
#include <cmath>
#include <iostream> 

#include <string>
#include <sstream>
#include <unordered_set>

#define NCPF 50
#define NIPF 100
#define NSV 5
#define NPART 6
#define PI 3.141592654
#define MINPT 0.01

using namespace baconhep;

struct JetHistory {
  int user_idx;
  int child_idx;
};

PerJetLoader::PerJetLoader(TTree *iTree,std::string iJet,std::string iAddJet,std::string iJetCHS,std::string iAddJetCHS,int iN, bool iData,std::string iLabel) { 
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
  fYear=iLabel;

  fJEC = new JECLoader(iData,iLabel,"AK8PFPuppi");
  r = new TRandom3(1993);

  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";

  int activeAreaRepeats = 1;
  double ghostArea = MINPT;
  double ghostEtaMax = 7.0;
  activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);

  // jetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.8);
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8);
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
  for (auto &iter : fIPFArrs) 
    delete iter.second;
  for (auto &iter : fSVArrs) 
    delete iter.second;
  for (auto &iter : fPartArrs) 
    delete iter.second;
  delete activeArea;
  delete areaDef;
  delete jetDef;
  delete fJEC;
}

void PerJetLoader::reset() { 
  fNLooseVJets        = 0;
  fNTightVJets        = 0;  
  for(int i0 = 0; i0 < int(fisTightVJet.size()); i0++) fisTightVJet[i0] = -999;
  selectedVJets.clear();
  fLooseVJets.clear();
  x1List.clear();
  for (auto &iter : fSingletons) {
    fSingletons[iter.first] = 0;
  }
  fN_cpf = 0; fN_ipf = 0; fN_sv = 0;
  for (auto &iter : fCPFArrs) {
    for (unsigned i = 0; i != NCPF; ++i) 
      fCPFArrs[iter.first][i] = 0;
  }
  for (auto &iter : fIPFArrs) {
    for (unsigned i = 0; i != NIPF; ++i) 
      fIPFArrs[iter.first][i] = 0;
  }
  for (auto &iter : fSVArrs) {
    for (unsigned i = 0; i != NSV; ++i) 
      fSVArrs[iter.first][i] = 0;
  }
  for (auto &iter : fPartArrs) {
    for (unsigned i = 0; i != NPART; ++i) 
      fPartArrs[iter.first][i] = 0;
  }
}

void PerJetLoader::setupTreeQbert(TTree *iTree, std::string iJetLabel) { 
  reset();

  fSingletons.clear();
  fCPFArrs.clear();
  fIPFArrs.clear();
  fSVArrs.clear();
  fPartArrs.clear();

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
  fSingletons["tau42"] = 0;
  fSingletons["tau41"] = 0;
  fSingletons["msd"] = 0;
  fSingletons["rho"] = 0;
  fSingletons["deepdoubleb"] = 0;
  fSingletons["deepdoublec"] = 0;
  fSingletons["deepdoublecvb"] = 0;
  fSingletons["deepdoubleb_nomasssculptpen"] = 0;
  fSingletons["deepdoublec_nomasssculptpen"] = 0;
  fSingletons["deepdoublecvb_nomasssculptpen"] = 0;
  fSingletons["partonFlavor"] = 0;
  fSingletons["hadronFlavor"] = 0;
  fSingletons["nCharged"] = 0;
  fSingletons["nNeutrals"] = 0;
  fSingletons["nParticles"] = 0;
  fSingletons["vertexFlavor"] = 0;
  fSingletons["vertexFlavorInfo"] = 0;
  fSingletons["ptraw"] = 0;
  fSingletons["genpt"] = 0;
  fSingletons["e2_b1"] = 0; // Correlation function inputs beta=1
  fSingletons["e3_b1"] = 0;
  fSingletons["e3_v1_b1"] = 0;
  fSingletons["e3_v2_b1"] = 0;
  fSingletons["e4_v1_b1"] = 0;
  fSingletons["e4_v2_b1"] = 0;
  fSingletons["e2_sdb1"] = 0; // Correlation function inputs beta=1 soft-dropped 
  fSingletons["e3_sdb1"] = 0;
  fSingletons["e3_v1_sdb1"] = 0;
  fSingletons["e3_v2_sdb1"] = 0;
  fSingletons["e4_v1_sdb1"] = 0;
  fSingletons["e4_v2_sdb1"] = 0;
  fSingletons["N2sdb1"] = 0; // 2-prong ECFs observables
  fSingletons["M2sdb1"] = 0;
  fSingletons["D2sdb1"] = 0;
  fSingletons["N2b1"] = 0;
  fSingletons["M2b1"] = 0;
  fSingletons["D2b1"] = 0;
  fSingletons["pt_old"] = 0;
  fSingletons["pt_JESUp"] = 0;
  fSingletons["pt_JESDown"] = 0;
  fSingletons["pt_JERUp"] = 0;
  fSingletons["pt_JERDown"] = 0;
  fSingletons["nProngs"] = 0;
  fSingletons["nResonanceProngs"] = 0;
  fSingletons["resonanceType"] = -1;
  fSingletons["decayType"] = -1;
  fSingletons["nB"] = 0; 
  fSingletons["nC"] = 0;
  fSingletons["partonPt"] = 0;
  fSingletons["partonEta"] = 0;
  fSingletons["partonPhi"] = 0;
  fSingletons["partonM"] = 0;
  fSingletons["lepCPt"] = 0; // lep in jet
  fSingletons["lepCEta"] = 0;
  fSingletons["lepCPhi"] = 0;
  fSingletons["lepCId"] = 0;
  fSingletons["lsfCInc"] = 0;
  fSingletons["lsfC_2"] = 0;
  fSingletons["lsfC_3"] = 0;
  fSingletons["lsfC_4"] = 0;
  fSingletons["dR2lepC"] = 0;
  fSingletons["dR2lepCsj1"] = 0;
  fSingletons["dR2lepCsj2"] = 0;
  fSingletons["dR2lepCsj3"] = 0;
  fSingletons["glepPt"] = 0; //gen lep and met
  fSingletons["glepEta"] = 0;
  fSingletons["glepPhi"] = 0;
  fSingletons["glepId"] = 0;
  fSingletons["gmetPt"] = 0;
  fSingletons["gmetEta"] = 0;
  fSingletons["gmetPhi"] = 0;
  fSingletons["glepsizeR"] = 0; // some temp variables
  fSingletons["glepsizeJ"] = 0;
  fSingletons["glepdecay"] = 0;
  fSingletons["metlepphi"] = 0; // met
  fSingletons["metphi"] = 0;

  fCPFArrs["cpf_pt"] = new float[NCPF]; 
  fCPFArrs["cpf_eta"] = new float[NCPF]; 
  fCPFArrs["cpf_phi"] = new float[NCPF]; 
  fCPFArrs["cpf_m"] = new float[NCPF]; 
  fCPFArrs["cpf_e"] = new float[NCPF]; 
  fCPFArrs["cpf_q"] = new float[NCPF]; 
  fCPFArrs["cpf_pfType"] = new float[NCPF]; 
  fCPFArrs["cpf_vtxID"] = new float[NCPF]; 
  fCPFArrs["cpf_trkChi2"] = new float[NCPF]; 
  fCPFArrs["cpf_pup"] = new float[NCPF]; 
  fCPFArrs["cpf_vtxChi2"] = new float[NCPF]; 
  fCPFArrs["cpf_ecalE"] = new float[NCPF]; 
  fCPFArrs["cpf_hcalE"] = new float[NCPF]; 
  fCPFArrs["cpf_d0"] = new float[NCPF]; 
  fCPFArrs["cpf_dz"] = new float[NCPF]; 
  fCPFArrs["cpf_d0Err"] = new float[NCPF]; 
  fCPFArrs["cpf_dptdpt"] = new float[NCPF]; 
  fCPFArrs["cpf_detadeta"] = new float[NCPF]; 
  fCPFArrs["cpf_dphidphi"] = new float[NCPF]; 
  fCPFArrs["cpf_dxydxy"] = new float[NCPF]; 
  fCPFArrs["cpf_dzdz"] = new float[NCPF]; 
  fCPFArrs["cpf_dxydz"] = new float[NCPF]; 
  fCPFArrs["cpf_dphidxy"] = new float[NCPF]; 
  fCPFArrs["cpf_dlambdadz"] = new float[NCPF]; 

  fIPFArrs["ipf_pt"] = new float[NIPF]; 
  fIPFArrs["ipf_eta"] = new float[NIPF]; 
  fIPFArrs["ipf_phi"] = new float[NIPF]; 
  fIPFArrs["ipf_m"] = new float[NIPF]; 
  fIPFArrs["ipf_e"] = new float[NIPF]; 
  fIPFArrs["ipf_pfType"] = new float[NIPF]; 
  fIPFArrs["ipf_pup"] = new float[NIPF]; 
  fIPFArrs["ipf_ecalE"] = new float[NIPF]; 
  fIPFArrs["ipf_hcalE"] = new float[NIPF]; 
  fIPFArrs["ipf_d0"] = new float[NIPF]; 
  fIPFArrs["ipf_dz"] = new float[NIPF]; 
  fIPFArrs["ipf_d0Err"] = new float[NIPF]; 
  fIPFArrs["ipf_dptdpt"] = new float[NIPF]; 
  fIPFArrs["ipf_detadeta"] = new float[NIPF]; 
  fIPFArrs["ipf_dphidphi"] = new float[NIPF]; 
  fIPFArrs["ipf_dxydxy"] = new float[NIPF]; 
  fIPFArrs["ipf_dzdz"] = new float[NIPF]; 
  fIPFArrs["ipf_dxydz"] = new float[NIPF]; 
  fIPFArrs["ipf_dphidxy"] = new float[NIPF]; 
  fIPFArrs["ipf_dlambdadz"] = new float[NIPF]; 

  fSVArrs["sv_pt"] = new float[NSV];
  fSVArrs["sv_eta"] = new float[NSV];
  fSVArrs["sv_phi"] = new float[NSV];
  fSVArrs["sv_mass"] = new float[NSV];
  fSVArrs["sv_etarel"] = new float[NSV];
  fSVArrs["sv_phirel"] = new float[NSV];
  fSVArrs["sv_deltaR"] = new float[NSV];
  fSVArrs["sv_ntracks"] = new float[NSV];
  fSVArrs["sv_chi2"] = new float[NSV];
  fSVArrs["sv_ndf"] = new float[NSV];
  fSVArrs["sv_normchi2"] = new float[NSV];
  fSVArrs["sv_dxy"] = new float[NSV];
  fSVArrs["sv_dxyerr"] = new float[NSV];
  fSVArrs["sv_dxysig"] = new float[NSV];
  fSVArrs["sv_d3d"] = new float[NSV];
  fSVArrs["sv_d3derr"] = new float[NSV];
  fSVArrs["sv_d3dsig"] = new float[NSV];
  fSVArrs["sv_enratio"] = new float[NSV];

  fPartArrs["parton_pt"] = new float[NPART];
  fPartArrs["parton_eta"] = new float[NPART];
  fPartArrs["parton_phi"] = new float[NPART];
  fPartArrs["parton_m"] = new float[NPART];
  fPartArrs["parton_pdgid"] = new float[NPART];
  
  fTree = iTree;
  for (auto &iter : fSingletons) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    fTree->Branch(bname.str().c_str(), &(iter.second), (bname.str()+"/F").c_str());
  }
  fTree->Branch("n_cpf",&fN_cpf,"n_cpf/I");
  for (auto &iter : fCPFArrs) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname.str() << "[" << NCPF << "]/F";
    fTree->Branch(bname.str().c_str(), (iter.second), bname2.str().c_str());
  }
  fTree->Branch("n_ipf",&fN_ipf,"n_ipf/I");
  for (auto &iter : fIPFArrs) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname.str() << "[" << NIPF << "]/F";
    fTree->Branch(bname.str().c_str(), (iter.second), bname2.str().c_str());
  }
  fTree->Branch("n_sv",&fN_sv,"n_sv/I");
  for (auto &iter : fSVArrs) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname.str() << "[" << NSV << "]/F";
    fTree->Branch(bname.str().c_str(), (iter.second), bname2.str().c_str());
  }
  fTree->Branch("n_part",&fN_part,"n_part/I");
  for (auto &iter : fPartArrs) {
    std::stringstream bname;
    bname << iJetLabel << "_" << iter.first;
    std::stringstream bname2;
    bname2 << bname.str() << "[" << NPART << "]/F";
    fTree->Branch(bname.str().c_str(), (iter.second), bname2.str().c_str());
  }

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
			       double metphi,
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
    if(!passJetTightSel(pVJet,fYear))                                      continue;

    addJet(pVJet,fLooseVJets);
    lCount++;
    x1List.push_back(x1);
    lCountT++;
  }
  addVJet(fLooseVJets,selectedVJets);
  fNLooseVJets = lCount;
  fNTightVJets = lCountT;

  fillVJet(fN,fLooseVJets,dR,iRho,metphi,runNum);
}

double SignedDeltaPhi(double phi1, double phi2) {
    double dPhi = phi1-phi2;
    if (dPhi<-PI)
        dPhi = 2*PI+dPhi;
    else if (dPhi>PI)
        dPhi = -2*PI+dPhi;
    return dPhi;
}

void PerJetLoader::fillVJet(int iN,
                            std::vector<TJet*> &iObjects,
                            double dR, 
                            double iRho, 
			    double metphi,
                            unsigned int runNum)
{ 
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) {

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
    
    double jetCorrPtSmear = jetCorrPt*jetEnergySmearFactor;
    double jetPtJESUp = jetCorrPt*jetEnergySmearFactor*(1+unc);
    double jetPtJESDown = jetCorrPt*jetEnergySmearFactor/(1+unc);
    double jetPtJERUp = jetCorrPt*jetEnergySmearFactorUp;
    double jetPtJERDown = jetCorrPt*jetEnergySmearFactorDown;

    TAddJet *pAddJet = getAddJet(iObjects[i0]);

    fSingletons["pt"] = jetCorrPtSmear;
    fSingletons["eta"] = iObjects[i0]->eta;
    fSingletons["phi"] = iObjects[i0]->phi;
    fSingletons["mass"]  = JEC*jetEnergySmearFactor*(iObjects[i0]->mass);
    fSingletons["csv"]  = iObjects[i0]->csv;
    fSingletons["CHF"]  = iObjects[i0]->chHadFrac;
    fSingletons["NHF"]  = iObjects[i0]->neuHadFrac;
    fSingletons["NEMF"]  = iObjects[i0]->neuEmFrac;
    fSingletons["tau21"]  = (pAddJet->tau2/pAddJet->tau1);
    fSingletons["tau32"]  = (pAddJet->tau3/pAddJet->tau2);
    fSingletons["tau42"]  = (pAddJet->tau4/pAddJet->tau2);
    fSingletons["tau41"]  = (pAddJet->tau4/pAddJet->tau1);
    fSingletons["msd"]  = pAddJet->mass_sd0;
    fSingletons["rho"]  = 2*log(pAddJet->mass_sd0/jetCorrPtSmear);
    fSingletons["deepdoubleb"] = pAddJet->deepdoubleb;
    fSingletons["deepdoublec"] = pAddJet->deepdoublec;
    fSingletons["deepdoublecvb"] = pAddJet->deepdoublecvb;
    fSingletons["deepdoubleb_nomasssculptpen"] = pAddJet->deepdoubleb_nomasssculptpen;
    fSingletons["deepdoublec_nomasssculptpen"] = pAddJet->deepdoublec_nomasssculptpen;
    fSingletons["deepdoublecvb_nomasssculptpen"] = pAddJet->deepdoublecvb_nomasssculptpen;
    fSingletons["partonFlavor"] = iObjects[i0]->partonFlavor;
    fSingletons["hadronFlavor"] = iObjects[0]->hadronFlavor;
    fSingletons["nCharged"] = iObjects[0]->nCharged;
    fSingletons["nNeutrals"] = iObjects[0]->nNeutrals;
    fSingletons["nParticles"] = iObjects[0]->nParticles;
    fSingletons["vtxFlavor"] = iObjects[0]->vtxFlavor;
    fSingletons["vtxFlavInfo"] = iObjects[0]->vtxFlavInfo;
    fSingletons["ptraw"] = iObjects[i0]->ptRaw;
    fSingletons["genpt"] = iObjects[i0]->genpt;
    fSingletons["e2_b1"] = pAddJet->e2_b1;
    fSingletons["e3_b1"] = pAddJet->e3_b1;
    fSingletons["e3_v1_b1"] = pAddJet->e3_v1_b1;
    fSingletons["e3_v2_b1"] = pAddJet->e3_v2_b1;
    fSingletons["e4_v1_b1"] = pAddJet->e4_v1_b1;
    fSingletons["e4_v2_b1"] = pAddJet->e4_v2_b1;
    fSingletons["e2_sdb1"] = pAddJet->e2_sdb1;
    fSingletons["e3_sdb1"] = pAddJet->e3_sdb1;
    fSingletons["e3_v1_sdb1"] = pAddJet->e3_v1_sdb1;
    fSingletons["e3_v2_sdb1"] = pAddJet->e3_v2_sdb1;
    fSingletons["e4_v1_sdb1"] = pAddJet->e4_v1_sdb1;
    fSingletons["e4_v2_sdb1"] = pAddJet->e4_v2_sdb1;
    fSingletons["N2sdb1"] = pAddJet->e3_v2_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    fSingletons["M2sdb1"] = pAddJet->e3_v1_sdb1/(pAddJet->e2_sdb1);
    fSingletons["D2sdb1"] = pAddJet->e3_sdb1/(pAddJet->e2_sdb1*pAddJet->e2_sdb1*pAddJet->e2_sdb1);
    fSingletons["N2b1"] = pAddJet->e3_v2_b1/(pAddJet->e2_b1*pAddJet->e2_b1);
    fSingletons["M2b1"] = pAddJet->e3_v1_b1/(pAddJet->e2_b1);
    fSingletons["D2b1"] = pAddJet->e3_b1/(pAddJet->e2_b1*pAddJet->e2_b1*pAddJet->e2_b1);
    fSingletons["pt_old"] = iObjects[i0]->pt;
    fSingletons["jetPtJESUp"] = jetPtJESUp;
    fSingletons["jetPtJESDown"] = jetPtJESDown;
    fSingletons["jetPtJERUp"] = jetPtJERUp;
    fSingletons["jetPtJERDown"] = jetPtJERDown;
    fSingletons["nC"] = (iObjects[i0]->vtxFlavor % 1000) / 100;
    fSingletons["nB"] = (iObjects[i0]->vtxFlavor - (iObjects[i0]->vtxFlavor % 1000)) / 1000; 

    fSingletons["lepCPt"]= pAddJet->lepCPt;
    fSingletons["lepCEta"] = pAddJet->lepCEta;
    fSingletons["lepCPhi"] = pAddJet->lepCPhi;
    fSingletons["lepCId"] = pAddJet->lepCId;
    fSingletons["lsfCInc"] = pAddJet->lsfCInc;
    fSingletons["lsfC_2"]= pAddJet->lsfC_2;
    fSingletons["lsfC_3"] = pAddJet->lsfC_3;
    fSingletons["lsfC_4"] = pAddJet->lsfC_4;
    if(pAddJet->lepCPt > 0) {
      fSingletons["dR2lepC"] = deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, pAddJet->lepCEta, pAddJet->lepCPhi);
      fSingletons["dR2lepCsj1"] = deltaR2(pAddJet->lepCEta, pAddJet->lepCPhi, pAddJet->lsfC_3_sj1_eta, pAddJet->lsfC_3_sj1_phi);
      fSingletons["dR2lepCsj2"] = deltaR2(pAddJet->lepCEta, pAddJet->lepCPhi, pAddJet->lsfC_3_sj2_eta, pAddJet->lsfC_3_sj2_phi);
      fSingletons["dR2lepCsj3"] = deltaR2(pAddJet->lepCEta, pAddJet->lepCPhi, pAddJet->lsfC_3_sj3_eta, pAddJet->lsfC_3_sj3_phi);
      fSingletons["metlepphi"] = SignedDeltaPhi(metphi, pAddJet->lepCPhi);

    }
    else{
      fSingletons["dR2lepC"] = -1;
      fSingletons["dR2lepCsj1"] = -1;
      fSingletons["dR2lepCsj2"] = -1;
      fSingletons["dR2lepCsj3"] = -1;
      fSingletons["dR2metlep"] = -1;
    }
    fSingletons["metphi"] = SignedDeltaPhi(metphi, iObjects[i0]->phi);

    unsigned nG = fGens->GetEntriesFast();
    unsigned nP = 0;
    double dR2 = dR * dR;

    double threshold = 0.2 * iObjects[i0]->pt;
    auto matchJet = [jet = iObjects[i0], dR2](TGenParticle *p) -> bool {
      return deltaR2(jet->eta, jet->phi, p->eta, p->phi) < dR2;
    };

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

      if (part->pt < threshold)
        continue;

      if (!matchJet(part))
        continue;

      TGenParticle *parent = part;
      TGenParticle *foundParent = NULL;
      while (parent->parent > 0) {
        parent = (TGenParticle*)((*fGens)[parent->parent]);
        if (partons.find(parent) != partons.end()) {
          foundParent = parent;
          break;
        }
      }

      // check if the particle has a 1->2 splitting where the daughters satisfy
      // the z-cut condition and match the jet cone
      TGenParticle *dau1 = NULL, *dau2 = NULL;
      for (unsigned jG = 0; jG != nG; ++jG) {
        TGenParticle *child = (TGenParticle*)(*fGens)[jG];
        if (child->parent != (int)iG)
          continue; // only consider splittings from the current particle 

        if (dau1)
          dau2 = child; 
        else 
          dau1 = child; 
        if (dau1 && dau2)
          break; // ! assume it's not possible to have 1->N for N>2
                 // ... this is not a good assumption in the miniaod-compressed
                 // showering history
      }
      if (dau1 && dau1->pt > threshold && matchJet(dau1) && 
          dau2 && dau2->pt > threshold && matchJet(dau2)) {
        if (foundParent) {
          partons.erase(partons.find(foundParent)); // remove the ancestor
        }
        partons.insert(dau1);
        partons.insert(dau2);
      } else if (foundParent) {
        // this means we found an ancestor parton, but this isn't the vertex that gives
        // a large 1->2 splitting, so we can skip it as an intermediary
        continue; 
      } else {
        // it passes all the checks and doesn't have an ancestor - keep it!
        partons.insert(part);
      }

    }  
    nP = std::min((int)partons.size(), NPART);
    fSingletons["nProngs"] = nP;

    std::vector<TGenParticle*> vPartons; vPartons.reserve(nP);
    for (auto &iter : partons) 
      vPartons.push_back(iter);
    auto ptsort = [](TGenParticle *p1, TGenParticle *p2) -> bool {
      return p1->pt > p2->pt;
    };
    std::sort(vPartons.begin(), vPartons.end(), ptsort);
    TLorentzVector vPartonSum(0,0,0,0);
    TLorentzVector vTmp;
    for (unsigned iP = 0; iP != nP; ++iP) {
      TGenParticle *part = vPartons.at(iP);
      fPartArrs["parton_pt"][iP] = part->pt;
      fPartArrs["parton_eta"][iP] = part->eta;
      fPartArrs["parton_phi"][iP] = part->phi;
      fPartArrs["parton_m"][iP] = part->mass;
      fPartArrs["parton_pdgid"][iP] = part->pdgId;
      vTmp.SetPtEtaPhiM(part->pt, part->eta, part->phi, part->mass);
      vPartonSum += vTmp;
    }
    fSingletons["partonPt"] = vPartonSum.Pt();
    fSingletons["partonEta"] = vPartonSum.Eta();
    fSingletons["partonPhi"] = vPartonSum.Phi();
    fSingletons["partonM"] = vPartonSum.M();

    ///////
    ///look for resonances
    fSingletons["resonanceType"] = -1;
    // start with top 
    unsigned target = 6;
    for (unsigned iG = 0; iG != nG; ++iG) {
      TGenParticle *part = (TGenParticle*)((*fGens)[iG]);
      if (abs(part->pdgId) != target)
        continue; 
      if (deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, part->eta, part->phi) > dR2)
        continue;

      TGenParticle *pW = 0, *pB = 0;
      for (unsigned jG = 0; jG != nG; ++jG) {
        TGenParticle *child = (TGenParticle*)((*fGens)[jG]);
        if (child->parent != (int)iG)
          continue;
        switch (abs(child->pdgId)) {
          case 5:
            pB = child; 
            break;
          case 24:
            pW = child; 
            break;
        };
        if (pW && pB)
          break;
      }

      if (!pW || !pB)
        continue;
      TGenParticle *pQ1 = 0, *pQ2 = 0;
      for (unsigned jG = 0; jG != nG; ++jG) {
        TGenParticle *child = (TGenParticle*)((*fGens)[jG]);
        if (abs(child->pdgId) > 5)
          continue; 
        if (child->parent < 0)
          continue;

        bool foundW = false; 
        // direct parent must be a W, but not necessarily the W directly from the t 
        TGenParticle *parent = (TGenParticle*)((*fGens)[child->parent]);
        if (abs(parent->pdgId) != 24)
          continue; 

        while (!foundW) {
          if (parent == pW) {
            foundW = true; 
            if (!pQ1) 
              pQ1 = child; 
            else 
              pQ2 = child;
          } 
          if (parent->parent < 0)
            break; 
          parent = (TGenParticle*)((*fGens)[parent->parent]);
        }
        if (pQ1 && pQ2)
          break;
      }

      if (!(pB && pW && pQ1 && pQ2))
        continue; 

      // now calculate the top size 
      double size = 0; 
      for (auto *child : {pB, pQ1, pQ2}) {
        size = TMath::Max(size,
                          deltaR2(part->eta, part->phi, child->eta, child->phi));
      }
      if (size > 1.8 * dR2)
        continue; 

      fSingletons["nResonanceProngs"] = 3;
      fSingletons["resonanceType"] = 4; // 0=q/g, 1=Z, 2=W, 3=H, 4=top
      fSingletons["decayType"] = 0;
    }
  
    // first check for H
    if (fSingletons["resonanceType"] < 0) {
      unsigned target = 25;
      for (unsigned iG = 0; iG != nG; ++iG) {
        TGenParticle *part = (TGenParticle*)((*fGens)[iG]);

        if (abs(part->pdgId) != target)
          continue;

        if (deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, part->eta, part->phi) > 4*dR2){
          continue;
	}

        TGenParticle *pQ1 = 0, *pQ2 = 0;
        unsigned decay_id = 0;
        for (unsigned jG = 0; jG != nG; ++jG) {
          TGenParticle *child = (TGenParticle*)((*fGens)[jG]);
          if (child->parent < 0)
            continue;
          TGenParticle *parent = (TGenParticle*)((*fGens)[child->parent]);

          if (abs(parent->pdgId) != 25)
            continue;
	  if (abs(child->pdgId) == 25)
	    continue;
	 
          bool foundP = false;
          while (!foundP) {
            if (parent == part) {
              foundP = true;
              if (!pQ1)
                pQ1 = child;
              else
                pQ2 = child;
            }
            if (parent->parent < 0)
              break;
            parent = (TGenParticle*)((*fGens)[parent->parent]);
          }
          if (pQ1 && pQ2)
            break;
        }

        if (!(pQ1 && pQ2))
          continue;

	if(abs(pQ1->pdgId) <= 1 || abs(pQ2->pdgId) <= 2) decay_id = 1;
        if((abs(pQ1->pdgId) > 2 && abs(pQ1->pdgId) <= 4) || (abs(pQ2->pdgId) > 2 && abs(pQ2->pdgId) <= 4)) decay_id = 2;
        if((abs(pQ1->pdgId) > 4 && abs(pQ1->pdgId) <= 6) || (abs(pQ2->pdgId) > 4 && abs(pQ2->pdgId) <= 6)) decay_id = 3;
        if(abs(pQ1->pdgId) == 15 || abs(pQ2->pdgId) == 15) decay_id = 4;
        if(abs(pQ1->pdgId) == 21 || abs(pQ2->pdgId) == 21) decay_id = 5;
        if(abs(pQ1->pdgId) == 23 || abs(pQ2->pdgId) == 23) decay_id = 6;
        if(abs(pQ1->pdgId) == 24 || abs(pQ2->pdgId) == 24) decay_id = 7;


	if(abs(pQ1->pdgId) == 24 || abs(pQ1->pdgId) == 23 || abs(pQ1->pdgId) == 15) {
	  TGenParticle *pD1 = 0, *pD2 = 0, *pD3 = 0, *pD4 = 0;
	  for (unsigned kG = 0; kG != nG; ++kG) {
	    TGenParticle *child = (TGenParticle*)((*fGens)[kG]);
	    if (child->parent < 0)
	      continue;
	    TGenParticle *parent = (TGenParticle*)((*fGens)[child->parent]);
	    if (abs(parent->pdgId) != abs(pQ1->pdgId))
	      continue;
	    if (abs(child->pdgId) == abs(pQ1->pdgId) || abs(child->pdgId) == 22)
	      continue;
	    bool foundPD = false;
	    while (!foundPD) {
	      if (parent == part) {
		foundPD = true;
		if (!pD1)
		  pD1 = child;
		else if (!pD2)
		  pD2 = child;
		else if (!pD3)
		  pD3 = child;
		else
		  pD4 = child;
	      }
	      if (parent->parent < 0)
		break;
	      parent = (TGenParticle*)((*fGens)[parent->parent]);
	    }
	    if (pD1 && pD2 && pD3 && pD4)
	      break;
	  }
	  if (!(pD1 && pD2 && pD3 && pD4))
	    continue;
	  
	  // this is assuming certain order on decay which is not unreasonable thing to ask
	  double size = 0, sizej = 0;
	  double lepPt(0.), lepEta(0.), lepPhi(0.), lepId(0.);
	  double metPt(0.), metEta(0.), metPhi(0.);
	  if((abs(pD1->pdgId) <= 5 && abs(pD2->pdgId) <= 5) && (abs(pD3->pdgId) <= 5 && abs(pD4->pdgId) <= 5)) {
	    decay_id = 8; //qqqq
	    for (auto *child : {pD1, pD2, pD3, pD4}) {
	      size = TMath::Max(size,deltaR2(part->eta, part->phi, child->eta, child->phi));
	    }
	  }
          if((abs(pD1->pdgId) <= 5 && abs(pD2->pdgId) <= 5) && (abs(pD3->pdgId) >= 11 && abs(pD4->pdgId) >= 11)) {
	    decay_id = 9; //qqlv
	    if(abs(pD3->pdgId) == 11 || abs(pD3->pdgId) == 13 || abs(pD3->pdgId) == 15){
	      lepPt = pD3->pt; lepEta = pD3->eta; lepPhi = pD3->phi; lepId = abs(pD3->pdgId);
	      metPt = pD4->pt; metEta = pD4->eta; metPhi = pD4->phi;
	      for (auto *child : {pD1, pD2, pD3}) {
		size = TMath::Max(size,deltaR2(part->eta, part->phi, child->eta, child->phi));
                sizej =  TMath::Max(sizej,deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, child->eta, child->phi));
	      }
	    }
            if(abs(pD4->pdgId) == 11 || abs(pD4->pdgId) == 13 || abs(pD4->pdgId) == 15){
	      lepPt = pD4->pt; lepEta = pD4->eta; lepPhi = pD4->phi; lepId = abs(pD4->pdgId);
	      metPt = pD3->pt; metEta = pD3->eta; metPhi = pD3->phi;
              for (auto *child : {pD1, pD2, pD4}) {
                size = TMath::Max(size,deltaR2(part->eta, part->phi, child->eta, child->phi));
                sizej =  TMath::Max(sizej,deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, child->eta, child->phi));
              }
	    }
	  }
	  if((abs(pD3->pdgId) <= 5 && abs(pD4->pdgId) <= 5) && (abs(pD1->pdgId) >= 11 && abs(pD2->pdgId) >= 11)) {
            decay_id = 9; //lvqq
            if(abs(pD1->pdgId) == 11 || abs(pD1->pdgId) == 13 || abs(pD1->pdgId) == 15){
              lepPt = pD1->pt; lepEta = pD1->eta; lepPhi = pD1->phi; lepId = abs(pD1->pdgId);
              metPt = pD2->pt; metEta = pD2->eta; metPhi = pD2->phi;
              for (auto *child : {pD3, pD4, pD1}) {
                size = TMath::Max(size,deltaR2(part->eta, part->phi, child->eta, child->phi));
		sizej =  TMath::Max(sizej,deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, child->eta, child->phi));
              }
            }
            if(abs(pD2->pdgId) == 11 || abs(pD2->pdgId) == 13 || abs(pD2->pdgId) == 15){
              lepPt = pD2->pt; lepEta = pD2->eta; lepPhi = pD2->phi; lepId = abs(pD2->pdgId);
              metPt = pD1->pt; metEta = pD1->eta; metPhi = pD1->phi;
              for (auto *child : {pD3, pD4, pD2}) {
                size = TMath::Max(size,deltaR2(part->eta, part->phi, child->eta, child->phi));
                sizej =  TMath::Max(sizej,deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, child->eta, child->phi));
              }
            }
          }

	  fSingletons["glepPt"] = lepPt;
	  fSingletons["glepEta"] = lepEta;
	  fSingletons["glepPhi"] = lepPhi;
	  fSingletons["glepId"] = lepId;
	  fSingletons["gmetPt"] = metPt;
	  fSingletons["gmetEta"] = metEta;
	  fSingletons["gmetPhi"] = metPhi;
	  fSingletons["glepsizeR"] = size;
	  fSingletons["glepsizeJ"] = sizej;
	  fSingletons["glepdecay"] = decay_id;

          if ((size > 8 * dR2) && (sizej > 8 * dR2)){
            continue;
	  }
	}
	else {
	  double size = 0;
	  for (auto *child : {pQ1, pQ2}) {
	    size = TMath::Max(size,
			      deltaR2(part->eta, part->phi, child->eta, child->phi));
	  }
	  if (size > 1.8 * dR2)
	    continue;
	}
	fSingletons["nResonanceProngs"] = 2;
	fSingletons["resonanceType"] = 3; // 0=q/g, 1=Z, 2=W, 3=H, 4=top, 5=Z'                                                                                                             
	fSingletons["decayType"] = decay_id; // 1:u/d, 2:c/s, 3:b, 4:tautau, 5:gluglu, 6:ZZ, 7:WW (taus or non-mathed), 8:qqqq, 9:qqlv or lvqq
      }
    }

    if (fSingletons["resonanceType"] < 0) {
      std::vector<unsigned> targets = {23, 24, 55};
      for (unsigned iG = 0; iG != nG; ++iG) {
        TGenParticle *part = (TGenParticle*)((*fGens)[iG]);

        auto found_target = std::find(targets.begin(),targets.end(),abs(part->pdgId));
        if (found_target == targets.end())
          continue; 
        unsigned found_id = *found_target;

        if (deltaR2(iObjects[i0]->eta, iObjects[i0]->phi, part->eta, part->phi) > dR2)
          continue;

        TGenParticle *pQ1 = 0, *pQ2 = 0; 
	unsigned decay_id = 0;
        for (unsigned jG = 0; jG != nG; ++jG) {
          TGenParticle *child = (TGenParticle*)((*fGens)[jG]);
          if (abs(child->pdgId) > 5) 
	    continue; 
          if (child->parent < 0)
            continue;

          // direct parent must be a W, but not necessarily the W
          TGenParticle *parent = (TGenParticle*)((*fGens)[child->parent]);

          if (abs(parent->pdgId) != found_id)
            continue; 

          bool foundP = false;
          while (!foundP) {
            if (parent == part) {
              foundP = true; 
              if (!pQ1) 
                pQ1 = child; 
              else 
                pQ2 = child;
            } 
            if (parent->parent < 0)
              break; 
            parent = (TGenParticle*)((*fGens)[parent->parent]);
          }
          if (pQ1 && pQ2)
            break;
        }

        if (!(pQ1 && pQ2)){
          continue; 
	}

	if(abs(pQ1->pdgId) <= 1 || abs(pQ2->pdgId) <= 2) decay_id = 1;
	if((abs(pQ1->pdgId) > 2 && abs(pQ1->pdgId) <= 4) || (abs(pQ2->pdgId) > 2 && abs(pQ2->pdgId) <= 4)) decay_id = 2;
        if((abs(pQ1->pdgId) > 4 && abs(pQ1->pdgId) <= 6) || (abs(pQ2->pdgId) > 4 && abs(pQ2->pdgId) <= 6)) decay_id = 3;
	if(abs(pQ1->pdgId) == 17 || abs(pQ2->pdgId) == 17) decay_id = 4;
        if(abs(pQ1->pdgId) == 21 || abs(pQ2->pdgId) == 21) decay_id = 5;
	if(abs(pQ1->pdgId) == 23 || abs(pQ2->pdgId) == 23) decay_id = 6;
	if(abs(pQ1->pdgId) == 24 || abs(pQ2->pdgId) == 24) decay_id = 7;

        // now calculate the size 
        double size = 0; 
        for (auto *child : {pQ1, pQ2}) {
          size = TMath::Max(size,
                            deltaR2(part->eta, part->phi, child->eta, child->phi));
        }
        if (size > 1.8 * dR2){
          continue; 
	}

        fSingletons["nResonanceProngs"] = 2;
        fSingletons["resonanceType"] = found_id - 22; // 0=q/g, 1=Z, 2=W, 3=H, 4=top, 5=Z'
	if(found_id==55) fSingletons["resonanceType"] = 5;
	fSingletons["decayType"] = decay_id; // 1:u/d, 2:c/s, 3:b, 4:tautau, 5:gluglu, 6:ZZ, 7:WW
      }
    }

    if (fSingletons["resonanceType"] < 0) {
      fSingletons["resonanceType"] = 0;
      fSingletons["nResonanceProngs"] = 1;
      fSingletons["decayType"] = 0;
    }

    //////

    // fill neutral and charged PF candidates
    std::vector<TPFPart*> jetPFs;
    if (aktSort) {
      TLorentzVector vtmp;
      std::vector<fastjet::PseudoJet> pjs;
      for (auto idx : iObjects[i0]->pfCands) {
        TPFPart *part = ((TPFPart*)(fPFs->At(idx)));
        vtmp.SetPtEtaPhiM(part->pt,part->eta,part->phi,part->m);
        pjs.emplace_back(vtmp.Px(), vtmp.Py(), vtmp.Pz(), vtmp.E());
        pjs.back().set_user_index(idx);
      }
      fastjet::ClusterSequenceArea seq(pjs, *jetDef, *areaDef);

      auto &history = seq.history();
      auto &jets = seq.jets();

      std::vector<JetHistory> ordered_jets;
      for (auto &h : history) {
        // printf("%i %i %i %i ", h.parent1, h.parent2, h.child, h.jetp_index);
        if (h.jetp_index >= 0) {
          // printf("%.3f ", jets.at(h.jetp_index).perp());
          auto &j = jets.at(h.jetp_index);
          if (j.user_index() >= 0) {
            // printf("known particle");
            JetHistory jh;
            jh.user_idx = j.user_index();
            jh.child_idx = h.child;
            ordered_jets.push_back(jh);
          }
        }
        // printf("\n");
      }
      std::sort(ordered_jets.begin(),
                ordered_jets.end(),
                [](JetHistory x, JetHistory y) {return x.child_idx < y.child_idx;});
      for (auto &jh : ordered_jets) {
        jetPFs.push_back( (TPFPart*)(fPFs->At(jh.user_idx)) );
      }
    } else { // pT sort 
      for (auto idx : iObjects[i0]->pfCands) {
        jetPFs.push_back( (TPFPart*)(fPFs->At(idx)) );
      }
      std::sort(jetPFs.begin(),
                jetPFs.end(),
                [](TPFPart *x, TPFPart *y) {return x->pt * x->pup > y->pt * y->pup;});
    }
    int iCPF=0, iIPF=0;
    for (auto *pf : jetPFs) {
      if (pf->q && iCPF < NCPF) { // charged PF
        fCPFArrs["cpf_pt"][iCPF] = pf->pup * pf->pt / iObjects[i0]->pt; 
        fCPFArrs["cpf_eta"][iCPF] = pf->eta - iObjects[i0]->eta; 
        fCPFArrs["cpf_phi"][iCPF] =SignedDeltaPhi(pf->phi, iObjects[i0]->phi); 
        fCPFArrs["cpf_m"][iCPF] = pf->m; 
        fCPFArrs["cpf_e"][iCPF] = pf->e; 
        fCPFArrs["cpf_q"][iCPF] = pf->q; 
        fCPFArrs["cpf_pfType"][iCPF] = pf->pfType; 
        fCPFArrs["cpf_vtxID"][iCPF] = pf->vtxId; 
        fCPFArrs["cpf_trkChi2"][iCPF] = pf->trkChi2; 
        fCPFArrs["cpf_pup"][iCPF] = pf->pup; 
        fCPFArrs["cpf_vtxChi2"][iCPF] = pf->vtxChi2; 
        fCPFArrs["cpf_ecalE"][iCPF] = pf->ecalE; 
        fCPFArrs["cpf_hcalE"][iCPF] = pf->hcalE; 
        fCPFArrs["cpf_d0"][iCPF] = pf->d0; 
        fCPFArrs["cpf_dz"][iCPF] = pf->dz; 
        fCPFArrs["cpf_d0Err"][iCPF] = pf->d0Err; 
        fCPFArrs["cpf_dptdpt"][iCPF] = pf->dptdpt; 
        fCPFArrs["cpf_detadeta"][iCPF] = pf->detadeta; 
        fCPFArrs["cpf_dphidphi"][iCPF] = pf->dphidphi; 
        fCPFArrs["cpf_dxydxy"][iCPF] = pf->dxydxy; 
        fCPFArrs["cpf_dzdz"][iCPF] = pf->dzdz; 
        fCPFArrs["cpf_dxydz"][iCPF] = pf->dxydz; 
        fCPFArrs["cpf_dphidxy"][iCPF] = pf->dphidxy; 
        fCPFArrs["cpf_dlambdadz"][iCPF] = pf->dlambdadz; 
        iCPF++;
      } 
      if (iIPF >= NIPF)
        continue;
      fIPFArrs["ipf_pt"][iIPF] = pf->pup * pf->pt / iObjects[i0]->pt; 
      fIPFArrs["ipf_eta"][iIPF] = pf->eta - iObjects[i0]->eta; 
      fIPFArrs["ipf_phi"][iIPF] = SignedDeltaPhi(pf->phi, iObjects[i0]->phi); 
      fIPFArrs["ipf_m"][iIPF] = pf->m; 
      fIPFArrs["ipf_e"][iIPF] = pf->e; 
      fIPFArrs["ipf_pfType"][iIPF] = pf->pfType; 
      fIPFArrs["ipf_pup"][iIPF] = pf->pup; 
      fIPFArrs["ipf_ecalE"][iIPF] = pf->ecalE; 
      fIPFArrs["ipf_hcalE"][iIPF] = pf->hcalE; 
      iIPF++;
    }
    fN_cpf = std::min(iCPF, NCPF);
    fN_ipf = std::min(iIPF, NIPF);

    // fill PF 
    std::vector<TSVtx*> jetSVs;
    for (auto idx : pAddJet->svtx) {
      jetSVs.push_back( (TSVtx*)(fSVs->At(idx)) );
    }
    std::sort(jetSVs.begin(),
              jetSVs.end(),
              [](TSVtx *x, TSVtx *y) {return x->pt > y->pt;});
    int iSV=0;
    for (auto *sv : jetSVs) {
      if (iSV == NSV)
        break;
      fSVArrs["sv_pt"][iSV] = sv->pt / iObjects[i0]->pt;
      fSVArrs["sv_eta"][iSV] = sv->eta - iObjects[i0]->eta;
      fSVArrs["sv_phi"][iSV] = SignedDeltaPhi(sv->phi, iObjects[i0]->phi);
      fSVArrs["sv_mass"][iSV] = sv->mass;
      fSVArrs["sv_etarel"][iSV] = sv->etarel;
      fSVArrs["sv_phirel"][iSV] = sv->phirel;
      fSVArrs["sv_deltaR"][iSV] = sv->sv_deltaR;
      fSVArrs["sv_ntracks"][iSV] = sv->sv_ntracks;
      fSVArrs["sv_chi2"][iSV] = sv->sv_chi2;
      fSVArrs["sv_ndf"][iSV] = sv->sv_ndf;
      fSVArrs["sv_normchi2"][iSV] = sv->sv_normchi2;
      fSVArrs["sv_dxy"][iSV] = sv->sv_dxy;
      fSVArrs["sv_dxyerr"][iSV] = sv->sv_dxyerr;
      fSVArrs["sv_dxysig"][iSV] = sv->sv_dxysig;
      fSVArrs["sv_d3d"][iSV] = sv->sv_d3d;
      fSVArrs["sv_d3derr"][iSV] = sv->sv_d3derr;
      fSVArrs["sv_d3dsig"][iSV] = sv->sv_d3dsig;
      fSVArrs["sv_enratio"][iSV] = sv->sv_enratio;      
      iSV++;
    }
    fN_sv = std::min(iSV, NSV);

    fpartonFlavor   = iObjects[0]->partonFlavor;
    fhadronFlavor   = iObjects[0]->hadronFlavor;
    fnCharged       = iObjects[0]->nCharged;
    fnNeutrals      = iObjects[0]->nNeutrals;
    fnParticles     = iObjects[0]->nParticles;
    fnVtxFlavor     = iObjects[0]->vtxFlavor;
    fnVtxFlavInfo   = iObjects[0]->vtxFlavInfo;

    if(jetCorrPtSmear <=  400 && pAddJet->mass_sd0 <= 40) continue;

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
