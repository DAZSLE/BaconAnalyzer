#include "TFile.h"
#include "TMatrixD.h"
#include "../include/EvtLoader.hh"
#include <iostream>
#include <sstream>

using namespace baconhep;

EvtLoader::EvtLoader(TTree *iTree,std::string iName,std::string iHLTFile,std::string iPUWeight) { 
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
    iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  fEvt      = new TEventInfo();
  iTree->SetBranchAddress("Info",       &fEvt);
  fEvtBr    = iTree->GetBranch("Info");
  std::cout << "trigger file " << iHLTFile << std::endl;
  fTrigger  = new TTrigger(iHLTFile);
  
  fVertices = new TClonesArray("baconhep::TVertex");
  iTree->SetBranchAddress("PV",       &fVertices);
  fVertexBr     = iTree->GetBranch("PV");
  
  TFile *lFile = TFile::Open(iPUWeight.c_str()); 
  fPUWeightHist =  (TH1F*) lFile->Get("puw");
  fPUWeightHist->SetDirectory(0);
  lFile->Close();
  fSample = (char*) (iName.c_str());

  // Lepton SFs
  TFile *fMuSF = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/scaleFactor_muon_looseid_12p9.root");
  fhMuLoose =  (TH2D*) fMuSF->Get("scaleFactor_muon_looseid_RooCMSShape");
  fhMuLoose->SetDirectory(0);
  fMuSF->Close();
  TFile *fMuSFTight = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/scaleFactor_muon_tightid_12p9.root");
  fhMuTight =  (TH2D*) fMuSFTight->Get("scaleFactor_muon_tightid_RooCMSShape");
  fhMuTight->SetDirectory(0);
  fMuSFTight->Close();
  TFile *fMuSFTrack = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/scaleFactor_muon_track.root");
  fhMuTrack = (TH1D*) fMuSFTrack->Get("htrack2");
  fhMuTrack->SetDirectory(0);
  fMuSFTrack->Close();
  TFile *fEleSF = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/scaleFactor_electron_vetoid_12p9.root");
  fhEleVeto =  (TH2D*) fEleSF->Get("scaleFactor_electron_vetoid_RooCMSShape");
  fhEleVeto->SetDirectory(0);
  fEleSF->Close();
  TFile *fEleSFTight = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/scaleFactor_electron_tightid_12p9.root");
  fhEleTight =  (TH2D*) fEleSFTight->Get("scaleFactor_electron_tightid_RooCMSShape");
  fhEleTight->SetDirectory(0);
  fEleSFTight->Close();
  TFile *fEleSFTrack = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/scaleFactor_electron_track.root");
  fhEleTrack = (TH2D*) fEleSFTrack->Get("EGamma_SF2D");
  fhEleTrack->SetDirectory(0);
  fEleSFTrack->Close();

  // Trigger Eff
  TFile *fEleTrigB = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/ele_trig_lowpt_rebinned.root");
  hEleTrigB = (TH1D*) fEleTrigB->Get("h_num");
  hEleTrigB->SetDirectory(0);
  fEleTrigB->Close();
  TFile *fEleTrigE = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/ele_trig_lowpt_rebinned.root");
  hEleTrigE = (TH1D*) fEleTrigE->Get("h_num");
  hEleTrigE->SetDirectory(0);
  fEleTrigE->Close();
  TFile *fEleTrigLow  = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/ele_trig_lowpt.root");
  hEleTrigLow = (TH2D*) fEleTrigLow->Get("hEffEtaPt");
  hEleTrigLow->SetDirectory(0);
  fEleTrigLow->Close();
  TFile *fMetTrig  = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/met_trig.root");
  hMetTrig = (TH1D*) fMetTrig->Get("numer");
  hMetTrig->SetDirectory(0);
  fMetTrig->Close();
  TFile *fPhoTrig  = TFile::Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/pho_trig.root");
  hPhoTrig = (TH1D*) fPhoTrig->Get("h_num");
  hPhoTrig->SetDirectory(0);
  fPhoTrig->Close();

}
EvtLoader::~EvtLoader() { 
  delete  fEvt;
  delete  fEvtBr;
  delete  fVertices;
  delete  fVertexBr;
}
void EvtLoader::reset() { 
  fEvtV         = 0; 
  fITrigger     = 0; 
  fIMoreTrigger = 0;
  fsf_eleTrig   = 1;
  fsf_phoTrig   = 1;
  fsf_metTrig   = 1;
  fselectBits   = 0;
  fNVtx         = 0; 
  fNPU          = 0;
  fPUWeight     = 0; 
  fScale        = 1;
  fevtWeight    = 1;
  fRho          = 0;

  fsf_lep       = 1;
  fsf_lepTrack  = 1;

  fkfactor      = 1;
  fkFactor_CENT = 1;
  fEwkCorr_CENT = 1;
  fPDF          = 1;
  fPDF_UP       = 1;
  fPDF_DO       = 1;
  fRenScale_UP  = 1;
  fRenScale_DO  = 1;
  fFacScale_UP  = 1;
  fFacScale_DO  = 1;

  fMet          = 0; 
  fMetPhi       = 0; 
  fPuppEt       = 0; 
  fPuppEtPhi    = 0; 
  fCaloMet      = 0;
  fCaloMetPhi   = 0;

  fUW           = 0;
  fUWPhi        = 0;
  fUZ           = 0;
  fUZPhi        = 0;
  fUA           = 0;
  fUAPhi        = 0;
}
void EvtLoader::setupTree(TTree *iTree) {
  reset();
  fTree = iTree;
  fTree->Branch("runNum"          ,&fRun            ,"fRun/i");
  fTree->Branch("lumiSec"         ,&fLumi           ,"fLumi/i");
  fTree->Branch("evtNum"          ,&fEvtV           ,"fEvtV/i");
  fTree->Branch("passJson"        ,&fPassJson       ,"fPassJson/i");
  fTree->Branch("metfilter"       ,&fMetFilters     ,"fMetFilters/i");
  fTree->Branch("triggerBits"     ,&fITrigger       ,"fITrigger/i");
  fTree->Branch("moreTriggerBits" ,&fIMoreTrigger   ,"fIMoreTrigger/i");
  fTree->Branch("selectBits"      ,&fselectBits     ,"fselectBits/i");
 
  fTree->Branch("sf_eleTrig"      ,&fsf_eleTrig     ,"fsf_eleTrig/D");
  fTree->Branch("sf_metTrig"      ,&fsf_metTrig     ,"fsf_metTrig/D");
  fTree->Branch("sf_phoTrig"      ,&fsf_phoTrig     ,"fsf_phoTrig/D");

  fTree->Branch("sf_lep"          ,&fsf_lep         ,"fsf_lep/D");
  fTree->Branch("sf_lepTrack"     ,&fsf_lepTrack    ,"fsf_lepTrack/D");

  fTree->Branch("npu"             ,&fNPU            ,"fNPU/i");
  fTree->Branch("npv"             ,&fNVtx           ,"fNVtx/i");
  fTree->Branch("rho"             ,&fRho            ,"fRho/F");
  fTree->Branch("puWeight"        ,&fPUWeight       ,"fPUWeight/F");
  fTree->Branch("scale1fb"        ,&fScale          ,"fScale/F");  
  fTree->Branch("evtWeight"       ,&fevtWeight      ,"fevtWeight/F");

  fTree->Branch("kfactor"         ,&fkfactor        ,"fkfactor/F");
  fTree->Branch("kfactorNLO"      ,&fkFactor_CENT   ,"fkFactor_CENT/F");
  fTree->Branch("PDF"             ,&fPDF            ,"fPDF/F");
  fTree->Branch("PDF_UP"          ,&fPDF_UP         ,"fPDF_UP/F");
  fTree->Branch("PDF_DO"          ,&fPDF_DO         ,"fPDF_DO/F");
  fTree->Branch("RenScale_UP"     ,&fRenScale_UP    ,"fRenScale_UP/F");
  fTree->Branch("RenScale_DO"     ,&fRenScale_DO    ,"fRenScale_DO/F");
  fTree->Branch("FacScale_UP"     ,&fFacScale_UP    ,"fFacScale_UP/F");
  fTree->Branch("FacScale_DO"     ,&fFacScale_DO    ,"fFacScale_DO/F");

  fTree->Branch("pfmet"           ,&fMet            ,"fMet/F");
  fTree->Branch("pfmetphi"        ,&fMetPhi         ,"fMetPhi/F");
  fTree->Branch("puppet"          ,&fPuppEt         ,"fPuppEt/F");
  fTree->Branch("puppetphi"       ,&fPuppEtPhi      ,"fPuppEtPhi/F");

  //fTree->Branch("UWmag"           ,&fUW             ,"fUW/F");
  //fTree->Branch("UWphi"           ,&fUWPhi          ,"fUWPhi/F");
  //fTree->Branch("UZmag"           ,&fUZ             ,"fUZ/F");
  //fTree->Branch("UZphi"           ,&fUZPhi          ,"fUZPhi/F");
  //fTree->Branch("UAmag"           ,&fUA             ,"fUA/F");
  //fTree->Branch("UAphi"           ,&fUAPhi          ,"fUAPhi/F");

  fRun   = 0;
  fLumi  = 0;
  fPassJson = 0;
}
void EvtLoader::load(int iEvent) { 
  fVertices ->Clear();
  fEvtBr    ->GetEntry(iEvent);
  fVertexBr ->GetEntry(iEvent);
  fRun   = fEvt->runNum;
  fLumi  = fEvt->lumiSec;
}
void EvtLoader::fillEvent(unsigned int trigBit,float lWeight, unsigned int passJson) { 
  reset();
  fRun          = fEvt->runNum;
  fLumi         = fEvt->lumiSec;
  fNPU          = fEvt->nPUmean;
  fRho          = fEvt->rhoJet;
  fPassJson     = passJson;
  fNVtx         = nVtx();
  fITrigger     = trigBit;
  fIMoreTrigger = triggerBit();
  fselectBits   = 1;
  fPUWeight     = puWeight(float(fNPU)); 
  fScale        = lWeight;
  fevtWeight    = 1;
  fMetFilters   = fEvt->metFilterFailBits;
  fEvtV         = fEvt->evtNum;
  fkfactor      = 1;
  fsf_lep       = 1;
  fsf_lepTrack  = 1;
  fMet          = fEvt->pfMETC;
  fMetPhi       = fEvt->pfMETCphi;
  fCaloMet      = fEvt->caloMET;
  fCaloMetPhi   = fEvt->caloMETphi;
  fPuppEt       = fEvt->puppETC;
  fPuppEtPhi    = fEvt->puppETCphi;
  return;
}
// Recoil
void  EvtLoader::fillRecoil(std::vector<TLorentzVector> &iLooseLeptons, std::vector<TLorentzVector> iLoosePhotons) { 
  TLorentzVector vObj1, vObj2; // stuff recoils are wrt to
  TLorentzVector vUW, vUZ, vUA;
  TLorentzVector vPuppiMET;
  vPuppiMET.SetPtEtaPhiM(fPuppEt,0,fPuppEtPhi,0);

  if(iLooseLeptons.size() > 0){
    vObj1 = iLooseLeptons[0];
    vUW = vPuppiMET+vObj1;
    fUW = vUW.Pt(); fUWPhi = vUW.Phi();

    if (iLooseLeptons.size()>1 && flep1Id + flep2Id ==0) {
      // Z recoil: opposite sign to check
      vUZ = vUW+vObj2;
      fUZ = vUZ.Pt(); fUZPhi = vUZ.Phi();
    }

  }
  if (iLoosePhotons.size()>0) {
    vObj1 = iLoosePhotons[0];

    vUA = vPuppiMET+vObj1;
    fUA = vUA.Pt(); fUAPhi = vUA.Phi();
  }
}
// Met Filter
bool EvtLoader::passFilter() { 
  return (fEvt->metFilterFailBits == 0);
}
// Vtx
int EvtLoader::nVtx() {
  unsigned int lNVertex = 0;
  for(int i0 = 0; i0 < fVertices->GetEntries(); i0++) {
    TVertex *pVertex = (TVertex*) ((*fVertices)[i0]);
    if(fabs(pVertex->z) > 24) continue;
    if(pVertex->ndof    <  4) continue;
    float pX = pVertex->x;
    float pY = pVertex->y;
    float pRho  = sqrt(pX*pX+pY*pY);
    if(pRho             > 2.) continue;
    lNVertex++;
  }
  return lNVertex;
}
bool EvtLoader::PV(){
  return fEvt->hasGoodPV;
}
// Trigger Efficiency
void EvtLoader::triggerEff(std::vector<TLorentzVector> iElectrons, std::vector<TLorentzVector> iPhotons) {
  fsf_metTrig = getVal(hMetTrig,fMetNoMu);
  if(iElectrons.size() > 0){
    float eff1=0, eff2=0;
    if (iElectrons[0].Pt()<100)                   
      eff1 = getVal2D(hEleTrigLow,iElectrons[0].Eta(),iElectrons[0].Pt());
    if (TMath::Abs(iElectrons[0].Eta())<1.4442)
      eff1 = getVal(hEleTrigB,iElectrons[0].Pt());
    if (1.5660<TMath::Abs(iElectrons[0].Eta()) && TMath::Abs(iElectrons[0].Eta())<2.5)   
      eff1 = getVal(hEleTrigE,iElectrons[0].Pt());
    if (iElectrons.size()>1) {
      if (iElectrons[1].Pt()<100)
	eff2 = getVal2D(hEleTrigLow,iElectrons[1].Eta(),iElectrons[1].Pt());
      if (TMath::Abs(iElectrons[1].Eta())<1.4442)
	eff2 = getVal(hEleTrigB,iElectrons[1].Pt());
      if (1.5660<TMath::Abs(iElectrons[1].Eta()) && TMath::Abs(iElectrons[1].Eta())<2.5)
	eff2 = getVal(hEleTrigE,iElectrons[1].Pt());
    }
    fsf_eleTrig = 1 - (1-eff1)*(1-eff2);
  }
  if(iPhotons.size()   > 0){
    fsf_phoTrig = getVal(hPhoTrig,iPhotons[0].Pt());
  }
}
void EvtLoader::addTrigger(std::string iName) { 
  fTrigString.push_back(iName);
}
bool EvtLoader::passTrigger() {
  bool lPass = false;
  for(unsigned int i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits)) lPass = true;
  }
  return lPass;
}
bool EvtLoader::passTrigger(std::string iTrigger) {
  return fTrigger->pass(iTrigger,fEvt->triggerBits);
}
unsigned int EvtLoader::triggerBit() {
  unsigned int lBit = 0;
  for(unsigned int i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits))  lBit |= 1 << i0;
  }
  return lBit;
}
// puWeight
float EvtLoader::puWeight(float iNPU) { 
  float lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(iNPU)));
  if(iNPU > 50) lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(50)));
  if(iNPU <  1) lNPVW = Float_t(fPUWeightHist->GetBinContent(fPUWeightHist->FindBin(0)));
  return lNPVW;
}
// mT
float  EvtLoader::mT(float iMet,float iMetPhi,TLorentzVector &iVec) { 
  float lDPhi = fabs(iMetPhi-iVec.Phi());
  if(fabs(lDPhi) > TMath::Pi()*2.-lDPhi) lDPhi = TMath::Pi()*2.-lDPhi;
  float lMt = sqrt(2.0*(iVec.Pt()*iMet*(1.0-cos(lDPhi))));
  return lMt;
}
void  EvtLoader::fillmT(float iMet, float iMetPhi,float iFMet, float iFMetPhi, std::vector<TLorentzVector> &lCorr, float &fmT) {
  if(lCorr.size()>0){
    TLorentzVector lVecCorr;
    for(unsigned int i0 =0; i0 < lCorr.size(); i0++) lVecCorr += lCorr[i0];
    fmT = (lVecCorr.Pt()>0) ?  mT(iMet,     iMetPhi,     lVecCorr): -999;
    if(iFMet>0) fmT = (lVecCorr.Pt()>0) ? mT(iFMet,     iFMetPhi,     lVecCorr): -999;
  }
}
// Lepton SFs / Alternative calculation
void EvtLoader::fillLepSF(std::vector<TLorentzVector> iLeptons, double lepPdgId, int islep0Tight, int islep1Tight){
  double lsf_lep =1, lTrack =1;
  if (iLeptons.size()>0) {
    if (fabs(lepPdgId)==11) {
      if (islep0Tight==1) lsf_lep = getVal2D(fhEleTight,fabs(iLeptons[0].Eta()),iLeptons[0].Pt());
      else                lsf_lep = getVal2D(fhEleVeto,fabs(iLeptons[0].Eta()),iLeptons[0].Pt());
      lTrack = getVal2D(fhEleTrack,fabs(iLeptons[0].Eta()),fNVtx);
    } else if (fabs(lepPdgId)==13) {
      if (islep0Tight==1) lsf_lep = getVal2D(fhMuTight,fabs(iLeptons[0].Eta()),iLeptons[0].Pt());
      else                lsf_lep = getVal2D(fhMuLoose,fabs(iLeptons[0].Eta()),iLeptons[0].Pt());
      lTrack = getVal(fhMuTrack,fNVtx);
    }
  }
  if (iLeptons.size()>1) {
    if (fabs(lepPdgId)==11) {
      if (islep1Tight==1) lsf_lep *= getVal2D(fhEleTight,fabs(iLeptons[1].Eta()),iLeptons[1].Pt());
      else                lsf_lep *= getVal2D(fhEleVeto,fabs(iLeptons[1].Eta()),iLeptons[1].Pt());
      lTrack *= getVal2D(fhEleTrack,fabs(iLeptons[0].Eta()),fNVtx);
    } else if (fabs(lepPdgId)==13) {
      if (islep1Tight==1) lsf_lep *= getVal2D(fhMuTight,fabs(iLeptons[1].Eta()),iLeptons[1].Pt());
      else                lsf_lep *= getVal2D(fhMuLoose,fabs(iLeptons[1].Eta()),iLeptons[1].Pt());
      lTrack *= getVal(fhMuTrack,fNVtx);
    }
  }
  fsf_lep = lsf_lep;
  fsf_lepTrack = lTrack;
}
// kFactor and EWK
void EvtLoader::computeCorr(float iPt,std::string iHist0,std::string iHist1,std::string iHist2,std::string iNLO,std::string ikfactor){
  std::string isuffix = "";
  if(iNLO.find("G")!=std::string::npos) isuffix ="_G";
  std::stringstream pPDF,pRUP,pRDO,pFUP,pFDO;
  pPDF << iNLO << "/PDF";
  pRUP << iNLO << "/ren_up" << isuffix;
  pRDO << iNLO << "/ren_down" << isuffix;
  pFUP << iNLO << "/fact_up" << isuffix;
  pFDO << iNLO << "/fact_down" << isuffix;

  TFile *lFile = TFile::Open(ikfactor.c_str());
  fHist0 =  (TH1F*) lFile->Get(iHist0.c_str()); // NLO
  fHist0->SetDirectory(0);
  fHist1 =  (TH1F*) lFile->Get(iHist1.c_str()); // LO
  fHist1->SetDirectory(0);
  fHist2 =  (TH1F*) lFile->Get(iHist2.c_str()); // NLO*EWK
  fHist2->SetDirectory(0);
  fHistPDF = (TH1F*) lFile->Get(pPDF.str().c_str());
  fHistPDF->SetDirectory(0);
  fHistRUP = (TH1F*) lFile->Get(pRUP.str().c_str());
  fHistRUP->SetDirectory(0);
  fHistRDO = (TH1F*) lFile->Get(pRDO.str().c_str());
  fHistRDO->SetDirectory(0);
  fHistFUP = (TH1F*) lFile->Get(pFUP.str().c_str());
  fHistFUP->SetDirectory(0);
  fHistFDO = (TH1F*) lFile->Get(pFDO.str().c_str());
  fHistFDO->SetDirectory(0);
  lFile->Close();

  fHist2->Divide(fHist1);
  fHist0->Divide(fHist1);

  fkFactor_CENT = Float_t(fHist0->GetBinContent(fHist0->FindBin(iPt)));
  //if(iPt > 700) fkFactor_CENT = Float_t(fHist0->GetBinContent(fHist0->FindBin(700)));
  if(iPt < 100) fkFactor_CENT = Float_t(fHist0->GetBinContent(fHist0->FindBin(100)));
  
  fEwkCorr_CENT = Float_t(fHist2->GetBinContent(fHist2->FindBin(iPt)));
  //if(iPt > 700) fEwkCorr_CENT = Float_t(fHist2->GetBinContent(fHist2->FindBin(700)));
  if(iPt < 100) fEwkCorr_CENT = Float_t(fHist2->GetBinContent(fHist2->FindBin(100)));

  fkfactor = fEwkCorr_CENT; //(NLO*ewk/LO)

  fPDF = Float_t(fHistPDF->GetBinContent(fHistPDF->FindBin(iPt)));
  if(iPt > 700) fPDF = Float_t(fHistPDF->GetBinContent(fHistPDF->FindBin(700)));
  if(iPt < 100) fPDF = Float_t(fHistPDF->GetBinContent(fHistPDF->FindBin(100)));
  fPDF_UP = 1 + fPDF;
  fPDF_DO = 1 - fPDF;
  fRenScale_UP = Float_t(fHistRUP->GetBinContent(fHistRUP->FindBin(iPt)));
  if(iPt > 700) fRenScale_UP = Float_t(fHistRUP->GetBinContent(fHistRUP->FindBin(700)));
  if(iPt < 100) fRenScale_UP = Float_t(fHistRUP->GetBinContent(fHistRUP->FindBin(100)));
  fRenScale_DO = Float_t(fHistRDO->GetBinContent(fHistRDO->FindBin(iPt)));
  if(iPt > 700) fRenScale_DO = Float_t(fHistRDO->GetBinContent(fHistRDO->FindBin(700)));
  if(iPt < 100) fRenScale_DO = Float_t(fHistRDO->GetBinContent(fHistRDO->FindBin(100)));
  fFacScale_UP = Float_t(fHistFUP->GetBinContent(fHistFUP->FindBin(iPt)));
  if(iPt > 700) fFacScale_UP = Float_t(fHistFUP->GetBinContent(fHistFUP->FindBin(700)));
  if(iPt < 100) fFacScale_UP = Float_t(fHistFUP->GetBinContent(fHistFUP->FindBin(100)));
  fFacScale_DO = Float_t(fHistFDO->GetBinContent(fHistFDO->FindBin(iPt)));
  if(iPt > 700) fFacScale_DO = Float_t(fHistFDO->GetBinContent(fHistFDO->FindBin(700)));
  if(iPt < 100) fFacScale_DO = Float_t(fHistFDO->GetBinContent(fHistFDO->FindBin(100)));
}
