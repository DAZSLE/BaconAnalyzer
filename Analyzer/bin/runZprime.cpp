//================================================================================================
//
// Perform preselection for W/Zprime(->qq)+jets events and produce bacon bits 
//
// Input arguments
//   argv[1] => lName = input bacon file name
//   argv[2] => lOption = dataset type: "mc", "data"
//   argv[3] => lJSON = JSON file for run-lumi filtering of data, specify "none" for MC or no filtering
//   argv[4] => lXS = cross section (pb), ignored for data 
//   argv[5] => weight = total weight, ignored for data
//________________________________________________________________________________________________

#include "../include/GenLoader.hh"
#include "../include/EvtLoader.hh"
#include "../include/ElectronLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/PhotonLoader.hh"
#include "../include/TauLoader.hh"
#include "../include/JetLoader.hh"
#include "../include/VJetLoader.hh"
#include "../include/RunLumiRangeMap.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>

// Object Processors
GenLoader       *fGen        = 0; 
EvtLoader       *fEvt        = 0; 
MuonLoader      *fMuon       = 0; 
ElectronLoader  *fElectron   = 0; 
TauLoader       *fTau        = 0; 
PhotonLoader    *fPhoton     = 0; 
JetLoader       *fJet4       = 0;
VJetLoader      *fVJet8      = 0;
VJetLoader      *fVJet15     = 0;
RunLumiRangeMap *fRangeMap   = 0; 

TH1F *fHist                  = 0;


const int NUM_PDF_WEIGHTS = 60;

// Load tree and return infile
TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  fHist        = (TH1F* ) lFile->FindObjectAny("TotalEvents");
  return lTree;
}

// For Json 
bool passEvent(unsigned int iRun,unsigned int iLumi) { 
  RunLumiRangeMap::RunLumiPairType lRunLumi(iRun,iLumi);
  return fRangeMap->HasRunLumi(lRunLumi);
}

// === MAIN =======================================================================================================
int main( int argc, char **argv ) {
  gROOT->ProcessLine("#include <vector>");
  const std::string lName        = argv[1];
  const std::string lOption      = argv[2];
  const std::string lJSON        = argv[3];
  const double      lXS          = atof(argv[4]);
  //const double      weight       = atof(argv[5]);

  fRangeMap = new RunLumiRangeMap();
  if(lJSON.size() > 0) fRangeMap->AddJSONFile(lJSON.c_str());

  bool isData;
  if(lOption.compare("data")!=0) isData = false;
  else isData = true;
				   
  TTree *lTree = load(lName); 
  // Declare Readers 
  fEvt       = new EvtLoader     (lTree,lName);                                             // fEvt, fEvtBr, fVertices, fVertexBr
  fMuon      = new MuonLoader    (lTree);                                                   // fMuon and fMuonBr, fN = 2 - muonArr and muonBr
  fElectron  = new ElectronLoader(lTree);                                                   // fElectrons and fElectronBr, fN = 2
  fTau       = new TauLoader     (lTree);                                                   // fTaus and fTaurBr, fN = 1
  fPhoton    = new PhotonLoader  (lTree);                                                   // fPhotons and fPhotonBr, fN = 1
  fJet4      = new JetLoader     (lTree, isData);                                                   // fJets, fJetBr => AK4PUPPI, sorted by pT
  fVJet8     = new VJetLoader    (lTree,"AK8Puppi","AddAK8Puppi","AK8CHS","AddAK8CHS",3, isData);     // fVJets, fVJetBr => AK8PUPPI
  fVJet15    = new VJetLoader    (lTree,"CA15Puppi","AddCA15Puppi","CA15CHS","AddCA15CHS",3, isData);
  if(lOption.compare("data")!=0) fGen      = new GenLoader     (lTree);                     // fGenInfo, fGenInfoBr => GenEvtInfo, fGens and fGenBr => GenParticle

  TFile *lFile = TFile::Open("Output.root","RECREATE");
  TTree *lOut  = new TTree("Events","Events");

  //Setup histograms containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);
  
    
  // Setup Tree
  fEvt      ->setupTree      (lOut); 
  fVJet8    ->setupTree      (lOut,"AK8Puppijet"); 
  fVJet8    ->setupTreeZprime(lOut,"AK8Puppijet");
  fVJet8    ->setupTreeCHS   (lOut,"AK8CHSjet");
  fVJet15   ->setupTree      (lOut,"CA15Puppijet");
  fVJet15   ->setupTreeZprime(lOut,"CA15Puppijet");
  fVJet15   ->setupTreeCHS   (lOut,"CA15CHSjet");
  fJet4     ->setupTree      (lOut,"AK4Puppijet");
  fMuon     ->setupTree      (lOut); 
  fElectron ->setupTree      (lOut); 
  fTau      ->setupTree      (lOut); 
  fPhoton   ->setupTree      (lOut); 
  if(lOption.compare("data")!=0) fGen ->setupTree (lOut,float(lXS));

  // Loop over events i0 = iEvent
  int neventstest = 0;
  for(int i0 = 0; i0 < int(lTree->GetEntriesFast()); i0++) {
  //for(int i0 = 0; i0 < int(1000); i0++){ // for testing

    // Check GenInfo
    fEvt->load(i0);
    float lWeight = 1;
    unsigned int passJson = 0;
    if(lOption.compare("data")!=0){
      fGen->load(i0);
      //lWeight = (float(lXS)*1000.*fGen->fWeight)/weight;
      lWeight = fGen->fWeight;
      passJson = 1;
    }
    else{
      if(passEvent(fEvt->fRun,fEvt->fLumi)) passJson = 1;
    }

    
    NEvents->SetBinContent(1, NEvents->GetBinContent(1)+lWeight);
    SumWeights->Fill(1.0, lWeight);

    // Primary vertex requirement
    if(!fEvt->PV()) continue;
    
    // Triggerbits
    unsigned int trigbits=1;   
    if(lOption.find("data")!=std::string::npos){
      if(fEvt ->passTrigger("HLT_AK8PFJet360_TrimMass30_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*") ||
	 fEvt ->passTrigger("HLT_PFHT800_v*") || 
	 fEvt ->passTrigger("HLT_PFHT900_v*") || 
	 fEvt ->passTrigger("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*") ||
	 fEvt ->passTrigger("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*") ||
	 fEvt ->passTrigger("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v*") ||
	 fEvt ->passTrigger("HLT_PFJet450_v*")
	 )  trigbits = trigbits | 2;  // hadronic signal region
      if( fEvt ->passTrigger("HLT_Mu45_eta2p1_v*") ||
	  fEvt ->passTrigger("HLT_Mu50_v*")) trigbits = trigbits | 4; // single muon control region
      if( fEvt ->passTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*")) trigbits = trigbits | 8; // single electron control region 
      // if(trigbits==1) continue;
    }
    fEvt      ->fillEvent(trigbits,lWeight,passJson);
    
    // Objects
    std::vector<TLorentzVector> cleaningMuons, cleaningElectrons, cleaningPhotons; 
    fMuon     ->load(i0);
    fMuon     ->selectMuons(cleaningMuons,fEvt->fMet,fEvt->fMetPhi);
    fElectron ->load(i0);
    fElectron ->selectElectrons(fEvt->fRho,fEvt->fMet,cleaningElectrons);
    fTau      ->load(i0);
    fTau      ->selectTaus(cleaningElectrons, cleaningMuons);
    fPhoton   ->load(i0);
    fPhoton   ->selectPhotons(fEvt->fRho,cleaningElectrons,cleaningPhotons);
        
    // Lepton and Photon SF    
    if(lOption.find("data")==std::string::npos){
      fEvt->fillLepSF(fMuon->looseMuons,13,fMuon->fismu0Tight,fMuon->fismu1Tight);
      fEvt->fillLepSF(fElectron->looseElectrons,11,fElectron->fisele0Tight,fElectron->fisele1Tight);
    }
    fEvt      ->triggerEff(fElectron->looseElectrons, cleaningPhotons);


    // CA15Puppi Jets
    fVJet15   ->load(i0);
    fVJet15   ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,1.5,fEvt->fRho,fEvt->fRun);
    fVJet15   ->selectVJetsCHS(cleaningElectrons,cleaningMuons,cleaningPhotons,1.5,fEvt->fRho,fEvt->fRun);
    if(fVJet15->selectedVJets.size()>0) fEvt->fselectBits =  fEvt->fselectBits | 4;
    for (int i1=0; i1<int(fVJet15->selectedVJets.size()); i1++) {
      // Match CA15 Puppi jet with CA15 CHS jet within dR = 1.5
      fVJet15 ->matchJet(fVJet15->selectedVJetsCHS,fVJet15->selectedVJets[i1],1.5,i1);
    }
      
    // AK8Puppi Jets    
    fVJet8    ->load(i0);
    fVJet8    ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt->fRho,fEvt->fRun);
    fVJet8    ->selectVJetsCHS(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt->fRho,fEvt->fRun);
    if(fVJet8->selectedVJets.size()>0) fEvt->fselectBits =  fEvt->fselectBits | 2;
    for (int i1=0; i1<int(fVJet8->selectedVJets.size()); i1++){
      // Match AK8 Puppi jet with AK8 CHS jet within dR = 0.8
      fVJet8 ->matchJet(fVJet8->selectedVJetsCHS,fVJet8->selectedVJets[i1],0.8,i1);
    }
    // Match leading AK8 Puppi jet with CA15 Puppi jet within dR = 0.4 (to get pT ratio)
    if(fVJet8->selectedVJets.size()>0) fVJet8 ->matchJet15(fVJet15->selectedVJets,fVJet8->selectedVJets[0],0.4);
    
    // computing AK8Puppi Jets pT > 500, sorted by double b-tag -> hold off for now
    //fVJet8    ->selectVJetsByDoubleBCHS(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt->fRho);
      
    // AK4Puppi Jets
    fJet4     ->load(i0); 
    fJet4     ->selectJets(cleaningElectrons,cleaningMuons,cleaningPhotons,fVJet8->selectedVJets,fEvt->fRho,fEvt->fRun);

    // Select at least one AK8 or one CA15 jet
    if(!(fEvt->fselectBits & 2) || !(fEvt->fselectBits & 4)) continue;
    if((fVJet8->selectedVJets[0].Pt() < 350) || (fVJet15->selectedVJets[0].Pt() < 350)) continue;

    // TTbar, EWK and kFactor correction
    if(lName.find("ZJets")!=std::string::npos || lName.find("DYJets")!=std::string::npos){
      fGen->findBoson(23,0);
      if(fGen->fBosonPt>0)      fEvt->computeCorr(fGen->fBosonPt,"ZJets_012j_NLO/nominal","ZJets_LO/inv_pt","EWKcorr/Z","ZJets_012j_NLO");
      if(fVJet8->selectedVJets.size()>0)  fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,23);
      if(fVJet15->selectedVJets.size()>0)  fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,23);
    }
    if(lName.find("WJets")!=std::string::npos){
      fGen->findBoson(24,1);
      if(fGen->fBosonPt>0)      fEvt->computeCorr(fGen->fBosonPt,"WJets_012j_NLO/nominal","WJets_LO/inv_pt","EWKcorr/W","WJets_012j_NLO");
      if(fVJet8->selectedVJets.size()>0)  fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,24);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,24);
    }
    if(lName.find("ZPrime")!=std::string::npos || lName.find("VectorDiJet")!=std::string::npos){
      fGen->findBoson(10031,0);
      if(fGen->fBosonPt>0)      fEvt->computeCorr(fGen->fBosonPt,"ZJets_012j_NLO/nominal","ZJets_LO/inv_pt","EWKcorr/Z","ZJets_012j_NLO");
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,10031);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,10031);
    }
    if(lName.find("DMSpin0")!=std::string::npos){
      fGen->findBoson(9900032,1);
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,9900032);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,9900032);
    }
    if(lName.find("TTJets")!=std::string::npos || lName.find("TT_")!=std::string::npos){
      float ttbarPtWeight = fGen->computeTTbarCorr();
      fEvt->fevtWeight *= ttbarPtWeight;
      fGen->fWeight *= ttbarPtWeight;
      fGen->saveTTbarType();
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,6);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,6);
    }
    if(lName.find("HToBB")!=std::string::npos || lName.find("HTobb")!=std::string::npos){
      fGen->findBoson(25,1);
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,25);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,25);
    }

    lOut->Fill();
    neventstest++;
  }
  std::cout << neventstest << std::endl;
  std::cout << lTree->GetEntriesFast() << std::endl;
  lFile->cd();
  lOut->Write();  
  NEvents->Write();
  SumWeights->Write();
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  lFile->Close();
}
