//================================================================================================
//
// Perform preselection for W/Qbert(->qq)+jets events and produce bacon bits 
//
// Input arguments
//   argv[1] => lName = input bacon file name
//   argv[2] => lOption = dataset type: "mc", "data"
//   argv[3] => lJSON = JSON file for run-lumi filtering of data, specify "none" for MC or no filtering
//   argv[4] => lXS = cross section (pb), ignored for data 
//   argv[5] => weight = total weight, ignored for data
//________________________________________________________________________________________________

#include "../include/EvtLoader.hh"
#include "../include/ElectronLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/PhotonLoader.hh"
#include "../include/VJetLoader.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TError.h>
#include <string>
#include <iostream>

// Object Processors
EvtLoader       *fEvt        = 0; 
MuonLoader      *fMuon       = 0; 
ElectronLoader  *fElectron   = 0; 
PhotonLoader    *fPhoton     = 0; 
VJetLoader      *fVJet8      = 0;

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
  const int      iSplit          = atoi(argv[6]);
  const int      maxSplit        = atoi(argv[7]);

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
  fPhoton    = new PhotonLoader  (lTree);                                                   // fPhotons and fPhotonBr, fN = 1
  fVJet8     = new VJetLoader    (lTree,"AK8Puppi","AddAK8Puppi","AK8CHS","AddAK8CHS",3, isData);     // fVJets, fVJetBr => AK8PUPPI
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
  fVJet8    ->setupTreeQbert(lOut,"AK8Puppijet");
  fMuon     ->setupTree      (lOut); 
  fElectron ->setupTree      (lOut); 
  fPhoton   ->setupTree      (lOut); 
  if(lOption.compare("data")!=0) fGen ->setupTree (lOut,float(lXS));

  // Loop over events i0 = iEvent
  int neventstest = 0;
  int neventsTotal = int(lTree->GetEntriesFast());
  int minEventsPerJob = neventsTotal / maxSplit;
  int leftoverEvents = neventsTotal % maxSplit;
  int minEvent = iSplit * minEventsPerJob;
  int maxEvent = (iSplit+1) * minEventsPerJob;
  if (iSplit + 1 == maxSplit) maxEvent = neventsTotal;
  std::cout << neventsTotal << " total events" << std::endl;
  std::cout << iSplit << " iSplit " << std::endl;
  std::cout << maxSplit << " maxSplit " << std::endl;
  std::cout << minEvent << " min event" << std::endl;
  std::cout << maxEvent << " max event" << std::endl;  
  for(int i0 = minEvent; i0 < maxEvent; i0++) {
    //for(int i0 = 0; i0 < int(10000); i0++){ // for testing
    if (i0%1000 == 0) std::cout << i0 << " events processed " << std::endl;
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
      if( fEvt ->passTrigger("HLT_Mu50_v*") ||
	  fEvt ->passTrigger("HLT_TkMu50_v*")
	  ) trigbits = trigbits | 4; // single muon control region
      if( fEvt ->passTrigger("HLT_Ele45_WPLoose_v*") ||
	  fEvt ->passTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*")
	  ) trigbits = trigbits | 8; // single electron control region 
      // if(trigbits==1) continue;
    }
    fEvt      ->fillEvent(trigbits,lWeight,passJson);
    
    // Objects
    gErrorIgnoreLevel=kError;
    std::vector<TLorentzVector> cleaningMuons, cleaningElectrons, cleaningPhotons; 
    fMuon     ->load(i0);
    fMuon     ->selectMuons(cleaningMuons,fEvt->fMet,fEvt->fMetPhi);
    fElectron ->load(i0);
    fElectron ->selectElectrons(fEvt->fRho,fEvt->fMet,cleaningElectrons);
    fPhoton   ->load(i0);
    fPhoton   ->selectPhotons(fEvt->fRho,cleaningElectrons,cleaningPhotons);


    fVJet8    ->load(i0);
    fVJet8    ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt->fRho,fEvt->fRun);
    if(fVJet8->selectedVJets.size()>0) fEvt->fselectBits =  fEvt->fselectBits | 2;


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
