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
#include "../include/PerJetLoader.hh"
#include "../include/RunLumiRangeMap.h"

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
PerJetLoader    *fVJet8      = 0;
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
  const std::string lLabel       = argv[3];
  const std::string lJSON        = argv[4];
  const std::string lOutput      = argv[5];
  const int         iSplit       = atoi(argv[6]);
  const int         maxSplit     = atoi(argv[7]);

  fRangeMap = new RunLumiRangeMap();
  if(lJSON.size() > 0) fRangeMap->AddJSONFile(lJSON.c_str());
  std::cout << "args " << argv[1] << " " << argv[2] << " "<< argv[3] << " " << argv[4] << " " << argv[5] << std::endl;

  bool isData;
  if(lOption.compare("data")!=0) isData = false;
  else isData = true;
				   
  TTree *lTree = load(lName); 
  // Declare Readers 
  fEvt       = new EvtLoader     (lTree,lName);                                             // fEvt, fEvtBr, fVertices, fVertexBr
  fMuon      = new MuonLoader    (lTree);                                                   // fMuon and fMuonBr, fN = 2 - muonArr and muonBr
  fElectron  = new ElectronLoader(lTree);                                                   // fElectrons and fElectronBr, fN = 2
  fPhoton    = new PhotonLoader  (lTree);                                                   // fPhotons and fPhotonBr, fN = 1
  fVJet8     = new PerJetLoader  (lTree,"AK8Puppi","AddAK8Puppi","AK8CHS","AddAK8CHS",3, isData, lLabel);     // fVJets, fVJetBr => AK8PUPPI

  TFile *lFile = TFile::Open(lOutput.c_str(),"RECREATE");
  TTree *lOut  = new TTree("Events","Events");

  //Setup histograms containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
    
  // Setup Tree
  fEvt      ->setupTreeQbert (lOut); 
  fVJet8    ->setupTreeQbert (lOut,"fj"); 
  fMuon     ->setupTreeQbert (lOut); 
  fElectron ->setupTreeQbert (lOut); 

  // Loop over events i0 = iEvent
  int neventstest = 0;
  int neventsTotal = int(lTree->GetEntriesFast());
  int minEventsPerJob = neventsTotal / maxSplit;
  int minEvent = iSplit * minEventsPerJob;
  int maxEvent = (iSplit+1) * minEventsPerJob;
  if (iSplit + 1 == maxSplit) maxEvent = neventsTotal;
  for(int i0 = minEvent; i0 < maxEvent; i0++) {
    if (i0%1000 == 0) std::cout << i0 << " events processed " << std::endl;
    fEvt->load(i0);
    float lWeight = 1;
    unsigned int passJson = 0;
    unsigned int trigbits=1;   
    if(lOption.compare("data")!=0){
      // we don't care about weights for now
      // fGen->load(i0);
      // //lWeight = (float(lXS)*1000.*fGen->fWeight)/weight;
      // lWeight = fGen->fWeight;
      passJson = 1;
    }
    else{
      if(passEvent(fEvt->fRun,fEvt->fLumi)) passJson = 1;
    }

    // NEvts
    NEvents->SetBinContent(1, NEvents->GetBinContent(1)+lWeight);
    SumWeights->Fill(1.0, lWeight);

    // Primary vertex requirement
    if(!fEvt->PV()) {
      continue;
    }
    fEvt      ->fillEvent(trigbits,lWeight,passJson);
    
    // Objects
    gErrorIgnoreLevel=kError;
    std::vector<TLorentzVector> cleaningMuons, cleaningElectrons, cleaningPhotons; 
    std::vector<TLorentzVector> Muons, Electrons;
    fMuon     ->load(i0);
    fMuon     ->selectMuons(Muons,fEvt->fMet,fEvt->fMetPhi);
    fElectron ->load(i0);
    fElectron ->selectElectrons(fEvt->fRho,fEvt->fMet,Electrons);
    fPhoton   ->load(i0);
    fPhoton   ->selectPhotons(fEvt->fRho,Electrons,cleaningPhotons);

    // Do not clean leptons from jet - but only photons?
    fVJet8    ->load(i0);
    fVJet8    ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt->fRho,fEvt->fMetPhi,fEvt->fRun);

    neventstest++;
    //break;
  }
  std::cout << neventstest << std::endl;
  std::cout << lTree->GetEntriesFast() << std::endl;
  lFile->cd();
  lOut->Write();  
  NEvents->Write();
  SumWeights->Write();
  lFile->Close();
}
