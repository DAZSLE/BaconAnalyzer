//================================================================================================
//
//
// PU per file re-weighting
//________________________________________________________________________________________________
#include "../include/GenLoader.hh"
#include "../include/EvtLoader.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TError.h>
#include <string>
#include <iostream>

// Object Processors
GenLoader       *fGen        = 0; 
EvtLoader       *fEvt        = 0; 

TH1F *fHist                  = 0;

// Load tree and return infile
TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  fHist        = (TH1F* ) lFile->FindObjectAny("TotalEvents");
  return lTree;
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

  std::string lJson="${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/data/";
  lJson.append(lJSON);
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = lJson.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
    lJson.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  TTree *lTree = load(lName); 

  fEvt       = new EvtLoader     (lTree,lName);    
  if(lOption.compare("data")!=0) fGen      = new GenLoader     (lTree);    

  TFile *lFile = TFile::Open("Output.root","RECREATE");

  // Setup histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);

  // Setup PU histograms
  TH1F *Pu = new TH1F("Pu", "Pu", 100, 0, 100);

  // Loop over events
  int neventstest = 0;
  int neventsTotal = int(lTree->GetEntriesFast());
  std::cout << maxSplit << std::endl;
  int minEventsPerJob = neventsTotal / maxSplit;
  int minEvent = iSplit * minEventsPerJob;
  int maxEvent = (iSplit+1) * minEventsPerJob;
  if (iSplit + 1 == maxSplit) maxEvent = neventsTotal;
  std::cout << neventsTotal << " total events" << std::endl;
  std::cout << iSplit << " iSplit " << std::endl;
  std::cout << maxSplit << " maxSplit " << std::endl;
  std::cout << minEvent << " min event" << std::endl;
  std::cout << maxEvent << " max event" << std::endl;  
  for(int i0 = minEvent; i0 < maxEvent; i0++) {
    if (i0%1000 == 0) std::cout << i0 << " events processed " << std::endl;

    fEvt->load(i0);
    float lWeight = 1;
    if(lOption.compare("data")!=0){
      fGen->load(i0);
      lWeight = fGen->fWeight;
    }
    
    NEvents->SetBinContent(1, NEvents->GetBinContent(1)+lWeight);
    SumWeights->Fill(1.0, lWeight);
    Pu->Fill(fEvt->fPu); // no weight?
    
    neventstest++;
  }
  std::cout << neventstest << std::endl;
  std::cout << lTree->GetEntriesFast() << std::endl;
  lFile->cd();
  NEvents->Write();
  SumWeights->Write();
  Pu->Write();
  lFile->Close();
}
