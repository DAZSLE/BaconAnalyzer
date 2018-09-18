//================================================================================================
//
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
#include <TError.h>
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

  fRangeMap = new RunLumiRangeMap();
  std::cout << "json " << lJson << std::endl;
  if(lJSON.size() > 0) fRangeMap->AddJSONFile(lJson.c_str());

  bool isData;
  if(lOption.compare("data")!=0) isData = false;
  else isData = true;
				   
  TTree *lTree = load(lName); 
  // Declare Readers 
  fEvt       = new EvtLoader     (lTree,lName);    
  fMuon      = new MuonLoader    (lTree);          
  fElectron  = new ElectronLoader(lTree);          
  fTau       = new TauLoader     (lTree);          
  fPhoton    = new PhotonLoader  (lTree);          
  fJet4      = new JetLoader     (lTree, isData);  
  fVJet8     = new VJetLoader    (lTree,"AK8Puppi","AddAK8Puppi",3, isData);  
  fVJet15    = new VJetLoader    (lTree,"CA15Puppi","AddCA15Puppi",3, isData);
  if(lOption.compare("data")!=0) fGen      = new GenLoader     (lTree);                 

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
  fVJet15   ->setupTree      (lOut,"CA15Puppijet");
  fVJet15   ->setupTreeZprime(lOut,"CA15Puppijet");
  fJet4     ->setupTree      (lOut,"AK4Puppijet");
  fMuon     ->setupTree      (lOut); 
  fElectron ->setupTree      (lOut); 
  fTau      ->setupTree      (lOut); 
  fPhoton   ->setupTree      (lOut); 
  if(lOption.compare("data")!=0) fGen ->setupTree (lOut,float(lXS));

  // Loop over events i0 = iEvent
  int neventstest = 0;
  int neventsTotal = int(lTree->GetEntriesFast());
  std::cout << maxSplit << std::endl;
  int minEventsPerJob = neventsTotal / maxSplit;
  //int leftoverEvents = neventsTotal % maxSplit;
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
      if(passEvent(fEvt->fRun,fEvt->fLumi)) { passJson = 1;}
    }

    
    NEvents->SetBinContent(1, NEvents->GetBinContent(1)+lWeight);
    SumWeights->Fill(1.0, lWeight);

    // Primary vertex requirement
    if(!fEvt->PV()) continue;
    
    // Triggerbits
    unsigned int trigbits=1;   
    if(lOption.find("data")!=std::string::npos){
      if(
	 fEvt ->passTrigger("HLT_AK8PFJet400_TrimMass30_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet420_TrimMass30_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFHT800_TrimMass50_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFHT850_TrimMass50_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFHT900_TrimMass50_v*") ||
	 fEvt ->passTrigger("HLT_PFJet500_v*") ||
	 fEvt ->passTrigger("HLT_PFHT1050_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet500_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet550_v*") ||
	 fEvt ->passTrigger("HLT_CaloJet500_NoJetID_v*") ||
	 fEvt ->passTrigger("HLT_CaloJet550_NoJetId_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p17_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p1_v*")
	 )  trigbits = trigbits | 2;  // hadronic signal region
      if( fEvt ->passTrigger("HLT_Mu50_v*") ||
	  fEvt ->passTrigger("HLT_TkMu50_v*")
	  ) trigbits = trigbits | 4; // single muon control region
      if( fEvt ->passTrigger("HLT_Ele45_WPLoose_v*") ||
	  fEvt ->passTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*")
	  ) trigbits = trigbits | 8; // single electron control region 
      if(
	 fEvt ->passTrigger("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") ||
	 fEvt ->passTrigger("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") ||
	 fEvt ->passTrigger("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") 
	)   trigbits = trigbits | 16; // AK4maxDeta + doubleBTag  
    }
    // More trigger bits

    //Add the triggers we want
    fEvt ->addTrigger("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*");
    fEvt ->addTrigger("HLT_PFHT800_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_PFHT900_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_PFHT1050_v*")  ; 
    fEvt ->addTrigger("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*");
    fEvt ->addTrigger("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*");
    fEvt ->addTrigger("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v*");
    fEvt ->addTrigger("HLT_PFJet450_v*");
    fEvt ->addTrigger("HLT_PFJet500_v*");

    fEvt ->addTrigger("HLT_IsoMu24_v*"); // W(munu)H(bb)
    fEvt ->addTrigger("HLT_IsoTkMu24_v*"); // W(munu)H(bb)
    fEvt ->addTrigger("HLT_Mu50_v*");
    fEvt ->addTrigger("HLT_TkMu50_v*");
    
    fEvt ->addTrigger("HLT_Ele27_WPTight_Gsf_v*"); // W(enu)H(bb)
    fEvt ->addTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"); //Z(ee)H(bb)
    fEvt ->addTrigger("HLT_Ele45_WPLoose_v*");
    fEvt ->addTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"); 

    fEvt ->addTrigger("HLT_PFMET110_PFMHT110_IDTight_v*"); // Z(nunu)H(bb) 
    fEvt ->addTrigger("HLT_PFMET120_PFMHT120_IDTight_v*"); // Z(nunu)H(bb)
    fEvt ->addTrigger("HLT_PFMET170_NoiseCleaned_v*"); // Z(nunu)H(bb)
    fEvt ->addTrigger("HLT_PFMET170_HBHECleaned_v*"); // Z(nunu)H(bb)
    fEvt ->addTrigger("HLT_PFMET170_HBHE_BeamHaloCleaned_v*"); // Z(nunu)H(bb)

    //ggH(bb)
    fEvt ->addTrigger("HLT_AK8PFJet360_TrimMass30_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_AK8PFJet380_TrimMass30_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_AK8PFJet400_TrimMass30_v*");
    fEvt ->addTrigger("HLT_AK8PFJet420_TrimMass30_v*");
    fEvt ->addTrigger("HLT_AK8PFHT800_TrimMass50_v*");
    fEvt ->addTrigger("HLT_AK8PFHT850_TrimMass50_v*");
    fEvt ->addTrigger("HLT_AK8PFHT900_TrimMass50_v*");
    fEvt ->addTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p17_v*");
    fEvt ->addTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p1_v*");
    fEvt ->addTrigger("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
    fEvt ->addTrigger("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
    fEvt ->addTrigger("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
    fEvt ->addTrigger("HLT_CaloJet500_NoJetID_v*");
    fEvt ->addTrigger("HLT_CaloJet550_NoJetId_v*");
    fEvt ->addTrigger("HLT_AK8PFJet500_v*");
    fEvt ->addTrigger("HLT_AK8PFJet550_v*");

    fEvt      ->fillEvent(trigbits,lWeight,passJson);
    
    // Objects
    gErrorIgnoreLevel=kError;
    std::vector<TLorentzVector> cleaningMuons, cleaningElectrons, cleaningPhotons; 
    fMuon     ->load(i0);
    fMuon     ->selectMuons(cleaningMuons,fEvt->fMet,fEvt->fMetPhi);
    fElectron ->load(i0);
    fElectron ->selectElectrons(fEvt->fRho,fEvt->fMet,cleaningElectrons);
    fTau      ->load(i0);
    fTau      ->selectTaus(cleaningElectrons, cleaningMuons);
    fPhoton   ->load(i0);
    fPhoton   ->selectPhotons(fEvt->fRho,cleaningElectrons,cleaningPhotons);
        
    // CA15Puppi Jets
    fVJet15   ->load(i0);
    fVJet15   ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,1.5,fEvt->fRho,fEvt->fRun);
    if(fVJet15->selectedVJets.size()>0) fEvt->fselectBits =  fEvt->fselectBits | 4;
      
    // AK8Puppi Jets    
    fVJet8    ->load(i0);
    fVJet8    ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt->fRho,fEvt->fRun);
    if(fVJet8->selectedVJets.size()>0) fEvt->fselectBits =  fEvt->fselectBits | 2;

    // Match leading AK8 Puppi jet with CA15 Puppi jet within dR = 0.4 (to get pT ratio)
    if(fVJet8->selectedVJets.size()>0) fVJet8 ->matchJet15(fVJet15->selectedVJets,fVJet8->selectedVJets[0],0.4);
    
    // AK4Puppi Jets
    fJet4     ->load(i0); 
    fJet4     ->selectJets(cleaningElectrons,cleaningMuons,cleaningPhotons,fVJet8->selectedVJets,fEvt->fRho,fEvt->fRun);

    // Select at least one AK8 or one CA15 jet
    //if(!(fEvt->fselectBits & 2) || !(fEvt->fselectBits & 4)) continue;
    //if((fVJet8->selectedVJets[0].Pt() < 200) || (fVJet15->selectedVJets[0].Pt() < 200)) continue;

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
      // 55 for newer, public samples (to be consistent with DM samples) 
      //fGen->findBoson(10031,0);
      fGen->findBoson(55,1);
      if(fGen->fBosonPt>0)      fEvt->computeCorr(fGen->fBosonPt,"ZJets_012j_NLO/nominal","ZJets_LO/inv_pt","EWKcorr/Z","ZJets_012j_NLO");
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,10031);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,10031);
    }
    if(lName.find("Spin0")!=std::string::npos){
      // 9900032 for older, private samples
      // 55 for newer, public samples (to be consistent with DM samples)
      fGen->findBoson(55,1);
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,55);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,55);
    }
    if(lName.find("TTJets")!=std::string::npos || lName.find("TT_")!=std::string::npos || lName.find("TTTo")!=std::string::npos){
      float ttbarPtWeight = fGen->computeTTbarCorr();
      fEvt->fevtWeight *= ttbarPtWeight;
      fGen->fWeight *= ttbarPtWeight;
      fGen->saveTTbarType();
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,624);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,624);
    }    
    if(lName.find("ST_")!=std::string::npos){
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,624);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,624);
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
