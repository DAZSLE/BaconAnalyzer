//================================================================================================
//
// Skim
//
//________________________________________________________________________________________________
//

//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TChain.h>
#include <TList.h>
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 
#include <TKey.h> 
#include <TApplication.h>
#include "../include/EvtLoader.hh"

//#endif

using namespace std;

//=== MAIN MACRO ================================================================================================= 

int main(int argc, char **argv) {

  std::string inputList(argv[1]);
  std::string outputTag(argv[2]);

  // input
  ifstream filein(inputList.c_str());
  std::string curFilename;
  std::vector<std::string> inputLines;
  while(getline(filein, curFilename)){
    //std::cout << "reading " << curFilename << std::endl;
    if(curFilename.at(0) != '#') inputLines.push_back(curFilename); //'#' denotes a comment                                                                                                             
    else std::cout << "(Skipping commented line in input)" << std::endl;
  }

  //loop over all inputfiles
  map<pair<uint,uint>, bool > processedRunEvents;

  for(int i0 = 0; i0 < int(inputLines.size()); i0++) {
    std::string inputfile = inputLines.at(i0);
    std::cout << " Opening " << inputfile << " " << i0 << std::endl;

    // create output file
    std::string newfile = inputfile;
    newfile.replace(inputfile.length()-5,5,"_");
    //std::cout << "new file " << newfile << std::endl;
    std::stringstream outputfile; outputfile << newfile << outputTag << ".root";
    TFile *outputFile = new TFile(outputfile.str().c_str(), "RECREATE");

    // loop over all trees in input file
    TFile *inputFile = TFile::Open(inputfile.c_str(), "READ");
    assert(inputFile);
    inputFile->cd();
    inputFile->Purge(); //purge unwanted TTree cycles in file
    TIter nextkey(inputFile->GetListOfKeys());
    TKey *key;
    while((key = (TKey*)nextkey())){
      std::string className = key->GetClassName();
      //std::cout << "Getting key from file.  Class type: " << className << std::endl;
      if(className.compare("TTree") != 0){
	//cout << "Skipping key (not a TTree)" << endl;
	continue;
      }
      
      TTree *inputTree = (TTree*)key->ReadObj();
      //std::cout << "Processing tree " << inputTree->GetName() <<std::endl;

      const std::string lName = inputTree->GetName();
      if(lName.find("Events")==std::string::npos) continue;

      // create new tree
      outputFile->cd();
      //std::stringstream treename; treename << "Events_" << i0;
      TTree *outputTree = inputTree->CloneTree(0);
      //outputTree->SetName(treename.str().c_str());
      //std::cout << "Events in the ntuple: " << inputTree->GetEntries() << endl;

      uint run;
      uint event;
      /*
      inputTree->SetBranchAddress("runNum", &run);
      inputTree->SetBranchAddress("evtNum", &event);
      */

      TEventInfo   *fEvt;
      fEvt = new TEventInfo();
      inputTree->SetBranchAddress("Info", &fEvt);
      TBranch      *fEvtBr;
      fEvtBr = inputTree->GetBranch("Info");

      //map<pair<uint,uint>, bool > processedRunEvents;
      for (int n=0;n<inputTree->GetEntries();n++) { 
	//if (n%10000==0) std::cout << "Processed Event " << n << "\n";
	inputTree->GetEntry(n);
	fEvtBr->GetEntry(n);
	run = fEvt->runNum;
	event = fEvt->evtNum;
	if(processedRunEvents.find(make_pair(run, event)) == processedRunEvents.end()){ 
          //not duplicate
          processedRunEvents[make_pair(run, event)] = true;
          outputTree->Fill(); 
	}
	else{
	  std::cout << "Duplicate event " << run << " " << event <<  "  file " << inputfile <<  endl;
	}
      }
      
      std::cout << "Number of Input Events: " << inputTree->GetEntries() << "\n";
      std::cout << "Number of Output Events: " << outputTree->GetEntries() << "\n";
      
      //save
      //std::cout << "adding " << outputTree->GetName() << " " << outputTree->GetEntries() << std::endl;
      //outputTree->SetName(treename.str().c_str());
      //outputList->Add(outputTree);
      //outputList->Print();
      outputTree->Write();
      inputFile->cd();
    }
	
    inputFile->Close();
    outputFile->Close();
    delete outputFile;
    //delete inputFile;
  }

  //gROOT->cd();

  //outputList->Print();
  //TTree *newtree = TTree::MergeTrees(outputList); 
  //outputFile->cd();
  //newtree->SetName("Events"); 
  //newtree->Write();


  //std::cout << "Number of Total Output Events: " << newtree->GetEntries() << "\n";
  //std::cout << "Closing output file." << endl;
  //outputFile->Close();
  //delete outputList;
  //delete newtree;
  //delete outputFile;
  //delete outputList;
  //gApplication->Terminate();
  return 0;
}

