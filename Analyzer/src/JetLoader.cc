#include "../include/JetLoader.hh"
#include <cmath>
#include <iostream> 
#include <sstream>

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree) { 
  fJets  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("AK4Puppi",       &fJets);
  fJetBr = iTree->GetBranch("AK4Puppi");

  fJetsCHS  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("AK4CHS",       &fJetsCHS);
  fJetBrCHS = iTree->GetBranch("AK4CHS");

  fN = 4;
  fNV = 3; // max number of V jets to consider for dR anti-matching
  fNVars = 3; // pt, eta, phi
  fNOtherVars = 5; // Mass, b-tag, qgid, dR, dPhi 
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
  delete fJetsCHS;
  delete fJetBrCHS;
}
void JetLoader::reset() { 
  fNJetsPt30           = 0;
  fNFwdPt30            = 0;
  fNBTagsLPt30         = 0;
  fNBTagsMPt30         = 0;
  fNBTagsTPt30         = 0;
  fLooseJets.clear();
  fGoodJets.clear();
  selectedJets8.clear();
  selectedJets15.clear();  
  for(int i0 = 0; i0 < int(fNJetsPt30dR08.size()); i0++) {
    fNJetsPt30dR08[i0] = -999;
    fNBTagsLPt50dR08[i0] = -999;
    fNBTagsMPt50dR08[i0] = -999;
    fNBTagsTPt50dR08[i0] = -999;
    fNBTagsLPt100dR08[i0] = -999;
    fNBTagsMPt100dR08[i0] = -999;
    fNBTagsTPt100dR08[i0] = -999;
    fNBTagsLPt150dR08[i0] = -999;
    fNBTagsMPt150dR08[i0] = -999;
    fNBTagsTPt150dR08[i0] = -999;
  }  
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
}
void JetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();
  fTree = iTree;
  std::stringstream pSNPt30,pSfwdPt30,pSbPt30,pSbLPt30,pSbMPt30,pSbTPt30;
  pSNPt30     << "n" << iJetLabel << "sPt30";
  pSfwdPt30   << "n" << iJetLabel << "sfwdPt30";
  pSbPt30     << "n" << iJetLabel << "sbtagPt30";
  pSbLPt30    << "n" << iJetLabel << "sbtagLPt30";
  pSbMPt30    << "n" << iJetLabel << "sbtagMPt30";
  pSbTPt30    << "n" << iJetLabel << "sbtagTPt30";


  fTree->Branch(pSNPt30.str().c_str()           ,&fNJetsPt30           ,(pSNPt30.str()+"/I").c_str());  // jet multiplicity
  fTree->Branch(pSfwdPt30.str().c_str()         ,&fNFwdPt30            ,(pSfwdPt30.str()+"/I").c_str());
  fTree->Branch(pSbLPt30.str().c_str()          ,&fNBTagsLPt30         ,(pSbLPt30.str()+"/I").c_str()); // b tags
  fTree->Branch(pSbMPt30.str().c_str()          ,&fNBTagsMPt30         ,(pSbMPt30.str()+"/I").c_str());
  fTree->Branch(pSbTPt30.str().c_str()          ,&fNBTagsTPt30         ,(pSbTPt30.str()+"/I").c_str());


  fNJetsPt30dR08.clear();
  fNBTagsLPt50dR08.clear();
  fNBTagsMPt50dR08.clear();
  fNBTagsTPt50dR08.clear();
  fNBTagsLPt100dR08.clear();
  fNBTagsMPt100dR08.clear();
  fNBTagsTPt100dR08.clear();
  fNBTagsLPt150dR08.clear();
  fNBTagsMPt150dR08.clear();
  fNBTagsTPt150dR08.clear();
  for(int i0 = 0; i0 < fNV; i0++) {
    fNJetsPt30dR08.push_back(-999);
    fNBTagsLPt50dR08.push_back(-999);
    fNBTagsMPt50dR08.push_back(-999);
    fNBTagsTPt50dR08.push_back(-999);
    fNBTagsLPt100dR08.push_back(-999);
    fNBTagsMPt100dR08.push_back(-999);
    fNBTagsTPt100dR08.push_back(-999);
    fNBTagsLPt150dR08.push_back(-999);
    fNBTagsMPt150dR08.push_back(-999);
    fNBTagsTPt150dR08.push_back(-999);
  }  
  for(int i0 = 0; i0 < fNV; i0++) {
    std::stringstream pSNPt30dR08,pSbLPt50dR08,pSbMPt50dR08,pSbTPt50dR08,pSbLPt100dR08,pSbMPt100dR08,pSbTPt100dR08,pSbLPt150dR08,pSbMPt150dR08,pSbTPt150dR08; 
    pSNPt30dR08 << "n" << iJetLabel << "sPt30dR08_" << i0;
    pSbLPt50dR08<< "n" << iJetLabel << "sLPt50dR08_" << i0;
    pSbMPt50dR08<< "n" << iJetLabel << "sMPt50dR08_" << i0;
    pSbTPt50dR08<< "n" << iJetLabel << "sTPt50dR08_" << i0;
    pSbLPt100dR08<< "n" << iJetLabel << "sLPt100dR08_" << i0;
    pSbMPt100dR08<< "n" << iJetLabel << "sMPt100dR08_" << i0;
    pSbTPt100dR08<< "n" << iJetLabel << "sTPt100dR08_ " << i0;
    pSbLPt150dR08<< "n" << iJetLabel << "sLPt150dR08_" << i0;
    pSbMPt150dR08<< "n" << iJetLabel << "sMPt150dR08_" << i0;
    pSbTPt150dR08<< "n" << iJetLabel << "sTPt150dR08_" << i0;
    fTree->Branch(pSNPt30dR08.str().c_str()       ,&fNJetsPt30dR08[i0]       ,(pSNPt30dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt50dR08.str().c_str()      ,&fNBTagsLPt50dR08[i0]     ,(pSbLPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt50dR08.str().c_str()      ,&fNBTagsMPt50dR08[i0]     ,(pSbMPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt50dR08.str().c_str()      ,&fNBTagsTPt50dR08[i0]     ,(pSbTPt50dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt100dR08.str().c_str()      ,&fNBTagsLPt100dR08[i0]     ,(pSbLPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt100dR08.str().c_str()      ,&fNBTagsMPt100dR08[i0]     ,(pSbMPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt100dR08.str().c_str()      ,&fNBTagsTPt100dR08[i0]     ,(pSbTPt100dR08.str()+"/I").c_str());
    fTree->Branch(pSbLPt150dR08.str().c_str()      ,&fNBTagsLPt150dR08[i0]     ,(pSbLPt150dR08.str()+"/I").c_str());
    fTree->Branch(pSbMPt150dR08.str().c_str()      ,&fNBTagsMPt150dR08[i0]     ,(pSbMPt150dR08.str()+"/I").c_str());
    fTree->Branch(pSbTPt150dR08.str().c_str()      ,&fNBTagsTPt150dR08[i0]     ,(pSbTPt150dR08.str()+"/I").c_str());
  }

  for(int i0 = 0; i0 < fN*(10)+4; i0++) {double pVar = 0; fVars.push_back(pVar);}           
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars);                                            // from MonoXUtils.cc => fN=4 j*_pt,j*_eta,j*_phi for j1,j2,j3,j4 (3*4=12)
  addOthers  (iJetLabel.c_str(),iTree,fN,fVars);                                            // Mass, b-tag, qgid, dR, dPhi for j1,j2,j3,j4 (5*4=20)
}
void JetLoader::load(int iEvent) { 
  fJets   ->Clear();
  fJetBr ->GetEntry(iEvent);
}
void JetLoader::selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets){
  reset(); 
  int lCountPt30 = 0, lNFwdPt30 = 0, lNBTagLPt30 = 0,lNBTagMPt30 = 0, lNBTagTPt30 = 0;
  
  for (int i1 = 0; i1 < int(iVJets.size()); i1++) {
    fNJetsPt30dR08[i1] = 0;
    fNBTagsLPt50dR08[i1] = 0;
    fNBTagsMPt50dR08[i1] = 0;
    fNBTagsTPt50dR08[i1] = 0;
    fNBTagsLPt100dR08[i1] = 0;
    fNBTagsMPt100dR08[i1] = 0;
    fNBTagsTPt100dR08[i1] = 0;
    fNBTagsLPt150dR08[i1] = 0;
    fNBTagsMPt150dR08[i1] = 0;
    fNBTagsTPt150dR08[i1] = 0;
  }
    
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet->eta,pJet->phi,0.4,iElectrons))                      continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iMuons))                          continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iPhotons))                        continue;
    if(pJet->pt        <=  30)                                            continue;
    if(fabs(pJet->eta) > 2.5 && fabs(pJet->eta) < 4.5) lNFwdPt30++;
    if(fabs(pJet->eta) >= 2.5)                                            continue;
    if(!passJetLooseSel(pJet))                                            continue;
    lCountPt30++;
    addJet(pJet,fLooseJets);

    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    fGoodJets.push_back(pJet);
    
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVL){
      lNBTagLPt30++;
    }    
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVM){ 
      lNBTagMPt30++;
    }
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVT){ 
      lNBTagTPt30++;
    }
    
    // jet and b-tag multiplicity
    for (int i1 = 0; i1 < int(iVJets.size()); i1++) {
      
      if(iVJets[i1].Pt()>350 && vPJet.DeltaR(iVJets[i1])>0.8) {
	fNJetsPt30dR08[i1]++;
	if(pJet->pt>50 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVL) fNBTagsLPt50dR08[i1]++;
	if(pJet->pt>100 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVL) fNBTagsLPt100dR08[i1]++;
	if(pJet->pt>150 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVL) fNBTagsLPt150dR08[i1]++;
	
	if(pJet->pt>50 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVM) fNBTagsMPt50dR08[i1]++;
	if(pJet->pt>100 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVM) fNBTagsMPt100dR08[i1]++;
	if(pJet->pt>150 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVM) fNBTagsMPt150dR08[i1]++;
	
	if(pJet->pt>50 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVT) fNBTagsTPt50dR08[i1]++;
	if(pJet->pt>100 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVT) fNBTagsTPt100dR08[i1]++;
	if(pJet->pt>150 && fabs(pJet->eta) < 2.5 && pJet->csv > CSVT) fNBTagsTPt150dR08[i1]++;
      }
    }
  }
  addVJet(fLooseJets,selectedJets);
  fNJetsPt30           = lCountPt30;
  fNFwdPt30            = lNFwdPt30;
  fNBTagsLPt30         = lNBTagLPt30;
  fNBTagsMPt30         = lNBTagMPt30;
  fNBTagsTPt30         = lNBTagTPt30;

  fillJet(fN,fLooseJets,fVars);
  fillOthers(fN,fLooseJets,fVars,iVJets);
}
void JetLoader::addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals) { 
  for(int i0 = 0; i0 < iN; i0++) { 
    int lBase = iN*fNVars+i0*fNOtherVars;
    std::stringstream pSMass,pSCSV,pSQGID,pSdR,pSdP;
    pSMass  << iHeader << i0 << "_mass";
    pSCSV   << iHeader << i0 << "_csv";
    pSQGID  << iHeader << i0 << "_qgid";
    pSdR    << iHeader << i0 << "_dR08";
    pSdP    << iHeader << i0 << "_dPhi08";
    iTree->Branch(pSMass .str().c_str(),&iVals[lBase+0],(pSMass .str()+"/D").c_str());
    iTree->Branch(pSCSV .str().c_str() ,&iVals[lBase+1],(pSCSV  .str()+"/D").c_str());
    iTree->Branch(pSQGID.str().c_str() ,&iVals[lBase+2],(pSQGID .str()+"/D").c_str());
    iTree->Branch(pSdR  .str().c_str() ,&iVals[lBase+3],(pSdR   .str()+"/D").c_str());
    iTree->Branch(pSdP  .str().c_str() ,&iVals[lBase+4],(pSdP   .str()+"/D").c_str());
  }
}
void JetLoader::fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, std::vector<TLorentzVector> iVJets){ 
  int lBase = fNVars*fN;
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) {
    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(iObjects[i0]->pt, iObjects[i0]->eta, iObjects[i0]->phi, iObjects[i0]->mass); 
    iVals[lBase+i0*fNOtherVars+0] = iObjects[i0]->mass;
    iVals[lBase+i0*fNOtherVars+1] = iObjects[i0]->csv;
    iVals[lBase+i0*fNOtherVars+2] = iObjects[i0]->qgid;
    if(iVJets.size()>0) {
      iVals[lBase+i0*fNOtherVars+3] = vPJet.DeltaR(iVJets[0]);
      iVals[lBase+i0*fNOtherVars+4] = vPJet.DeltaPhi(iVJets[0]);
    }
  }
}
