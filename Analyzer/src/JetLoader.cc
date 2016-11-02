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
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
  delete fJetsCHS;
  delete fJetBrCHS;
}
void JetLoader::reset() { 
  fNJets           = 0;
  fNFwd            = 0;
  fNBTagsL         = 0;
  fNBTagsM         = 0;
  fNBTagsT         = 0;
  fLooseJets.clear();
  fGoodJets.clear();
  selectedJets8.clear();
  selectedJets15.clear();
  for(unsigned int i0 = 0; i0 < fVars.size(); i0++) fVars[i0] = 0;
  fNJetsdR08       = 0;
  fNBTagsLdR08     = 0;
  fNBTagsMdR08     = 0;
  fNBTagsTdR08     = 0;
}
void JetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();
  fTree = iTree;
  std::stringstream pSN,pSfwd,pSb,pSbL,pSbM,pSbT,pSNdR08,pSbLdR08,pSbMdR08,pSbTdR08;   
  pSN     << "n" << iJetLabel << "s";
  pSfwd   << "n" << iJetLabel << "sfwd";
  pSb     << "n" << iJetLabel << "sbtag";
  pSbL    << "n" << iJetLabel << "sbtagL";
  pSbM    << "n" << iJetLabel << "sbtagM";
  pSbT    << "n" << iJetLabel << "sbtagT";
  pSNdR08 << "n" << iJetLabel << "sdR08";
  pSbLdR08<< "n" << iJetLabel << "sLdR08";
  pSbMdR08<< "n" << iJetLabel << "sMdR08";
  pSbTdR08<< "n" << iJetLabel << "sTdR08";

  fTree->Branch(pSN.str().c_str()           ,&fNJets           ,(pSN.str()+"/I").c_str());  // jet multiplicity
  fTree->Branch(pSfwd.str().c_str()         ,&fNFwd            ,(pSfwd.str()+"/F").c_str());
  fTree->Branch(pSbL.str().c_str()          ,&fNBTagsL         ,(pSbL.str()+"/I").c_str()); // b tags
  fTree->Branch(pSbM.str().c_str()          ,&fNBTagsM         ,(pSbM.str()+"/I").c_str());
  fTree->Branch(pSbT.str().c_str()          ,&fNBTagsT         ,(pSbT.str()+"/I").c_str());
  fTree->Branch(pSNdR08.str().c_str()       ,&fNJetsdR08       ,(pSNdR08.str()+"/I").c_str());
  fTree->Branch(pSbLdR08.str().c_str()      ,&fNBTagsLdR08     ,(pSbLdR08.str()+"/I").c_str());
  fTree->Branch(pSbMdR08.str().c_str()      ,&fNBTagsMdR08     ,(pSbMdR08.str()+"/I").c_str());
  fTree->Branch(pSbTdR08.str().c_str()      ,&fNBTagsTdR08     ,(pSbTdR08.str()+"/I").c_str());

  for(int i0 = 0; i0 < fN*(10)+4; i0++) {double pVar = 0; fVars.push_back(pVar);}           
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars);                                            // from MonoXUtils.cc => fN =4 j*_pt,j*_eta,j*_phi for j1,j2,j3,j4 (3*4=12)
  addOthers  (iJetLabel.c_str(),iTree,fN,fVars);                                            // Mass + b-tag +qgid for j1,j2,j3,j4 (3*4=12)
}
void JetLoader::load(int iEvent) { 
  fJets   ->Clear();
  fJetBr ->GetEntry(iEvent);
}
void JetLoader::selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets){
  reset(); 
  int lCount = 0, lNFwd = 0, lNBTagL = 0,lNBTagM = 0,lNBTagT = 0, lCountdR08 = 0,lNBTagLdR08 = 0,lNBTagMdR08 = 0,lNBTagTdR08 = 0;
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet->eta,pJet->phi,0.4,iElectrons))                      continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iMuons))                          continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iPhotons))                        continue;
    if(pJet->pt        <=  30)                                            continue;
    if(fabs(pJet->eta) > 3.0) lNFwd++;
    if(fabs(pJet->eta) >= 4.5)                                            continue;
    if(!passJetLooseSel(pJet))                                            continue;
    lCount++;
    addJet(pJet,fLooseJets);

    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    fGoodJets.push_back(pJet);

    // jet and b-tag multiplicity
    if(iVJets.size()>0 && iVJets[0].Pt()>0 && vPJet.DeltaR(iVJets[0])>0.8) lCountdR08++;
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVL){
      lNBTagL++;
      if(iVJets.size()>0) {if(pJet->pt>0 && vPJet.DeltaR(iVJets[0])>0.8) lNBTagLdR08++;}
    }
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVM){ 
      lNBTagM++;
      if(iVJets.size()>0) {if(pJet->pt>0 && vPJet.DeltaR(iVJets[0])>0.8) lNBTagMdR08++;}
    }
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVT){
      lNBTagT++;
      if(iVJets.size()>0) {if(pJet->pt>0 && vPJet.DeltaR(iVJets[0])>0.8) lNBTagTdR08++;}
    }
  }
  addVJet(fLooseJets,selectedJets);
  fNJets           = lCount;
  fNFwd            = lNFwd;
  fNBTagsL         = lNBTagL;
  fNBTagsM         = lNBTagM;
  fNBTagsT         = lNBTagT;
  fNJetsdR08       = lCountdR08;
  fNBTagsLdR08     = lNBTagLdR08;
  fNBTagsMdR08     = lNBTagMdR08;
  fNBTagsTdR08     = lNBTagTdR08;

  fillJet(fN,fLooseJets,fVars);
  fillOthers(fN,fLooseJets,fVars,iVJets[0]);
}
void JetLoader::addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals) { 
  for(int i0 = 0; i0 < iN; i0++) { 
    int lBase = iN*3.+i0*2.;
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
void JetLoader::fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals, TLorentzVector iVJet){ 
  int lBase = 3.*fN;
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) {
    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(iObjects[i0]->pt, iObjects[i0]->eta, iObjects[i0]->phi, iObjects[i0]->mass); 
    iVals[lBase+i0*6+0] = iObjects[i0]->mass;
    iVals[lBase+i0*6+1] = iObjects[i0]->csv;
    iVals[lBase+i0*6+2] = iObjects[i0]->qgid;
    iVals[lBase+i0*6+3] = vPJet.DeltaR(iVJet);
    iVals[lBase+i0*6+4] = vPJet.DeltaPhi(iVJet);

  }
}
