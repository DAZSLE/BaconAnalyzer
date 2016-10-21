#include "../include/JetLoader.hh"
#include <cmath>
#include <iostream> 
#include <sstream>

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree, std::string iJet) { 
  fJets  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress(iJet.c_str(),       &fJets);
  fJetBr = iTree->GetBranch(iJet.c_str());
  fN = 4;
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
}
void JetLoader::reset() { 
  fNJets           = 0;
  fMT              = 0;
  fMinDPhi         = 1000;
  fMinDFPhi        = 1000; 
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
  fNJetsdR15       = 0;
  fNBTagsLdR15     = 0;
  fNBTagsMdR15     = 0;
  fNBTagsTdR15     = 0;
}
void JetLoader::setupTree(TTree *iTree, std::string iJetLabel) { 
  reset();
  fTree = iTree;
  std::stringstream pSN,pSdPhi,pSdFPhi,pSMT,pSb,pSbL,pSbM,pSbT,pSbLdR,pSbMdR,pSbTdR,pSNdR08,pSbLdR08,pSbMdR08,pSbTdR08;   
  pSN     << "n" << iJetLabel << "s";
  pSdPhi  << iJetLabel << "mindPhi";
  pSdFPhi << iJetLabel << "mindFPhi";
  pSMT    << iJetLabel << "mT";
  pSb     << iJetLabel << "sbtag";
  pSbL    << iJetLabel << "sbtagL";
  pSbM    << iJetLabel << "sbtagM";
  pSbT    << iJetLabel << "sbtagT";
  pSNdR08 << iJetLabel << "sdR08";
  pSbLdR08<< iJetLabel << "sLdR08";
  pSbMdR08<< iJetLabel << "sMdR08";
  pSbTdR08<< iJetLabel << "sTdR08";

  fTree->Branch(pSN.str().c_str()           ,&fNJets           ,(pSN.str()+"/I").c_str());  // jet multiplicity
  fTree->Branch(pSdPhi.str().c_str()        ,&fMinDPhi         ,(pSdPhi.str()+"/D").c_str());
  fTree->Branch(pSdFPhi.str().c_str()       ,&fMinDFPhi        ,(pSdFPhi.str()+"/D").c_str());
  fTree->Branch(pSMT.str().c_str()          ,&fMT              ,(pSMT.str()+"/F").c_str());
  fTree->Branch(pSbL.str().c_str()          ,&fNBTagsL         ,(pSbL.str()+"/I").c_str()); // b tags
  fTree->Branch(pSbM.str().c_str()          ,&fNBTagsM         ,(pSbM.str()+"/I").c_str());
  fTree->Branch(pSbT.str().c_str()          ,&fNBTagsT         ,(pSbT.str()+"/I").c_str());
  fTree->Branch(pSNdR08.str().c_str()       ,&fNJetsdR08       ,(pSNdR08.str()+"/I").c_str());
  fTree->Branch(pSbLdR08.str().c_str()      ,&fNBTagsLdR08     ,(pSbLdR08.str()+"/I").c_str());
  fTree->Branch(pSbMdR08.str().c_str()      ,&fNBTagsMdR08     ,(pSbMdR08.str()+"/I").c_str());
  fTree->Branch(pSbTdR08.str().c_str()      ,&fNBTagsTdR08     ,(pSbTdR08.str()+"/I").c_str());
  
  for(int i0 = 0; i0 < fN*(10)+3; i0++) {double pVar = 0; fVars.push_back(pVar);}           // declare array of 43 vars
  setupNtuple(iJetLabel.c_str(),iTree,fN,fVars);                                            // from MonoXUtils.cc => fN =4 j*_pt,j*_eta,j*_phi for j1,j2,j3,j4 (3*4=12)
  addOthers  (iJetLabel.c_str(),iTree,fN,fVars);                                            // Mass + b-tag for j1,j2,j3,j4 (2*4=8)
}
void JetLoader::load(int iEvent) { 
  fJets   ->Clear();
  fJetBr ->GetEntry(iEvent);
}
void JetLoader::selectJets(std::vector<TLorentzVector> &iElectrons, std::vector<TLorentzVector> &iMuons, std::vector<TLorentzVector> &iPhotons, std::vector<TLorentzVector> &iVJets,float iMet, float iMetPhi,float iFMet,float iFMetPhi){ 
  reset(); 
  int lCount = 0,lCountdR2 = 0,lNBTagL = 0,lNBTagM = 0,lNBTagT = 0,lNBTagLdR2 = 0,lNBTagMdR2 = 0,lNBTagTdR2 = 0;
  double pDPhi(999),pDFPhi(999);
  float  pMhtX(0.), pMhtY(0.);

  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet->eta,pJet->phi,0.4,iElectrons))                      continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iMuons))                          continue;
    if(passVeto(pJet->eta,pJet->phi,0.4,iPhotons))                        continue;
    if(pJet->pt        <=  30)                                            continue;
    if(fabs(pJet->eta) >= 4.5)                                            continue;
    if(!passJetLooseSel(pJet))                                            continue;
    lCount++;
    addJet(pJet,fLooseJets);

    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    fHT += vPJet.Pt(); pMhtX += vPJet.Px(); pMhtY += vPJet.Py();
    fGoodJets.push_back(pJet);

    if(acos(cos(iMetPhi-pJet->phi))    < pDPhi)    pDPhi  = acos(cos(iMetPhi-pJet->phi));
    if(iFMet > 0){
      if(acos(cos(iFMetPhi-pJet->phi)) < pDFPhi)   pDFPhi = acos(cos(iFMetPhi-pJet->phi)); 
    }

    // jet and b-tag multiplicity
    if(iVJets.size()>0 && iVJets[0].Pt()>0 && vPJet.DeltaR(iVJets[0])>2) lCountdR2++;
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVL){
      lNBTagL++;
      if(iVJets.size()>0) {if(pJet->pt>0 && vPJet.DeltaR(iVJets[0])>2) lNBTagLdR2++;}
    }
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVM){ 
      lNBTagM++;
      if(iVJets.size()>0) {if(pJet->pt>0 && vPJet.DeltaR(iVJets[0])>2) lNBTagMdR2++;}
    }
    if(fabs(pJet->eta) < 2.5 && pJet->csv > CSVT){
      lNBTagT++;
      if(iVJets.size()>0) {if(pJet->pt>0 && vPJet.DeltaR(iVJets[0])>2) lNBTagTdR2++;}
    }
  }
  addVJet(fLooseJets,selectedJets);
  fNJets           = lCount;
  fNJetsdR2        = lCountdR2;
  fMinDPhi         = pDPhi;
  fMinDFPhi        = pDFPhi;
  fNBTagsL         = lNBTagL;
  fNBTagsM         = lNBTagM;
  fNBTagsT         = lNBTagT;
  fNBTagsLdR2      = lNBTagLdR2;
  fNBTagsMdR2      = lNBTagMdR2;
  fNBTagsTdR2      = lNBTagTdR2;
  fillJet(fN,fLooseJets,fVars);
  fillOthers(fN,fLooseJets,fVars);
  TLorentzVector fMyMHT; fMyMHT.SetPxPyPzE(-pMhtX, -pMhtY, 0, sqrt(pow(pMhtX,2) + pow(pMhtY,2)));
  fMHT             = fMyMHT.Pt();
}
void JetLoader::fillGoodJets(std::vector<TLorentzVector> iVJets,double dR, std::vector<const TJet*> &iGoodJets){
  iGoodJets.clear();
  int lCountdR2 = 0,lNBTagLdR2 = 0,lNBTagMdR2 = 0,lNBTagTdR2 = 0;
  for(unsigned int i0 = 0; i0 < fGoodJets.size(); i0++) {
    TLorentzVector vPJet; vPJet.SetPtEtaPhiM(fGoodJets[i0]->pt, fGoodJets[i0]->eta, fGoodJets[i0]->phi, fGoodJets[i0]->mass);
    if(iVJets.size()>0 && iVJets[0].Pt()>0){
      if(vPJet.DeltaR(iVJets[0])>dR){
	iGoodJets.push_back(fGoodJets[i0]);
	lCountdR2++;
	if(fGoodJets[i0]->pt>0 && fabs(fGoodJets[i0]->eta) < 2.5 && fGoodJets[i0]->csv > CSVL) lNBTagLdR2++;
	if(fGoodJets[i0]->pt>0 && fabs(fGoodJets[i0]->eta) < 2.5 && fGoodJets[i0]->csv > CSVM) lNBTagMdR2++;
	if(fGoodJets[i0]->pt>0 && fabs(fGoodJets[i0]->eta) < 2.5 && fGoodJets[i0]->csv > CSVT) lNBTagTdR2++;
      }
    }
  }
  if(dR == 0.8){
    fNJetsdR08        = lCountdR2;
    fNBTagsLdR08      = lNBTagLdR2;
    fNBTagsMdR08      = lNBTagMdR2;
    fNBTagsTdR08      = lNBTagTdR2;  
  }
 if(dR == 1.5){
    fNJetsdR15        = lCountdR2;
    fNBTagsLdR15      = lNBTagLdR2;
    fNBTagsMdR15      = lNBTagMdR2;
    fNBTagsTdR15      = lNBTagTdR2;
  }
}
void JetLoader::addOthers(std::string iHeader,TTree *iTree,int iN,std::vector<double> &iVals) { 
  for(int i0 = 0; i0 < iN; i0++) { 
    int lBase = iN*3.+i0*2.;
    std::stringstream pSMass,pSCSV;
    pSMass  << iHeader << i0 << "_mass";
    pSCSV   << iHeader << i0 << "_csv";
    iTree->Branch(pSMass .str().c_str(),&iVals[lBase+0],(pSMass .str()+"/D").c_str());
    iTree->Branch(pSCSV .str().c_str() ,&iVals[lBase+1],(pSCSV  .str()+"/D").c_str());
  }
}
void JetLoader::fillOthers(int iN,std::vector<TJet*> &iObjects,std::vector<double> &iVals){ 
  int lBase = 3.*fN;
  int lMin = iObjects.size();
  if(iN < lMin) lMin = iN;
  for(int i0 = 0; i0 < lMin; i0++) { 
    iVals[lBase+i0*6+0] = iObjects[i0]->mass;
    iVals[lBase+i0*6+1] = iObjects[i0]->csv;
    iVals[lBase+i0*6+2] = iObjects[i0]->qgid;
    iVals[lBase+i0*6+3] = iObjects[i0]->chHadFrac;
    iVals[lBase+i0*6+4] = iObjects[i0]->neuHadFrac;
    iVals[lBase+i0*6+5] = iObjects[i0]->neuEmFrac;
  }
}
