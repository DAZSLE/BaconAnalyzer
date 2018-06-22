#include <iostream>
#include <assert.h> 
#include <string> 
#include <sstream>

#include "../include/GenLoader.hh"

using namespace baconhep;

GenLoader::GenLoader(TTree *iTree) { 
  fGenInfo  = new TGenEventInfo();
  iTree->SetBranchAddress("GenEvtInfo",       &fGenInfo);
  fGenInfoBr  = iTree->GetBranch("GenEvtInfo");

  fGens  = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("GenParticle",       &fGens);
  fGenBr  = iTree->GetBranch("GenParticle");

}
GenLoader::~GenLoader() { 
  delete fGenInfo;
  delete fGenInfoBr;

  delete fGens;
  delete fGenBr;
}
void GenLoader::reset() { 
  fBosonPt  = -1;
  fBosonPhi = -999;
  fBosonEta = -999;
  fBosonMass = -1;
  fBosonPdgId = -1;
  genEleFromW = -1;
  genMuFromW = -1;
  genTauFromW = -1;
  fTopPt = -1;
  fAntitopPt = -1;
  fTopPtWeight = -1;
}
void GenLoader::setupTree(TTree *iTree,float iXSIn) { 
  reset();
  fTree = iTree;
  fTree->Branch("genVPt"     ,&fBosonPt   ,"fBosonPt/F");
  fTree->Branch("genVPhi"    ,&fBosonPhi  ,"fBosonPhi/F");
  fTree->Branch("genVMass"    ,&fBosonMass  ,"fBosonMass/F");
  fTree->Branch("genVEta"    ,&fBosonEta  ,"fBosonEta/F");
  fTree->Branch("genVPdfId"    ,&fBosonPdgId  ,"fBosonPdgId/I");
  fTree->Branch("genEleFromW"    ,&genEleFromW  ,"genEleFromW/I");
  fTree->Branch("genMuFromW"    ,&genMuFromW  ,"genMuFromW/I");
  fTree->Branch("genTauFromW"    ,&genTauFromW  ,"genTauFromW/I");
  fTree->Branch("topPt"     ,&fTopPt   ,"fTopPt/F");
  fTree->Branch("antitopPt"     ,&fAntitopPt   ,"fAntitopPt/F");
  fTree->Branch("topPtWeight"     ,&fTopPtWeight   ,"fTopPtWeight/F");
}
void GenLoader::resetHiggs() {
  for(int i0 = 0; i0 < int(fgenHPt.size()); i0++) {
    fgenHPt[i0] = -99;
    fgenHEta[i0] = -99;
    fgenHPhi[i0] = -99;
    fgenHMass[i0] = -99;
    for(int i1 = 1; i1 < int(fgenHDauPt[i0].size()); i1++) {
      fgenHDauPt[i0][i1] = -99;
      fgenHDauEta[i0][i1] = -99;
      fgenHDauPhi[i0][i1] = -99;
      fgenHDauM[i0][i1] = -99;
      fgenHDauId[i0][i1] = -99;
      fgenHDauDecay[i0][i1] = -99;
    }
  }
}
void GenLoader::setupTreeHiggs(TTree *iTree) {
  resetHiggs();
  fTree = iTree;
  fgenHPt.clear(); fgenHEta.clear(); fgenHPhi.clear(); fgenHMass.clear(); 
  fgenHDauPt.clear(); fgenHDauEta.clear(); fgenHDauPhi.clear(); fgenHDauM.clear(); fgenHDauId.clear(); fgenHDauDecay.clear();
  for(int i0 = 0; i0 < 2; i0++) {
    fgenHPt.push_back(-999);
    fgenHEta.push_back(-999);
    fgenHPhi.push_back(-999);
    fgenHMass.push_back(-999);
    std::vector<float> lPt,lEta,lPhi,lM,lId,lDecay;
    for(int i1 = 0; i1 < 4; i1++) {
      lPt.push_back(-999);
      lEta.push_back(-999);
      lPhi.push_back(-999);
      lM.push_back(-999);
      lId.push_back(-999);
      lDecay.push_back(-999);
    }
    fgenHDauPt.push_back(lPt);
    fgenHDauEta.push_back(lEta);
    fgenHDauPhi.push_back(lPhi);
    fgenHDauM.push_back(lM);
    fgenHDauId.push_back(lId);
    fgenHDauDecay.push_back(lDecay);
  }
  for(int i0 = 0; i0 <  int(fgenHPt.size()); i0++) {
    std::stringstream pSPt;    pSPt << "genHPt" << i0;
    std::stringstream pSEta;   pSEta << "genHEta" << i0;
    std::stringstream pSPhi;   pSPhi << "genHPhi" << i0;
    std::stringstream pSMass;  pSMass << "genHMass" << i0;
    fTree->Branch(pSPt.str().c_str()   ,&fgenHPt[i0]   ,(pSPt.str()+"/F").c_str());
    fTree->Branch(pSEta.str().c_str()  ,&fgenHEta[i0]  ,(pSEta.str()+"/F").c_str());
    fTree->Branch(pSMass.str().c_str() ,&fgenHMass[i0] ,(pSMass.str()+"/F").c_str());
    fTree->Branch(pSPhi.str().c_str()  ,&fgenHPhi[i0]  ,(pSPhi.str()+"/F").c_str());
    for(int i1 = 0; i1 < int(fgenHDauPt[i0].size()); i1++) {
      std::stringstream pSPt;   pSPt << "genHDau" << i0 << "Pt" << i1;
      std::stringstream pSEta;  pSEta << "genHDau" << i0 << "Eta" << i1;
      std::stringstream pSPhi;  pSPhi << "genHDau" << i0 << "Phi" << i1;
      std::stringstream pSM;    pSM << "genHDau" << i0 << "Mass" << i1;
      std::stringstream pSId;   pSId << "genHDau" << i0 << "Id" << i1;
      std::stringstream pSDecay;   pSDecay << "genHDau" << i0 << "Decay" << i1;
      fTree->Branch(pSPt.str().c_str()   ,&fgenHDauPt[i0][i1]   ,(pSPt.str()+"/F").c_str());
      fTree->Branch(pSEta.str().c_str()  ,&fgenHDauEta[i0][i1]  ,(pSEta.str()+"/F").c_str());
      fTree->Branch(pSPhi.str().c_str()  ,&fgenHDauPhi[i0][i1]  ,(pSPhi.str()+"/F").c_str());
      fTree->Branch(pSM.str().c_str()    ,&fgenHDauM[i0][i1]    ,(pSM.str()+"/F").c_str());
      fTree->Branch(pSId.str().c_str()   ,&fgenHDauId[i0][i1]   ,(pSId.str()+"/F").c_str());
      fTree->Branch(pSDecay.str().c_str()  ,&fgenHDauDecay[i0][i1]     ,(pSDecay.str()+"/F").c_str());
    }
  }
}
void GenLoader::load(int iEvent) { 
  reset();
  fGens     ->Clear();
  fGenBr    ->GetEntry(iEvent);
  fGenInfoBr->GetEntry(iEvent);

  fWeight = fGenInfo->weight;
}
bool GenLoader::isType(std::string boson,std::string mode)
{
  int iPDGID,iId;
  if (boson.find("Z")==0) iPDGID = 23;
  if (boson.find("W")==0) iPDGID = 24;
  if (boson.find("H")==0) iPDGID = 25;
  if (boson.find("Zprime")==0) iPDGID = 10031;
  if (boson.find("DMSpin0")==0) iPDGID = 9900032;  
  if (mode.find("bb")==0) iId = 5;
  if (mode.find("cc")==0) iId = 4;

  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if (mode.find("bb")==0 || mode.find("cc")==0) {
      if(abs(pGen->pdgId)==iId) {
	if(pGen->parent<0) continue;
	TGenParticle *parent = (TGenParticle*)((*fGens)[pGen->parent]);
	if(abs(parent->pdgId)==iPDGID) return true;
      }
    }
    if (mode.find("cs")==0) {
      if(abs(pGen->pdgId)==4 || abs(pGen->pdgId)==3) {
	if(pGen->parent<0) continue;
	TGenParticle *parent = (TGenParticle*)((*fGens)[pGen->parent]);
	if(abs(parent->pdgId)==iPDGID) return true;
      }
    }
  }
  
  return false;
}

bool GenLoader::hard(int &iP)
{
  TGenParticle *p = (TGenParticle*)((*fGens)[iP]);
  // if particle itself from hard process                                                                                                                                                                                                                     
  if(p->status>20 && p->status<30) return true;
  return false;
}
bool GenLoader::hasChild(int &iparent, bool beHard)
{
  TGenParticle *parent = (TGenParticle*)((*fGens)[iparent]);
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *child = (TGenParticle*)((*fGens)[i0]);
    if (child->pdgId != parent->pdgId)
      continue;
    if (beHard && !hard(i0))
      continue;
    if (child->parent !=-2 &&
        child->parent == iparent) {
      return true;
    }
  }
  return false;
}
TGenParticle* GenLoader::findDaughter(int iparent, int dauId)
{
  for(int k = iparent+1; k < fGens->GetEntriesFast(); k++) {
    TGenParticle *genp = (TGenParticle*)((*fGens)[k]);
    if(genp->parent == iparent) {
      if(abs(genp->pdgId) == dauId) {
        return genp;
      }
    }
  }
  return 0;
}
int GenLoader::findDaughterId(int iparent, int dauId)
{
  for(int k = iparent+1; k < fGens->GetEntriesFast(); k++) {
    const baconhep::TGenParticle *genp = (TGenParticle*)((*fGens)[k]);;
    if(genp->parent == iparent) {
      if(abs(genp->pdgId) == dauId) {
        return k;
      }
    }
  }
  return -1;
}
int GenLoader::findLastParent(int iparent,int iId){
  Bool_t foundLast = kFALSE;
  int iLast = iparent;
  while (!foundLast) {
    int tmpId = findDaughterId(iLast,iId);
    if (tmpId>=0) iLast = tmpId;
    else foundLast = kTRUE;
  }
  return iLast;
}
void GenLoader::findBoson(int iId, int lOption){
  reset();
  float pbosonPt(-1),pbosonPhi(-999),pbosonMass(-1),pbosonEta(-999);
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *genp0 = (TGenParticle*)((*fGens)[i0]);

    // find highest Pt boson G(22)                                                                                                                                                                                                                            
    if(lOption == 0){
      if(fabs(genp0->pdgId)==iId && genp0->pt > pbosonPt){
        pbosonPt = genp0->pt;
        pbosonPhi = genp0->phi;
        pbosonEta = genp0->eta;
        pbosonMass = genp0->mass;
      }
    }

    // find last boson Z(23),W(24),Z'(10031),H(25),DMSpin(9900032)                                                                                                                                                                                            
    if(lOption == 1){
      if(fabs(genp0->pdgId)==iId){
        int iL0 = findLastParent(i0,iId);
        for(int k0 = 0; k0 < fGens->GetEntriesFast(); k0++) {
          TGenParticle *genp1 = (TGenParticle*)((*fGens)[k0]);
          if(k0==iL0){
            pbosonPt = genp1->pt;
            pbosonPhi = genp1->phi;
            pbosonEta = genp1->eta;
            pbosonMass = genp1->mass;
            break;
          }
        }
      }
    }

    // find W(24) for ttbar semilep                                                                                                                                                                                                                           
    if(lOption == 2){
      if(fabs(genp0->pdgId)==iId) {
        TGenParticle *dau1 = findDaughter(i0, 11);
        TGenParticle *dau2 = findDaughter(i0, 13);
        if(dau1 || dau2){
          pbosonPt = genp0->pt;
          pbosonPhi = genp0->phi;
          pbosonEta = genp0->eta;
          pbosonMass = genp0->mass;
        }
      }
    }

    // find W for ttbar dileptonic (6)                                                                                                                                                                                                                        
    if(lOption == 3){
      if(fabs(genp0->pdgId)==iId) {
        int iW0 = findLastParent(i0,24);
        for(int k0 = 0; k0 < fGens->GetEntriesFast(); k0++) {
          TGenParticle *dau0 = (TGenParticle*)((*fGens)[k0]);
          TGenParticle *ele0 = findDaughter(iW0, 11); TGenParticle *muo0 = findDaughter(iW0, 13);
          if(k0==iW0 && (ele0 || muo0)){
            for(int i1=iW0+1; i0 < fGens->GetEntriesFast(); i1++) {
              TGenParticle *genp1 = (TGenParticle*)((*fGens)[i1]);
              if(genp1->pdgId == iId) {
                int iW1 = findLastParent(i1,24);
                for(int k1 = 0; k1 < fGens->GetEntriesFast(); k1++) {
                  TGenParticle *dau1 = (TGenParticle*)((*fGens)[k1]);
                  TGenParticle *ele1 = findDaughter(iW1, 11); TGenParticle *muo1 = findDaughter(iW1, 13);
                  if(k1==iW1 && (ele1 || muo1)){
                    TLorentzVector vDau0; vDau0.SetPtEtaPhiM(dau0->pt, dau0->eta, dau0->phi, dau0->mass);
                    TLorentzVector vDau1; vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
                    pbosonPt = (vDau0 + vDau1).Pt();
                    pbosonPhi = (vDau0 + vDau1).Phi();
                  }
                }
              }
            }
          }
        }
      }
    }

    // find hadronic W(24) for ttbar semilep                                                                                                                                                                                                                  
    if(lOption == 4){
      if(fabs(genp0->pdgId)==iId) {
        TGenParticle *dau1 = findDaughter(i0, 2);
        TGenParticle *dau2 = findDaughter(i0, 3);
        if(dau1 || dau2){
          pbosonPt = genp0->pt;
          pbosonPhi = genp0->phi;
          pbosonEta = genp0->eta;
          pbosonMass = genp0->mass;
        }
      }
    }


  }
  fBosonPt = pbosonPt;
  fBosonPhi = pbosonPhi;
  fBosonMass = pbosonMass;
  fBosonEta = pbosonEta;
  fBosonPdgId = iId;
}
// is W(qq) in Top->Wb?                                                                                                                                                                                                                                       
int GenLoader::isHadronicWInTop(TGenParticle *genp,int j,TLorentzVector jet,double dR,double &wMatching, double &wSize)
{
  TLorentzVector vW,vDau1,vDau2;
  wMatching = -999.; wSize = -999.;
  double tmpWMatching(0), tmpWSize(0);
  if(abs(genp->pdgId)==24) {
    int iW = j;
    vW.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
    int iQ=0, jQ=0;
    for (; iQ<fGens->GetEntriesFast(); ++iQ) {
      TGenParticle *dau1 = (TGenParticle*)((*fGens)[iQ]);
      if(dau1->parent==iW && abs(dau1->pdgId)<6) {
        vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
        // if (vDau1.DeltaR(jet) > dR) return false;                                                                                                                                                                                                          
        tmpWMatching = TMath::Max(tmpWMatching,jet.DeltaR(vDau1));
        tmpWSize     = TMath::Max(tmpWSize,vW.DeltaR(vDau1));
        //std::cout << "found the first quark: tmpWMatching = " << jet.DeltaR(vDau1) <<  std::endl;                                                                                                                                                           
        //std::cout << "found the first quark: tmpWSize = " << vW.DeltaR(vDau1) <<  std::endl;                                                                                                                                                                
        break; // found the first quark                                                                                                                                                                                                                       
      }
    }
    for (jQ=iQ+1; jQ<fGens->GetEntriesFast(); ++jQ) {
      TGenParticle *dau2 = (TGenParticle*)((*fGens)[jQ]);
      if(dau2->parent==iW && abs(dau2->pdgId)<6) {
        vDau2.SetPtEtaPhiM(dau2->pt, dau2->eta, dau2->phi, dau2->mass);
        // if (vDau2.DeltaR(jet) > dR) return false;                                                                                                                                                                                                          
        tmpWMatching = TMath::Max(tmpWMatching,jet.DeltaR(vDau2));
        tmpWSize     = TMath::Max(tmpWSize,vW.DeltaR(vDau2));
        wMatching    = tmpWMatching;
        wSize        = tmpWSize;
        //std::cout << "found the second quark: tmpWMatching = " << jet.DeltaR(vDau2) <<  std::endl;                                                                                                                                                          
        //std::cout << "found the second quark: tmpWSize = " << vW.DeltaR(vDau2) <<  std::endl;                                                                                                                                                               
        //std::cout << "best quark: wMatching = " << wMatching <<  std::endl;                                                                                                                                                                                 
        //std::cout << "best quark: tmpWSize = " << wSize <<  std::endl;                                                                                                                                                                                      
        return 1;
      }
    }
  }
  return 0;
}
// t->W(qq)b                                                                                                                                                                                                                                                  
int GenLoader::isHadronicTop(TGenParticle *genp,int j,TLorentzVector jet,double dR,double &topMatching, double &topSize)
{
  TLorentzVector vTop,vB,vDau1,vDau2;
  topMatching = -999.; topSize = -999.;
  double tmpTopMatching(0), tmpTopSize(0);
  if(abs(genp->pdgId)==6) {
    vTop.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
    TGenParticle *mcB = findDaughter(j,5); //                                                                                                                                                                                                                 
    if(mcB){
      vB.SetPtEtaPhiM(mcB->pt, mcB->eta, mcB->phi, mcB->mass);
    }
    TGenParticle *mcW = findDaughter(j,24); //                                                                                                                                                                                                                
    if (!mcW || !mcB) return 0;     // this shouldn't happen                                                                                                                                                                                                  
    // if (vB.DeltaR(jet) > dR) return false; // all decay products fall into jet cone                                                                                                                                                                        
    tmpTopMatching = TMath::Max(tmpTopMatching,jet.DeltaR(vB));
    tmpTopSize     = TMath::Max(tmpTopSize,vTop.DeltaR(vB));

    int iW = findLastParent(j,24);

    int iQ=0, jQ=0;
    for (; iQ<fGens->GetEntriesFast(); ++iQ) {
      TGenParticle *dau1 = (TGenParticle*)((*fGens)[iQ]);
      if(dau1->parent==iW && abs(dau1->pdgId)<6) {
        vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
        // if (vDau1.DeltaR(jet) > dR) return false;                                                                                                                                                                                                          
        tmpTopMatching = TMath::Max(tmpTopMatching,jet.DeltaR(vDau1));
        tmpTopSize     = TMath::Max(tmpTopSize,vTop.DeltaR(vDau1));
        break; // found the first quark                                                                                                                                                                                                                       
      }
    }
    for (jQ=iQ+1; jQ<fGens->GetEntriesFast(); ++jQ) {
      TGenParticle *dau2 = (TGenParticle*)((*fGens)[jQ]);
      if(dau2->parent==iW && abs(dau2->pdgId)<6) {
        vDau2.SetPtEtaPhiM(dau2->pt, dau2->eta, dau2->phi, dau2->mass);
        // if (vDau2.DeltaR(jet) > dR) return false;                                                                                                                                                                                                          
        tmpTopMatching = TMath::Max(tmpTopMatching,jet.DeltaR(vDau2));
        tmpTopSize     = TMath::Max(tmpTopSize,vTop.DeltaR(vDau2));
        topMatching    = tmpTopMatching;
        topSize        = tmpTopSize;
        return 1;
      }
    }
  }
  return 0;
}
// jet matched to V(qq)                                                                                                                                                                                                                                       
int GenLoader::isHadronicV(TGenParticle *genp,int j,int iId, TLorentzVector jet,double dR,double &vMatching, double &vSize)
{
  TLorentzVector vV,vDau1,vDau2;
  vMatching = -999.; vSize = -999.;
  double tmpVMatching(0), tmpVSize(0);
  if(abs(genp->pdgId)==iId){
    vV.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
    int iV = findLastParent(j,iId);

    int iQ=0, jQ=0;
    for (; iQ<fGens->GetEntriesFast(); ++iQ) {
      TGenParticle *dau1 = (TGenParticle*)((*fGens)[iQ]);
      if(dau1->parent==iV && abs(dau1->pdgId)<6) {
        vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
        tmpVMatching = TMath::Max(tmpVMatching,jet.DeltaR(vDau1));
        tmpVSize     = TMath::Max(tmpVSize,vV.DeltaR(vDau1));
        break;
      }
    }
    for (jQ=iQ+1; jQ<fGens->GetEntriesFast(); ++jQ) {
      TGenParticle *dau2 = (TGenParticle*)((*fGens)[jQ]);
      if(dau2->parent==iV && abs(dau2->pdgId)<6) {
        vDau2.SetPtEtaPhiM(dau2->pt, dau2->eta, dau2->phi, dau2->mass);
        tmpVMatching = TMath::Max(tmpVMatching,jet.DeltaR(vDau2));
        tmpVSize     = TMath::Max(tmpVSize,vV.DeltaR(vDau2));
        vMatching    = tmpVMatching;
        vSize        = tmpVSize;
        return 1;
      }
    }
  }
  return 0;
}

// jet matched to V(qq) with quark flavor?                                                                                                                                                                                                                    
int GenLoader::isHadronicVflav(TGenParticle *genp,int j,int iId, TLorentzVector jet,double dR,double &vMatching, double &vSize, int dauId)
{
  TLorentzVector vV,vDau1,vDau2;
  vMatching = -999.; vSize = -999.;
  double tmpVMatching(0), tmpVSize(0);
  if(abs(genp->pdgId)==iId){
    vV.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
    int iV = findLastParent(j,iId);

    int iQ=0, jQ=0;
    for (; iQ<fGens->GetEntriesFast(); ++iQ) {
      TGenParticle *dau1 = (TGenParticle*)((*fGens)[iQ]);
      if(dau1->parent==iV && abs(dau1->pdgId)==dauId) {
        vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
        tmpVMatching = TMath::Max(tmpVMatching,jet.DeltaR(vDau1));
        tmpVSize     = TMath::Max(tmpVSize,vV.DeltaR(vDau1));
        break;
      }
    }
    for (jQ=iQ+1; jQ<fGens->GetEntriesFast(); ++jQ) {
      TGenParticle *dau2 = (TGenParticle*)((*fGens)[jQ]);
      if(dau2->parent==iV && abs(dau2->pdgId)<=dauId) {  //22 or 21; 43 or 44 we don't check for 33, 11, 55                                                                                                                                                   
        vDau2.SetPtEtaPhiM(dau2->pt, dau2->eta, dau2->phi, dau2->mass);
        tmpVMatching = TMath::Max(tmpVMatching,jet.DeltaR(vDau2));
        tmpVSize     = TMath::Max(tmpVSize,vV.DeltaR(vDau2));
        vMatching    = tmpVMatching;
        vSize        = tmpVSize;
        return 1;
      }
    }
  }
  return 0;
}

// jet matched to something,                                                                                                                                                                                                                                  
// For topcitos(6, or 624 in 2016 sample?) result = 1                                                                                                                                                                                                         
// For bosons(24,23,10031,25,55) if result:                                                                                                                                                                                                                   
// 1: V(qq)                                                                                                                                                                                                                                                   
// 2: W->ud, Z->dd, H->dd                                                                                                                                                                                                                                     
// 3: W->cs, Z->cc, H->cc                                                                                                                                                                                                                                     
// 4: Z->bb; H->bb                                                                                                                                                                                                                                            
int GenLoader::ismatchedJet(TLorentzVector jet0, double dR,double &matching, double &size, int iId){
  //std::cout << "new event" << std::endl;                                                                                                                                                                                                                    
  int result=0;
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *genp0 = (TGenParticle*)((*fGens)[i0]);
    TLorentzVector mcMom; mcMom.SetPtEtaPhiM(genp0->pt,genp0->eta,genp0->phi,genp0->mass);
    if (mcMom.DeltaR(jet0) < dR) {
      if(iId == 624 && isHadronicWInTop(genp0,i0,jet0,dR,matching,size)==1) {
        result= 1;
        //std::cout<< result << " isMatched "<<iId << std::endl;                                                                                                                                                                                              
        //std::cout<< " matching "<< matching << std::endl;                                                                                                                                                                                                   
        //std::cout<< " size "<< size << std::endl;                                                                 if(iId == 24 || iId == 23 || iId == 10031 || iId == 25 || iId == 55){
        if (isHadronicV(genp0,i0,iId,jet0,dR,matching,size)==1) result= 1;
        if (isHadronicVflav(genp0,i0,iId,jet0,dR,matching,size,2)==1) result= 2; //W->ud, Z->dd, H->dd                                                                                                                                                        
        if (isHadronicVflav(genp0,i0,iId,jet0,dR,matching,size,4)==1) result= 3; //W->cs, Z->cc, H->cc                                                                                                                                                        
        if (isHadronicVflav(genp0,i0,iId,jet0,dR,matching,size,5) ==1) result= 4;//Z->bb; H->bb                                                                                                                                                               
      }
    }
  }
  return result;
}
int GenLoader::ismatchedSubJet(TLorentzVector subjet0){
  int lOption =0;
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *genp0 = (TGenParticle*)((*fGens)[i0]);
    TLorentzVector vq;
    if(abs(genp0->pdgId)==5) {
      vq.SetPtEtaPhiM(genp0->pt, genp0->eta, genp0->phi, 5);
      if(vq.DeltaR(subjet0) < 0.4) lOption = 1; //isB                                                                                                                                                                                                         
    }
    else if(abs(genp0->pdgId)==4) {
      vq.SetPtEtaPhiM(genp0->pt, genp0->eta, genp0->phi, 1.29);
      if(vq.DeltaR(subjet0) < 0.4) lOption = 2; //isC                                                                                                                                                                                                         
    }
    else lOption = 3; //isLF                                                                                                                                                                                                                                  
  }
  return lOption;
}
int GenLoader::isHadronicBoson(int iV,int iId, float &genSize)
{
  TGenParticle *genV = (TGenParticle*)((*fGens)[iV]);
  TLorentzVector vV,vDau1,vDau2;
  if(abs(genV->pdgId)==iId) {
    vV.SetPtEtaPhiM(genV->pt, genV->eta, genV->phi, genV->mass);
    int iQ=0, jQ=0;
    for (; iQ<fGens->GetEntriesFast(); ++iQ) {
      TGenParticle *dau1 = (TGenParticle*)((*fGens)[iQ]);
      if(dau1->parent==iV && abs(dau1->pdgId)<6) {
        vDau1.SetPtEtaPhiM(dau1->pt, dau1->eta, dau1->phi, dau1->mass);
        break; // found the first quark                                                                                                                                                                                                                       
      }
    }
    for (jQ=iQ+1; jQ<fGens->GetEntriesFast(); ++jQ) {
      TGenParticle *dau2 = (TGenParticle*)((*fGens)[jQ]);
      if(dau2->parent==iV && abs(dau2->pdgId)<6) { // found second                                                                                                                                                                                            
        vDau2.SetPtEtaPhiM(dau2->pt, dau2->eta, dau2->phi, dau2->mass);
        genSize = TMath::Max(vV.DeltaR(vDau1),vV.DeltaR(vDau2));
        return 1;
      }
    }
  }
  return 0;
}
int GenLoader::isHDau(int iId, int iDauId, TLorentzVector jet, int iHiggs){
  int iH(-1),  iR(-1);
  std::vector<TLorentzVector> lParts;
  std::vector<int> lHiggs;
  for(int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(fabs(pGen->pdgId == iId)) {
      lHiggs.push_back(findLastParent(i0,iId));
    }
  }
  unsigned int lMax =2;
  float lMaxR = 0.8;
  if(lHiggs.size()<lMax) lMax = lHiggs.size();
  for(unsigned int i0 = 0; i0 < lMax; i0++)
    {
      int iD(-1);
      iH = lHiggs[i0];
      TGenParticle *pH = (TGenParticle*)((*fGens)[iH]);
      TLorentzVector ivH; ivH.SetPtEtaPhiM(pH->pt, pH->eta, pH->phi, pH->mass);
      fgenHPt[i0] = pH->pt;
      fgenHEta[i0] = pH->eta;
      fgenHPhi[i0] = pH->phi;
      fgenHMass[i0] = pH->mass;
      for(int i1 = 0; i1 < fGens->GetEntriesFast(); i1++) {
        TGenParticle *pGen = (TGenParticle*)((*fGens)[i1]);
        if(pGen->parent == iH) {
          //if(fabs(pGen->pdgId) == iDauId) {
	  iD += 1;
	  if(iD > 3) break;
	  int iP = findLastParent(i1,iDauId);
	  TGenParticle *pP = (TGenParticle*)((*fGens)[iP]);
	  //std::cout << i1 << " iD  "<<  iD << std::endl;
	  //std::cout << pP->pt << " id " << pP->pdgId << std::endl;
	  fgenHDauPt[i0][iD] = pP->pt;
	  fgenHDauEta[i0][iD] = pP->eta;
	  fgenHDauPhi[i0][iD] = pP->phi;
	  fgenHDauM[i0][iD] = pP->mass;
	  fgenHDauId[i0][iD] = pP->pdgId;
	  for(int i2 = 0; i2 < fGens->GetEntriesFast(); i2++) {
	    TGenParticle *pPart = (TGenParticle*)((*fGens)[i2]);
	    TLorentzVector vPart; vPart.SetPtEtaPhiM(pPart->pt, pPart->eta, pPart->phi, pPart->mass);
	    if(fabs(pPart->parent) ==iP) {
	      fgenHDauDecay[i0][iD] = pPart->pdgId;
	      break;
	    }
	  }
	  if(fabs(pGen->pdgId) == iDauId) {
	    TLorentzVector ivP; ivP.SetPtEtaPhiM(pP->pt, pP->eta, pP->phi, pP->mass);
	    lParts.push_back(ivP);
	  }
	}
      }

      if(jet.DeltaR(ivH)<lMaxR && iD > -1) {
	lMaxR = jet.DeltaR(ivH);
        for(unsigned int i1 = 0; i1<lParts.size(); i1++) {
          if(jet.DeltaR(lParts[i1])<lMaxR){ iR = i0;}
        }
      }
    }
  return iR;
}
int GenLoader::isttbarType(int lepPdgId)
{
  assert(fGens);
  int nlep=0;
  for(int i0=0; i0<fGens->GetEntriesFast(); i0++) {
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(abs(pGen->pdgId)==abs(lepPdgId)) {
      if(pGen->parent<0) continue;
      TGenParticle *lparent = (TGenParticle*)((*fGens)[pGen->parent]);
      if(abs(lparent->pdgId)==24) nlep++;
    }
  }
  return nlep;
}

void GenLoader::saveTTbarType() {
  genEleFromW = isttbarType(11);
  genMuFromW = isttbarType(13);
  genTauFromW = isttbarType(15);
}

float GenLoader::computeTTbarCorr() {
  
  //                                                                                                                                                                                                                            
  // compute ttbar MC pT correction                                                                                                                                                                                                                         
  //                                                                                                                                                                            
  const int TOP_PDGID = 6;
  double pt1=0, pt2=0;
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *p = (TGenParticle*)((*fGens)[i0]);
    if(p->pdgId ==  TOP_PDGID) {
      pt1 = p->pt;
      fTopPt = pt1;
    }
    if(p->pdgId == -TOP_PDGID) {
      pt2 = p->pt;
      fAntitopPt = pt2;
    }
  }

  // Reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#MC_SFs_Reweighting
  // 8 TeV values:
  //double w1 = exp(0.156 - 0.00137*pt1);
  //double w2 = exp(0.156 - 0.00137*pt2);
  // 13 TeV values:
  double w1 = exp(0.0615 - 0.0005*pt1);
  double w2 = exp(0.0615 - 0.0005*pt2);

  fTopPtWeight = sqrt(w1*w2);

  // 8 TeV:
  //return 1.001*sqrt(w1*w2);
  // 13 TeV:
  return sqrt(w1*w2);
}
