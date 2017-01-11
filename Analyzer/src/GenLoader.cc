#include <iostream>
#include <assert.h> 
#include <string> 
#include "../include/GenLoader.hh"

using namespace baconhep;

GenLoader::GenLoader(TTree *iTree) { 
  fGenInfo  = new TGenEventInfo();
  iTree->SetBranchAddress("GenEvtInfo",       &fGenInfo);
  fGenInfoBr  = iTree->GetBranch("GenEvtInfo");

  fGens  = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("GenParticle",       &fGens);
  fGenBr  = iTree->GetBranch("GenParticle");
  fBoson = 0;
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
float GenLoader::lepmatched(int iId, std::vector<TLorentzVector> vec, double dR){
  if(vec.size() > 0){
    for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
      TGenParticle *genp0 = (TGenParticle*)((*fGens)[i0]);
      if(abs(genp0->pdgId) == iId) {
	TLorentzVector vGenp0; vGenp0.SetPtEtaPhiM(genp0->pt, genp0->eta, genp0->phi, genp0->mass);
	if(vec[0].DeltaR(vGenp0)<dR){
	  if(vec.size() == 1) return 1;
	  if(vec.size() == 2){
	    for(int i1=i0+1; i1<fGens->GetEntriesFast(); i1++) {
	      TGenParticle *genp1 = (TGenParticle*)((*fGens)[i1]);
	      if(abs(genp1->pdgId) == iId) {
		TLorentzVector vGenp1; vGenp1.SetPtEtaPhiM(genp1->pt, genp1->eta, genp1->phi, genp1->mass);
		if(vec[1].DeltaR(vGenp1)<dR) return 1;
	      }
	    }
	    return -1;
	  }
	}
      }
    }
  }
  return 0;
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
int GenLoader::findLastBoson(int iparent,int iId){
  Bool_t foundLast = kFALSE;
  int iLast = iparent;
  while (!foundLast) {
    int tmpId = findDaughterId(iLast,iId);
    if (tmpId>=0) iLast = tmpId;
    else foundLast = kTRUE;
  }
  return iLast;
}
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

    int iW = findLastBoson(j,24);
 
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
int GenLoader::isHadronicV(TGenParticle *genp,int j,int iId, TLorentzVector jet,double dR,double &vMatching, double &vSize)
{
  TLorentzVector vV,vDau1,vDau2;
  vMatching = -999.; vSize = -999.;
  double tmpVMatching(0), tmpVSize(0);
  if(abs(genp->pdgId)==iId){
    vV.SetPtEtaPhiM(genp->pt, genp->eta, genp->phi, genp->mass);
    int iV = findLastBoson(j,iId);

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
int GenLoader::ismatchedJet(TLorentzVector jet0, double dR,double &matching, double &size, int iId){
  for(int i0=0; i0 < fGens->GetEntriesFast(); i0++) {
    TGenParticle *genp0 = (TGenParticle*)((*fGens)[i0]);
    TLorentzVector mcMom; mcMom.SetPtEtaPhiM(genp0->pt,genp0->eta,genp0->phi,genp0->mass);
    if (mcMom.DeltaR(jet0) < dR) {
      if(iId == 6 && isHadronicTop(genp0,i0,jet0,dR,matching,size)==1) return 1;
      if(iId == 24 || iId == 23 || iId == 10031 || iId == 25){
        if (isHadronicV(genp0,i0,iId,jet0,dR,matching,size)==1) return 1;
      }
    }
  }
  return 0;
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
        int iL0 = findLastBoson(i0,iId);
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
	int iW0 = findLastBoson(i0,24);
	for(int k0 = 0; k0 < fGens->GetEntriesFast(); k0++) {
	  TGenParticle *dau0 = (TGenParticle*)((*fGens)[k0]);
	  TGenParticle *ele0 = findDaughter(iW0, 11); TGenParticle *muo0 = findDaughter(iW0, 13);
	  if(k0==iW0 && (ele0 || muo0)){
	    for(int i1=iW0+1; i0 < fGens->GetEntriesFast(); i1++) {
	      TGenParticle *genp1 = (TGenParticle*)((*fGens)[i1]);
	      if(genp1->pdgId == iId) {
		int iW1 = findLastBoson(i1,24);
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
    
  }
  fBosonPt = pbosonPt;
  fBosonPhi = pbosonPhi;
  fBosonMass = pbosonMass;
  fBosonEta = pbosonEta;
  fBosonPdgId = iId;
}
