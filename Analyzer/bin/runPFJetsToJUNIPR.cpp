// run:
// ./runPFJetsToJUNIPR input_directory input_file recluster_def label isdata

#include "../include/PerJetLoader.hh"
#include "../include/JuniprJetLoader.hh"

#include "fastjet/ClusterSequence.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace fastjet;
using namespace std;

#define PI 3.141592654

PerJetLoader      *fVJet8      = 0;
JuniprJetLoader   *fJuniprJet  = 0;


TTree* load(std::string iName) {
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}

double SignedDeltaPhi(double phi1, double phi2) {
  double dPhi = phi1-phi2;
  if (dPhi<-PI)
    dPhi = 2*PI+dPhi;
  else if (dPhi>PI)
    dPhi = -2*PI+dPhi;
  return dPhi;
}

int main(int argc, char *argv[]) {

  if (argc != 7) {    
    cout << "Wrong number of arguments" << endl;
    return 0;
  }
  
  // User input
  string input_file = argv[1];
  string output_file = argv[2];
  double recluster_def = atof(argv[3]);
  int label = atoi(argv[4]);
  string isdata = argv[5];
  const std::string lLabel = argv[6];

  bool isData;
  if(isdata.compare("data")!=0) isData = false;
  else isData = true;

  // Input, output
  stringstream s_infile;
  s_infile << input_file;

  stringstream s_outfile;
  s_outfile << "JUNIPR_format_" << output_file ;
  ofstream outfile;
  outfile.open( s_outfile.str().c_str());
  outfile << "{ \"JuniprJets\": [\n";

  // Read Tree
  const std::string lName = s_infile.str().c_str();
  TTree *lTree = load( lName );
  fVJet8     = new PerJetLoader  (lTree,"AK8Puppi","AddAK8Puppi","AK8CHS","AddAK8CHS",3, isData, lLabel);

  // Select jets from Tree
  int nsubjets;
  int jet_counter = 0;
  for(int i0 = 0; i0 < int(lTree->GetEntriesFast()); i0++) {
    fVJet8->load(i0);
    //int ifill = 0;
    for  (int i1 = 0; i1 < fVJet8->fVJets->GetEntriesFast(); i1++) {
      TJet *pVJet = (TJet*)((*fVJet8->fVJets)[i1]);
      TAddJet *pAddJet = fVJet8->getAddJet(pVJet);
      if(pVJet->pt   <=  200) continue;
      if(fabs(pVJet->eta) >  2.4) continue;
      // select pt bin
      if(pVJet->pt   <  500 || pVJet->pt > 600) continue;
      // get jet variables
      std::vector<double> vars;
      double pJet_Pt; pJet_Pt = pVJet->pt; vars.push_back(pJet_Pt);
      double pJet_Eta; pJet_Eta = pVJet->eta; vars.push_back(pJet_Eta);
      double pJet_Phi; pJet_Phi = pVJet->phi; vars.push_back(pJet_Phi);
      double pJet_Msd; pJet_Msd = pAddJet->mass_sd0; vars.push_back(pJet_Msd);
      double pJet_Tau21; pJet_Tau21 = (pAddJet->tau2/pAddJet->tau1); vars.push_back(pJet_Tau21);
      //if(ifill > 1) continue;
      //ifill += 1;
      std::vector<TPFPart*> jetPFs;
      for (auto idx : pVJet->pfCands) {
        jetPFs.push_back( (TPFPart*)(fVJet8->fPFs->At(idx)) );
      }
      vector<PseudoJet> particles;
      nsubjets = jetPFs.size();
      // Particle loop
      for (auto *pf : jetPFs) {
	double pPt, pEta, pPhi, pM;
	pPt = pf->pup * pf->pt / pVJet->pt;
	pEta = pf->eta - pVJet->eta;
	pPhi = SignedDeltaPhi(pf->phi, pVJet->phi);
	pM = pf->m;
	TLorentzVector jetPF; jetPF.SetPtEtaPhiM(pPt,pEta,pPhi,pM);
	PseudoJet particle(jetPF.Px(),jetPF.Py(),jetPF.Pz(),jetPF.E());
	particles.push_back(particle);
      } // end loop (particles)
      // Recluster jet constituents
      JetDefinition reclust_def(ee_genkt_algorithm, 7, recluster_def); // radius unused
      ClusterSequence reclust_seq(particles, reclust_def);
      vector<PseudoJet> reclust_jets = reclust_seq.exclusive_jets(1);
      if (jet_counter == 0)
	cout << "Reclustered with " << reclust_def.description() << endl;
      if(int(reclust_jets[0].constituents().size()) != nsubjets){
	cout << "Mismatch between nsubjets and reclustered jet" << endl;
      }

      if (jet_counter>0){
	outfile << ",\n";
      }

      fJuniprJet = new JuniprJetLoader(reclust_seq, label);
      fJuniprJet->write_juniprjet(outfile,vars);

      jet_counter++;
    } // end loop (jets)
  } // end loop (events)
  outfile << "\n]\n}";

  outfile.close();
  return 0;

} // end function (main)

