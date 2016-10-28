//===============================================================================================
// 
//
//________________________________________________________________________________________________

//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting styles
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TH1D.h>                     // 1D histogram class
#include <TLorentzVector.h>           // 4-vector class
#include <vector>                     // STL vector class
#include <iostream>                   // standard I/O
#include <iomanip>                    // functions to format standard I/O
#include <fstream>                    // functions for file I/O
#include <string>                     // C++ string class
#include <cmath>                      // C++ math library
#include <cassert>

#include "CPlot.hh"                   // helper class for plots
#include "KStyle.hh"
#include "CSample.hh"
#include "ZprimeBitsLoader.hh"

//#endif

using namespace std;

//Object Processors                                                                                                                                                                                      
ZprimeBitsLoader       *fBits      = 0;

//=== FUNCTION DECLARATIONS ======================================================================================

// make "standard" plot
void makePlot(TCanvas *c, const string outname, const string xlabel, const string ylabel,
              const vector<TH1D*>& histv, const vector<CSample*>& samplev, TH1D* hExp, TH1D* hPull,
              const bool doBlind, const double lumi, const bool doLogy=false, const double legdx=0, const double legdy=0,
              const double ymin=-1, const double ymax=-1);
TH1D* makePullHist(TH1D* hData, TH1D* hMC, const string name, const bool doBlind);
float CalcSig1(TH1D*sig1, TH1D*bkg1);

//=== MAIN MACRO =================================================================================================

void plotggHbb(const string selection, const string algo, const string jet, float cut, float csv)
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================

  const bool doBlind = false;

  // Create output directory 
  const string outputDir("ggHbbPlots/"+selection+"_"+algo+"_"+jet);
  gSystem->mkdir(outputDir.c_str(), true);
  CPlot::sOutDir = outputDir;

  //
  // Samples
  // Note: macro assumes samplev[0] is data
  //
  vector<CSample*> samplev;

  samplev.push_back(new CSample("data",0,0));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/JetHTRun2016B.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/JetHTRun2016C.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/JetHTRun2016D.root");
  samplev.push_back(new CSample("QCD", kMagenta - 10, kMagenta - 10));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/QCD_HT100to200.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/QCD_HT200to300.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/QCD_HT300to500.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/QCD_HT500to700.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cvernier/zprimebits/QCD_HT700to1000_13TeV.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/QCD_HT1000to1500.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/QCD_HT1500to2000.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/QCD_HT2000toInf.root");
  samplev.push_back(new CSample("W+jets",kGreen - 10,kGreen - 10));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/W.root");
  samplev.push_back(new CSample("Z+jets", kCyan - 9, kCyan - 9));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/DY.root");
  samplev.push_back(new CSample("Single Top",kRed - 9,kRed - 9));
  samplev.back()->fnamev.push_back("/eos/user/c/cvernier/zprimebits/ST_tW_antitop_5f_inclusiveDecays_13TeV.root");
  samplev.back()->fnamev.push_back("/eos/user/c/cvernier/zprimebits/ST_tW_top_5f_inclusiveDecays_13TeV.root");
  samplev.push_back(new CSample("t#bar{t}",kOrange - 3,kOrange - 3));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv4/TTbar_madgraphMLM.root");
  samplev.push_back(new CSample("m_{Z'}=100 GeV",2,2));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv3/VectorDiJet1Jet_M100.root");
  samplev.push_back(new CSample("m_{Z'}=125 GeV",4,4));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv3/VectorDiJet1Jet_M125.root");
  samplev.push_back(new CSample("m_{Z'}=150 GeV",6,6));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv3/VectorDiJet1Jet_M150.root");
  samplev.push_back(new CSample("m_{Z'}=200 GeV",7,7));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv3/VectorDiJet1Jet_M200.root");
  samplev.push_back(new CSample("m_{Z'}=250 GeV",8,8));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv3/VectorDiJet1Jet_M250.root");
  samplev.push_back(new CSample("m_{Z'}=300 GeV",3,3));
  samplev.back()->fnamev.push_back("/eos/user/c/cmantill/VectorDiJet1Jetv3/VectorDiJet1Jet_M300.root");


  // integrated luminosity to scale MC
  const double LUMI = 12.89;

  // histograms for various corrections
  const string cmssw_base = getenv("CMSSW_BASE");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  //
  // Declare histograms
  //
  char hname[100];
  vector<TH1D*> hFatJetMassv, hFatJetTau21v, hFatJetTau21DDTv, hFatJetPtv;     
  vector<double> neventsv;
  
  for(unsigned int isam=0; isam<12; isam++) {
    sprintf(hname,"hFatJetMass_%i",isam);     hFatJetMassv.push_back(new TH1D(hname,"",40,0,400));         hFatJetMassv[isam]->Sumw2();
    sprintf(hname,"hFatJetPt_%i",isam);       hFatJetPtv.push_back(new TH1D(hname,"",20,250,1000));        hFatJetPtv[isam]->Sumw2();
    sprintf(hname,"hFatJetTau21_%i",isam);    hFatJetTau21v.push_back(new TH1D(hname,"",25,0,1));          hFatJetTau21v[isam]->Sumw2();
    sprintf(hname,"hFatJetTau21DDT_%i",isam); hFatJetTau21DDTv.push_back(new TH1D(hname,"",25,0,1));       hFatJetTau21DDTv[isam]->Sumw2();
    neventsv.push_back(0);
  }

  TH1D *hFatJetMassMC       = (TH1D*)hFatJetMassv[0]      ->Clone("hFatJetMassMC");
  TH1D *hFatJetPtMC         = (TH1D*)hFatJetPtv[0]      ->Clone("hFatJetPtMC");
  TH1D *hFatJetTau21MC      = (TH1D*)hFatJetTau21v[0]     ->Clone("hFatJetTau21MC");
  TH1D *hFatJetTau21DDTMC   = (TH1D*)hFatJetTau21DDTv[0]  ->Clone("hFatJetTau21DDTMC");

  double neventsMC=0;

  TFile *infile=0;
  TTree *intree=0;
  
  // Loop over samples
 
  for(unsigned int isam=0; isam<12; isam++) {
    CSample *sample  = samplev[isam];
    cout << "Sample: " << sample->label << endl;
    bool isData    = (isam==0);

    for(unsigned int ifile=0; ifile<sample->fnamev.size(); ifile++) {
      string infilename = sample->fnamev[ifile];
      cout << " ==> Processing " << infilename << "... "; cout.flush();
      infile = new TFile(infilename.c_str()); assert(infile);
      intree = (TTree*)infile->Get("Events"); assert(intree);
      fBits  = new ZprimeBitsLoader(intree,algo,jet);

      double nevts=0;
      int noweight=0;

      for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
      //for(unsigned int ientry=0; ientry<10000; ientry++) {

        intree->GetEntry(ientry);
	// BLINDING POLICY
        //if(isData && ientry % 5 != 0) continue;

	if(!fBits->selectJetAlgoAndSize(algo))      continue;
	if(!isData && fBits->metfilter!=0)          continue;
	if(!fBits->passBoostedZprimePreSelection()) continue;

	// Apply weigths
        double wgt = 1;
	wgt *= fBits->getWgt(isData,algo,LUMI);
	
	nevts += wgt;
	noweight++;
        neventsv[isam]+=wgt;

	hFatJetMassv[isam]     ->Fill(fBits->bst_jet0_msd, wgt);
        hFatJetPtv[isam]       ->Fill(fBits->bst_jet0_pt, wgt);
	hFatJetTau21v[isam]    ->Fill(fBits->bst_jet0_tau21, wgt);
        hFatJetTau21DDTv[isam] ->Fill(fBits->tau21DDT(), wgt);
      }
      if(isData && doBlind) {
        cout << endl;
      } else {
        cout << nevts << " " << noweight  <<  " " << fBits->scale1fb << endl;
      }
      delete infile;
      infile=0;
      intree=0;
    }
   
  }


  // 
  // QCD SF
  //
  double QCDSF = 1.0;
  QCDSF = (neventsv[0]-(neventsv[2]+neventsv[3]+neventsv[4]+neventsv[5]))/neventsv[1];
  std::cout << "QCDSF" << QCDSF << std::endl;

  hFatJetMassv[1]     ->Scale(QCDSF);
  hFatJetPtv[1]       ->Scale(QCDSF);
  hFatJetTau21v[1]    ->Scale(QCDSF);
  hFatJetTau21DDTv[1] ->Scale(QCDSF);

  for(unsigned int isam=1; isam<6; isam++) {
    cout << "Adding " << samplev[isam]->label <<" to MC"<<endl;
    hFatJetMassMC     ->Add(hFatJetMassv[isam]);
    hFatJetPtMC       ->Add(hFatJetPtv[isam]);
    hFatJetTau21MC    ->Add(hFatJetTau21v[isam]);
    hFatJetTau21DDTMC ->Add(hFatJetTau21DDTv[isam]);
  }
  
  neventsMC = neventsv[1]*QCDSF+neventsv[2]+neventsv[3]+neventsv[4]+neventsv[5];
  std::cout << neventsMC << std::endl;
  //
  // Make pull histograms
  //
  
  TH1D *hFatJetMassPull     = makePullHist(hFatJetMassv[0],     hFatJetMassMC,     "hFatJetMassPull",      doBlind);
  TH1D *hFatJetPtPull       = makePullHist(hFatJetPtv[0],       hFatJetPtMC,       "hFatJetPtPull",        doBlind);
  TH1D *hFatJetTau21Pull    = makePullHist(hFatJetTau21v[0],    hFatJetTau21MC,    "hFatJetTau21Pull",     doBlind);
  TH1D *hFatJetTau21DDTPull = makePullHist(hFatJetTau21DDTv[0], hFatJetTau21DDTMC, "hFatJetTau21DDTPull",  doBlind);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  ofstream txtfile;
  char txtfname[200];
  sprintf(txtfname,"%s/summary.txt",outputDir.c_str());
  txtfile.open(txtfname,std::ios_base::app);
  txtfile << setprecision(6) << fixed;
  float max = samplev.size();
  for(unsigned int isam=1; isam<8; isam++) {
    txtfile << setw(35) << samplev[isam]->label;
    txtfile << setw(15) << neventsv[isam] << endl;
  }
  txtfile << "---------------------------------------------"  << endl;
  txtfile << setw(35) << "SM Expected:" << setw(15) << neventsMC << endl;
  if(!doBlind) { txtfile << setw(35) << "Observed:" << setw(15) << neventsv[0] << endl; }

  txtfile << "QCD Scale Factor:" << QCDSF << endl;
  txtfile << "---------------------------------------------"  << endl;
  txtfile.close();

   //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0);
  c->cd(1)->SetLeftMargin(0.15);
  c->cd(1)->SetRightMargin(0.07);
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.1);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);

  char ylabel[100];

  sprintf(ylabel,"Events");
  makePlot(c, "msd", "Soft Drop Mass [GeV]", ylabel, hFatJetMassv, samplev, hFatJetMassMC, hFatJetMassPull, doBlind, LUMI, true, -0.03, -0.08,
           2e-5*(hFatJetMassMC->GetBinContent(hFatJetMassMC->GetMaximumBin()))/(hFatJetMassMC->GetBinWidth(hFatJetMassMC->GetMaximumBin())),
           4e2*(hFatJetMassMC->GetBinContent(hFatJetMassMC->GetMaximumBin()))/(hFatJetMassMC->GetBinWidth(hFatJetMassMC->GetMaximumBin())));

  sprintf(ylabel,"Events");
  makePlot(c, "pt", "p_{T} [GeV]", ylabel, hFatJetPtv, samplev, hFatJetPtMC, hFatJetPtPull, doBlind, LUMI, false, -0.03, -0.08,
           0.1, 2.1*(hFatJetPtMC->GetBinContent(hFatJetPtMC->GetMaximumBin()))/(hFatJetPtMC->GetBinWidth(hFatJetPtMC->GetMaximumBin())));

  sprintf(ylabel,"Events");
  makePlot(c, "tau21DDT", "#tau_{21}^{DDT}", ylabel, hFatJetTau21DDTv, samplev, hFatJetTau21DDTMC, hFatJetTau21DDTPull, doBlind, LUMI, false, -0.03, -0.08,
           0.1, 2.1*(hFatJetTau21DDTMC->GetBinContent(hFatJetTau21DDTMC->GetMaximumBin()))/(hFatJetTau21DDTMC->GetBinWidth(hFatJetTau21DDTMC->GetMaximumBin())));

  cout << endl;
  cout << " <> Output saved in " << outputDir << endl;
  cout << endl;


}

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makePlot(TCanvas *c, const string outname, const string xlabel, const string ylabel,
              const vector<TH1D*>& histv, const vector<CSample*>& samplev, TH1D* hExp, TH1D* hPull,
              bool doBlind, const double lumi, const bool doLogy, const double legdx, const double legdy,
              const double ymin, const double ymax)
{
  const int uncColor = kGray+3;

  // Should divide by bin width                                                                                                                                                                             
  for (int iB=1; iB<hExp->GetNbinsX()+1; ++iB) {
    float currentVal = hExp->GetBinContent(iB);
    float currentErr = hExp->GetBinError(iB);
    float binWidth = hExp->GetBinWidth(iB);
    hExp->SetBinContent(iB,currentVal/binWidth);
    hExp->SetBinError(iB,currentErr/binWidth);
  }
  for(unsigned int i=0; i<histv.size(); i++) {
    for (int iB=1; iB<histv[i]->GetNbinsX()+1; ++iB) {
      float currentVal = histv[i]->GetBinContent(iB);
      float currentErr = histv[i]->GetBinError(iB);
      float binWidth = histv[i]->GetBinWidth(iB);
      histv[i]->SetBinContent(iB,currentVal/binWidth);
      histv[i]->SetBinError(iB,currentErr/binWidth);
    }
  }

  histv[0]->SetMarkerSize(0.9);

  CPlot plot(outname.c_str(),"",xlabel.c_str(),ylabel.c_str());
  plot.AddHist1D(hExp,"E2",uncColor,1,3004);
  if(!doBlind) { plot.AddHist1D(histv[0],samplev[0]->label,"E"); }
  float max = samplev.size()-6;
  for(unsigned int i=max-1; i>=1; i--) {
  //for(unsigned int i=1; i<max; i++) {
    histv[i]->SetMarkerSize(0.7);
    plot.AddToStack(histv[i],samplev[i]->label,samplev[i]->fillcolor,samplev[i]->linecolor);
  }
  
  for(unsigned int i=max; i<histv.size(); i++) {
    plot.AddHist1D(histv[i],samplev[i]->label,"hist",samplev[i]->fillcolor,samplev[i]->linecolor);
  }

  TLegend *leg = plot.GetLegend();
  leg->SetTextSize(0.02);
  
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1} (13 TeV)",lumi);
  plot.AddTextBox(lumitext,0.66,0.99,0.95,0.925,0,kBlack);
  plot.AddTextBox("CMS",0.18,0.88,0.30,0.82,0,kBlack,62);
  plot.AddTextBox("Preliminary",0.18,0.82,0.37,0.77,0,kBlack,52);
  plot.TransLegend(legdx, legdy);

  if(doLogy) plot.SetLogy();
  if(ymin!=ymax) plot.SetYRange(ymin,ymax);

  hPull->SetMarkerSize(0.8);
  const double xmin = histv[0]->GetXaxis()->GetBinLowEdge(1);
  const double xmax = histv[0]->GetXaxis()->GetBinUpEdge(histv[0]->GetNbinsX());
  TH1D *hExpPull = (TH1D*)hExp->Clone("hExpPull");
  for (int iB=1; iB<hExpPull->GetNbinsX()+1; ++iB) {
    float currentVal = hExpPull->GetBinContent(iB);
    float currentErr = hExpPull->GetBinError(iB);
    hExpPull->SetBinContent(iB,1.);
    hExpPull->SetBinError(iB,currentErr/currentVal);
  }

  CPlot plotPull(outname.c_str(),"",xlabel.c_str(),"Pull");
  plotPull.AddHist1D(hPull,"EX0",kBlack);
  plotPull.AddHist1D(hExpPull,"E2",uncColor,1,3004);
  plotPull.SetYRange(0.,2.);
  plotPull.AddLine(xmin,1,xmax,1,kBlack,3);

  plot.Draw(c,false,"png",1);
  plot.Draw(c,false,"pdf",1);
  plotPull.Draw(c,true,"png",2);
  plotPull.Draw(c,true,"pdf",2);
}

//--------------------------------------------------------------------------------------------------
TH1D* makePullHist(TH1D* hData, TH1D* hMC, const string name, const bool doBlind)
{
  const Int_t NBINS = 5;
  Double_t edges[NBINS + 1] = {250,300,350,400,500,1000};
  TH1D *hPull = new TH1D(name.c_str(),"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  if (name=="hMETPull" || name=="hMETLogPull")
    hPull = new TH1D(name.c_str(),"",NBINS,edges);
  for(int ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    double numer = hData->GetBinContent(ibin);
    double denom = hMC->GetBinContent(ibin);
    double pull  = (denom>0) ? numer/denom : 0;
    double err   = (denom>0) ? sqrt(hData->GetBinContent(ibin))/hMC->GetBinContent(ibin) : 0;

    if(doBlind) {
      pull = 1;
      err  = 1;
    }

    hPull->SetBinContent(ibin,pull);
    hPull->SetBinError(ibin,err);
  }
  hPull->GetYaxis()->SetTitleOffset(0.42);
  hPull->GetYaxis()->SetTitleSize(0.13);
  hPull->GetYaxis()->SetLabelSize(0.10);
  hPull->GetYaxis()->SetNdivisions(5);
  hPull->GetYaxis()->CenterTitle();
  hPull->GetXaxis()->SetTitleOffset(1.2);
  hPull->GetXaxis()->SetTitleSize(0.13);
  hPull->GetXaxis()->SetLabelSize(0.12);

  return hPull;
}












