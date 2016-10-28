#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include "KStyle.hh"

void KStyle() {
  printf("\n");
  printf(" Use: MakeCanvas(name,title)\n");
  printf(" InitSubPad(pad,nPad)\n");
  printf(" InitHist(hist,xTitle,yTitle,color)\n");
  printf("\n");
  SetStyle();
}

TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor (0);
  canvas->SetBorderMode (0);
  canvas->SetBorderSize (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin (0.18);
  canvas->SetRightMargin (0.05);
  canvas->SetTopMargin (0.08);
  canvas->SetBottomMargin (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void SetStyle()
{
  TStyle *KStyle = new TStyle("K-Style","Cuz Kevin can't think of a good name");
  gStyle = KStyle;

  // Canvas
  KStyle->SetCanvasColor (0);
  KStyle->SetCanvasBorderSize(10);
  KStyle->SetCanvasBorderMode(0);
  KStyle->SetCanvasDefH (700);
  KStyle->SetCanvasDefW (700);
  KStyle->SetCanvasDefX (100);
  KStyle->SetCanvasDefY (100);

  // Pads
  KStyle->SetPadColor (0);
  KStyle->SetPadBorderSize (10);
  KStyle->SetPadBorderMode (0);
  KStyle->SetPadBottomMargin(0.13);
  KStyle->SetPadTopMargin (0.08);
  KStyle->SetPadLeftMargin (0.15);
  KStyle->SetPadRightMargin (0.05);
  KStyle->SetPadGridX (0);
  KStyle->SetPadGridY (0);
  KStyle->SetPadTickX (0);
  KStyle->SetPadTickY (0);

  // Frames
  KStyle->SetFrameFillStyle ( 0);
  KStyle->SetFrameFillColor ( 0);
  KStyle->SetFrameLineColor ( 1);
  KStyle->SetFrameLineStyle ( 0);
  KStyle->SetFrameLineWidth ( 1);
  KStyle->SetFrameBorderSize(10);
  KStyle->SetFrameBorderMode( 0);

  // Histograms
  KStyle->SetHistFillColor(2);
  KStyle->SetHistFillStyle(0);
  KStyle->SetHistLineColor(1);
  KStyle->SetHistLineStyle(0);
  KStyle->SetHistLineWidth(2);
  KStyle->SetNdivisions(505);

  // Functions
  KStyle->SetFuncColor(1);
  KStyle->SetFuncStyle(0);
  KStyle->SetFuncWidth(2);

  // Various
  KStyle->SetMarkerStyle(20);
  KStyle->SetMarkerColor(kBlack);
  KStyle->SetMarkerSize (1.2);

  KStyle->SetTitleBorderSize(0);
  KStyle->SetTitleFillColor (0);
  KStyle->SetTitleX (0.2);

  KStyle->SetTitleSize (0.055,"X");
  KStyle->SetTitleOffset(1.200,"X");
  KStyle->SetLabelOffset(0.005,"X");
  KStyle->SetLabelSize (0.050,"X");
  KStyle->SetLabelFont (42 ,"X");

  KStyle->SetStripDecimals(kFALSE);

  KStyle->SetTitleSize (0.055,"Y");
  KStyle->SetTitleOffset(1.400,"Y");
  KStyle->SetLabelOffset(0.010,"Y");
  KStyle->SetLabelSize (0.050,"Y");
  KStyle->SetLabelFont (42 ,"Y");

  KStyle->SetTextSize (0.055);
  KStyle->SetTextFont (42);

  KStyle->SetStatFont (42);
  KStyle->SetTitleFont (42);
  KStyle->SetTitleFont (42,"X");
  KStyle->SetTitleFont (42,"Y");

  KStyle->SetOptStat (0);
  return;
}
