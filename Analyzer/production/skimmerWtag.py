#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
from optparse import OptionParser
from submitZprime import samplesDict
from skimmer import getFilesRecursively

import ROOT
import json
CMSSW = os.environ['CMSSW_BASE']
PT_CUT = 200.
#94X
fLCSV = 0.5803
fMCSV = 0.8838
fTCSV = 0.9693
#Pu
fpuDir = "root://cmsxrootd.fnal.gov//store/group/lpcbacon/dazsle/zprimebits-v12.07-Pu/hadd/"
fDataDir = CMSSW+"/src/BaconAnalyzer/Analyzer/data/"
fpuData = fDataDir+"pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root"

# msd correction
def setupMassCorrection():
    fcorrGEN = ROOT.TF1("corrGEN","[0]+[1]*pow(x*[2],-[3])",200,3500);
    fcorrRECO_cen = ROOT.TF1("corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200,3500);
    fcorrRECO_for = ROOT.TF1("corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200,3500);
    fcorrGEN.SetParameter(0,1.00626)
    fcorrGEN.SetParameter(1, -1.06161)
    fcorrGEN.SetParameter(2,0.0799900)
    fcorrGEN.SetParameter(3,1.20454)
    fcorrRECO_cen.SetParameter(0,1.09302)
    fcorrRECO_cen.SetParameter(1,-0.000150068)
    fcorrRECO_cen.SetParameter(2,3.44866e-07)
    fcorrRECO_cen.SetParameter(3,-2.68100e-10)
    fcorrRECO_cen.SetParameter(4,8.67440e-14)
    fcorrRECO_cen.SetParameter(5,-1.00114e-17)
    fcorrRECO_for.SetParameter(0,1.27212)
    fcorrRECO_for.SetParameter(1,-0.000571640)
    fcorrRECO_for.SetParameter(2,8.37289e-07)
    fcorrRECO_for.SetParameter(3,-5.20433e-10)
    fcorrRECO_for.SetParameter(4,1.45375e-13)
    fcorrRECO_for.SetParameter(5,-1.50389e-17)
    return fcorrGEN,fcorrRECO_cen,fcorrRECO_for

# 2018 k-factors for now
def setupkFactors(iPt,iType,iFilename=fDataDir+"kfactors.root"):
    iPtMin=150; iPtMax=1000;
    f_kfactors = ROOT.TFile.Open(iFilename)
    hQCD_Z = f_kfactors.Get('ZJets_012j_NLO/nominal')
    hQCD_W = f_kfactors.Get('WJets_012j_NLO/nominal')
    hLO_Z = f_kfactors.Get('ZJets_LO/inv_pt')
    hLO_W = f_kfactors.Get('WJets_LO/inv_pt')
    hEWK_Z = f_kfactors.Get('EWKcorr/Z')
    hEWK_W = f_kfactors.Get('EWKcorr/W')
    hQCD_Z.SetDirectory(0)
    hQCD_W.SetDirectory(0)
    hLO_Z.SetDirectory(0)
    hLO_W.SetDirectory(0)
    hEWK_Z.SetDirectory(0)
    hEWK_W.SetDirectory(0)
    f_kfactors.Close()
    iPtMin = hQCD_Z.GetBinCenter(1);
    iPtMax = hQCD_Z.GetBinCenter(hQCD_Z.GetNbinsX())
    if iPt < iPtMin: iPt = iPtMin
    if iPt > iPtMax: iPt = iPtMax
    hEWK_Z.Divide(hQCD_Z);
    hEWK_W.Divide(hQCD_W);
    hQCD_Z.Divide(hLO_Z);
    hQCD_W.Divide(hLO_W);
    if iType == 0: #VectorDiJet                                                                                                                                             
        iQCDKF = hQCD_Z.GetBinContent(hQCD_Z.FindBin(bosonpt));
        ivjetsKF = DY_SF*iQCDKF;
    elif iType == 1: #Wjets                                                                                                                                                       
        iQCDKF = hQCD_W.GetBinContent(hQCD_W.FindBin(bosonpt));
        iEWKKF = hEWK_W.GetBinContent(hEWK_W.FindBin(bosonpt));
        ivjetsKF = W_SF*iEWKKF*iQCDKF;
        wscale=[1.0,1.0,1.0,1.20,1.25,1.25,1.0];
        ptscale=[0, 500, 600, 700, 800, 900, 1000,3000];
        ptKF=1.
        for i in range(0, len(ptscale)):
            if iPt > ptscale[i] and iPt<ptscale[i+1]:  ptKF=wscale[i]
        ivjetsKF = W_SF*iEWKKF*iQCDKF*ptKF;
    elif iType == 2: #DYJets                                                                                                                                                    
        iEWKKF = hEWK_Z.GetBinContent(hEWK_Z.FindBin(bosonpt));
        ivjetsKF = DY_SF*iEWKKF;
    else:
        ivjetsKF = 1;
    return ivjetsKF;

def setuph2ddt(ifilename=fDataDir+"GridOutput_v13.root",iddt="Rho2D"):
    f_h2ddt = ROOT.TFile.Open(ifilename)
    if "Rho2D" in iddt:
        ltrans_h2ddt = f_h2ddt.Get("Rho2D");
    else:
        lHTmp= f_h2ddt.Get(iddt);
        ltrans_h2ddt = getDDT(lHTmp,0.05)
    ltrans_h2ddt.SetDirectory(0)
    f_h2ddt.Close()
    return ltrans_h2ddt

def getN2DDT(iMass,iPt,ih2ddt):
    iRho = 2.*math.log(iMass/iPt);
    lRho = ih2ddt.GetXaxis().FindBin(iRho);
    lPt = ih2ddt.GetYaxis().FindBin(iPt);
    if iRho > ih2ddt.GetXaxis().GetBinUpEdge( ih2ddt.GetXaxis().GetNbins() ): lRho = ih2ddt.GetXaxis().GetNbins();
    if iRho < ih2ddt.GetXaxis().GetBinLowEdge( 1 ): lRho = 1;
    if iPt > ih2ddt.GetYaxis().GetBinUpEdge( ih2ddt.GetYaxis().GetNbins() ): lPt = ih2ddt.GetYaxis().GetNbins();
    if iPt < ih2ddt.GetYaxis().GetBinLowEdge( 1 ): lRho = 1;
    return ih2ddt.GetBinContent(lRho,lPt);

# 2016 puweights
def setuppuw(ifilename=fDataDir+"puWeights_All.root"):
    f_pu = ROOT.TFile.Open(ifilename)
    lpuw = f_pu.Get("puw")
    lpuw_up = f_pu.Get("puw_p")
    lpuw_down = f_pu.Get("puw_m")
    lpuw.SetDirectory(0)
    lpuw_up.SetDirectory(0)
    lpuw_down.SetDirectory(0)
    f_pu.Close()
    return lpuw,lpuw_up,lpuw_down

# 2017 puweights
def setuppuw2017(iSample):
    print fpuDir+'/'+iSample
    f_puMC = ROOT.TFile.Open(fpuDir+'/'+iSample)
    lpuMC= f_puMC.Get("Pu")
    lpuMC.Scale(1/lpuMC.Integral())
    lpuMC.SetDirectory(0)
    f_puMC.Close()
    f_puData = ROOT.TFile.Open(fpuData)
    lpuData= f_puData.Get("pileup")
    lpuData.Scale(1/lpuData.Integral())
    lpuData.SetDirectory(0)
    f_puData.Close()
    lpuData.Divide(lpuMC)
    f_puData = ROOT.TFile.Open(fpuData)
    lpuData_up= f_puData.Get("pileup_plus")
    lpuData_up.Scale(1/lpuData_up.Integral())
    lpuData_up.SetDirectory(0)
    lpuData_down= f_puData.Get("pileup_minus")
    lpuData_down.Scale(1/lpuData_down.Integral())
    lpuData_down.SetDirectory(0)
    f_puData.Close()
    lpuData_up.Divide(lpuMC)
    lpuData_down.Divide(lpuMC)
    return lpuData,lpuData_up,lpuData_down

def correctEff(iEff,iX,iY,iType=1,iName=""):
    if iType == 1:
        lweight = iEff.GetBinContent(iEff.FindBin(iX,iY))
        lweightUp = lweight + iEff.GetBinError(iEff.FindBin(iX,iY))
        lweightDown =lweight - iEff.GetBinError(iEff.FindBin(iX,iY))
    else:
        for etaKey, values in sorted(iEff[iName]["abseta_pt"].iteritems()):
            if float(etaKey[8:12]) < iY and float(etaKey[13:17]) > iY:
                for ptKey, result in sorted(values.iteritems()) :
                    if float(ptKey[4:9]) < iX and float(ptKey[10:15]) > iX:
                        lweight = result["value"]
                        lweightUp = result["value"] + result["error"]
                        lweightDown = result["value"] - result["error"]
                    else:
                        lweight = 1; lweightDown = 1; lweightUp = 1;
            else: lweight = 1; lweightDown = 1; lweightUp = 1;
    if lweight <= 0 or lweightDown <= 0 or lweightUp <= 0:
        lweight = 1
        lweightDown = 1
        lweightUp = 1
    return lweight,lweightUp,lweightDown

class miniTreeProducer:
    def __init__(self, isMc, isPu, ofile, otree, ifile, itree, ijet = 'AK8'):
        self.isMc = isMc
        self.ibase = os.path.basename( ifile)
        self.Pu = False
        if self.isMc is True:
            self.cutFormula = "1==1"
            if isPu:
                self.Pu = True
                print 'loading puweight from file '
                self.puw,self.puw_up,self.puw_down = setuppuw2017(self.ibase.replace('_1000pb_weighted',''))
            else:
                self.puw, self.puw_up, self.puw_down = setuppuw()
        else:
            self.cutFormula = "(triggerBits&4)&&passJson"
        self.ofile = ofile
        self.otree = otree
        self.ifile = ifile
        self.itree = itree
        self.jet = ijet
        self.corrGEN,self.corrRECO_cen,self.corrRECO_for = setupMassCorrection()

    def correct(self,iEta,iPt,iMass):
        genCorr  = 1.
        recoCorr = 1.
        genCorr =  self.corrGEN.Eval( iPt )
        if( abs(iEta)  < 1.3 ): recoCorr = self.corrRECO_cen.Eval( iPt )
        else: recoCorr = self.corrRECO_for.Eval( iPt )
        #print recoCorr*genCorr
        return iMass*recoCorr*genCorr

    def runProducer(self,ih2ddt,fmutrig_eff,fmuid_eff,fmuiso_eff):

        self.Puppijet0_N2 = array('f', [-100.0])
        self.Puppijet0_Tau21 = array('f', [-100.0])
        self.Puppijet0_N2DDT = array('f', [-100.0])
        self.Puppijet0_doublecsv = array('f', [-100.0])
        self.Puppijet0_pt = array('f', [-100.0])
        self.Puppijet0_msd = array('f', [-100.0])
        self.Puppijet0_vMatching = array('f', [-100.0])
        self.Puppijet0_isHadronicV = array('f', [-100.0])

        self.AK4Puppijet0_pt = array('f', [-100.0])
        self.AK4Puppijet0_eta = array('f', [-100.0])
        self.AK4Puppijet0_phi = array('f', [-100.0])
        self.AK4Puppijet0_csv = array('f', [-100.0])
        self.AK4Puppijet0_mass = array('f', [-100.0])
        self.AK4Puppijet1_pt = array('f', [-100.0])
        self.AK4Puppijet1_eta = array('f', [-100.0])
        self.AK4Puppijet1_phi = array('f', [-100.0])
        self.AK4Puppijet1_csv = array('f', [-100.0])
        self.AK4Puppijet1_mass = array('f', [-100.0])
        self.AK4Puppijet2_pt = array('f', [-100.0])
        self.AK4Puppijet2_eta = array('f', [-100.0])
        self.AK4Puppijet2_phi = array('f', [-100.0])
        self.AK4Puppijet2_csv = array('f', [-100.0])
        self.AK4Puppijet2_mass = array('f', [-100.0])
        self.AK4Puppijet3_pt = array('f', [-100.0])
        self.AK4Puppijet3_eta = array('f', [-100.0])
        self.AK4Puppijet3_phi = array('f', [-100.0])
        self.AK4Puppijet3_csv = array('f', [-100.0])
        self.AK4Puppijet3_mass = array('f', [-100.0])

        self.vmuoLoose0_pt = array('f', [-100.0])
        self.vmuoLoose0_eta = array('f', [-100.0])
        self.vmuoLoose0_phi = array('f', [-100.0])
        self.nLooseMu = array('f', [-100.0])
        self.nTightMu = array('f', [-100.0])

        self.pfmet = array('f', [-100.0])
        self.triggerBits = array('f', [-100.0])
        self.nEvents = array('f', [-100.0])
        self.LeadingJet_MatchedHadW = array('f', [-100.0])
        self.triggerpassbb = array('f', [-100.0])
        self.scale1fb = array('f',[-100.0])
        self.weight = array('f',[-100.0])

        # branches we need
        self.otree.Branch('Puppijet0_N2', self.Puppijet0_N2, 'Puppijet0_N2/F')
        self.otree.Branch('Puppijet0_Tau21', self.Puppijet0_Tau21, 'Puppijet0_Tau21/F')
        self.otree.Branch('Puppijet0_N2DDT', self.Puppijet0_N2DDT, 'Puppijet0_N2DDT/F')
        self.otree.Branch('Puppijet0_doublecsv', self.Puppijet0_doublecsv, 'Puppijet0_doublecsv/F')
        self.otree.Branch('Puppijet0_pt', self.Puppijet0_pt, 'Puppijet0_pt/F')
        self.otree.Branch('Puppijet0_msd', self.Puppijet0_msd, 'Puppijet0_msd/F')
        self.otree.Branch('Puppijet0_vMatching', self.Puppijet0_vMatching, 'Puppijet0_vMatching/F')
        self.otree.Branch('Puppijet0_isHadronicV', self.Puppijet0_isHadronicV, 'Puppijet0_isHadronicV/F')
        self.otree.Branch('nLooseMu', self.nLooseMu, 'nLooseMu/F')
        self.otree.Branch('nTightMu', self.nTightMu, 'nTightMu/F')
        self.otree.Branch('nEvents', self.nEvents, 'nEvents/F')
        self.otree.Branch('LeadingJet_MatchedHadW', self.LeadingJet_MatchedHadW, 'LeadingJet_MatchedHadW/F')
        self.otree.Branch('scale1fb', self.scale1fb, 'scale1fb/F')
        self.otree.Branch('triggerBits', self.triggerBits, 'triggerBits/F')
        self.otree.Branch('weight', self.weight, 'weight/F')

        self.bb = ROOT.TH1F("bb", "No Cuts", 3, -0.5, 1.5)
        self.bb0 = ROOT.TH1F("bb0", "After MET", 3, -0.5, 1.5)
        self.bb1 = ROOT.TH1F("bb1", "After Tight Muon", 3, -0.5, 1.5)
        self.bb2 = ROOT.TH1F("bb2", "After Loose Muon", 3, -0.5, 1.5)
        self.bb3 = ROOT.TH1F("bb3", "After Hadronic  Jet", 3, -0.5, 1.5)
        self.bb4 = ROOT.TH1F("bb4", "After b-tagged AK4 Jet", 3, -0.5, 1.5)
        self.bb5 = ROOT.TH1F("bb5", "After Leptonic  Jet", 3, -0.5, 1.5)
        self.bb6 = ROOT.TH1F("bb6", "After VHadronic Jet", 3, -0.5, 1.5)
        self.bb7 = ROOT.TH1F("bb7", "After VMatching Jet", 3, -0.5, 1.5)
        self.bb8 = ROOT.TH1F("bb8", "After VSize Jet", 3, -0.5, 1.5)

        self.f1 = ROOT.TFile.Open(self.ifile,'read')
        self.treeMine = self.f1.Get(self.itree)
        try:
            if not self.treeMine.InheritsFrom("TTree"):
                return -1
        except:
            return -1
        
        nent = self.treeMine.GetEntries()
        fcutFormula = ROOT.TTreeFormula("cutFormula",self.cutFormula,self.treeMine)
        for i in range(0,nent):

            self.treeMine.LoadTree(i)
            selected = False
            for j in range(fcutFormula.GetNdata()):
                if (fcutFormula.EvalInstance(j)):
                    selected = True
                    break
            if not selected: continue

            if( nent/100 > 0 and i % (1 * nent/100) == 0):
                sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done here")
                sys.stdout.flush()
            self.treeMine.GetEntry(i)
            
            # jets
            lMsd = getattr(self.treeMine,"%sPuppijet0_msd"%(self.jet));
            lPt = getattr(self.treeMine,"%sPuppijet0_pt"%(self.jet));
            lEta =  getattr(self.treeMine,"%sPuppijet0_eta"%(self.jet));
            lPhi =  getattr(self.treeMine,"%sPuppijet0_phi"%(self.jet));
            lTight = getattr(self.treeMine,"%sPuppijet0_isTightVJet"%(self.jet));
            lN2 = getattr(self.treeMine,"%sPuppijet0_N2sdb1"%(self.jet));
            lDcsv = getattr(self.treeMine,"%sPuppijet0_doublecsv"%(self.jet));
            lVMatching = getattr(self.treeMine,"%sPuppijet0_vMatching"%(self.jet));
            lVHadronic = getattr(self.treeMine,"%sPuppijet0_isHadronicV"%(self.jet));
            lVSize = getattr(self.treeMine,"%sPuppijet0_vSize"%(self.jet));
            lTau21 = getattr(self.treeMine,"%sPuppijet0_tau21"%(self.jet));
            lMass = self.correct(lEta,lPt,lMsd);
            if lMass <=0: lMass = 0.01
            pRho = 2.*math.log(lMass/lPt);
            lN2DDT = lN2 - getN2DDT(lMass,lPt,ih2ddt)
            lConeSize = 0.8;
            if 'CA15' in self.jet: lConeSize = 1.5;

            self.Puppijet0_N2[0] = lN2
            self.Puppijet0_Tau21[0] = lTau21
            self.Puppijet0_N2DDT[0] = lN2DDT
            self.Puppijet0_doublecsv[0] = lDcsv
            self.Puppijet0_pt[0] = lPt
            self.Puppijet0_msd[0] = lMass
            self.Puppijet0_vMatching[0] = lVMatching
            self.Puppijet0_isHadronicV[0] = lVHadronic
                
            # little jets and muons
            self.AK4Puppijet0_pt[0] = self.treeMine.AK4Puppijet0_pt
            self.AK4Puppijet0_eta[0] = self.treeMine.AK4Puppijet0_eta
            self.AK4Puppijet0_phi[0] = self.treeMine.AK4Puppijet0_phi
            self.AK4Puppijet0_csv[0] = self.treeMine.AK4Puppijet0_csv
            self.AK4Puppijet1_pt[0] = self.treeMine.AK4Puppijet1_pt
            self.AK4Puppijet1_eta[0] = self.treeMine.AK4Puppijet1_eta
            self.AK4Puppijet1_phi[0] = self.treeMine.AK4Puppijet1_phi
            self.AK4Puppijet1_csv[0] = self.treeMine.AK4Puppijet1_csv
            self.AK4Puppijet2_pt[0] = self.treeMine.AK4Puppijet2_pt
            self.AK4Puppijet2_eta[0] = self.treeMine.AK4Puppijet2_eta
            self.AK4Puppijet2_phi[0] = self.treeMine.AK4Puppijet2_phi
            self.AK4Puppijet2_csv[0] = self.treeMine.AK4Puppijet2_csv
            self.AK4Puppijet3_pt[0] = self.treeMine.AK4Puppijet3_pt
            self.AK4Puppijet3_eta[0] = self.treeMine.AK4Puppijet3_eta
            self.AK4Puppijet3_phi[0] = self.treeMine.AK4Puppijet3_phi
            self.AK4Puppijet3_csv[0] = self.treeMine.AK4Puppijet3_csv
            
            self.pfmet[0] =self.treeMine.pfmet
            self.triggerBits[0] = self.treeMine.triggerBits&4
            self.nEvents[0] = self.f1.NEvents.GetBinContent(1)
            self.scale1fb[0] = self.treeMine.scale1fb

            # muon weights
            if self.isMc is True:
                mutrigweight = 1; mutrigweightUp = 1; mutrigweightDown = 1;
                muidweight = 1;muidweightUp = 1; muidweightDown = 1;
                muisoweight = 1; muisoweightUp = 1; muisoweightDown = 1;
                if self.treeMine.nmuLoose > 0:
                    mutrigweight,mutrigweightUp,mutrigweightDown = correctEff(fmutrig_eff,self.treeMine.vmuoLoose0_pt,abs(self.treeMine.vmuoLoose0_eta))
                    muidweight,muidweightUp,muidweightDown = correctEff(fmuid_eff,self.treeMine.vmuoLoose0_pt,abs(self.treeMine.vmuoLoose0_eta),2,"NUM_SoftID_DEN_genTracks")
                    muisoweight,muisoweightUp,muisoweightDown = correctEff(fmuiso_eff,self.treeMine.vmuoLoose0_pt,abs(self.treeMine.vmuoLoose0_eta),2,"NUM_LooseRelIso_DEN_LooseID")
            
            #puweight                
            if self.isMc is True:
                nPuForWeight = min(self.treeMine.npu,49.5)
                puweight = self.treeMine.puWeight
                if self.Pu:
                    puweight = self.puw.GetBinContent(self.puw.FindBin(self.treeMine.npu));
                    if puweight<= 0: puweight =1;
                puweight_up = self.puw_up.GetBinContent(self.puw_up.FindBin(nPuForWeight))
                puweight_down = self.puw_down.GetBinContent(self.puw_down.FindBin(nPuForWeight))

            #kfactor
            if 'WJets' in self.ifile:
                vjetsKF = setupkFactors(2,self.treeMine.fBosonPt)
            elif 'DYJets' in self.ifile:
                vjetsKF = setupkFactors(3,self.treeMine.fBosonPt)
            else:
                vjetsKF = 1;

            # final weight
            if self.isMc is False:
                self.weight[0] = 1
            if self.isMc is True:
                self.weight[0] = puweight*self.treeMine.scale1fb*vjetsKF*mutrigweight*muidweight*muisoweight
                #print 'puweight scale1fb vjets ltrig lid liso ',puweight,self.treeMine.scale1fb,vjetsKF,mutrigweight,muidweight,muisoweight
                #self.weight[0] = puweight*vjetsKF*mutrigweight # scale1fb is useless if sample not Norm

            # evt cuts
            self.triggerpassbb[0] = 1.
            self.bb.Fill(self.triggerpassbb[0])

            # hadronic W matching
            if lVHadronic == 1 and lVMatching < lConeSize and lVMatching > 0.0:
                lMatchedHadW = 1.
            else:
                lMatchedHadW = 0.
            self.LeadingJet_MatchedHadW[0]  = lMatchedHadW;

            if lVHadronic == 1:
                self.bb6.Fill(self.triggerpassbb[0])
            if lVMatching < lConeSize:
                self.bb7.Fill(self.triggerpassbb[0])
            if lVSize < lConeSize:
                self.bb8.Fill(self.triggerpassbb[0])

            # met
            if self.treeMine.pfmet > 40:
                self.bb0.Fill(self.triggerpassbb[0])
            else:
                continue

            # number of loose muons
            lnLooseMuons = 0
            if self.treeMine.nmuLoose == 1:
                if self.treeMine.vmuoLoose0_pt > 20 and abs(self.treeMine.vmuoLoose0_eta) < 2.4:
                    lnLooseMuons+= 1
            if self.treeMine.nmuLoose == 2:
                if self.treeMine.vmuoLoose0_pt > 20 and abs(self.treeMine.vmuoLoose0_eta) < 2.4:
                    lnLooseMuons+= 1
                if self.treeMine.vmuoLoose1_pt > 20 and abs(self.treeMine.vmuoLoose1_eta) < 2.4:
                    lnLooseMuons+= 1
            self.nLooseMu[0] = lnLooseMuons

            # number of tight muons
            lnTightMuons = 0
            lmuPt = 0; lmuEta = 0; lmuPhi = 0;
            if self.treeMine.nmuTight == 1:
                if self.treeMine.vmuoLoose0_pt > 53 and abs(self.treeMine.vmuoLoose0_eta) < 2.1:
                    lnTightMuons+= 1
                    lmuPt = self.treeMine.vmuoLoose0_pt
                    lmuEta = self.treeMine.vmuoLoose0_eta
                    lmuPhi = self.treeMine.vmuoLoose0_phi
            if self.treeMine.nmuTight == 2:
                if self.treeMine.vmuoLoose0_pt > 53 and abs(self.treeMine.vmuoLoose0_eta) < 2.1:
                    lnTightMuons+= 1
                    lmuPt = self.treeMine.vmuoLoose0_pt
                    lmuEta = self.treeMine.vmuoLoose0_eta
                    lmuPhi = self.treeMine.vmuoLoose0_phi
                if self.treeMine.vmuoLoose1_pt > 53 and abs(self.treeMine.vmuoLoose1_eta) < 2.1:
                    lnTightMuons+= 1
                    lmuPt = self.treeMine.vmuoLoose1_pt
                    lmuEta = self.treeMine.vmuoLoose1_eta
                    lmuPhi = self.treeMine.vmuoLoose1_phi
            self.nTightMu[0] = lnTightMuons
            self.vmuoLoose0_pt[0] = lmuPt
            self.vmuoLoose0_eta[0] = lmuEta
            self.vmuoLoose0_phi[0] = lmuPhi

            if lnTightMuons == 1:
                self.bb1.Fill(self.triggerpassbb[0])
                if lnLooseMuons == 1:
                    self.bb2.Fill(self.triggerpassbb[0])
            if lnTightMuons != 1 or lnLooseMuons != 1:
                continue

            # hadronic jet
            lHadronicJet = 0.
            lDPhi_Jet_TightMuon = lPhi - lmuPhi
            if lDPhi_Jet_TightMuon >= math.pi:
                lDPhi_Jet_TightMuon -= 2.*math.pi
            elif lDPhi_Jet_TightMuon < -math.pi:
                lDPhi_Jet_TightMuon += 2.*math.pi
            lDR_Jet_TightMuon = math.sqrt((lEta - lmuEta)*(lEta - lmuEta)+lDPhi_Jet_TightMuon*lDPhi_Jet_TightMuon)
            if abs(lEta) < 2.4 and lDR_Jet_TightMuon > 1.0:
                lHadronicJet = 1
            if lHadronicJet != 1:
                continue
            else:
                self.bb3.Fill(self.triggerpassbb[0])

            # at least one b-tagged AK4 jet, isolated from jet and muon
            lnBTaggedAK4Jet = 0
            for i0 in range(0,4):
                lDR_TightMuon_Matched = 0
                lDR_Jet_Matched = 0
                lAK4Pt = getattr(self.treeMine,"AK4Puppijet%i_pt"%i0);    
                lAK4Phi = getattr(self.treeMine,"AK4Puppijet%i_phi"%i0);
                lAK4Eta = getattr(self.treeMine,"AK4Puppijet%i_eta"%i0);
                lAK4Csv = getattr(self.treeMine,"AK4Puppijet%i_csv"%i0);
                
                # (AK4,tight muon): DPhi,DR      
                lDPhi_AK4Jet_TightMuon = lAK4Phi - lmuPhi
                if lDPhi_AK4Jet_TightMuon >= math.pi:
                    lDPhi_AK4Jet_TightMuon -= 2.*math.pi
                elif lDPhi_AK4Jet_TightMuon < -math.pi:
                    lDPhi_AK4Jet_TightMuon += 2.*math.pi
                lDR_AK4Jet_TightMuon = math.sqrt((lAK4Eta - lmuEta)*(lAK4Eta - lmuEta)+lDPhi_AK4Jet_TightMuon*lDPhi_AK4Jet_TightMuon)
                if lDR_AK4Jet_TightMuon < 0.3:
                    lDR_TightMuon_Matched = 1                                     

                # (AK4,jet): DPhi,DR
                lDPhi_AK4Jet_Jet = lAK4Phi - lPhi
                if lDPhi_AK4Jet_Jet >= math.pi:
                    lDPhi_AK4Jet_Jet -= 2.*math.pi
                elif lDPhi_AK4Jet_Jet < -math.pi:
                    lDPhi_AK4Jet_Jet += 2.*math.pi
                lDR_AK4Jet_Jet = math.sqrt((lAK4Eta - lEta)*(lAK4Eta - lEta)+lDPhi_AK4Jet_Jet*lDPhi_AK4Jet_Jet)
                if lDR_AK4Jet_Jet < lConeSize:
                    lDR_Jet_Matched = 1

                if abs(lAK4Eta) < 2.4 and lAK4Csv > fMCSV and lDR_TightMuon_Matched ==0 and lDR_Jet_Matched == 0 and lAK4Pt > 30.:
                    lnBTaggedAK4Jet += 1

            if lnBTaggedAK4Jet > 0:
                self.bb4.Fill(self.triggerpassbb[0])
            else:
                continue

            # leptonic W
            lPassingLeptonicW = 0.
            lpxMet = self.treeMine.pfmet*math.cos(self.treeMine.pfmetphi)
            lpyMet = self.treeMine.pfmet*math.sin(self.treeMine.pfmetphi)
            lpxMuon = lmuPt*math.cos(lmuPhi)
            lpyMuon = lmuPt*math.sin(lmuPhi)
            if math.sqrt((lpxMet+lpxMuon)*(lpxMet+lpxMuon) + (lpyMet+lpyMuon)*(lpyMet+lpyMuon)) > 200.:
                lPassingLeptonicW = 1.
            if lPassingLeptonicW > 0:
                self.bb5.Fill(self.triggerpassbb[0])
            else: 
                continue

            # final selection
            #print 'MET %f, tightMu %i, looseMu %i, hadronicjet %i, btagAK4 %i, leptonicW %f, json %i, triggerbits %i'%(self.pfmet[0],self.nTightMu[0],self.nLooseMu[0],HadronicJet,nBTaggedAK4Jet,PassingLeptonicW,self.treeMine.passJson,self.treeMine.triggerBits&4)
            if self.treeMine.pfmet > 40 and lnTightMuons == 1 and lnLooseMuons == 1 and lHadronicJet == 1 and lnBTaggedAK4Jet > 0 and lPassingLeptonicW > 0:
                self.otree.Fill()

        self.f1.Close()

def main(options,args):

    DataDir = options.idir
    OutDir = options.odir

    try:
        samples = samplesDict[options.sample]
        tags = []
        for sample in samples:
            tags.append([sample, 0])
    except KeyError:
        tags = [[options.sample,0]]

    ftrans_h2ddt = setuph2ddt(fDataDir+options.ddt,options.iddt)

    f_mutrig = ROOT.TFile.Open(fDataDir+"EfficienciesAndSF_RunBtoF_Nov17Nov2017.root", "read")
    fmutrig_eff = f_mutrig.Get("Mu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
    fmutrig_eff.Sumw2()
    fmutrig_eff.SetDirectory(0)
    f_mutrig.Close()
    
    with open(fDataDir+"RunBCDEF_data_ID.json") as ID_input_file:
        fmuid_eff = json.load(ID_input_file)
        
    with open(fDataDir+"RunBCDEF_data_ISO.json") as ISO_input_file:
        fmuiso_eff = json.load(ISO_input_file)

    for i in range(len(tags)):
        filesToConvert = getFilesRecursively(DataDir,tags[i][0],None,None)
        print "files To Convert = ",filesToConvert

        for iFile in filesToConvert:
            basename = os.path.basename( iFile )
            oFile =  ROOT.TFile.Open(OutDir+'/'+basename, 'recreate')
            oFile.cd()
            oTree =  ROOT.TTree('otree2', 'otree2')
            prod = miniTreeProducer(options.isMc, options.isPu,oFile,oTree, iFile, options.itree, options.jet)
            prod.runProducer(ftrans_h2ddt,fmutrig_eff,fmuid_eff,fmuiso_eff)
            oFile.cd()
            oFile.Write()
            oFile.Close()
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-f", "--pathIn", dest="inputFile", help="inputFile path")
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with bacon bits', metavar='idir')
    parser.add_option('-o','--odir', dest='odir', default = 'skim/',help='directory to write skimmed backon bits', metavar='odir')
    parser.add_option("--isMc", dest="isMc", default=False, action='store_true',help="MC or data")
    parser.add_option("--isPu", dest="isPu", default=False, action='store_true',help="2017 pu or old")
    parser.add_option("--ddt", type=str, default='GridOutput_v13.root', help="ddt")
    parser.add_option("--iddt", type=str, default='Rho2D', help="iddt")
    parser.add_option('--jet', dest='jet', default='AK8', help='jet type')
    parser.add_option('--itree', dest='itree', default='otree', help='itree name')
    parser.add_option('-s','--sample',dest="sample", default="All",type='string', help="samples to produce")

    (options, args) = parser.parse_args()

    main(options,args)
