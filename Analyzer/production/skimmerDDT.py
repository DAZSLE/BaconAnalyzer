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

PT_CUT = 200.

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

    # make a tmp dir
    #####
    postfix = ''
    for i in range(len(tags)):
        filesToConvert = getFilesRecursively(DataDir,tags[i][0],None,None)
        print "files To Convert = ",filesToConvert

        for f in filesToConvert:
            status = sklimAdd(f,OutDir,tags[i][1])
            print status


def sklimAdd(fn,odir,mass=0):

    basename = os.path.basename( fn )

    f1 = ROOT.TFile.Open(fn,'read')
    tree = f1.Get("Events")
    try:
        if not tree.InheritsFrom("TTree"):
            return -1
    except:
        return -1
    
    ofile = ROOT.TFile.Open(odir+'/'+basename,'RECREATE')
    ofile.cd()
    f1.cd()	
    obj = ROOT.TObject
    for key in ROOT.gDirectory.GetListOfKeys():
        f1.cd()
        obj = key.ReadObj()
        print obj.GetName()
        if obj.GetName() == 'Events':
            continue
        ofile.cd()
        print key.GetName()
        obj.Write(key.GetName())

        
    
    otree2 = tree.CloneTree(0)
    otree2.SetName("otree")

    otree2.SetBranchStatus("AK8Puppijet1*",0)
    otree2.SetBranchStatus("AK8Puppijet2*",0)
    otree2.SetBranchStatus("AK8Puppijet0_ratioCA15_04",0)
    otree2.SetBranchStatus("*Puppijet0_e3*",0)
    otree2.SetBranchStatus("*Puppijet0_e4*",0)
    otree2.SetBranchStatus("*Puppijet0_e2*",0)
    otree2.SetBranchStatus("AK8Puppijet0_M2sdb1",0)
    otree2.SetBranchStatus("AK8Puppijet0_M2sdb2",0)
    otree2.SetBranchStatus("AK8Puppijet0_D2sdb1",0)
    otree2.SetBranchStatus("AK8Puppijet0_D2sdb2",0)
    otree2.SetBranchStatus("AK8Puppijet0_N2b1",0)
    otree2.SetBranchStatus("AK8Puppijet0_N2b2",0)
    otree2.SetBranchStatus("AK8Puppijet0_M2b1",0)
    otree2.SetBranchStatus("AK8Puppijet0_M2b2",0)
    otree2.SetBranchStatus("AK8Puppijet0_D2b1",0)
    otree2.SetBranchStatus("AK8Puppijet0_D2b2",0)
    otree2.SetBranchStatus("AK8Puppijet0_phi",0)
    otree2.SetBranchStatus("AK8Puppijet0_mass",0)
    otree2.SetBranchStatus("CA15Puppijet1*",0)
    otree2.SetBranchStatus("CA15Puppijet2*",0)  
    otree2.SetBranchStatus("CA15Puppijet0_ratioCA15_04",0)
    otree2.SetBranchStatus("CA15Puppijet0_M2sdb1",0)
    otree2.SetBranchStatus("CA15Puppijet0_M2sdb2",0)
    otree2.SetBranchStatus("CA15Puppijet0_D2sdb1",0)
    otree2.SetBranchStatus("CA15Puppijet0_D2sdb2",0)
    otree2.SetBranchStatus("CA15Puppijet0_N2b1",0)
    otree2.SetBranchStatus("CA15Puppijet0_N2b2",0)
    otree2.SetBranchStatus("CA15Puppijet0_M2b1",0)
    otree2.SetBranchStatus("CA15Puppijet0_M2b2",0)
    otree2.SetBranchStatus("CA15Puppijet0_D2b1",0)
    otree2.SetBranchStatus("CA15Puppijet0_D2b2",0)
    otree2.SetBranchStatus("CA15Puppijet0_N2sdb2",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_old",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JESUp",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JESDown",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JERUp",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JERDown",0)
    otree2.SetBranchStatus("nCA15Puppijets",0)
    otree2.SetBranchStatus("CA15Puppijet0_isHadronicV",0)
    otree2.SetBranchStatus("CA15Puppijet0_isTightVJet",0)
    otree2.SetBranchStatus("CA15Puppijet0_vMatching",0)
    otree2.SetBranchStatus("CA15Puppijet0_vSize",0)
    otree2.SetBranchStatus("CA15Puppijet0_partonFlavor",0)
    otree2.SetBranchStatus("CA15Puppijet0_hadronFlavor",0)
    otree2.SetBranchStatus("CA15Puppijet0_nCharged",0)
    otree2.SetBranchStatus("CA15Puppijet0_nNeutrals",0)
    otree2.SetBranchStatus("CA15Puppijet0_nParticles",0)
    otree2.SetBranchStatus("CA15Puppijet0_mass",0)
    otree2.SetBranchStatus("CA15Puppijet0_rho",0)
    otree2.SetBranchStatus("CA15Puppijet0_ptraw",0)
    otree2.SetBranchStatus("CA15Puppijet0_genpt",0)
    otree2.SetBranchStatus("CA15Puppijet0_phi",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JESUp",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JESDown",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JERUp",0)
    otree2.SetBranchStatus("CA15Puppijet0_pt_JERDown",0)
    otree2.SetBranchStatus("*Met*Corr*",0)
    otree2.SetBranchStatus("*AK4*",0)
    otree2.SetBranchStatus("selectBits",0)
    otree2.SetBranchStatus("vmu*",0)
    otree2.SetBranchStatus("vele*",0)
    otree2.SetBranchStatus("vpho*",0)
    otree2.SetBranchStatus("*tau21*",0)
    otree2.SetBranchStatus("*tau32*",0)
    otree2.SetBranchStatus("*csv*",0)
    otree2.SetBranchStatus("*doublesub*",0)
    otree2.SetBranchStatus("*CHF*",0)
    otree2.SetBranchStatus("*NHF*",0)
    otree2.SetBranchStatus("*NEMF*",0)
    otree2.SetBranchStatus("*CHS*",0)
    otree2.SetBranchStatus("runNum",0)
    otree2.SetBranchStatus("lumiSec",0)
    otree2.SetBranchStatus("metfilter",0)
    otree2.SetBranchStatus("*trigger*",0)
    otree2.SetBranchStatus("moreTriggerBits",0)
    otree2.SetBranchStatus("neleHEEP",0)
    otree2.SetBranchStatus("neleTight",0)
    otree2.SetBranchStatus("nmuTight",0)
    otree2.SetBranchStatus("nmuHighPt",0)
    otree2.SetBranchStatus("sf_*",0)
    otree2.SetBranchStatus("rho",0)
    otree2.SetBranchStatus("*Scale_*",0)
    otree2.SetBranchStatus("puppet",0)
    otree2.SetBranchStatus("puppetphi",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JESUp",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JESDown",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JERUp",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JERDown",0)
    otree2.SetBranchStatus("nAK8Puppijets",0)
    otree2.SetBranchStatus("AK8Puppijet0_isHadronicV",0)
    otree2.SetBranchStatus("AK8Puppijet0_isTightVJet",0)
    otree2.SetBranchStatus("AK8Puppijet0_vMatching",0)
    otree2.SetBranchStatus("AK8Puppijet0_vSize",0)
    otree2.SetBranchStatus("AK8Puppijet0_partonFlavor",0)
    otree2.SetBranchStatus("AK8Puppijet0_hadronFlavor",0)
    otree2.SetBranchStatus("AK8Puppijet0_nCharged",0)
    otree2.SetBranchStatus("AK8Puppijet0_nNeutrals",0)
    otree2.SetBranchStatus("AK8Puppijet0_nParticles",0)
    otree2.SetBranchStatus("AK8Puppijet0_N2sdb2",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_old",0)
    otree2.SetBranchStatus("nmu*",0)
    otree2.SetBranchStatus("nele*",0)
    otree2.SetBranchStatus("ntau*",0)
    otree2.SetBranchStatus("ispho0Tight",0)
    otree2.SetBranchStatus("npho*",0)
    otree2.SetBranchStatus("genVPt",0)
    otree2.SetBranchStatus("genVPhi",0)
    otree2.SetBranchStatus("genVMass",0)
    otree2.SetBranchStatus("genVEta",0)
    otree2.SetBranchStatus("genVPdfId",0)
    otree2.SetBranchStatus("genEleFromW",0)
    otree2.SetBranchStatus("genMuFromW",0)
    otree2.SetBranchStatus("genTauFromW",0)
    otree2.SetBranchStatus("topPt",0)
    otree2.SetBranchStatus("antitopPt",0)
    otree2.SetBranchStatus("topPtWeight",0)
    otree2.SetBranchStatus("*Met*Corr*",0)
    otree2.SetBranchStatus("*AK4*",0)
    otree2.SetBranchStatus("selectBits",0)
    otree2.SetBranchStatus("vmu*",0)
    otree2.SetBranchStatus("vele*",0)
    otree2.SetBranchStatus("vpho*",0)
    otree2.SetBranchStatus("*tau21*",0)
    otree2.SetBranchStatus("*tau32*",0)
    otree2.SetBranchStatus("*csv*",0)
    otree2.SetBranchStatus("*doublesub*",0)
    otree2.SetBranchStatus("runNum",0)
    otree2.SetBranchStatus("lumiSec",0)
    otree2.SetBranchStatus("metfilter",0)
    otree2.SetBranchStatus("triggerBits",0)
    otree2.SetBranchStatus("neleHEEP",0)
    otree2.SetBranchStatus("neleTight",0)
    otree2.SetBranchStatus("nmuTight",0)
    otree2.SetBranchStatus("nmuHighPt",0)
    otree2.SetBranchStatus("sf_*",0)
    otree2.SetBranchStatus("rho",0)
    otree2.SetBranchStatus("*Scale_*",0)
    otree2.SetBranchStatus("puppet",0)
    otree2.SetBranchStatus("puppetphi",0)
    otree2.SetBranchStatus("PDF*",0)
    otree2.SetBranchStatus("evtNum",0)
    otree2.SetBranchStatus("passJson",0)
    otree2.SetBranchStatus("npu",0)
    otree2.SetBranchStatus("npv",0)
    otree2.SetBranchStatus("evtWeight",0)
    otree2.SetBranchStatus("kfactor*",0)
    otree2.SetBranchStatus("pfmet*",0)
    otree2.SetBranchStatus("AK8Puppijet0_mass",0)
    otree2.SetBranchStatus("AK8Puppijet0_rho",0)
    otree2.SetBranchStatus("AK8Puppijet0_ptraw",0)
    otree2.SetBranchStatus("AK8Puppijet0_genpt",0)
    otree2.SetBranchStatus("AK8Puppijet0_phi",0)

    nent = tree.GetEntries()
    # print nent
    # finfo = ROOT.TFile.Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/signalXS/sig_vectordijet_xspt.root")
    # fvbf = ROOT.TFile.Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/signalXS/vbf_ptH_n3lo.root")
    # fr = ROOT.TFile.Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/signalXS/Higgs_v2.root")
    # h_ggh_num = fr.Get('gghpt_amcnlo012jmt')
    # h_ggh_den = fr.Get('ggh_hpt')
    # h_ggh_den.Scale(28.45024/h_ggh_den.Integral())
    # h_vbf_num = fvbf.Get('h_nnnlo_ptH')
    # h_vbf_den = fvbf.Get('h_lo_ptH')

    # # # h_rw = ROOT.TH1F()
    # h_rw = None
    # if 'VectorDiJet' in fn and mass > 0: 	
    #     hname = "med_"+str(mass)+"_0.1_proc_800"
    #     if '75' in fn: hname = "med_"+str(mass)+"_0.1_proc_801"
    #     hinfo = finfo.Get(hname)
    #     hinfo.Scale(100*1000.) # 100. for coupling, 1000. for conversion to pb is the cross-section
    #     hinfo_nbins = hinfo.GetNbinsX()
    #     hinfo_xlo = hinfo.GetXaxis().GetBinLowEdge(1)
    #     hinfo_xhi = hinfo.GetXaxis().GetBinUpEdge(hinfo_nbins)
    #     htmp = ROOT.TH1F("htmp","htmp",hinfo_nbins,hinfo_xlo,hinfo_xhi)
    #     for i in range(nent):
    #         tree.GetEntry(i)
    #         htmp.Fill(tree.genVPt,tree.scale1fb) 

    #     h_rw = ROOT.TH1F( hinfo.Clone() )
    #     h_rw.Divide(htmp)

    #newscale1fb = array( 'f', [ 0. ] ) #rewriting this guy
    # newkfactor  = array( 'f', [ 0. ] ) #rewriting this guy
    #tree.SetBranchAddress("scale1fb",newscale1fb)
    # tree.SetBranchAddress("kfactor",newkfactor)

    for i in range(nent):

        if( nent/100 > 0 and i % (1 * nent/100) == 0):
            sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done here")
            sys.stdout.flush()

        tree.GetEntry(i)

        #if (tree.AK8Puppijet0_pt > PT_CUT or tree.AK8Puppijet0_pt_JESUp > PT_CUT or tree.AK8Puppijet0_pt_JERUp > PT_CUT or tree.AK8Puppijet0_pt_JESDown > PT_CUT or tree.AK8Puppijet0_pt_JERDown > PT_CUT  ):
	    # if 'GluGluHToBB_M125_13TeV_powheg' in fn:  		
	    #     newscale1fb[0] =  tree.scale1fb*NLOcorr(tree.genVPt )
            # if 'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix' in fn and tree.genVPt<1000. : newscale1fb[0] =  tree.scale1fb*h_vbf_num.GetBinContent( h_vbf_num.FindBin(tree.genVPt) )/h_vbf_den.GetBinContent( h_vbf_den.FindBin(tree.genVPt) )
	    # if 'VectorDiJet' in fn and mass > 0:
            #     ptToWeightFrom = tree.genVPt
            #     if ptToWeightFrom < PT_CUT: ptToWeightFrom = PT_CUT # protection
            #     newscale1fb[0] = tree.scale1fb*h_rw.GetBinContent( h_rw.FindBin(ptToWeightFrom) )

        #newscale1fb[0] = tree.scale1fb
        otree2.Fill()   


    print "\n"
    otree2.AutoSave()
    ofile.cd()
    otree2.Write()
    ofile.Close()
    return 0
 	

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('--train', action='store_true', dest='train', default=False, help='train')
    parser.add_option("--lumi", dest="lumi", default = 30,type='float',help="luminosity", metavar="lumi")
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with bacon bits', metavar='idir')
    parser.add_option('-o','--odir', dest='odir', default = 'skim/',help='directory to write skimmed backon bits', metavar='odir')
    parser.add_option('-s','--sample',dest="sample", default="All",type='string',
                      #choices=['All','Hbb','QCD','JetHT','SingleMuon','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produces")


    (options, args) = parser.parse_args()

    main(options,args)
