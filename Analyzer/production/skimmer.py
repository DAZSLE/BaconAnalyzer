#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
from optparse import OptionParser
from submitZprime import samplesDict


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

        
    
    otree = tree.CloneTree(0)
    otree.SetName("otree")

    otree.SetBranchStatus("*Puppijet0_e2*",0)
    otree.SetBranchStatus("*Puppijet0_e3*",0)
    otree.SetBranchStatus("*Puppijet0_e4*",0)
    otree.SetBranchStatus("CA15Puppi*",0)	
    #otree.SetBranchStatus("bst8_PUPPIjet0_pt",1)

    nent = tree.GetEntries()
    print nent
    finfo = ROOT.TFile.Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/signalXS/sig_vectordijet_xspt.root")
    fvbf = ROOT.TFile.Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/signalXS/vbf_ptH_n3lo.root")
    fr = ROOT.TFile.Open("${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/signalXS/Higgs_v2.root")
    h_ggh_num = fr.Get('gghpt_amcnlo012jmt')
    h_ggh_den = fr.Get('ggh_hpt')
    h_ggh_den.Scale(28.45024/h_ggh_den.Integral())
    h_vbf_num = fvbf.Get('h_nnnlo_ptH')
    h_vbf_den = fvbf.Get('h_lo_ptH')

    # # h_rw = ROOT.TH1F()
    h_rw = None
    if 'VectorDiJet' in fn and mass > 0: 	
        hname = "med_"+str(mass)+"_0.1_proc_800"
        if '75' in fn: hname = "med_"+str(mass)+"_0.1_proc_801"
        hinfo = finfo.Get(hname)
        hinfo.Scale(100*1000.) # 100. for coupling, 1000. for conversion to pb is the cross-section
        hinfo_nbins = hinfo.GetNbinsX()
        hinfo_xlo = hinfo.GetXaxis().GetBinLowEdge(1)
        hinfo_xhi = hinfo.GetXaxis().GetBinUpEdge(hinfo_nbins)
        htmp = ROOT.TH1F("htmp","htmp",hinfo_nbins,hinfo_xlo,hinfo_xhi)
        for i in range(nent):
            tree.GetEntry(i)
            htmp.Fill(tree.genVPt,tree.scale1fb) 

        h_rw = ROOT.TH1F( hinfo.Clone() )
        h_rw.Divide(htmp)

    newscale1fb = array( 'f', [ 0. ] ) #rewriting this guy
    # newkfactor  = array( 'f', [ 0. ] ) #rewriting this guy
    tree.SetBranchAddress("scale1fb",newscale1fb)
    # tree.SetBranchAddress("kfactor",newkfactor)

    for i in range(nent):

        if( nent/100 > 0 and i % (1 * nent/100) == 0):
            sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done")
            sys.stdout.flush()

        tree.GetEntry(i)

        if (tree.AK8Puppijet0_pt > PT_CUT or tree.AK8Puppijet0_pt_JESUp > PT_CUT or tree.AK8Puppijet0_pt_JERUp > PT_CUT or tree.AK8Puppijet0_pt_JESDown > PT_CUT or tree.AK8Puppijet0_pt_JERDown > PT_CUT  ):
	    #if 'GluGluHToBB_M125_13TeV_powheg' in fn:  newscale1fb[0] =  h_ggh_num.GetBinContent( h_ggh_num.FindBin(tree.genVPt) )/h_ggh_den.GetBinContent( h_ggh_den.FindBin(tree.genVPt) )
            if 'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix' in fn and tree.genVPt<1000. : newscale1fb[0] =  tree.scale1fb*h_vbf_num.GetBinContent( h_vbf_num.FindBin(tree.genVPt) )/h_vbf_den.GetBinContent( h_vbf_den.FindBin(tree.genVPt) )
	    if 'VectorDiJet' in fn and mass > 0:
                ptToWeightFrom = tree.genVPt
                if ptToWeightFrom < PT_CUT: ptToWeightFrom = PT_CUT # protection
                newscale1fb[0] = tree.scale1fb*h_rw.GetBinContent( h_rw.FindBin(ptToWeightFrom) )

            else: newscale1fb[0] = tree.scale1fb
            otree.Fill()   


    print "\n"
    otree.AutoSave()
    ofile.cd()
    otree.Write()
    ofile.Close()
    return 0

def getFilesRecursively(dir,searchstring,additionalstring = None, skipString = None):
	
    # thesearchstring = "_"+searchstring+"_"
    thesearchstring = searchstring

    theadditionalstring = None
    if not additionalstring == None: 
        theadditionalstring = additionalstring

    cfiles = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            # print file	
            if thesearchstring in file:
                if skipString != None and (skipString in file or skipString in dir or skipString in root):
                    print "already skimmed"
                    return []
                if theadditionalstring == None or theadditionalstring in file:
                    cfiles.append(os.path.join(root, file))
    return cfiles

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
