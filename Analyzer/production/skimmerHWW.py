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

PT_CUT = 300.

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
            print f
            status = sklimAdd(f,OutDir,options.sel1a,options.sel1b,tags[i][1])
            print status


def sklimAdd(fn,odir,sel1a,sel1b,mass=0):

    basename = os.path.basename( fn )

    print fn
    f1 = ROOT.TFile.Open(fn,'read')
    tree = f1.Get("Events")
    try:
        if not tree.InheritsFrom("TTree"):
            print 'not Events'
            return -1
    except:
        print 'exception'
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

    #otree.SetBranchStatus("*Puppijet0_e2*",0)
    #otree.SetBranchStatus("*Puppijet0_e3*",0)
    #otree.SetBranchStatus("*Puppijet0_e4*",0)

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

    newscale1fb = array( 'f', [ 0. ] ) #rewriting this guy
    # newkfactor  = array( 'f', [ 0. ] ) #rewriting this guy
    tree.SetBranchAddress("scale1fb",newscale1fb)
    # tree.SetBranchAddress("kfactor",newkfactor)

    for i in range(nent):

        if( nent/100 > 0 and i % (1 * nent/100) == 0):
            sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done here")
            sys.stdout.flush()

        tree.GetEntry(i)

        if (tree.AK8Puppijet0_pt > PT_CUT or tree.AK8Puppijet0_pt_JESUp > PT_CUT or tree.AK8Puppijet0_pt_JERUp > PT_CUT or tree.AK8Puppijet0_pt_JESDown > PT_CUT or tree.AK8Puppijet0_pt_JERDown > PT_CUT):

            # Take 3 highest pt jets:
            # Take the one with the largest lsf3: hww jet
            lmaxlsf3 = -99
            lhww = 1
            for i0 in range(0,3):
                lLepClsf3 = getattr(tree,"AK8Puppijet%s_lsfC_3"%str(i0))
                if lLepClsf3 > lmaxlsf3:
                    lhww = i0
            # Pick one of the other two:
            # Take the one with either the larger double-b or the larger msd
            lmaxdoublecsv = -99
            lmaxmsd = -99
            lhbb = 0
            for i0 in range(0,3):
                if i0 == lhww: continue
                if sel1a:
                    lDoublecsv = getattr(tree,"AK8Puppijet%s_doublecsv"%str(i0))
                    if lDoublecsv > lmaxdoublecsv:
                        lhbb = i0
                elif sel1b:
                    lMsd = getattr(tree,"AK8Puppijet%s_msd"%str(i0))
                    if lMsd > lmaxmsd:
                        lhbb = i0
                else:
                  lhbb = 0
                  lhww = 1
  
            # cuts on hww side (hww.lsf3 > 0.4,hww.tau31 < 0.6)
            lhwwlsf3 = getattr(tree,"AK8Puppijet%s_lsfC_3"%str(lhww))
            #lhwwmsd =  getattr(tree,"AK8Puppijet%s_msd"%str(lhww))
            lhwwtau31 = getattr(tree,"AK8Puppijet%s_tau31"%str(lhww))
            if lhwwlsf3 < 0.4 or lhwwtau31 > 0.6: #or lhwwmsd < 20
                continue

            # cuts on hbb side (hbb.doubleb > 0.2 and hbb.msd > 50
            lhbbdoublecsv = getattr(tree,"AK8Puppijet%s_doublecsv"%str(lhbb))
            lhbbmsd = getattr(tree,"AK8Puppijet%s_msd"%str(lhbb))
            if lhbbdoublecsv < 0.2 or lhbbmsd < 50:
                continue

            newscale1fb[0] = tree.scale1fb
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
            print file	
            if thesearchstring in file:
                if skipString != None and (skipString in file or skipString in dir or skipString in root):
                    print "already skimmed"
                    return []
                if theadditionalstring == None or theadditionalstring in file:
                    cfiles.append(os.path.join(root, file))
    return cfiles

def NLOcorr(Hpt=200.):
        NLO_= ROOT.TF1("NLO_", "pol2", 200, 1200)
	NLO_.SetParameter(0, 2.70299e+00)
	NLO_.SetParameter(1, -2.18233e-03)
	NLO_.SetParameter(2,5.22287e-07 )

        Weight = 1.
	if (Hpt>200) :
        	Weight = NLO_.Eval(Hpt)
	print Weight
        return Weight
 	

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('--train', action='store_true', dest='train', default=False, help='train')
    parser.add_option("--lumi", dest="lumi", default = 30,type='float',help="luminosity", metavar="lumi")
    parser.add_option('--sel1a', action='store_true', dest='sel1a', default=True, help='selection 1a')
    parser.add_option('--sel1b', action='store_true', dest='sel1b', default=False, help='selection 1b')
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with bacon bits', metavar='idir')
    parser.add_option('-o','--odir', dest='odir', default = 'skim/',help='directory to write skimmed backon bits', metavar='odir')
    parser.add_option('-s','--sample',dest="sample", default="All",type='string',
                      #choices=['All','Hbb','QCD','JetHT','SingleMuon','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produces")


    (options, args) = parser.parse_args()

    main(options,args)
