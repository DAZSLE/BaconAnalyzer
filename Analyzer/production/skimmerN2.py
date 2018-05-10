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
    otree2.SetBranchStatus("AK8Puppijet0_mass",0)
    otree2.SetBranchStatus("nAK8Puppijets",0)
    otree2.SetBranchStatus("AK8Puppijet0_partonFlavor",0)
    otree2.SetBranchStatus("AK8Puppijet0_hadronFlavor",0)
    otree2.SetBranchStatus("AK8Puppijet0_nCharged",0)
    otree2.SetBranchStatus("AK8Puppijet0_nNeutrals",0)
    otree2.SetBranchStatus("AK8Puppijet0_nParticles",0)
    otree2.SetBranchStatus("AK8Puppijet0_N2sdb2",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_old",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JESUp",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JESDown",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JERUp",0)
    otree2.SetBranchStatus("AK8Puppijet0_pt_JERDown",0)
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
    otree2.SetBranchStatus("CA15Puppijet0_partonFlavor",0)
    otree2.SetBranchStatus("CA15Puppijet0_hadronFlavor",0)
    otree2.SetBranchStatus("CA15Puppijet0_nCharged",0)
    otree2.SetBranchStatus("CA15Puppijet0_nNeutrals",0)
    otree2.SetBranchStatus("CA15Puppijet0_nParticles",0)
    otree2.SetBranchStatus("CA15Puppijet0_mass",0)
    otree2.SetBranchStatus("CA15Puppijet0_rho",0)
    otree2.SetBranchStatus("CA15Puppijet0_ptraw",0)
    otree2.SetBranchStatus("CA15Puppijet0_genpt",0)
    otree2.SetBranchStatus("*Met*Corr*",0)
    otree2.SetBranchStatus("*tau21*",0)
    otree2.SetBranchStatus("*tau32*",0)
    otree2.SetBranchStatus("*doublesub*",0)
    otree2.SetBranchStatus("*CHF*",0)
    otree2.SetBranchStatus("*NHF*",0)
    otree2.SetBranchStatus("*NEMF*",0)
    otree2.SetBranchStatus("*CHS*",0)
    otree2.SetBranchStatus("sf_*",0)
    otree2.SetBranchStatus("rho",0)
    otree2.SetBranchStatus("*Scale_*",0)
    otree2.SetBranchStatus("PDF*",0)
    otree2.SetBranchStatus("puppet",0)
    otree2.SetBranchStatus("puppetphi",0)


    nent = tree.GetEntries()
    newscale1fb = array( 'f', [ 0. ] ) #rewriting this guy
    tree.SetBranchAddress("scale1fb",newscale1fb)

    for i in range(nent):

        if( nent/100 > 0 and i % (1 * nent/100) == 0):
            sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done here")
            sys.stdout.flush()

        tree.GetEntry(i)

        if (tree.AK8Puppijet0_pt > PT_CUT or tree.AK8Puppijet0_pt_JESUp > PT_CUT or tree.AK8Puppijet0_pt_JERUp > PT_CUT or tree.AK8Puppijet0_pt_JESDown > PT_CUT or tree.AK8Puppijet0_pt_JERDown > PT_CUT  ):

            newscale1fb[0] = tree.scale1fb
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
                      help="samples to produces")


    (options, args) = parser.parse_args()

    main(options,args)
