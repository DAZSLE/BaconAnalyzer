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

PT_CUT = 350.

def main(options,args):

    DataDir = options.idir
    OutDir = options.odir

    if options.ifile is None:
        try:
            samples = samplesDict[options.sample]
            tags = []
            for sample in samples:
                tags.append([sample, 0])
        except KeyError:
            tags = [[options.sample,0]]

        for i in range(len(tags)):
            lFiles = getFilesRecursively(DataDir,tags[i][0],None,None)
            print "files To Convert = ",lFiles

            for iFile in lFiles:
                print iFile
                status = sklimAdd(iFile,OutDir,tags[i][1])
                print status
    else:
        iFile = options.ifile
        print iFile
        tags = [[options.sample,0]]
        status = sklimAdd(iFile,OutDir,tags[0][1])
        print status

def sklimAdd(fn,odir,mass=0):

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
    nent = tree.GetEntries()

    newscale1fb = array( 'f', [ 0. ] )
    tree.SetBranchAddress("scale1fb",newscale1fb)

    for i in range(nent):

        if( nent/100 > 0 and i % (1 * nent/100) == 0):
            sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done here")
            sys.stdout.flush()

        tree.GetEntry(i)

        if (tree.AK8Puppijet0_pt > PT_CUT or 
            tree.AK8Puppijet0_pt_JESUp > PT_CUT or 
            tree.AK8Puppijet0_pt_JERUp > PT_CUT or 
            tree.AK8Puppijet0_pt_JESDown > PT_CUT or 
            tree.AK8Puppijet0_pt_JERDown > PT_CUT or 
            tree.CA15Puppijet0_pt > PT_CUT or 
            tree.CA15Puppijet0_pt_JESUp > PT_CUT or 
            tree.CA15Puppijet0_pt_JERUp > PT_CUT or 
            tree.CA15Puppijet0_pt_JESDown > PT_CUT or 
            tree.CA15Puppijet0_pt_JERDown > PT_CUT ) :
	    if 'GluGluHToBB_M125_13TeV_powheg' in fn:  		
	        newscale1fb[0] = tree.scale1fb*NLOcorr(tree.genVPt )
            elif 'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix' in fn and tree.genVPt<1000. : 
                newscale1fb[0] = tree.scale1fb*h_vbf_num.GetBinContent( h_vbf_num.FindBin(tree.genVPt) )/h_vbf_den.GetBinContent( h_vbf_den.FindBin(tree.genVPt) )
            else:
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
    parser.add_option('--ifile', dest='ifile', default = None, help='file to skim')
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with bacon bits', metavar='idir')
    parser.add_option('-o','--odir', dest='odir', default = 'skim/',help='directory to write skimmed backon bits', metavar='odir')
    parser.add_option('-s','--sample',dest="sample", default="All",type='string',
                      help="samples to produces")


    (options, args) = parser.parse_args()

    main(options,args)
