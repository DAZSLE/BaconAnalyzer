#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
from optparse import OptionParser
from submitZprime import samplesDict, exec_me

import ROOT
normDict = {'DYJetsToQQ_HT180_13TeV': 'DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8',
            'WJetsToQQ_HT180_13TeV': 'WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8',
            'QCD_HT100to200_13TeV': 'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT100to200_13TeV_ext': 'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT200to300_13TeV': 'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT200to300_13TeV_ext': 'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT300to500_13TeV': 'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT300to500_13TeV_ext': 'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT500to700_13TeV': 'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT500to700_13TeV_ext': 'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_13TeV': 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_13TeV_ext': 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1000to1500_13TeV': 'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1000to1500_13TeV_ext': 'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1500to2000_13TeV': 'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT1500to2000_13TeV_ext': 'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT2000toInf_13TeV_ext': 'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'QCD_HT2000toInf_13TeV': 'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin': 'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
            'ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV_powhegV2_madspin': 'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
            'ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4': 'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
            'ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M2T4': 'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
            'TTJets_13TeV': 'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'TT_powheg': 'TT_TuneCUETP8M1_13TeV-powheg-pythia8',
            'WJetsToQQ_HT_600ToInf_13TeV': 'WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'ZJetsToQQ_HT_600ToInf_13TeV': 'ZJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'WWTo4Q_13TeV_amcatnlo': 'WWTo4Q_4f_13TeV_amcatnloFXFX_madspin_pythia8',
            'WWTo4Q_13TeV_powheg': 'WWTo4Q_13TeV_powheg',
            'WZ_13TeV_pythia8': 'WZ_TuneCUETP8M1_13TeV-pythia8',
            'ZZTo4Q_13TeV_amcatnlo': 'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8': 'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8',
            'GluGluHToBB_M125_13TeV_powheg_herwigpp': 'GluGluHToBB_M125_13TeV_powheg_herwigpp',
            'GluGluHToBB_M125_13TeV_powheg_pythia8': 'GluGluHToBB_M125_13TeV_powheg_pythia8',
            'VBFHToBB_M125_13TeV_amcatnlo_pythia8': 'VBFHToBB_M125_13TeV_amcatnlo_pythia8',
            'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix': 'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix',
            'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8',
            'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8': 'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
            'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
            'Spin0_ggPhi12j_g1_100_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_100',
            'Spin0_ggPhi12j_g1_100_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_100_Pseudo',
            'Spin0_ggPhi12j_g1_125_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_125_Pseudo',
            'Spin0_ggPhi12j_g1_125_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_125',
            'Spin0_ggPhi12j_g1_150_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_150_Pseudo',
            'Spin0_ggPhi12j_g1_150_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_150',
            'Spin0_ggPhi12j_g1_200_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_200_Pseudo',
            'Spin0_ggPhi12j_g1_200_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_200',
            'Spin0_ggPhi12j_g1_250_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_250_Pseudo',
            'Spin0_ggPhi12j_g1_250_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_250',
            'Spin0_ggPhi12j_g1_25_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_25_Pseudo',
            'Spin0_ggPhi12j_g1_25_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_25',
            'Spin0_ggPhi12j_g1_300_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_300_Pseudo',
            'Spin0_ggPhi12j_g1_300_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_300',
            'Spin0_ggPhi12j_g1_350_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_350_Pseudo',
            'Spin0_ggPhi12j_g1_350_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_350',
            'Spin0_ggPhi12j_g1_400_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_400',
            'Spin0_ggPhi12j_g1_400_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_400_Pseudo',
            'Spin0_ggPhi12j_g1_500_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_500_Pseudo',
            'Spin0_ggPhi12j_g1_500_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_500',
            'Spin0_ggPhi12j_g1_600_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_600_Pseudo',
            'Spin0_ggPhi12j_g1_600_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_600',
            'Spin0_ggPhi12j_g1_800_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_800_Pseudo',
            'Spin0_ggPhi12j_g1_800_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_800',
            'Spin0_ggPhi12j_g1_1000_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_1000_Pseudo',
            'Spin0_ggPhi12j_g1_1000_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_1000',
            'Spin0_ggPhi12j_g1_5_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_5',
            'Spin0_ggPhi12j_g1_5_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_5_Pseudo',
            'Spin0_ggPhi12j_g1_75_PseudoScalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_75_Pseudo',
            'Spin0_ggPhi12j_g1_75_Scalar_13TeV_madgraph': 'DMSpin0_ggPhibb1j_75',
            'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8',
            'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8': 'WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8',
            'WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8': 'WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8',
            'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8': 'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8',
            'ttHTobb_M125_13TeV_powheg_pythia8': 'ttHTobb_M125_13TeV_powheg_pythia8',
            'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV_powheg_pythia8': 'ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV_powheg_pythia8',
            'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8': 'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8',
            'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8': 'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8',
            'bbHToBB_M_125_4FS_yb2_13TeV_amcatnlo': 'bbHToBB_M_125_4FS_yb2_13TeV_amcatnlo',
            'bbHToBB_M_125_4FS_ybyt_13TeV_amcatnlo': 'bbHToBB_M_125_4FS_ybyt_13TeV_amcatnlo',
            'VectorDiJet1Jet_100_13TeV_madgraph':'VectorDiJet1Jet_M100',
            'VectorDiJet1Jet_125_13TeV_madgraph':'VectorDiJet1Jet_M125',
            'VectorDiJet1Jet_150_13TeV_madgraph':'VectorDiJet1Jet_M150',
            'VectorDiJet1Jet_200_13TeV_madgraph':'VectorDiJet1Jet_M200',
            'VectorDiJet1Jet_25_13TeV_madgraph':'VectorDiJet1Jet_M25',
            'VectorDiJet1Jet_300_13TeV_madgraph':'VectorDiJet1Jet_M300',
            'VectorDiJet1Jet_400_13TeV_madgraph':'VectorDiJet1Jet_M400',
            'VectorDiJet1Jet_500_13TeV_madgraph':'VectorDiJet1Jet_M500',
            'VectorDiJet1Jet_50_13TeV_madgraph':'VectorDiJet1Jet_M50',
            'VectorDiJet1Jet_600_13TeV_madgraph':'VectorDiJet1Jet_M600',
            'VectorDiJet1Jet_800_13TeV_madgraph':'VectorDiJet1Jet_M800',
            }
    
def main(options,args):

    DataDir = options.idir
    OutDir = options.idir
    
    samples = samplesDict[options.sample]
    
    EOS = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'
    postfix = ''
    exec_me('%s mkdir -p /%s/hadd'%(EOS,OutDir),options.dryRun)
    exec_me('%s mkdir -p /%s/sklim'%(EOS,OutDir),options.dryRun)
    exec_me('%s mkdir -p /%s/norm'%(EOS,OutDir),options.dryRun)
    for label, isMc in samples.iteritems():
        basename = label + '.root'
        filesToConvert, badFiles = getFilesRecursively(DataDir,label+'/',None,None)
        print "files To Convert = ",filesToConvert
        print "bad files = ", badFiles
        cwd = os.getcwd()
        for i in range(0,len(filesToConvert)/500+1):         
            haddCommand = '#!/bin/bash\n'
            haddCommand += 'pwd\n'
            haddCommand += 'cd %s\n'%cwd   
            haddCommand += 'pwd\n'
            haddCommand += 'eval `scramv1 runtime -sh`\n'
            haddCommand += 'cd -\n'
            haddCommand += 'pwd\n'
            haddCommand += 'mkdir -p $PWD/hadd\n'
            haddCommand += 'hadd -f hadd/%s %s\n'%(basename.replace('.root','_%i.root'%i),(' '.join(filesToConvert[i*500:(i+1)*500])).replace('eos','root://eoscms.cern.ch//eos'))
        haddCommand += 'hadd -f $PWD/hadd/%s $PWD/hadd/%s\n'%(basename,basename.replace('.root','_*.root'))
        haddCommand += 'rm $PWD/hadd/%s\n'%(basename.replace('.root','_*.root'))     
        haddCommand += '%s cp $PWD/hadd/%s /%s/hadd/%s\n'%(EOS,basename,OutDir,basename)        
        haddCommand += 'mkdir -p $PWD/sklim\n'
        haddCommand += 'export PYTHONPATH=${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/:${PYTHONPATH}\n'
        haddCommand += 'python ${CMSSW_BASE}/src/BaconAnalyzer/Analyzer/production/skimmer.py -i $PWD/hadd/ -o $PWD/sklim/ -s %s\n'%(basename.replace('.root',''))
        haddCommand += '%s cp $PWD/sklim/%s /%s/sklim/%s\n'%(EOS,basename,OutDir,basename)        
        if isMc=='mc':
            haddCommand += 'echo "%s\t${PWD}/sklim/%s" > normlist.txt\n'%(normDict[basename.replace('.root','')],basename)
            haddCommand += 'NormalizeNtuple normlist.txt\n'
            haddCommand += '%s cp $PWD/sklim/%s /%s/norm/%s\n'%(EOS,basename.replace('.root','_1000pb_weighted.root'),OutDir,basename.replace('.root','_1000pb_weighted.root'))
        haddCommand += 'rm -r $PWD/hadd\n'
        haddCommand += 'rm -r $PWD/sklim\n'
            
            
        with open('hadd_command_%s.sh'%(basename),'w') as f:
            f.write(haddCommand)

            exec_me('bsub -q 8nh -o $PWD/hadd_command_%s.log source $PWD/hadd_command_%s.sh'%(basename,basename),options.dryRun)

        with open('bad_files_%s.txt'%basename,'w') as f:
            for badFile in badFiles:
                f.write(badFile+'\n')


def getFilesRecursively(dir,searchstring,additionalstring = None, skipString = None):
    
    # thesearchstring = "_"+searchstring+"_"
    thesearchstring = searchstring

    theadditionalstring = None
    if not additionalstring == None: 
        theadditionalstring = additionalstring

    cfiles = []
    badfiles = []
    for root, dirs, files in os.walk(dir+'/'+thesearchstring):
        nfiles = len(files)
        for ifile, file in enumerate(files):
            
            if ifile%100==0:
                print '%i/%i files checked in %s'%(ifile,nfiles,dir+'/'+thesearchstring)
            try:
                #f = ROOT.TFile.Open((os.path.join(root, file)).replace('eos','root://eoscms.cern.ch//eos'))
                f = ROOT.TFile.Open((os.path.join(root, file)))
                if f.IsZombie():
                    print 'file is zombie'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                elif not f.Get('Events'):
                    print 'tree is false'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                elif not f.Get('Events').InheritsFrom('TTree'):
                    print 'tree is not a tree'
                    f.Close()
                    badfiles.append(os.path.join(root, file))                    
                    continue
                else:
                    f.Close()
                    cfiles.append(os.path.join(root, file))                    
            except:
                print 'could not open file or tree'
                badfiles.append(os.path.join(root, file))                    
                continue
                
    return cfiles, badfiles



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('--train', action='store_true', dest='train', default=False, help='train')
    parser.add_option("--lumi", dest="lumi", default = 30,type='float',help="luminosity", metavar="lumi")
    parser.add_option('-i','--idir', dest='idir', default = 'data/',help='directory with bacon bits', metavar='idir')
    #parser.add_option('-o','--odir', dest='odir', default = 'skim/',help='directory to write hadded bits', metavar='odir')
    parser.add_option('-s','--sample',dest="sample", default="All",type='string',
                      #choices=['All','Hbb','QCD','JetHT','SingleMuon','DMSpin0','TT','DY','W','Diboson','Triboson','SingleTop','VectorDiJet1Jet','VectorDiJet1Gamma','MC','Data'],
                      help="samples to produces")
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")
    (options, args) = parser.parse_args()

    main(options,args)
