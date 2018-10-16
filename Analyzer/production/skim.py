#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
from optparse import OptionParser
from submitZprime import samplesDict, exec_me

EOS = 'eos root://cmseos.fnal.gov'
cmssw = os.getenv('CMSSW_VERSION', 'CMSSW_9_4_7')
cmssw_base = os.getenv('CMSSW_BASE', 'CMSSW_9_4_7')

import ROOT

filesToTransfer = "{0}.tgz, {0}/bin/slc6_amd64_gcc630/NormalizeNtuple, {0}/src/BaconAnalyzer/Analyzer/data.tgz".format(cmssw_base)
filesToTransfer += ",{0}/src/BaconAnalyzer/Analyzer/production/skimmer.py, {0}/src/BaconAnalyzer/Analyzer/production/skimmerDDT.py, {0}/src/BaconAnalyzer/Analyzer/production/skimmerN2.py, {0}/src/BaconAnalyzer/Analyzer/production/skimmerWtag.py, {0}/src/BaconAnalyzer/Analyzer/production/skimmerHWW.py, {0}/src/BaconAnalyzer/Analyzer/production/submitZprime.py".format(cmssw_base)

def addCommand(iBasename,isMc,options,iFile=None):
    if iFile is None:
        pCommand = 'python %s -i $PWD/hadd/ -o $PWD/skim/ -s %s '%(options.skimmer,iBasename.replace('.root',''))
    else:
        pCommand = 'python %s --ifile %s -o $PWD/skim/ -s %s '%(options.skimmer,iFile,iBasename.replace('.root',''))
    if 'Wtag' in options.skimmer:
        pCommand+= '--jet %s --ddt %s --iddt %s '%(options.jet,options.ddt,options.iddt)
        if isMc=='mc':
            pCommand+= '--isMc '
            if not options.is2016:
                pCommand+= '--isPu '
    pCommand += '\n'
    return pCommand

def main(options,args):
    global iDataDir,iOutDir,iJobDir;
    iDataDir = options.idir
    iOutDir = options.idir
    iJobDir = 'skim_jobs'
    nFiles = int(options.filestohadd)
    lSamples = samplesDict[options.sample]
    
    exec_me('mkdir -p $PWD/%s/'%iJobDir,options.dryRun)
    exec_me('%s mkdir -p /%s/%s'%(EOS,iOutDir,skimDir),options.dryRun)

    for iLabel, isMc in lSamples.iteritems():
        pBasename = iLabel + '.root'
        if options.job>-1:
            pBasename = iLabel + '_%i.root'%options.job

        if os.path.isfile(iOutDir+'/'+skimDir+'/'+pBasename):
            print 'skim existing skip %s'%iLabel
            continue

        lFiles = []
        lBadFiles = []
        lFiles, lBadFiles = getFiles(iDataDir,iLabel+'/',None,None)
        #print "files To Convert = ",lFiles
        print "bad files = ", lBadFiles

        cwd = os.getcwd()
        pCommand = '#!/bin/bash\n'
        pCommand += 'source /cvmfs/cms.cern.ch/cmsset_default.sh\n'
        pCommand += 'pwd\n'
        pCommand += 'tar -xf %s.tgz\n'% (cmssw)
        pCommand += 'rm %s.tgz\n'% (cmssw)
        pCommand += 'export SCRAM_ARCH=slc6_amd64_gcc630\n'
        pCommand += 'mkdir -p %s/src\n'% (cmssw)
        pCommand += 'scramv1 project CMSSW %s\n'%cmssw
        pCommand += 'cd %s/src\n'%cmssw
        pCommand += 'scram b ProjectRename\n'
        pCommand += 'eval `scramv1 runtime -sh`\n'
        pCommand += 'pwd\n'
        pCommand += 'cp ../../data.tgz .\n'
        pCommand += 'mkdir -p ${PWD}/BaconAnalyzer/Analyzer/\n'
        pCommand += 'tar -xvzf data.tgz -C ${PWD}/BaconAnalyzer/Analyzer/\n'
        pCommand += 'mkdir -p $PWD/skim\n'
        pCommand += 'cp ../../submitZprime.py .\n'
        pCommand += 'cp ../../%s .\n'%options.skimmer
        pCommand += 'eval `scramv1 runtime -sh`\n'
        pCommand += 'rm -r $PWD/hadd\n'
        pCommand += 'rm -r $PWD/skim\n'
        pCommand += 'mkdir -p $PWD/hadd \n'
        pCommand += 'mkdir -p $PWD/skim \n'
        for i0 in range(0,len(lFiles)/nFiles+1):
            pBasename_i = pBasename.replace('.root','_%i.root'%i0)
            pInit = i0*nFiles;
            pFin = i0*nFiles+len(lFiles[i0*nFiles:(i0+1)*nFiles])
            for i1 in range(pInit,pFin):
                pFile = lFiles[i1]
                pCommand += addCommand(pBasename,isMc,options,pFile)
            pCommand += 'hadd -f $PWD/hadd/%s $PWD/skim/* \n'%(pBasename_i)
            pCommand += 'xrdcp -s $PWD/hadd/%s root://cmseos.fnal.gov//%s/%s/%s\n'%(pBasename_i,iOutDir,skimDir,pBasename_i)
            pCommand += 'rm -r -f $PWD/hadd/*\n'
        pCommand += 'rm -r $PWD/hadd\n'
        pCommand += 'rm -r $PWD/skim\n'
            
        with open('%s/skim_command_%s.sh'%(iJobDir,pBasename),'w') as f:
            f.write(pCommand)

        os.system('rm -f %s/skim_command_%s.stdout' % (iJobDir,pBasename))
        os.system('rm -f %s/skim_command_%s.stderr' % (iJobDir,pBasename))
        os.system('rm -f %s/skim_command_%s.log' % (iJobDir,pBasename))
        os.system('rm -f %s/skim_command_%s.jdl'% (iJobDir,pBasename))
        condor_file = open('%s/skim_command_%s.jdl' % (iJobDir,pBasename), 'w')
        condor_file.write('universe = vanilla\n')
        condor_file.write('Executable = %s/skim_command_%s.sh\n'% (iJobDir,pBasename))
        condor_file.write('Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )\n')
        condor_file.write('request_disk = 3000000\n')
        condor_file.write('request_memory = 10000\n')
        condor_file.write('Should_Transfer_Files = YES\n')
        condor_file.write('WhenToTransferOutput = ON_EXIT\n')
        condor_file.write('Transfer_Input_Files = %s\n'%filesToTransfer)
        condor_file.write('use_x509userproxy = true\n')
        condor_file.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
        condor_file.write('Output = %s.stdout\n' % os.path.abspath(condor_file.name))
        condor_file.write('Error = %s.stderr\n' % os.path.abspath(condor_file.name))
        condor_file.write('Log = %s.log\n' % os.path.abspath(condor_file.name))
        condor_file.write('Queue 1\n')
        condor_file.close()
        os.system('chmod +x %s'% os.path.abspath(condor_file.name))
        exec_me('condor_submit %s'%(os.path.abspath(condor_file.name)),options.dryRun)
                
        with open('%s/bad_files_%s.txt'%(iJobDir,pBasename),'w') as f:
            for iBadFile in lBadFiles:
                f.write(iBadFile+'\n')

def getFiles(dir,searchstring,additionalstring = None, skipString = None):
    
    thesearchstring = searchstring
    theadditionalstring = None
    if not additionalstring == None: 
        theadditionalstring = additionalstring

    cfiles = []
    badfiles = []
    files = []
    os.system('%s ls %s/%s > tmp.txt'%(EOS,dir,thesearchstring))
    with open("tmp.txt", 'r') as mylist:
        files = [(myfile.replace('\n', ''), True) for myfile in mylist.readlines()]
    nfiles = len(files)
    for ifile, fi in enumerate(files):
        if ifile%100==0:
            print '%i/%i files checked in %s'%(ifile,nfiles,dir+'/'+thesearchstring)
        if 'runPu' not in options.executable:
            try:
                filename = 'root://cmseos.fnal.gov//%s/%s/%s'%(dir,thesearchstring,fi[0])
                f = ROOT.TFile.Open(filename)
                if f.IsZombie():
                    print 'file is zombie'
                    f.Close()
                    badfiles.append(filename)
                    continue
                elif not f.Get('Events'):
                    print 'tree is false'
                    f.Close()
                    badfiles.append(filename)
                    continue
                elif not f.Get('Events').InheritsFrom('TTree'):
                    print 'tree is not a tree'
                    f.Close()
                    badfiles.append(filename)
                    continue
                else:
                    f.Close()
                    cfiles.append(filename)
            except:
                print 'could not open file or tree'
                badfiles.append(filename)
                continue
        else:
            filename = 'root://cmseos.fnal.gov//%s/%s/%s'%(dir,thesearchstring,fi[0])
            cfiles.append(filename)
                
    return cfiles, badfiles


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, 
                      help='no X11 windows')
    parser.add_option('-i','--idir', dest='idir', default = 'data/',
                      help='directory with bacon bits', metavar='idir')
    parser.add_option('-j','--job', dest='job', default =-1,type='int',
                      help='just do part of it', metavar='idir')
    parser.add_option('-s','--sample',dest="sample", default="All",type='string',
                      help="samples to produce")
    parser.add_option('-e','--executable',dest="executable", default="runZprime", 
                      help = "executable name")
    parser.add_option("--tagddt", type=str, default='', 
                      help="tagddt")
    parser.add_option("--ddt", type=str, default='GridOutput_v13.root', 
                      help="ddt")
    parser.add_option("--iddt", type=str, default='Rho2D', 
                      help="iddt")
    parser.add_option('--jet', dest='jet', default='AK8', 
                      help='jet type')
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                      help="Just print out commands to run")
    parser.add_option('--skimmer',dest="skimmer",default="skimmer.py", 
                      help = "skimmer")
    parser.add_option('--files-to-hadd',dest="filestohadd",default=50,type=int,
                      help = "number of files to hadd together after skimming")
    (options, args) = parser.parse_args()

    global skimDir;
    skimDir = "skim"
    if 'DDT' in options.skimmer: skimDir +='DDT'
    if 'N2' in options.skimmer: skimDir +='N2'
    if 'Wtag' in options.skimmer: skimDir +='Wtag%s%s'%(options.jet,options.tagddt)
    if 'HWW' in  options.skimmer: skimDir +='HWW1a'

    if '2016' in options.sample and not 'Wtag' in options.skimmer: skimDir+='2016'

    main(options,args)
    
