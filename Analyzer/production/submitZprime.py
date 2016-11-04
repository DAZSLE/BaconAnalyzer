#!/usr/bin/env python

import sys, commands, os, fnmatch
from optparse import OptionParser

def exec_me(command, dryRun=False):
    print command
    if not dryRun:
        os.sytem(command)
        
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")    
    parser.add_option("--monitor",default='',help="Monitor mode (sub/resub/check directory of jobs)")
    
    (options,args) = parser.parse_args()

    monitorOption = ''
    if options.monitor is not '':
        monitorOption = '--monitor %s'%options.monitor
    
    json = "$PWD/../data/Cert_271036-282037_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"
    xsec = 1

    optionsDataMc = {
        'mc': "Output.root --passSumEntries 5:Events  -a 2:mc -a 3:none  -n 7000 -q 2nw4cores",
        'data': "Output.root -a 5:1  -a 2:data -a 3:%s -n 7000 -q 1nh"%(json)
    }
        
    analysisdir = "zprimebits"
    executable = "runZprime"
    
    samplesDict = {}
    samplesDict['Hbb'] = {
        'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8': 'mc',  
        #'GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8_ext': 'mc',        
        #'GluGluHToBB_M125_13TeV_powheg_herwigpp': 'mc',
        'GluGluHToBB_M125_13TeV_powheg_pythia8': 'mc',
        #'VBFHToBB_M125_13TeV_amcatnlo_pythia8': 'mc',
        #'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix': 'mc',
        #'VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix_ext': 'mc',
        #'ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8': 'mc'
        }

    samples = samplesDict['Hbb']

    for label, ismc in samples.iteritems():    
        exec_me("python baconBatch.py %s %s -a 4:%f --list 1:../lists/production11/%s.txt -o $PWD/../%s/%s_%s %s"%(executable,optionsDataMc[ismc],xsec,label,analysisdir,label,ismc,monitorOption),options.dryRun)
