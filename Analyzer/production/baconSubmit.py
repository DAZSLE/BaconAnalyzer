#!/usr/bin/env python                                                                                                                                                                
import glob
import sys, commands, os, fnmatch, json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--list', action="store", dest="list")
parser.add_argument('--index', action="store", type=int, dest="index")
parser.add_argument('--nfiles-per-job', dest="nfiles_per_job", type=int, default=1)
parser.add_argument("--njobs-per-file", dest="njobs_per_file", type=int, default=1)
parser.add_argument("--eosoutdir", default='')
parser.add_argument('arguments', nargs='*')
args = parser.parse_args()

exec_line = ''
for arg in args.arguments:
    exec_line+= arg + ' '

with open(args.list, 'r') as mylist:
    allfiles = [(myfile.replace('\n', ''), True) for myfile in mylist.readlines()]

job_i = args.index
files=[]
lIFile = job_i*args.nfiles_per_job
for subjob_i in range(args.nfiles_per_job):
    files.append(allfiles[job_i*args.nfiles_per_job+subjob_i])
    lFFile = lIFile+subjob_i
print 'files ',files

cwd = os.getcwd()
for subjob_i in range(args.njobs_per_file):
    lDir = 'tmp_%s_job%i_subjob%i'%(args.list,job_i,subjob_i)
    os.system('mkdir -p %s'%lDir)
    os.chdir(lDir)
    nfiles_i = 0
    outfile = 'Output_job%d'%job_i
    job_hadd = 'hadd -f %s.root '%outfile

    for fil_i, fil in enumerate(files):
        if not fil[1]: continue
        new_line = exec_line.replace('fileinput'," "+fil[0]+" ")
        new_line = new_line.replace('subjob_i','%d'%subjob_i)
        new_line = new_line.replace('Output.root','Output_%d.root'%fil_i)
        time_line = "timeout 5m "
        print time_line+new_line
        os.system(time_line+new_line)
        job_hadd += 'Output_%d.root '%fil_i
        nfiles_i += 1
        
    print 'to hadd'
    os.listdir(os.getcwd())
    # hadd and copy
    print job_hadd
    os.system(job_hadd)
    lOut = ''
    if args.njobs_per_file > 1:                    
        lOut = 'root://cmseos.fnal.gov/%s/%s_file%dto%d_subjob%d.root'%(
            args.eosoutdir, outfile, lIFile, lFFile, subjob_i)
    else:
        lOut = 'root://cmseos.fnal.gov/%s/%s_file%dto%d.root'%(
            args.eosoutdir, outfile, lIFile, lFFile)
    job_copy = 'xrdcp -s %s.root %s; \n' %( outfile, lOut)
    os.system(job_copy)
    os.listdir(os.getcwd())
    os.system('rm *.root')
    os.chdir(cwd)
