# to create sub files
python submitZprime.py -s JetHTrereco_9X -t zprimebits-v16.04 --nfiles-per-job 10 --njobs-per-file 10 

# to submit
python submitZprime.py -s JetHTrereco_9X -t zprimebits-v16.04 --nfiles-per-job 10 --njobs-per-file 10  --monitor sub

# to submit skim (when raw bits are ready)
#python skim.py -s JetHTrereco_9X -i /eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v16.04/ --nsplit 5 --files-to-hadd 100 --nsubmit 2



