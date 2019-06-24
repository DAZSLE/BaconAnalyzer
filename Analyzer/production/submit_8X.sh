python submitZprime.py -s TT_9X -t zprimebits-v16.03 --nfiles-per-job 20  --production 16 
python submitZprime.py -s W_9X -t zprimebits-v16.03 --nfiles-per-job 15   --production 16
python submitZprime.py -s DY_9X -t zprimebits-v16.03 --nfiles-per-job 10  --production 16
python submitZprime.py -s Hbb_9X -t zprimebits-v16.03 --nfiles-per-job 10 --monitor sub

python submitZprime.py -s TT_9X -t zprimebits-v16.03 --nfiles-per-job 20 --monitor sub --production 16
python submitZprime.py -s W_9X -t zprimebits-v16.03 --nfiles-per-job 15 --monitor sub --production 16
python submitZprime.py -s DY_9X -t zprimebits-v16.03 --nfiles-per-job 10 --monitor sub  --production 16
#python submitZprime.py -s QCD_9X -t zprimebits-v15.04-vbf --nfiles-per-job 10 --monitor sub --production 16
python submitZprime.py -s Hbb_9X -t zprimebits-v16.03 --nfiles-per-job 10 --monitor sub --production 16

# python skim.py -s TT_9X -i /eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v16.03/ --nsubmit 2
# python skim.py -s W_9X -i /eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v16.03/ --nsubmit 2
# python skim.py -s DY_9X -i /eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v16.03/ --nsubmit 2
# python skim.py -s Hbb_9X -i /eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v16.03/ --nsubmit 2


