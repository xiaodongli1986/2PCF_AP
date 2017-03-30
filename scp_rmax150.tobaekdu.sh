
#scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/DR12v4-CMASS-N/xyzw.finesplitted/data*.2pcf ./
#dataran.xyzw.1of3.cosmo-converted.om0.2100_w-0.7000.rmax51.51rbins.120mubins.rr
#for om in 0.3100 #0.1100 0.2100 0.2600 0.3100 0.4100
#do
#for w in  0.0000 -0.5000 -1.0000 -1.5000 -2.5000
#do
for ibin in 1 2 3
do
for catname in DR12v4-CMASS DR12v4-LOWZ
do
	nowlocaldir=~/SparseFilaments/data/input/boss2pcf/data/$catname/xyzw.binsplitted/
	mkdir -p $nowlocaldir
	sshpass -p lx821118 scp -r $nowlocaldir/data.xyzw.${ibin}of3.rmax150.150rbins.120mubins.2pcf xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/$catname/xyzw.binsplitted/
done
done
#done
#done
