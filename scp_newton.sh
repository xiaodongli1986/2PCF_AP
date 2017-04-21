
#scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/DR12v4-CMASS-N/xyzw.finesplitted/data*.2pcf ./
#dataran.xyzw.1of3.cosmo-converted.om0.2100_w-0.7000.rmax51.51rbins.120mubins.rr
#for om in 0.1100 0.2100 0.2600 0.3100 0.4100
#do
#for w in  0.0000 -0.5000 -1.0000 -1.5000 -2.5000
#do
#for ibin in 1 2 3
#do
for kw in om0.3100* om0.4100* om0.2600* #om0.1100_w-2.0000.rmax150.750rbins.600mubins om0.1100_w-1.0000.rmax150.750rbins.600mubins
do
for catname in DR12v4-CMASS DR12v4-LOWZ
do
	echo $kw : $catname
	mkdir -p $catname/xyzw.binsplitted/
	#sshpass -p lx821118 scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/$catname/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.om${om}_w${w}.rmax51.51rbins.120mubins.2pcf /home/xiaodongli/SparseFilaments/data/local_input/boss2pcf/data/$catname/xyzw.binsplitted/
	#sshpass -p lx821118 scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/$catname/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.om${om}_w${w}.rmax*.*rbins.*mubins.2pcf /home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/${catname}_Files/xyzw.binsplitted/
	sshpass -p lx821118 scp -r xiaodongli@newton.kias.re.kr:~/$catname/data*${kw}*.2pcf /home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/${catname}/xyzw.binsplitted/
done
done
#done
#done
#done
