
#scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/DR12v4-CMASS-N/xyzw.finesplitted/data*.2pcf ./
#dataran.xyzw.1of3.cosmo-converted.om0.2100_w-0.7000.rmax51.51rbins.120mubins.rr
for om in 0.26
do
for w in  -1.0
do
for wa in -2.0000 -1.0000
do
for ibin in 1 2 3
do
for catname in DR12v4-CMASS DR12v4-LOWZ
do
	mkdir -p $catname/xyzw.binsplitted/
	#sshpass -p lx821118 scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/$catname/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.om${om}_w${w}.rmax51.51rbins.120mubins.2pcf /home/xiaodongli/SparseFilaments/data/local_input/boss2pcf/data/$catname/xyzw.binsplitted/
	#sshpass -p lx821118 scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/$catname/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.om${om}_w${w}.rmax*.*rbins.*mubins.2pcf /home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/${catname}_Files/xyzw.binsplitted/
	sshpass -p lx821118 scp -r xiaodongli@newton.kias.re.kr:${catname}//xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.om0.2600_w-1.0000_wa${wa}.rmax51.51rbins.120mubins.2pcf /home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/${catname}/xyzw.binsplitted/
done
done
done
done
done
