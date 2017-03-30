
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
	mkdir -p $catname/xyzw.binsplitted
	for catname2 in J08 LC93 M12 V13 B08 HR3 HR4PSB J08.dat.z_0 J08.dat.z_0.5 PatchyV6C
	do 
		sshpass -p lx821118 scp -r xiaodongli@baekdu.kias.re.kr:~/boss2pcf/data/$catname/xyzw.binsplitted/$catname2.*.2pcf ./$catname/xyzw.binsplitted
	done
done
done
#done
#done
