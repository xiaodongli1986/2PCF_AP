
mkdir -p ../APCPL/covmats
cp -r src/*.f90 ../APCPL/src/
cp -r src/Makefile ../APCPL/src/

#exit

for mubin in 20 21 22 23 24 25
do
 cp covmats/${mubin}mubins.mumax0.97.iz6.CovMock_2000.s6.0to40.0.covmat ../APCPL/covmats
done

mkdir -p ../APCPL/2pCFs/



for catname in DR12v4-CMASS DR12v4-LOWZ
do
mkdir -p ../APCPL/2pCFs/$catname/xyzw.binsplitted
for ibin in 1 2 3
do
cp /home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/$catname/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf ../APCPL/2pCFs/$catname/xyzw.binsplitted
for imock in 0 1 2 3 
do
filename=J08.RSD.00$imock.xyzw.${ibin}of3.rmax150.150rbins.120mubins.2pcf
cp /home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/$catname/xyzw.binsplitted/$filename ../APCPL/2pCFs/$catname/xyzw.binsplitted
done
done
done
