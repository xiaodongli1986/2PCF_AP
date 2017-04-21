

program main
use LSS_ximu_tests
USE de_model_init
!USE de_chisqs_JLA
implicit none
  
  integer :: i,j,k,i1,i2,iz
  integer, parameter :: numomwstds = 2
  real(rt) :: omegam, omstds(1000), wstds(1000), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    smutabstds(nbins_database,mubins_database,3,nz,numomwstds)
  character(charlen) :: inputMCMCfile, outputMCMCfile, mcmcdir, nowchisqstr
  type(omwpar) :: nowpar
  logical :: smutabstds_inited 

  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'

  inputMCMCfile = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA_post_lensing'

  omstds(1)  = 0.26_rt;  wstds(1)  = -1.00_rt
  omstds(2)  = 0.31_rt;  wstds(2)  = -1.00_rt
!  omstds(3)  = 0.26_rt;  wstds(3)  = -1.40_rt
!  omstds(4)  = 0.26_rt;  wstds(4)  = -0.60_rt
!  omstds(5)  = 0.31_rt;  wstds(5)  = -1.40_rt
!  omstds(6)  = 0.31_rt;  wstds(6)  = -0.60_rt
!  omstds(7)  = 0.41_rt;  wstds(7)  = -1.00_rt

  outputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(inputMCMCfile))//'___'//&
    trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
  print *, trim(adjustl(outputMCMCfile))
  
!  stop
  
  de_model_lab  = de_wcdm_lab
	
  de_CP%Ob0hsq  =  0.02253
  omegam 	=  0.264936E+00
  de_CP%h	=  0.711833E+00
  de_CP%alpha	=  0.142125E+01
  de_CP%beta	=  0.325121E+01  
  de_CP%Odm0 	=  omegam - de_CP%Ob0hsq/de_CP%h**2.0
  de_CP%Ok0	= 0
	
  call de_init()

  print *, 'Load in covmats...'
  call load_covmats('')
  print *, 'Invert covmats...'
  call invert_covmats()
  print *, 'Compute systematic correction...'
  call calc_syscor()

  
  smutabstds_inited = .false.
  nowpar%omegam = 0.06_rt; nowpar%w = -1.5_rt
!  omstds(1)  = 0.26_rt; wstds(1)  = -1.00_rt
  do iz = 1, nz
    DAs(iz) = DAz_wcdm(nowpar, zeffs(iz))
    Hs(iz) = Hz_wcdm(nowpar, zeffs(iz))
  enddo
  if(.true.) &
   call smu_ximu_CalcDAHChisqs(& 
    DAs, Hs, & ! Values of DA, H in six cosmologies
    omstds(1:numomwstds), wstds(1:numomwstds), numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
    chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
    weightedstds = .true. &
    ) 
   do i = 1, N1
   do j = 1, N2
     if ( mubins(i) .eq. 25 .and. j.eq.1) then
       print *, '* mubin / mucut = ', mubins(i),mucuts(j)
       print *, '  chisqs (no cor) = ', real(chisqs_nosyscor(i,j,1:nz-1))
       print *, '  chisqs (corred) = ', real(chisqs_syscor(i,j,1:nz-1))
     endif
   enddo
   enddo
end program
