! Fiducial values:
!   z  =   0.215409830      0.315654516      0.386915922      0.478598893      0.540820956      0.618656099    
!   DA =    458.461761       584.131958       652.222473       720.240295       756.434937       792.691040    
!   H  =    133.994934       150.910355       163.336349       179.798065       191.266510       205.941605

program main
use LSS_ximu_tests
USE de_model_init
!USE de_chisqs_JLA
implicit none
  
  integer :: i,j,k,i1,i2,iz,num_MCMCchain,numDAratwstds,nlines
  
  
  real(rt) :: omegam, omstds(1000), wstds(1000), DAfids(nz), Hfids(nz), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    nowlnlike,nowweight,nowDArat,nowHrat,nowH0,nowwa,nowomk,APlnlikemin,&
    t0,t1,t2,dt, DAratmin, DAratmax, Hratmin, Hratmax, wamin, wamax
  integer :: iline, numDArat, numHrat0, numHrata, iDArat, iHrat, iwa
  real(rt), allocatable :: APlnlikes(:), smutabstds(:,:,:,:,:)
  character(charlen) :: outputallinfofile, outputMCMCfile, mcmcdir, nowchisqstr, fileindexstr, MCMCfilestr, &
    suffixstr='', sbinsuffix, tmpstr1, tmpstr2, tmpstr3, tmpstr4
  type(omwpar) :: nowpar
  integer,parameter :: model_wcdm=3, model_cpl=4, model_owcdm=5, model_ocpl=6, model_lcdm=7, model_olcdm=8
  integer :: nowmodel
  logical :: smutabstds_inited, debug=.false., print_allinfo=.true.

!!! Next step: write this for wbinned model!!!


  DAratmin=0.4;  DAratmax=1.2; 
  Hratmin=0.4; Hratmax=1.2; 
  wamin=0.0; wamax=0.1;
  numDArat=16; numHrat0=21; numHrata= 1

  nowmodel = model_wcdm; !DAratmin=0.31; DAratmax=0.31; numDArat=1


  sbinsuffix = ''
  write(tmpstr1, '(f6.1)') ints1
  write(tmpstr2, '(f6.1)') ints2
  write(tmpstr3, *) mubins(1)
  write(tmpstr4, *) mubins(N1)
  suffixstr = 'base1omws_om0.2600_w-1.0000_nbin'//trim(adjustl(tmpstr3))//'to'//trim(adjustl(tmpstr4))//&
   '_s'//trim(adjustl(tmpstr1))//'to'//trim(adjustl(tmpstr2))//'_mucut0.7'//trim(adjustl(sbinsuffix))!_fixom0.31'
! sbinsuffix = '_2sbin.s20'
! suffixstr = 'base1omws_om0.2600_w-1.0000_nbin10to15_'//trim(adjustl(sbinsuffix))!_fixom0.31'
  

  
!---------------------------------------------------------
  !--------------------------------
  ! Settings of the model

  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/'
  MCMCfilestr = 'DAH_scan_'
  suffixstr = 'scan1'

!---------------------------------------------------------
  !--------------------------------
  ! Preparation for the compute of AP likelihood
  
  numDAratwstds = 1
  omstds(1)  = 0.26_rt;  wstds(1)  = -1.00_rt
!  omstds(2)  = 0.26_rt;  stds(2)  = -0.60_rt
!  omstds(3)  = 0.26_rt;  wstds(3)  = -1.40_rt
!  omstds(4)  = 0.26_rt;  wstds(4)  = -0.60_rt
!  omstds(5)  = 0.31_rt;  wstds(5)  = -1.40_rt
!  omstds(6)  = 0.31_rt;  wstds(6)  = -0.60_rt

  print *, '(Begin) Load in necessary files.'
!  call system('sleep 100'); print *, 'Compute/output covmats...';call calc_covmats();call output_covmats(sbinsuffix)
  print *, '* Load in covmats:'
  call load_covmats(sbinsuffix)
  print *, '* Invert covmats:'
  call invert_covmats()
  print *, '* Compute systematic correction:'
  call calc_syscor()
  smutabstds_inited = .false.
  allocate(smutabstds(nbins_database,mubins_database,3,nz,numDAratwstds))
  ! End 
  !--------------------------------

!---------------------------------------------------------  
  !--------------------------------
  ! Scan of grid
  if(.true.) then

    fileindexstr = '_1.txt'
    if(trim(adjustl(suffixstr)).eq.'') &
      suffixstr = trim(adjustl( AP_MCMCstr(numDAratwstds, omstds(1:numDAratwstds), wstds(1:numDAratwstds)) ))
    outputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//'___'//&
      trim(adjustl(suffixstr))//trim(adjustl(fileindexstr))
    outputallinfofile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//'___'//&
      trim(adjustl(suffixstr))//'_allinfo'

    print *
    print *, '###################################################'
    print *, '** Compuate AP chisqs for: '
    print *, '   DAratmin, DAratmax = ', DAratmin, DAratmax
    print *, '   Hratmin, Hratmax = ', Hratmin, Hratmax
    print *, '   wamin, wamax = ', wamin, wamax
    print *, '** Key-word: '
    print *, '   ', trim(adjustl( AP_MCMCstr(numDAratwstds, omstds(1:numDAratwstds), wstds(1:numDAratwstds)) ))
    print *, '** outputfile name: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    
    nlines = numDArat*numHrat0*numHrata
    allocate(APlnlikes(nlines))
    
    print *, '** Computing ', nlines, 'chisqs...'
    iline = 1
    call cpu_time(t0); t1=t0; dt = 60.0;

    ! Values of parameters
    ! Using Planck+ext, Planck 2015 result (om = 0.3089, H0=67.74)
    
    nowH0=67.74d0;nowomk=0.0
    de_CP%Ob0hsq  =  0.02230!0.02253
    de_CP%h	    =  nowH0 / 100.0
    de_CP%alpha   =  0.142125E+01
    de_CP%beta    =  0.325121E+01  
    de_CP%Odm0    =  0.3089d0 & ! Omega_m, Planck 
       - de_CP%Ob0hsq/de_CP%h**2.0
    de_CP%Ok0	    =  0.0
    de_CP%wcdm%w  =  nowHrat !! model dependent
    de_CP%CPL%w0  =  nowHrat !! model dependent
    de_CP%CPL%wa  =  0.0d0 !! model dependent

    call de_init()
    
    ! DAs & Hzs, fiducial values
    do iz = 1, nz
        DAfids(iz) = de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
        Hfids(iz) = 100.0 / de_inv_e(zeffs(iz)) 
        if(debug) then
          nowpar%omegam = 0.06_rt; nowpar%w=-1.5_rt
          DAfids(iz) = DAz_wcdm(nowpar,zeffs(iz))
          Hfids(iz)  = Hz_wcdm(nowpar,zeffs(iz))
        endif
        !print *, 'Check DA: ', DAz_wcdm(nowpar,zeffs(iz)), DAs(iz)
        !print *, 'Check H:  ', Hz_wcdm(nowpar,zeffs(iz)), Hs(iz)
        !Hs(iz) = Hz_wcdm(nowpar, zeffs(iz))
    enddo
    print *, 'At six redshifts, fiducial values are: '
    print *, '  z  = ', real(zeffs)
    print *, '  DA = ', real(DAfids)
    print *, '  H  = ', real(Hfids)
    
    open(unit=2319084,file=outputallinfofile,action='write')
    do iDArat = 1, numDArat
    do iHrat = 1, numHrat0
    do iwa = 1, numHrata
      
      ! Begin model dependent
      nowDArat=DAratmin + (DAratmax-DAratmin)/max(dble(numDArat-1),1.0d0)*(iDArat-1);
      nowHrat=Hratmin + (Hratmax-Hratmin)/max(dble(numHrat0-1),1.0d0)*(iHrat-1); 
      nowwa=wamin + (wamax-wamin)/max(dble(numHrata-1),1.0d0)*(iwa-1);
      
      DAs = DAfids * nowDArat      
      Hs = Hfids * nowHrat
      if (print_allinfo) then
       print *, 'Value of DAs = ', real(DAs)
       print *, 'Value of Hs  = ', real(Hs)
      endif
      ! AP likelihood
      if(.true.) then
        call smu_ximu_CalcDAHChisqs(& 
          DAs, Hs, & ! Values of DA, H in six cosmologies
          omstds(1:numDAratwstds), wstds(1:numDAratwstds), numDAratwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
          smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
          chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
          chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
          weightedstds = .false., avg_counts = .false. &
          ) 
       APlnlikes(iline) = sum(chisqs_syscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
!        APlnlikes(iline) = sum(chisqs_nosyscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
      else
        APlnlikes(iline) = 0.0d0
      endif

      if(debug) then
        print *, 'chisqs_nosyscor(1,1):', real(chisqs_nosyscor(1,1,:)), real(sum(chisqs_nosyscor(1,1,:)))
        print *, 'chisqs_syscor(1,1):  ', real(chisqs_syscor(1,1,:)), real(sum(chisqs_syscor(1,1,:)))
        print *, 'chisqs_nosyscor_all:', real(chisqs_nosyscor_all)
        print *, 'chisqs_syscor_all:  ', real(chisqs_syscor_all)
        stop
      endif
      
      call cpu_time(t2)
      if (t2-t1.gt.dt.or.print_allinfo) then
        write(*,'(f10.1,A,i5,A,f4.1,A)') (t2-t0)/dt, ' minutes passed.   #-parameters = ', &
           iline, ' (',100*float(iline)/float(nlines),'%)'
        write(*,'(A,1x,6(f9.4))') '             Current set of wei / chi2 / APchi2 / par:  ', &
          APlnlikes(iline), nowDArat, nowH0/100.0, nowHrat, nowwa, nowomk
        t1=t2
      endif
      write(2319084,'(2f10.5, 3x, 2f10.3, 5x,  5f10.3, 5x, 5f10.3)') nowDArat, nowHrat, &
         sum(chisqs_nosyscor_all(1:nz-1)) * (4.0/5.0), sum(chisqs_syscor_all(1:nz-1)) * (4.0/5.0), &
         chisqs_nosyscor_all(1:nz-1)* (4.0/5.0), chisqs_syscor_all(1:nz-1)* (4.0/5.0)
      if (print_allinfo) &
         write(*,'(2f10.5, 3x, 2f10.3, 5x,  5f10.3, 5x, 5f10.3)') nowDArat, nowHrat, &
           sum(chisqs_nosyscor_all(1:nz-1)) * (4.0/5.0), sum(chisqs_syscor_all(1:nz-1)) * (4.0/5.0), &
           chisqs_nosyscor_all(1:nz-1)* (4.0/5.0), chisqs_syscor_all(1:nz-1)* (4.0/5.0)

      iline = iline +1
    enddo
    enddo
    enddo
    close(2319084)
    
    APlnlikemin = minval(APlnlikes(1:nlines))
    
    ! write the values to new file
    print *, '  Write   AP chisqs to file: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    open(unit=3294,file=outputMCMCfile,action='write')
    iline = 1
    do iDArat = 1, numDArat
    do iHrat = 1, numHrat0
    do iwa = 1, numHrata
      ! Begin model dependent
      nowDArat=DAratmin + (DAratmax-DAratmin)/max(dble(numDArat-1),1.0d0)*(iDArat-1);
      nowHrat=Hratmin + (Hratmax-Hratmin)/max(dble(numHrat0-1),1.0d0)*(iHrat-1); 
      nowwa=wamin + (wamax-wamin)/max(dble(numHrata-1),1.0d0)*(iwa-1);
      nowH0=70.0;nowomk=0.0
      
      nowlnlike = APlnlikes(iline)
      nowweight = exp(APlnlikemin - nowlnlike)
      write(3294,'(7(e14.7,1x))') nowweight, nowlnlike, &
        nowDArat, nowH0/100.0, nowHrat, nowwa, nowomk !model dependent
      iline = iline+1
    enddo 
    enddo
    enddo     
    close(3294)
    deallocate(APlnlikes)
  endif     
  
 
!  nowpar%omegam = 0.06_rt; nowpar%w = -1.5_rt

!  if(.true.) &
!   call smu_ximu_CalcDAHChisqs(& 
!    DAs, Hs, & ! Values of DA, H in six cosmologies
!    omstds(1:numDAratwstds), wstds(1:numDAratwstds), numDAratwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
!    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
!    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
!    chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
!    weightedstds = .true. &
!    ) 
!   do i = 1, N1
!   do j = 1, N2
!     if ( mubins(i) .eq. 25 .and. j.eq.1) then
!       print *, '* mubin / mucut = ', mubins(i),mucuts(j)
!       print *, '  chisqs (no cor) = ', real(chisqs_nosyscor(i,j,1:nz-1)), '; ', real(sum(chisqs_nosyscor(i,j,1:nz-1)))
!       print *, '  chisqs (corred) = ', real(chisqs_syscor(i,j,1:nz-1)), '; ', real(sum(chisqs_syscor(i,j,1:nz-1)))
!     endif
!   enddo
!   enddo
end program
