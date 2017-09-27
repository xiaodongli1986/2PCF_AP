

!Rhct
! Sys-corred:      16.0076275       26.95636 56       12.7971764       13.1067686       19.3084583       88.176397148907327     
! Sys-uncorred:    16.5017948       23.0967922       12.0079880       11.4178925       14.0837126       77.108180501000902

! LCDM model. om =   0.27   &nbs p;

! Sys-corred:      14.0279408       21.7550316       11.8282537       10.6070719       12.7000446       70.918342580960200     
! param  mean           sddev          lower1         upper1         lower2         upper2         lower3         upper3
!    1  0.2774861E+00  0.2212044E-01  0.2600098E+00  0.2961426E+00  0.2395020E+00  0.3320312E+00  0.2160645E+00  0.3835449E+00   \Omega_m

! Sys-uncorred:    14.5382814       18.5175323       12.0481081       9.10813332       11.9272099       66.1392657 22333093  

program main
use LSS_ximu_tests
USE de_model_init
!USE de_chisqs_JLA
implicit none
  
  integer :: i,j,k,i1,i2,iz,num_MCMCchain,numomwstds,nlines
  
  
  real(rt) :: omegam, omstds(1000), wstds(1000), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    nowlnlike,nowweight,nowom,noww0,nowH0,nowwa,nowomk,APlnlikemin,&
    t0,t1,t2,dt, ommin, ommax, w0min, w0max, wamin, wamax
  integer :: iline, numom, numw0, numwa, iom, iw0, iwa
  real(rt), allocatable :: APlnlikes(:), smutabstds(:,:,:,:,:)
  character(charlen) :: outputMCMCfile, mcmcdir, nowchisqstr, fileindexstr, MCMCfilestr, suffixstr=''
  type(omwpar) :: nowpar
  integer,parameter :: model_wcdm=3, model_cpl=4, model_owcdm=5, model_ocpl=6, model_lcdm=7, model_olcdm=8
  integer :: nowmodel
  logical :: smutabstds_inited, debug=.false., print_allinfo=.false., do_lcdm_constraint=.false.

!!! Next step: write this for wbinned model!!!

  de_model_lab  = de_wcdm_lab
  ommin=0.0;  ommax=1.0; 
   numw0=1; numwa=1;
!  print *, log(2.732); stop
  nowmodel = model_wcdm;   wamin=0.0;wamax=0.0;numwa=1
!  nowmodel = model_wcdm; !ommin=0.31; ommax=0.31; numom=1
!  nowmodel = model_olcdm

 
  !ommin=0.1; ommax=0.2; numom=34; suffixstr = 'base1omws_om0.2600_w-1.0000_1'!_fixom0.31'
  !ommin=0.20303030303d0; ommax=0.3; numom=33; suffixstr = 'base1omws_om0.2600_w-1.0000_2'!_fixom0.31'
  !ommin=0.30303030303d0; ommax=0.4; numom=33; suffixstr = 'base1omws_om0.2600_w-1.0000_3'!_fixom0.31'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_ExcludeLastThreeBins_B'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_nbins35to40'
!  suffixstr = 'base1omws_om0.2600_w-1.0000_smax40'
!  suffixstr = 'base1omws_om0.3100_w-1.0000_smax40'
!  suffixstr = 'base1omws_om0.1100_w-2.0000'
!  print_allinfo = .true.
  

  
!---------------------------------------------------------
  !--------------------------------
  ! Settings of the model

!! LCDM comments
  if(do_lcdm_constraint) then
    mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/'
    MCMCfilestr = 'base_lcdm_'
    numom = 2001
  else
    mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/'
    MCMCfilestr = 'base_Rhct_'
    numom = 31; ommin=0.1; ommax=0.4;
  endif
  
  
!---------------------------------------------------------
  !--------------------------------
  ! Preparation for the compute of AP likelihood
  
  numomwstds = 1
!  omstds(1)  = 0.11_rt;  wstds(1)  = -2.0_rt
  omstds(1)  = 0.26_rt;  wstds(1)  = -1.00_rt
  omstds(2)  = 0.26_rt;  wstds(2)  = -0.60_rt
!  omstds(3)  = 0.26_rt;  wstds(3)  = -1.40_rt
!  omstds(4)  = 0.26_rt;  wstds(4)  = -0.60_rt
!  omstds(5)  = 0.31_rt;  wstds(5)  = -1.40_rt
!  omstds(6)  = 0.31_rt;  wstds(6)  = -0.60_rt

  print *, '(Begin) Load in necessary files.'
!  call system('sleep 0'); print *, 'Compute/output covmats...';call calc_covmats();call output_covmats()
  print *, '* Load in covmats:'
  call load_covmats()
  print *, '* Invert covmats:'
  call invert_covmats()
  print *, '* Compute systematic correction:'
  call calc_syscor()
  smutabstds_inited = .false.
  allocate(smutabstds(nbins_database,mubins_database,3,nz,numomwstds))
  ! End 
  !--------------------------------

!---------------------------------------------------------  
  !--------------------------------
  ! Scan of grid
  if(.true.) then

    fileindexstr = '_1.txt'
    if(trim(adjustl(suffixstr)).eq.'') &
      suffixstr = trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    outputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//'___'//&
      trim(adjustl(suffixstr))//trim(adjustl(fileindexstr))

    print *
    print *, '###################################################'
    print *, '** Compuate AP chisqs for: '
    print *, '   ommin, ommax = ', ommin, ommax
    print *, '   w0min, w0max = ', w0min, w0max
    print *, '   wamin, wamax = ', wamin, wamax
    print *, '** Key-word: '
    print *, '   ', trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    print *, '** outputfile name: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    
    nlines = numom*numw0*numwa
    allocate(APlnlikes(nlines))
    
    print *, '** Computing ', nlines, 'chisqs...'
    iline = 1
    call cpu_time(t0); t1=t0; dt = 60.0;

    do iom = 0, numom
    do iw0 = 1, numw0
    do iwa = 1, numwa
      
      ! Begin model dependent
      nowom=ommin + (ommax-ommin)/max(dble(numom-1),1.0d0)*(iom-1);
      noww0=-1.0
      nowwa=0.0
      nowH0=70.0;nowomk=0.0
      
      ! Values of parameters
      de_CP%Ob0hsq  =  0.02253
      de_CP%h	    =  nowH0 / 100.0
      de_CP%alpha   =  0.142125E+01
      de_CP%beta    =  0.325121E+01  
      de_CP%Odm0    =  nowom - de_CP%Ob0hsq/de_CP%h**2.0
      de_CP%Ok0	    =  nowomk
      de_CP%wcdm%w  =  noww0 !! model dependent
      de_CP%CPL%w0  =  noww0 !! model dependent
      de_CP%CPL%wa  =  nowwa !! model dependent
      call de_init()
      ! End model dependent
      
      ! DAs & Hzs
      do iz = 1, nz
      
      !! LCDM comments
        if(do_lcdm_constraint) then
          DAs(iz) = de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
          Hs(iz) = 100.0 / de_inv_e(zeffs(iz)) 
        else
          if(iom.eq.0) then
            print *, 'DA of LCDM / Rhct', de_Inte(zeffs(iz)) / (1.0+zeffs(iz)), log(1.0+zeffs(iz)) / (1.0+zeffs(iz))
            print *, 'Hz of LCDM / Rhct', 100.0 / de_inv_e(zeffs(iz)) , 100.0 * (1.0+zeffs(iz))
	    DAs(iz) =  log(1.0+zeffs(iz)) / (1.0+zeffs(iz))*CONST_C/100.d0
	    Hs(iz) = 100.0 * (1.0+zeffs(iz))
	    if(iz.eq.1) print *, 'Rhct model. '
	  else
	    DAs(iz) = de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
            Hs(iz) = 100.0 / de_inv_e(zeffs(iz)) 
	    if(iz.eq.1) print *, 'LCDM model. om = ', nowom
	  endif
	endif
        if(debug) then
          nowpar%omegam = 0.06_rt; nowpar%w=-1.5_rt
          DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
          Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
        endif
!        print *, 'Check DA: ', DAz_wcdm(nowpar,zeffs(iz)), DAs(iz)
!        print *, 'Check H:  ', Hz_wcdm(nowpar,zeffs(iz)), Hs(iz)
        !Hs(iz) = Hz_wcdm(nowpar, zeffs(iz))
      enddo
!      stop
      
      ! AP likelihood
      if(.true.) then
        call smu_ximu_CalcDAHChisqs(& 
          DAs, Hs, & ! Values of DA, H in six cosmologies
          omstds(1:numomwstds), wstds(1:numomwstds), numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
          smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
          chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
          chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
          weightedstds = .false., avg_counts = .false. &
          ) 
        APlnlikes(iline) = sum(chisqs_syscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
        print *, 'Sys-corred:   ', real(chisqs_syscor_all(1:nz-1))* (4.0/5.0),  sum(chisqs_syscor_all(1:nz-1))* (4.0/5.0)
        print *, 'Sys-uncorred: ', real(chisqs_nosyscor_all(1:nz-1))* (4.0/5.0),  sum(chisqs_nosyscor_all(1:nz-1))* (4.0/5.0)
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
        write(*,'(A,1x,6(f9.4))') '             Current set of APchi2 / par:  ', &
          APlnlikes(iline), nowom, nowH0/100.0, noww0, nowwa, nowomk
        t1=t2
      endif
      iline = iline +1
    enddo
    enddo
    enddo
    
    APlnlikemin = minval(APlnlikes(1:nlines))
    
    ! write the values to new file
    print *, '  Write   AP chisqs to file: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    open(unit=3294,file=outputMCMCfile,action='write')
    iline = 1
    do iom = 1, numom
    do iw0 = 1, numw0
    do iwa = 1, numwa
      ! Begin model dependent
      nowom=ommin + (ommax-ommin)/max(dble(numom-1),1.0d0)*(iom-1);
      noww0=-1.0
      nowwa=0.0
      nowH0=70.0;nowomk=0.0
      
      nowlnlike = APlnlikes(iline)
      nowweight = exp(APlnlikemin - nowlnlike)
      write(3294,'(7(e14.7,1x))') nowweight, nowlnlike, &
        nowom, nowH0/100.0, noww0, nowwa, nowomk !model dependent
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
!    omstds(1:numomwstds), wstds(1:numomwstds), numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
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
