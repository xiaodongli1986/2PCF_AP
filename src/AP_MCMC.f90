

program main
use LSS_ximu_tests
!use redbin_weights_bossdr12_6bins
USE de_model_init !!! Xiao-Dong: This is an outer package which compute DA, H of theoretical models.
                  !!!            In your case, you may comment this line, and just write your own
                  !!!            subroutines to compute DA, H in your model 
!USE de_chisqs_JLA
implicit none
  
  integer :: i,j,k,i1,i2,iz,num_MCMCchain,numomwstds,nlines
  
  
  integer :: iomcol=1,iw0col=1,iwacol=1,iomkcol=1,iH0col=1,maxcol, ifile,iline
  
  real(rt) :: omegam, omstds(1000), wstds(1000), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    tmpX(1000),nowlnlike,nowAPlnlike,nowweight,nowom,noww0,nowH0,nowwa,nowomk,APlnlikemin,&
     DAvalues(numgal_nbin),Hvalues(numgal_nbin),DAarray(numgal_nbin),Harray(numgal_nbin),&
    t0,t1,t2,dt
  real(rt), allocatable :: APlnlikes(:), smutabstds(:,:,:,:,:)
  character(charlen) :: inputMCMCfile, outputMCMCfile, mcmcdir, nowchisqstr, fileindexstr, MCMCfilestr, suffixstr=''
  type(omwpar) :: nowpar
  integer,parameter :: model_wcdm=3, model_cpl=4, model_owcdm=5, model_ocpl=6, model_lcdm=7, model_olcdm=8
  integer :: nowmodel
  logical :: smutabstds_inited, debug=.false., avg_counts = .false., print_allinfo=.false., VolumeWeightedDAH=.false.

  if(.not. numgal_inited) call numgal_init()
  !stop

!  nowmodel = model_wcdm
  nowmodel = model_cpl
!  nowmodel = model_olcdm
  suffixstr = 'base1omws_om0.2600_w-1.0000_Check'
!  print_allinfo = .true.
  
!---------------------------------------------------------
  !--------------------------------
  ! Settings of the model
  
  if(nowmodel .eq. model_wcdm) then! .or. nowmodel .eq. model_owcdm) then
	  de_model_lab  = de_wcdm_lab
	  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
          MCMCfilestr = 'base_w_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA'	  
          iw0col=5; iH0col=37; iomcol=39; 
  elseif(nowmodel .eq. model_cpl ) then
  	  de_model_lab = de_CPL_lab
  	  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
  	  MCMCfilestr = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA'
  	  iw0col=5; iwacol=6; iH0col=38; iomcol=40; 
  elseif(nowmodel .eq. model_olcdm ) then
  	  de_model_lab = de_LCDM_lab
  	  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
  	  MCMCfilestr = 'base_omegak_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA'
  	  iw0col=1; iwacol=1; iH0col=37; iomcol=39; iomkcol=5;
  else
          print *, 'Wrong model : ', nowmodel
          stop
  endif

!  mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA/'    
  
  ! Special setting for CMB+BAO chain
  !if(.false.) then
  ! mcmcdir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/plikHM_TTTEEE_lowTEB_BAO/'
  ! MCMCfilestr = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO'; iw0col=5; iwacol=6; iH0col=36; iomcol=38; 
  !endif
!   MCMCfilestr = 'base_w_wa_plikHM_TTTEEE_lowTEB_BAO_H070p6_JLA_post_lensing'

  maxcol=max(iw0col,iwacol,iH0col,iomcol,iomkcol)+2
  num_MCMCchain = 4

  ! End of settings
  !--------------------------------

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
!  call system('sleep 0'); print *, 'Compute/output covmats...';call calc_covmats();call output_covmats('_2sbin.s20')
!  print *, '* Load in covmats:';  call load_covmats('_2sbin.s20')
  print *, '* Load in covmats:'; call load_covmats('')
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
  ! Loop of all MCMC files
  do ifile = 1, num_MCMCchain

    ! File names
    write(fileindexstr,*) ifile
    fileindexstr = '_'//trim(adjustl(fileindexstr))//'.txt'
    inputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//trim(adjustl(fileindexstr))
    if(avg_counts) fileindexstr = '_avg_counts'//trim(adjustl(fileindexstr)) 
    if(trim(adjustl(suffixstr)).eq.'') &
      suffixstr = trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    outputMCMCfile = trim(adjustl(mcmcdir))//'/'//trim(adjustl(MCMCfilestr))//'___'//&
      trim(adjustl(suffixstr))//trim(adjustl(fileindexstr))

    print *
    print *, '###################################################'
    print *, '** Compuate AP chisqs from file: '
    print *, '   ', trim(adjustl(inputMCMCfile))
    print *, '** Key-word: '
    print *, '   ', trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    print *, '** outputfile name: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    
    ! Open file and compute likelihoods...
    call de_count_line_number (trim(adjustl(inputMCMCfile)), nlines); allocate(APlnlikes(nlines))
    print *, '** Computing ', nlines, 'chisqs...'
    open(unit=3293,file=inputMCMCfile,action='read')
    iline = 1
    call cpu_time(t0); t1=t0; dt = 60.0;
    ! Loop of all files
    do while(.true.)
      read(3293,*,end=100) tmpX(1:maxcol)
      
      ! Begin model dependent
      nowweight=tmpX(1); nowlnlike=tmpX(2); nowom=tmpX(iomcol+2); noww0=tmpX(iw0col+2); nowH0=tmpX(iH0col+2)
      nowwa=tmpX(iwacol+2); nowomk=tmpX(iomkcol+2)
      
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_owcdm .or. &
      	 nowmodel.eq.model_lcdm .or. nowmodel.eq.model_olcdm) &
      	    nowwa = 0.0; 
      	    
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_cpl   .or. &
         nowmodel.eq.model_lcdm)  &
            nowomk = 0.0; 
            
      if(nowmodel.eq. model_lcdm .or. nowmodel.eq.model_olcdm) &
            noww0 = 0
            
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
      if(.not.VolumeWeightedDAH) then
        do iz = 1, nz
         DAs(iz) = de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
         Hs(iz) = 100.0 / de_inv_e(zeffs(iz)) 
         if(debug) then
           nowpar%omegam = 0.06_rt; nowpar%w=-1.5_rt
           DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
           Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
         endif
         !print *, 'Check DA: ', DAz_wcdm(nowpar,zeffs(iz)), DAs(iz)
         !print *, 'Check H:  ', Hz_wcdm(nowpar,zeffs(iz)), Hs(iz)
         !Hs(iz) = Hz_wcdm(nowpar, zeffs(iz))
        enddo
      else
        do i = 1, numgal_nbin
            DAvalues(i) = de_Inte(numgal_zcenters(i))*CONST_C/100.d0 / (1.0+numgal_zcenters(i))
            Hvalues(i) = 100.0 / de_inv_e(numgal_zcenters(i))
        enddo
        do iz = 1, nz
          DAarray = DAvalues(1:numgal_nbin)*redbin_weights(iz,1:numgal_nbin)
          Harray = Hvalues(1:numgal_nbin)*redbin_weights(iz,1:numgal_nbin)
          DAs(iz) = sum(DAarray) / sum(redbin_weights(iz,1:numgal_nbin))
          Hs(iz) = sum(Harray) / sum(redbin_weights(iz,1:numgal_nbin))
          !print *, 'weighted redshift: ', iz
          !print *, DAs(iz), de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
          !print *, Hs(iz), 100.0 / de_inv_e(zeffs(iz))
        enddo
      endif
      !stop !!! We stop here. 
           !!! If we use "volume weighted DA/H", we shall use this also for the fiducial model.
           !!!  (for this, too complicated. We can add an option "input DAs, Hs", and cmpute wcdm DAs, Hs and input
      
      ! AP likelihood
      if(.true.) then
        call smu_ximu_CalcDAHChisqs(& 
          DAs, Hs, & ! Values of DA, H in six cosmologies
          omstds(1:numomwstds), wstds(1:numomwstds), numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
          smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
          chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
          chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
          weightedstds = .false., avg_counts = avg_counts, VolumeWeightedDAH=VolumeWeightedDAH &
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
        write(*,'(A,e12.4,1x,f10.3,A,f10.3,A,f8.3,A,5x,5(f9.4))') '             Current set of wei / chi2 /  par:  ', &
          nowweight, nowlnlike+APlnlikes(iline), '  (',nowlnlike, '  +', APlnlikes(iline), ' )', &
          nowom, nowH0/100.0, noww0, nowwa, nowomk
          
!       0.2 minutes passed.   #-parameters =     1 ( 0.0%)
!             Current set of wei / chi2 /  par:    0.1000E+01   6870.000  (  6833.867+    36.133)        0.2965   0.6927  -1.0644   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     2 ( 0.0%)
!             Current set of wei / chi2 /  par:    0.3000E+01   6869.146  (  6833.193+    35.953)        0.2970   0.6908  -1.0478   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     3 ( 0.0%)
!             Current set of wei / chi2 /  par:    0.1000E+01   6872.053  (  6835.938+    36.115)        0.2978   0.6918  -1.0667   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     4 ( 0.0%)
!             Current set of wei / chi2 /  par:    0.1000E+01   6869.462  (  6833.537+    35.925)        0.2996   0.6895  -1.0601   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     5 ( 0.0%)
!             Current set of wei / chi2 /  par:    0.2000E+01   6867.157  (  6830.227+    36.930)        0.3065   0.6820  -1.0364   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     6 ( 0.0%)
!             Current set of wei / chi2 /  par:    0.1000E+01   6875.092  (  6838.224+    36.868)        0.3070   0.6818  -1.0398   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     7 ( 0.1%)
!             Current set of wei / chi2 /  par:    0.4000E+01   6867.282  (  6830.157+    37.125)        0.3108   0.6778  -1.0218   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     8 ( 0.1%)
!             Current set of wei / chi2 /  par:    0.1000E+01   6866.756  (  6829.730+    37.026)        0.3087   0.6799  -1.0294   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =     9 ( 0.1%)
!             Current set of wei / chi2 /  par:    0.4000E+01   6867.033  (  6830.743+    36.290)        0.3031   0.6902  -1.0890   0.0000   0.0000
!       0.2 minutes passed.   #-parameters =    10 ( 0.1%)
!             Current set of wei / chi2 /  par:    0.3000E+01   6864.832  (  6828.445+    36.387)        0.3027   0.6909  -1.0907   0.0000   0.0000

        t1=t2
      endif
      iline = iline +1
      cycle
100   exit 
    enddo
    close(3293)
    
    APlnlikemin = minval(APlnlikes(1:nlines))
    
    ! write the values to new file
    print *, '  Write   AP chisqs to file: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    open(unit=3293,file=inputMCMCfile,action='read')
    open(unit=3294,file=outputMCMCfile,action='write')
    iline = 1
    do while(.true.)
      read(3293,*,end=101) tmpX(1:maxcol)
      nowweight=tmpX(1); nowlnlike=tmpX(2); nowom=tmpX(iomcol+2); noww0=tmpX(iw0col+2); nowH0=tmpX(iH0col+2) ! model dependent
      nowwa=tmpX(iwacol+2); nowomk=tmpX(iomkcol+2)
      
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_owcdm .or. &
      	 nowmodel.eq.model_lcdm .or. nowmodel.eq.model_olcdm) &
      	    nowwa = 0.0; 
      	    
      if(nowmodel.eq.model_wcdm .or. nowmodel.eq.model_cpl   .or. &
         nowmodel.eq.model_lcdm)  &
            nowomk = 0.0; 
            
      if(nowmodel.eq. model_lcdm .or. nowmodel.eq.model_olcdm) &
            noww0 = 0
            
      nowAPlnlike = APlnlikes(iline)
      nowweight = nowweight * exp(APlnlikemin - nowAPlnlike)
      write(3294,'(7(e14.7,1x))') nowweight, nowlnlike+nowAPlnlike, &
        nowom, nowH0/100.0, noww0, nowwa, nowomk !model dependent
      iline = iline+1
      cycle
101   exit
    enddo      
    close(3293); close(3294)
    deallocate(APlnlikes)
  enddo     
  
 
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
