

program main
use LSS_ximu_tests
USE de_model_init
!USE de_chisqs_JLA
implicit none
  
  integer :: i,j,k,i1,i2,iz,num_MCMCchain,numomwstds,nlines
  
  real(rt) :: omegam, omstds(1000), wstds(1000), DAs(nz), Hs(nz), & 
    chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1), &
    nowlnlike,nowweight,nowH0,nowomk,APlnlikemin,&
    t0,t1,t2,dt, &
    nowom, nowwde, nowxi, ommin, ommax, wdemin, wdemax, ximin, ximax !!! model dept
  integer :: iline, &
    numom, numwde, numxi, iom, iwde, ixi !!! model dept
  real(rt), allocatable :: APlnlikes(:), smutabstds(:,:,:,:,:)
  character(charlen) :: outputMCMCfile, mcmcdir='./', nowchisqstr, fileindexstr='', MCMCfilestr, suffixstr=''
  type(omwpar) :: nowpar
  logical :: smutabstds_inited, debug=.true., print_allinfo=.true.

!!! Next step: write this for wbinned model!!!

  de_model_lab = de_coupled_de_lab
  de_CP%coupled_de%use_xi1 = .true.
  MCMCfilestr = 'coupled_de__xi1'
!   numom=1;   ommin=0.27;  ommax=0.28;
!   numwde=1;  wdemin=-1.0; wdemax=-0.9;
!   numxi=1; ximin=0.0;  ximax=0.1;

!   numom=1;   ommin=(0.02224+0.08725)/0.6845**2.0;  ommax=0.28;  
!   numwde=1;  wdemin=-0.9434; wdemax=-0.9;
!   numxi=1; ximin=-0.0929;  ximax=0.1;

   numom=1;   ommin=(0.02232+0.121)/0.6793**2.0;  ommax=0.28;  
   numwde=1;  wdemin=-1.06; wdemax=-0.9;
   numxi=1; ximin=0.0007127;  ximax=0.1;

   
   !  suffixstr = 'base1omws_om0.2600_w-1.0000_ExcludeLastThreeBins_B'
  
!---------------------------------------------------------
  !--------------------------------
  ! Preparation for the compute of AP likelihood
  
  numomwstds = 1
  omstds(1)  = 0.26_rt;  wstds(1)  = -1.00_rt

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
    print *, '   wdemin, wdemax = ', wdemin, wdemax
    print *, '   ximin, ximax = ', ximin, ximax
    print *, '** Key-word: '
    print *, '   ', trim(adjustl( AP_MCMCstr(numomwstds, omstds(1:numomwstds), wstds(1:numomwstds)) ))
    print *, '** outputfile name: '
    print *, '   ', trim(adjustl(outputMCMCfile))
    
    nlines = numom*numwde*numxi
    allocate(APlnlikes(nlines))

    print *, '** Computing ', nlines, 'chisqs...'
    iline = 1
    call cpu_time(t0); t1=t0; dt = 60.0;

    do iom  = 1, numom
    do iwde = 1, numwde
    do ixi  = 1, numxi
      
      ! Begin model dependent
      nowom=ommin + (ommax-ommin)/max(dble(numom-1),1.0d0)*(iom-1);
      nowwde=wdemin + (wdemax-wdemin)/max(dble(numwde-1),1.0d0)*(iwde-1);
      nowxi=ximin + (ximax-ximin)/max(dble(numxi-1),1.0d0)*(ixi-1);
      nowH0=70.0;nowomk=0.0
      
      ! Values of parameters
!      de_CP%Ob0hsq  =  0.02253
      de_CP%Ob0hsq  =  0.02232
      de_CP%h	    =  nowH0 / 100.0
      de_CP%alpha   =  0.142125E+01
      de_CP%beta    =  0.325121E+01  
      de_CP%Odm0    =  nowom - de_CP%Ob0hsq/de_CP%h**2.0
      de_CP%Ok0	    =  nowomk
      de_CP%wcdm%w  =  nowwde !! model dependent
      de_CP%coupled_de%wde = nowwde
      de_CP%coupled_de%xi1 = nowxi      
      de_CP%coupled_de%xi2 = nowxi
      call de_init()
      ! End model dependent
      
      ! DAs & Hzs
      do iz = 1, nz
      
	DAs(iz) = de_Inte(zeffs(iz))*CONST_C/100.d0 / (1.0+zeffs(iz))
        Hs(iz) = 100.0 / de_inv_e(zeffs(iz)) 

        if(debug) then
          print *, 'DA: ', DAs(iz)
          print *, 'Hz: ', Hs(iz)
          nowpar%omegam = 0.26_rt; nowpar%w=-1.0_rt
          print *, 'Check DA of wcdm / our model: ', DAz_wcdm(nowpar,zeffs(iz)), DAs(iz)
          print *, 'Check H  of wcdm / our model ', Hz_wcdm(nowpar,zeffs(iz)), Hs(iz)
        endif
      enddo
      
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
!        APlnlikes(iline) = sum(chisqs_nosyscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
        APlnlikes(iline) = sum(chisqs_syscor_all(1:nz-1)) / 2.0 * (4.0/5.0)
        print *, 'Sys-corred:   ', real(chisqs_syscor_all(1:nz-1))* (4.0/5.0),  sum(chisqs_syscor_all(1:nz-1))* (4.0/5.0)
        print *, 'Sys-uncorred: ', real(chisqs_nosyscor_all(1:nz-1))* (4.0/5.0),  sum(chisqs_nosyscor_all(1:nz-1))* (4.0/5.0)
      else
        APlnlikes(iline) = 0.0d0
      endif

      if(debug.and..false.) then
        print *, 'chisqs_nosyscor(1,1):', real(chisqs_nosyscor(1,1,:)), real(sum(chisqs_nosyscor(1,1,:)))
        print *, 'chisqs_syscor(1,1):  ', real(chisqs_syscor(1,1,:)), real(sum(chisqs_syscor(1,1,:)))
        print *, 'chisqs_nosyscor_all:', real(chisqs_nosyscor_all)
        print *, 'chisqs_syscor_all:  ', real(chisqs_syscor_all)
      endif
      
      call cpu_time(t2)
      if (t2-t1.gt.dt.or.print_allinfo) then
        write(*,'(f10.1,A,i5,A,f4.1,A)') (t2-t0)/dt, ' minutes passed.   #-parameters = ', &
           iline, ' (',100*float(iline)/float(nlines),'%)'
        write(*,'(A,1x,6(f9.4))') '             Current set of APchi2 / par:  ', &
          APlnlikes(iline)*2, nowom, nowH0/100.0, nowwde, nowxi, nowomk
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
    do iom  = 0, numom
    do iwde = 1, numwde
    do ixi  = 1, numxi
      
      nowom=ommin + (ommax-ommin)/max(dble(numom-1),1.0d0)*(iom-1);
      nowwde=wdemin + (wdemax-wdemin)/max(dble(numwde-1),1.0d0)*(iwde-1);
      nowxi=ximin + (ximax-ximin)/max(dble(numxi-1),1.0d0)*(ixi-1);
      nowH0=70.0;nowomk=0.0
      
      nowlnlike = APlnlikes(iline)
      nowweight = exp(APlnlikemin - nowlnlike)
      write(3294,'(7(e14.7,1x))') nowweight, nowlnlike, &
        nowom, nowH0/100.0, nowwde, nowxi, nowomk !model dependent
      iline = iline+1
    enddo 
    enddo
    enddo     
    close(3294)
    deallocate(APlnlikes)
  endif     
  
 

end program

