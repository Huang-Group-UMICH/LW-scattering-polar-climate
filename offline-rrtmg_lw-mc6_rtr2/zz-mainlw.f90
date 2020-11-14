      use parrrtm
      use rrtmg_lw_init, only : rrtmg_lw_ini
      use rrtmg_lw_rad, only : rrtmg_lw
      use shr_kind_mod, only: r8 => shr_kind_r8
      use rrtmg_lw_cldprmc, only: cldprmc

      implicit none
      integer, parameter :: ncol = 1
      integer, parameter :: nlay = 49
      integer :: icld 
      integer :: inflglw 
      integer :: iceflglw 
      integer :: liqflglw 

      real(kind=r8) :: play(ncol,nlay)
      real(kind=r8) :: plev(ncol,nlay+1)
      real(kind=r8) :: tlay(ncol,nlay)
      real(kind=r8) :: tlev(ncol,nlay+1)
      real(kind=r8) :: tsfc(ncol)
      real(kind=r8) :: h2ovmr(ncol,nlay)
      real(kind=r8) :: o3vmr(ncol,nlay)
      real(kind=r8) :: co2vmr(ncol,nlay)
      real(kind=r8) :: ch4vmr(ncol,nlay)
      real(kind=r8) :: o2vmr(ncol,nlay)
      real(kind=r8) :: n2ovmr(ncol,nlay)
      real(kind=r8) :: cfc11vmr(ncol,nlay)
      real(kind=r8) :: cfc12vmr(ncol,nlay)
      real(kind=r8) :: cfc22vmr(ncol,nlay)
      real(kind=r8) :: ccl4vmr(ncol,nlay)
      real(kind=r8) :: emis(ncol,nbndlw)

      real(kind=r8) :: cldfmcl(ngptlw,ncol,nlay)
      real(kind=r8) :: taucmcl(ngptlw,ncol,nlay)
      real(kind=r8) :: abscmcl(ngptlw,ncol,nlay)
      real(kind=r8) :: ssacmcl(ngptlw,ncol,nlay)
      real(kind=r8) :: xmomcmcl(0:16,ngptlw,ncol,nlay)
      real(kind=r8) :: ciwpmcl(ngptlw,ncol,nlay)
      real(kind=r8) :: clwpmcl(ngptlw,ncol,nlay)
      real(kind=r8) :: reicmcl(ncol,nlay)
      real(kind=r8) :: relqmcl(ncol,nlay)
      real(kind=r8) :: tauaer(ncol,nlay,nbndlw)

      real(kind=r8) :: rei(ncol,nlay)
      real(kind=r8) :: iciwp(ncol,nlay)

!      real(kind=r8) :: taucmcl(ngptlw,ncol,nlay-1)
!      real(kind=r8) :: tauaer(ncol,nlay,nbndlw-1)

      real(kind=r8) :: uflx(ncol,nlay+1)
      real(kind=r8) :: dflx(ncol,nlay+1)
      real(kind=r8) :: hr(ncol,nlay)
      real(kind=r8) :: uflxc(ncol,nlay+1)
      real(kind=r8) :: dflxc(ncol,nlay+1)
      real(kind=r8) :: hrc(ncol,nlay)
      real(kind=r8) :: uflxs(nbndlw,ncol,nlay+1)
      real(kind=r8) :: dflxs(nbndlw,ncol,nlay+1)
      real(kind=r8) :: netflx(ncol,nlay+1)
      real(kind=r8) :: netflxc(ncol,nlay+1)
      real(kind=r8) :: h2ogg(ncol,nlay)

      real(kind=r8) :: netflxs(nbndlw,ncol,nlay+1)
      real(kind=r8) :: hrs(nbndlw,ncol,nlay)
      real(kind=r8), parameter :: heatfac  = 8.43339130434783


      !logical, parameter :: lw_cloud_scat = .True.
      logical lw_cloud_scat

! local varibles
      character*10 pp
      character*10 choice_file
      character*100 filename_input
      character*100 filename_flux
      character*100 filename_spectral
      character*100 filename_vars1
      character*100 filename_out
      character*50 flag_scat
     
      integer i,k,k_cloud
      real(kind=r8) :: tt1 
      real(kind=r8) :: tmq
      real(kind=r8) :: rei0,iciwp0,q_factor,taucavg

      real(r8) emis_ice(16), emis_snow(16)

      data emis_snow / 0.9917, 0.9887, 0.9771, 0.9680, 0.9597,  0.9796,  0.9898, 0.9845,  &
         0.9787, 0.9747, 0.9742, 0.9682, 0.9654, 0.9654, 0.9654, 0.9654/ ! coarse snow
      data emis_ice / 0.8819, 0.9517,  0.9308,  0.9197, 0.9107, 0.9534,  0.9768, 0.9690, &
         0.9636, 0.9614, 0.9639, 0.9625, 0.9609, 0.9609, 0.9609, 0.9609/ ! ice
  
!----------------------------

!!!!!!!!!!!!!!!!!
! NEXT
!   check tauc, ssa with Chia-Pang's paper
!!!!!!!!!!!!!!!!!

!------------------------
! read from input files
!------------------------
!       !choice_file = "saw"
!       choice_file = "sas"
!
!       if (choice_file.eq."saw") then
!         filename_input = "zz-input-saw.txt"
!       else if (choice_file.eq."sas") then
!         filename_input = "zz-input-sas.txt"
!       else
!         print*,"ERROR: unsupported file"
!       end if
!
!      k_cloud = 45
!      cldfmcl(:,1,k_cloud) = 1._r8
!      rei(:,k_cloud) = 40.   ! effective radius, microns
!      iciwp(:,k_cloud) = 10. ! in-cloud ice water path, kg/m2
!
!      !flag_scat = "Scat"
!      flag_scat = "noScat"
!
!      q_factor = 1._r8
!
!      filename_flux = "./aa-dd1.txt"
!      filename_spectral = "./aa-dd2.txt"
!      filename_out = "./aa-dd3.txt"

      open(12,file='zz-input-parameter.txt',form='formatted')
        read(12,*) choice_file
        read(12,*) k_cloud
        read(12,*) rei0
        read(12,*) iciwp0
        read(12,*) flag_scat
        read(12,*) q_factor
        read(12,*) filename_flux
        read(12,*) filename_spectral
        read(12,*) filename_out
        read(12,*) filename_vars1
      close(12)

       if (choice_file.eq."saw") then
         filename_input = "zz-input-saw.txt"
       else if (choice_file.eq."sas") then
         filename_input = "zz-input-sas.txt"
       else
         print*,"ERROR: unsupported file"
       end if

      ! because rrtmg will add an additional layer, so move the index of k_cloud up a layer
      cldfmcl(:,1,k_cloud-1) = 1._r8
      rei(:,k_cloud-1) = rei0  ! effective radius, microns
      iciwp(:,k_cloud-1) = iciwp0 ! in-cloud ice water path, kg/m2

!        write(*,*) choice_file
!        write(*,*) k_cloud
!        write(*,*) rei0
!        write(*,*) iciwp0
!        write(*,*) flag_scat
!        write(*,*) q_factor
!        write(*,*) filename_flux
!        write(*,*) filename_spectral
!        write(*,*) filename_out


      !goto 9999

!------------------------
! read from input files
!------------------------

       print*,'rrtmg_lw program'
       !print*,'input file  : ./zz-input-scam-test.txt'
       !print*,'output file : ./zz-output-now.txt' 
       !print*,'output file : ./zz-output-spectral.txt' 
       print*,'---------------------------------------------'

!       print*,'call rrtmg_lw_ini'
!       call rrtmg_lw_ini

       !icld = 2
       !inflglw = 2  ! use cldprop
       !iceflglw = 4 ! 4 -> use MC6 
       !liqflglw = 1

       icld = 2 ! cloud overlap
       inflglw = 0
       iceflglw = 0
       liqflglw = 0

       cldfmcl(:,:,:) = 0._r8
       taucmcl(:,:,:) = 0._r8
       ciwpmcl(:,:,:) = 0._r8
       clwpmcl(:,:,:) = 0._r8
       tauaer(:,:,:)  = 0._r8
       reicmcl(:,:)   = 0._r8
       relqmcl(:,:)   = 0._r8
       ssacmcl(:,:,:)   = 0._r8
       xmomcmcl(:,:,:,:)  = 0._r8

      !print*,'rei',rei
      !print*,'iciwp',iciwp
      !print*,'taucmcl',taucmcl(:,1,k)

      open(12,file=filename_input,form='formatted')
      !open(12,file='./zz-input-scam-arm-02.txt',form='formatted')
      !open(12,file='./zz-input-chou-mls.txt',form='formatted')
        read(12,"(49(F11.6,2X))") play(ncol,1:nlay)
        read(12,"(50(F11.6,2X))") plev(ncol,1:nlay+1)
        read(12,"(49(F11.6,2X))") tlay(ncol,1:nlay)
        read(12,"(50(F11.6,2X))") tlev(ncol,1:nlay+1)
        read(12,*)                tsfc(ncol)
        read(12,"(49(E15.7,2X))")  h2ovmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  o3vmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  co2vmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  ch4vmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  o2vmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  n2ovmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  cfc11vmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  cfc12vmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  cfc22vmr(ncol,1:nlay)
        read(12,"(49(E15.7,2X))")  ccl4vmr(ncol,1:nlay)
        read(12,"(16(F11.6,2X))") emis(ncol,1:nbndlw)
      close(12)

!      o3vmr = 0._r8
!      h2ovmr = 0._r8
!      co2vmr = 0._r8
!      ch4vmr = 0._r8
!      o2vmr = 0._r8
!      n2ovmr = 0._r8
!      cfc11vmr = 0._r8
!      cfc12vmr = 0._r8
!      cfc22vmr = 0._r8
!      ccl4vmr = 0._r8

      call rrtmg_lw_ini

      !subroutine tamu_ice_get_rad_props_lw_OFFLINE(rei, iciwpth, ext_od, abs_od, ssa_od, xmomc_od)
      call tamu_ice_get_rad_props_lw_OFFLINE(ncol, nlay, ngptlw, rei, iciwp, taucmcl,abscmcl, ssacmcl, xmomcmcl)

      !print*,taucmcl(0,1,:)
      !taucmcl(:,1,k_cloud-1) = 10._r8
      !abscmcl(:,1,k_cloud-1) = 8._r8

!flag11
  !print*,taucmcl(:,1,k_cloud)
  !print*,abscmcl(:,1,k_cloud)
  !print*,ssacmcl(:,1,k_cloud)
  !taucmcl(:,1,k_cloud) = 4._r8
  !abscmcl(:,1,k_cloud) = 3.7_r8
!flag11

      taucavg = sum(taucmcl(:,1,k_cloud-1)) / ngptlw
      !print*,taucavg

      if (flag_scat .eq. "noScat") then
        ! yihsuan setup
        taucmcl(:,:,:)  = abscmcl(:,:,:)
        ssacmcl  = 0._r8
        xmomcmcl = 0._r8

        ! Chia-Pang setup
        !taucmcl(:,:,:)  = taucmcl(:,:,:)
        !ssacmcl  = 0._r8
        !xmomcmcl = 0._r8

      else if (flag_scat .eq. "Scat") then
        i=0

      else if (flag_scat .eq. "Emis_noScat") then
        ! yihsuan setup
        taucmcl(:,:,:)  = abscmcl(:,:,:)
        ssacmcl  = 0._r8
        xmomcmcl = 0._r8
        emis(1,:) = emis_ice(:)       

      else if (flag_scat .eq. "Emis_Scat") then
        i=0
        emis(1,:) = emis_ice(:)       

      else if (flag_scat .eq. "noScatEmis_ice") then
        taucmcl(:,:,:)  = abscmcl(:,:,:)
        ssacmcl  = 0._r8
        xmomcmcl = 0._r8
        emis(1,:) = emis_ice(:)       

      else if (flag_scat .eq. "noScatEmis_black") then
        taucmcl(:,:,:)  = abscmcl(:,:,:)
        ssacmcl  = 0._r8
        xmomcmcl = 0._r8
        !emis(1,:) = emis_ice(:)       
        emis(1,:) = 0.5_r8

      else if (flag_scat .eq. "noScatEmis_test") then
        taucmcl(:,:,:)  = abscmcl(:,:,:)
        ssacmcl  = 0._r8
        xmomcmcl = 0._r8
        emis(1,:) = 0.5_r8

      else if (flag_scat .eq. "ScatEmis_ice") then
        !tsfc = 287.8424
        !tlev(1,nlay+1) = tsfc(1)
        emis(1,:) = emis_ice(:)       

      else if (flag_scat .eq. "ScatEmis_black") then
        emis(1,:) = emis_ice(:)       
        i=1

      else 
        print*,"EREOR: unsupported flag_scat"      
        goto 9999
      end if

!================
! TAMU RRTMG
!================

      lw_cloud_scat = .True.
      h2ovmr = h2ovmr * q_factor

      call rrtmg_lw &
            (2   ,ncol    ,nlay    ,icld    , lw_cloud_scat,        &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,&
             cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             !cldfmcl ,taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             cldfmcl ,taucmcl , ssacmcl, xmomcmcl, ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, uflxs, dflxs )

      if (flag_scat .eq. "ScatEmis_black" .or. flag_scat .eq. "noScatEmis_black") then
        k=1
        tsfc(1) = (uflx(1,k)/5.67/1.E-8)**0.25

        !tlev(1,nlay+1) = tsfc(1)
        emis(1,:) = 1._r8
        !print*,'uflx,tsfc',uflx(1,k),tsfc(1)

      call rrtmg_lw &
            (2   ,ncol    ,nlay    ,icld    , lw_cloud_scat,        &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,&
             cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             !cldfmcl ,taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             cldfmcl ,taucmcl , ssacmcl, xmomcmcl, ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, uflxs, dflxs )

      end if

!----------------------
! write out to a flie
!----------------------
    do i=1,ncol
    do k=1,nlay
      h2ogg(i,k) = h2ovmr(i,k) * (18._r8/28.8_r8)
    end do
    end do

    tmq = 0._r8
    do k=1,nlay
      tmq = tmq + h2ogg(1,k) * (plev(1,k+1)-plev(1,k))*100._r8/9.8_r8
    end do

   !print*,'tmq (mm)=',tmq
   !print*,k_cloud,play(1,k_cloud),rei(1,k_cloud),iciwp(1,k_cloud),dflx(1,1),tmq
 
   open(12,file=filename_out,form='formatted')
     write(12,'(8(F7.3,2X))') play(1,k_cloud),rei(1,k_cloud-1),iciwp(1,k_cloud),tmq,taucavg,dflx(1,1),uflx(1,nlay+1)
   close(12)

!----------------------
! write out to a flie
!----------------------


      do k=1,nlay+1
        netflx (ncol,k) = uflx (ncol,k)-dflx (ncol,k)
        netflxc(ncol,k) = uflxc(ncol,k)-dflxc(ncol,k)
      enddo

! compute net spectral flux
      do i = 1 , nbndlw
        do k = 1,nlay+1
          netflxs(i,1,k) = uflxs(i,1,k) - dflxs(i,1,k)
        enddo
      enddo

! compute spectral heating rate
      do i = 1, nbndlw
        do k = 1,nlay
          hrs (i,1,k) = -heatfac*(netflxs(i,1,k+1)-netflxs(i,1,k)) &
                         / (plev(1,nlay-k+2)-plev(1,nlay-k+1))
        enddo
      enddo

      open(11,file=filename_flux,form='formatted')

        write(11,*) "-------------------------------------"
        write(11,*) "rrtmg longwave calculation result"
        write(11,*) "-------------------------------------"
        write(11,"(2X,A50)") " "
        write(11,"(2X,A20)") "*** input data ****"
        write(11,"(2X,A50)") " "
        write(11,"(2X,A12,E9.3)") "co2vmr = ",co2vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "n2ovmr = ",n2ovmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "ch4vmr = ",ch4vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "cfc11vmr = ",cfc11vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "cfc12vmr = ",cfc12vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "cfc22vmr = ",cfc22vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "o2vmr = ",o2vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "ccl4vmr = ",ccl4vmr(ncol,1)

        write(11,*) " "
        write(11,"(1X,A1,2(6X,A4),2(4X,A4),6X,A6,4X,A5)") &
              "k","plev","play","tlev","tlay","h2ovmr","o3vmr"
        do k=1,nlay
          write(11,"(I2,2X,F10.5,10X,F7.3)") k,plev(ncol,k),tlev(ncol,k)
          write(11,"(13X,F10.5,10X,F7.3,3(2X,E9.3))") &
                      play(ncol,k),tlay(ncol,k), &
                      h2ovmr(ncol,k), o3vmr(ncol,k) 
        enddo
        write(11,"(I2,2X,F10.5,10X,F7.3)") &
               k,plev(ncol,nlay+1),tlev(ncol,nlay+1)
        
        write(11,*) " "
        write(11,*) " "
        write(11,"(2X,A20)") "*** output data ***"
        write(11,"(2X,A50)") " "
        write(11,"(1X,A1,6X,A4,2(6X,A5),4X,A7,6X,A3,2X,2(6X,A5),4X,A7,6X,A3)") &
               "k","plev", & 
               "uflx","dflx","netflx","hr", &
               "uflxc","dflxc","netflxc","hrc"
! write flux
        do k=1,nlay
          write(11,"(I2,4(2X,F9.4),11X,3(2X,F9.4))") &
            k,plev(1,k),& 
              uflx(1,nlay-k+2),dflx(1,nlay-k+2), netflx(1,nlay-k+2), &
              uflxc(1,nlay-k+2),dflxc(1,nlay-k+2), netflxc(1,nlay-k+2)
! write heating rate
          write(11,"(48X,F8.3,35X,F8.3)") hr(1,nlay-k+1), hrc(1,nlay-k+1)
        enddo

! write surface flux
        k=nlay+1
          write(11,"(I2,4(2X,F9.4),11X,3(2X,F9.4))") &
            k,plev(1,k),& 
              uflx(1,1),dflx(1,1), netflx(1,1), &
              uflxc(1,1),dflxc(1,1), netflxc(1,1)
      close(11)

! write out spectal flux & heating rate
      open(11,file=filename_spectral,form='formatted')

        write(11,*) "-------------------------------------------------"
        write(11,*) "rrtmg longwave calculation result (each band)"
        write(11,*) "-------------------------------------------------"
        write(11,"(2X,A50)") " "
        write(11,"(2X,A20)") "*** input data ****"
        write(11,"(2X,A50)") " "
        write(11,"(2X,A12,E9.3)") "co2vmr = ",co2vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "n2ovmr = ",n2ovmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "ch4vmr = ",ch4vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "cfc11vmr = ",cfc11vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "cfc12vmr = ",cfc12vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "cfc22vmr = ",cfc22vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "o2vmr = ",o2vmr(ncol,1)
        write(11,"(2X,A12,E9.3)") "ccl4vmr = ",ccl4vmr(ncol,1)

        write(11,*) " "
        write(11,"(1X,A1,2(6X,A4),2(4X,A4),6X,A6,4X,A5)") &
              "k","plev","play","tlev","tlay","h2ovmr","o3vmr"
        do k=1,nlay
          write(11,"(I2,2X,F10.5,10X,F7.3)") k,plev(ncol,k),tlev(ncol,k)
          write(11,"(13X,F10.5,10X,F7.3,3(2X,E9.3))") &
                      play(ncol,k),tlay(ncol,k), &
                      h2ovmr(ncol,k), o3vmr(ncol,k)
        enddo
        write(11,"(I2,2X,F10.5,10X,F7.3)") &
               k,plev(ncol,nlay+1),tlev(ncol,nlay+1)

        write(11,*) " "
        write(11,*) " "
        write(11,"(2X,A20)") "*** output data ***"

        do i = 1,nbndlw
          write(11,7000)
          write(11,7001) i
          write(11,"(2X,A50)") " "
          write(11,"(1X,A1,6X,A4,2(6X,A5),5X,A7,6X,A3)") &
               "k","plev","uflx","dflx","netflx","hr"
          write(11,7002)

          do k=1,nlay
            write(11,"(I2,4(2X,F9.4))") &
              k,plev(1,k),uflxs(i,1,nlay-k+2),dflxs(i,1,nlay-k+2), &
                netflxs(i,1,nlay-k+2)
            write(11,"(48X,F8.3)") hrs(i,1,nlay-k+1)
          enddo

          k=nlay+1
          write(11,"(I2,4(2X,F9.4))") &
              k,plev(1,k),uflxs(i,1,nlay-k+2),dflxs(i,1,nlay-k+2), &
                netflxs(i,1,nlay-k+2)
          write(11,"(2X,A50)") " "
        enddo
      close(11)

      open(11,file=filename_vars1,form='formatted')
        write(11,*) '; ice effective radius (micron) = ',rei0

        k=nlay+1
        write(11,"(A20,16(2X,F9.4),A5)") 'FLUS( , :) = (/',uflxs(:,1,nlay-k+2),'/)'
        write(11,"(A20,16(2X,F9.4),A5)") 'FLDS( , :) = (/',dflxs(:,1,nlay-k+2),'/)'
        write(11,"(A20,16(2X,F9.4),A5)") 'FLUT( , :) = (/',uflxs(:,1,nlay),'/)'
        write(11,"(2X,A50)") " "  
      close(11)


7000  format('------------------------------------------')
7001  format('band ',I2)
7002  format(8x,'(mb)',5x,'(W/m2)',6x,'(W/m2)',4x,'(W/m2)',5x,'(K/day)')
7003  format("(A20,16(2X,F9.4),A5)")

9999  continue
      stop
      end

      subroutine tamu_ice_get_rad_props_lw_OFFLINE(pcols, pver, nsubcol, rei, iciwp, ext_od, abs_od, ssa_od, xmomc_od)
       
       !-----------------------------------------------------
       ! yihsuan@UMich
       !
       ! Description:
       !   TAMU ice cloud optics scheme
       !
       ! History:
       !   2016/05/20  ver 1.0 , ADD
       !   2017/08/03  ver 1.1 , ADD extinction optical depth and absorption optical depth
       !   2017/09/07  ver 2.1 , DEBUG use REI, ice effective radius, directly from CESM becasue it's already times a factor of 1.5 
       !   2017/09/19  ver 2.2 , DEBUG iciwp, previous version use inconsistent iciwp with the Mitchell one, and the unit is also wrong
       !
       ! Author:
       !   Yi-Hsuan Chen
       !   yihsuan@umich.edu
       !-----------------------------------------------------
       
        use shr_kind_mod, only: r8 => shr_kind_r8
        use parrrtm
          use rrlw_cld, only: abscld1, absliq0, absliq1, &
                              absice0, absice1, absice2, absice3, &
                              absice4, extice4, ssaice4, asyice4 ! CPKuo@TAMU
        use rrlw_wvn, only: ngb
       
       !--- Input arguments ---!
          integer, intent(in) :: pcols, pver, nsubcol
          real(r8), intent(in) :: rei(pcols,pver)     ! ice effective radius (microns)
          real(r8), intent(in) :: iciwp(pcols,pver) ! in-cloud ice water path (kg/m2)

       !--- Output arguments ---!
          real(r8), intent(out) :: ext_od  (nsubcol,pcols,pver)        ! cloud ice extinction optical depth, i.e. including absorption and scattering
          real(r8), intent(out) :: abs_od  (nsubcol,pcols,pver)        ! cloud ice absorption optical depth
          real(r8), intent(out) :: ssa_od  (nsubcol,pcols,pver)        ! single scattering albedo
          real(r8), intent(out) :: xmomc_od(0:16,nsubcol,pcols,pver)   ! phase function
       
       !--- Output arguments ---!
          real(r8) :: ext_od0  (nbndlw,pcols,pver)        ! cloud ice extinction optical depth, i.e. including absorption and scattering
          real(r8) :: abs_od0  (nbndlw,pcols,pver)        ! cloud ice absorption optical depth
          real(r8) :: ssa_od0  (nbndlw,pcols,pver)        ! single scattering albedo
          real(r8) :: xmomc_od0(0:16,nbndlw,pcols,pver)   ! phase function
       
       !--- Local ---!
          real(r8) :: diaice                        ! cloud ice effective diameter (microns)
          real(r8) :: invrad 
          real(r8) :: extcoice(nbndlw)            ! ice mass-extinction coefficients of TAMU scheme (m2/g)
          real(r8) :: ssacoice(nbndlw)            ! ice single scattering albedo of TAMU scheme (unitless)
          real(r8) :: asycoice(nbndlw)            ! ice asymmetric factor of TAMU scheme (unitless)
          real(r8) :: ext_tamu(nbndlw,pcols,pver) ! ice cloud extinction optical thickness (unitless)
          real(r8) :: abs_tamu(nbndlw,pcols,pver) ! ice cloud absorption optical thickness (unitless)
          real(r8) :: ssa_tamu(nbndlw,pcols,pver) ! ice cloud single scattering albedo (unitless)
          real(r8) :: xmomc_tamu(0:16,nbndlw,pcols,pver) ! ice cloud asymmetric factor (unitless)
          !real(r8) :: iciwp(pcols,pver)             ! work array of in-cloud ice water path (kg/m2)
          real(r8) :: taucloud, ssacloud, asycloud, xmomcloud(0:16), fpeak
       
          real(r8) :: gicewp(pcols,pver)            ! work array of grid-box ice water path (kg/m2)
          real(r8) :: cicewp(pcols,pver)            ! work array of in-cloud ice water path (kg/m2)    
       
          integer  :: icb(nbndlw,0:2)
          data icb /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
                    1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
                    1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/
       
          logical :: option_pbuf
       
          !integer :: rei_idx, iciwp_idx, cld_idx, err, itim
          integer :: rei_idx, iciwp_idx, err, itim
          integer :: ib,i,k,ncol,lchnk,iceind
          integer :: ilev, isubcol,nlay
       !----------------------------
       ! initialize return arrays
       !----------------------------
          ext_tamu   = 0._r8
          abs_tamu   = 0._r8
          ssa_tamu   = 0._r8
          xmomc_tamu = 0._r8
       
          ext_od   = 0._r8
          abs_od   = 0._r8
          ssa_od   = 0._r8
          xmomc_od = 0._r8

          ncol = pcols
          nlay = pver
       !--------------------------------------------------------
       ! compute absorption coefficient and cloud optical depth
       !--------------------------------------------------------
!print*,ncol,pver,nbndlw
          do i = 1,ncol
             do k = 1,pver

!print*,'aa1.0'
               diaice = 2._r8*rei(i,k)  ! effective diameter is defined by that particle-volume divided by particle-projected-area (i.e. effective diameter) times 1.5.
                                        ! Chia-Pang Kuo @ TAMU, personal communication
                                        ! because the REI in CESM is already times 1.5, so don't need to times 1.5 here.
                                        ! (ref: eq. 4.154 in CAM5 scientific description)
!print*,'aa1.1'
               extcoice(:) = 0._r8
               ssacoice(:) = 0._r8
               asycoice(:) = 0._r8
!print*,'aa1.2'
       

               !*** if there is no ice cloud ***
               if( iciwp(i,k) < 1.e-80_r8 .or. diaice .eq. 0._r8) then
                 ext_tamu (:,i,k) = 0._r8
                 abs_tamu (:,i,k) = 0._r8
                 ssa_tamu (:,i,k) = 0._r8
                 xmomc_tamu (:,:,i,k) = 0._r8
                 iceind = 0
       
               !*** if there is an ice cloud layer ***
               else
!print*,'iciwp=',iciwp(i,k),', rei=',rei(i,k),', diaice=',diaice

                 diaice = min(diaice, 370._r8)  ! upper bound of ice effective diameter of TAMU scheme is 370 microns
                                                ! because mass absorption coefficient larger than 370um vary slightly with diameter,
                                                ! apply mass absorption coefficient of 370um to those larger ones is not a bad approximation
                                                ! Chia-Pang Kuo @ TAMU, personal communication
       
                 ! if ice effective diameter in approciate range (3-370 microns)
                 if (diaice .ge. 3.0_r8 .and. diaice .le. 370.0_r8) then
                    invrad = 1.0/diaice - 0.04_r8
         
                    do ib=1,nbndlw           
                      if (diaice .ge. 25.0_r8) then
                         extcoice(ib) = ((extice4(1,1,ib) * invrad + &
                                          extice4(1,2,ib)) * invrad + &
                                          extice4(1,3,ib)) * invrad + &
                                          extice4(1,4,ib)
                         ssacoice(ib) = ((ssaice4(1,1,ib) * invrad + &
                                          ssaice4(1,2,ib)) * invrad + &
                                          ssaice4(1,3,ib)) * invrad + &
                                          ssaice4(1,4,ib)
                         asycoice(ib) = ((asyice4(1,1,ib) * invrad + &
                                          asyice4(1,2,ib)) * invrad + &
                                          asyice4(1,3,ib)) * invrad + &
                                          asyice4(1,4,ib)
       
                      else
                         extcoice(ib) = ((((((extice4(2,1,ib) * invrad + &
                                              extice4(2,2,ib)) * invrad + &
                                              extice4(2,3,ib)) * invrad + &
                                              extice4(2,4,ib)) * invrad + &
                                              extice4(2,5,ib)) * invrad + &
                                              extice4(2,6,ib)) * invrad + &
                                              extice4(2,7,ib)) * invrad + &
                                              extice4(2,8,ib)
                         ssacoice(ib) = ((((((ssaice4(2,1,ib) * invrad + &
                                              ssaice4(2,2,ib)) * invrad + &
                                              ssaice4(2,3,ib)) * invrad + &
                                              ssaice4(2,4,ib)) * invrad + &
                                              ssaice4(2,5,ib)) * invrad + &
                                              ssaice4(2,6,ib)) * invrad + &
                                              ssaice4(2,7,ib)) * invrad + &
                                              ssaice4(2,8,ib)
                         asycoice(ib) = ((((((asyice4(2,1,ib) * invrad + &
                                              asyice4(2,2,ib)) * invrad + &
                                              asyice4(2,3,ib)) * invrad + &
                                              asyice4(2,4,ib)) * invrad + &
                                              asyice4(2,5,ib)) * invrad + &
                                              asyice4(2,6,ib)) * invrad + &
                                              asyice4(2,7,ib)) * invrad + &
                                              asyice4(2,8,ib)
                      endif  ! end if of diaice .ge. 25.
                    enddo    ! end do of lwbands for cloud radiative coefficients
                    iceind = 2
       
                    do ib=1,nbndlw
                       taucloud = iciwp(i,k) * 1000._r8 * extcoice(icb(ib,iceind))  ! change iwp from kg/m2 to g/m2 then compute ice cloud optical thickness
                       ssacloud = ssacoice(icb(ib,iceind))
                       ! delta-scaling (Liou, 2002, 313p)
                       asycloud = asycoice(icb(ib,iceind))
                       fpeak = asycloud * asycloud
       
                       abs_tamu  (ib,i,k)    = taucloud * (1._r8-ssacloud)  ! absorption optical depth
       
                       taucloud = (1._r8-ssacloud*fpeak) * taucloud   ! delta-scaling technique, ref: Joseph, Wiscombe and Weinman (1976, JAS)
                       ssacloud = (1._r8-fpeak)*ssacloud / &
                                          (1._r8-ssacloud*fpeak)
                       asycloud = (asycloud-fpeak) / (1.0_r8-fpeak)
                       xmomcloud(0) = 1._r8
                       xmomcloud(1) = asycloud
       
                       ext_tamu  (ib,i,k)    = taucloud       ! cloud extinction optical depth, i.e. including absorption and scattering
                       ssa_tamu  (ib,i,k)    = ssacloud       ! single scattering albedo
                       xmomc_tamu(0:,ib,i,k) = xmomcloud(0:)  ! asymmetric factor
                    enddo    ! end do of lwbands for cloud radiative coefficients
         
                 ! if ice effective diameter is less than 3 microns
                 else
                   !call endrun('ice effective diameter is less than lower boundl 3 microns) in TAMU cloud ice optics')
         
                 endif ! end if of ice effective diameter range
               endif   ! end if of ice cloud layer
             enddo     ! end do of k
          enddo        ! end do of i
       
       
       !---------
       ! return
       !---------
          ext_od0  (:,:,:) = ext_tamu  (:,:,:) 
          abs_od0  (:,:,:) = abs_tamu  (:,:,:) ! return abs_tamu to return variable, abs_od
          ssa_od0  (:,:,:) = ssa_tamu  (:,:,:)
          xmomc_od0(:,:,:,:) = xmomc_tamu(:,:,:,:)

!print*,'aa1.4'
!print*,nlay,ncol,nsubcol

      do ilev = 1,nlay
         do i = 1,ncol
            do isubcol = 1,nsubcol
                  n = ngb(isubcol)
                  ext_od(isubcol,i,ilev) = ext_od0(n,i,ilev)
                  abs_od(isubcol,i,ilev) = abs_od0(n,i,ilev)
                  ssa_od(isubcol,i,ilev) = ssa_od0(n,i,ilev)
                  xmomc_od(:,isubcol,i,ilev) = xmomc_od0(:,n,i,ilev)
            enddo
         enddo
      enddo

       
       !  open(11,file="/scratch/climate_flux/yihsuan/b0-E2000_rrtmg_emis-step5/SourceMods/src.cam/xx.tamu_mc6",form="formatted",position="append")
       !    write(11,*) 'yaya tamu'
       !    write(11,*) 'i,k,ib,iciwp(i,k),abs_od(ib,i,k)'
       !    do i=1,ncol
       !    do ib=1,nbndlw
       !        !write(11,*) k,ib,abs_od(ib,i,k),abs_tamu(ib,i,k)
       !        write(11,*) 'aaa ext_od,ncol,ib',i,ib
       !        write(11,*) ext_od(ib,i,:)
       !        write(11,*) 'aaa abs_od,ncol,ib',i,ib
       !        write(11,*) abs_od(ib,i,:)
       !        write(11,*) 'aaa ssa_od,ncol,ib',i,ib
       !        write(11,*) ssa_od(ib,i,:)
       !        write(11,*) 'aaa xmomc_od(0),ncol,ib',i,ib
       !        write(11,*) xmomc_od(0,ib,i,:)
       !        write(11,*) 'aaa xmomc_od(1),ncol,ib',i,ib
       !        write(11,*) xmomc_od(1,ib,i,:)
       !    enddo
       !    enddo
       !  close(11)
       
       !open(11,file='/scratch/climate_flux/yihsuan/b0-E2000_rrtmg_emis-step5/SourceMods/src.cam/xx.lwice.mc6', &
       !        form='formatted',position='append')
       !  write(11,*) 'bbb mc6, abs_od',i_cld
       !  do i=1,ncol
       !    !write(11,*) 'ncol,iciwp(kg/m2)',i
       !    !write(11,*) iciwp(i,:)
       !    !write(11,*) 'ncol,cldn(kg/m2)',i
       !    !write(11,*) cldn(i,:)
       !  do k=1,nbndlw
       !    write(11,*) 'ncol=',i,', nlwband=',k
       !    write(11,*) abs_od(k,i,:)
       !  enddo
       !  enddo
       !close(11)
       
         return
       end subroutine tamu_ice_get_rad_props_lw_OFFLINE
       
