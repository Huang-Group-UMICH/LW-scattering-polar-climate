      use parrrtm
      use rrtmg_lw_init, only : rrtmg_lw_ini
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none
      integer, parameter :: ncol = 1
      integer, parameter :: nlay = 1
     
      integer i,j,k,k_cloud,idei
      real(kind=r8) :: dei, iciwp 
      real(kind=r8) :: taucmcl(nbndlw,ncol,nlay)
      real(kind=r8) :: abscmcl(nbndlw,ncol,nlay)
      real(kind=r8) :: extcmcl(nbndlw,ncol,nlay)
      real(kind=r8) :: ssacmcl(nbndlw,ncol,nlay)
      real(kind=r8) :: xmomcmcl(0:16,nbndlw,ncol,nlay)
      real(kind=r8) :: taucavg
      character*30  :: band_range(16)
      data band_range/"10-350" , "350-500", "500-630" , "630-700"  , &
                      "700-820", "820-980", "980-1080", "1080-1180", &
                      "1180-1390", "1390-1480", "1480-1800", "1800-2080",  &
                      "2080-2250", "2250-2380", "2380-2600", "2600-3250"/
      real(kind=r8) :: deis(5)
      data deis/10._r8, 20._r8, 40._r8, 100._r8, 200._r8/ 
      integer band

!----------------------------
      iciwp = 0.1_r8 ! in-cloud ice water path (kg/m2)
      !iciwp = 0.05_r8 ! in-cloud ice water path (kg/m2)
      band = 6        ! print this band 

      do j=1,5
        dei = deis(j)

        !dei = 5._r8    ! ice effectice diameter (microns)
        !dei = 5._r8    ! ice effectice diameter (microns)
        !dei = 100._r8    ! ice effectice diameter (microns)
      
        call rrtmg_lw_ini
        call mc6_ice_get_rad_props_lw_OFFLINE(ncol, nlay, ngptlw, dei, iciwp, taucmcl,abscmcl, ssacmcl, xmomcmcl)

        taucavg = sum(taucmcl(:,1,1))/nbndlw
        !print*,'dei=',dei,'microns, iciwp=',iciwp,'kg/m2, band-avg taucloud=',taucavg
        print*,'dei=',dei,'microns, iciwp=',iciwp,'kg/m2, taucloud=',taucmcl(band,1,1),'at band ',band,' ',band_range(band)
        print*,''
      !do i=1,16
      !  print*,'band=',i,band_range(i),', tau=',taucmcl(i,1,1)
      !enddo
      enddo ! j

      extcmcl = taucmcl/iciwp/1000. ! unit: m2/g

      !print*,'taucmcl',taucmcl
      !print*,'extcmcl',extcmcl
      !print*,'ssacmcl',ssacmcl
      !print*,'xmomcmcl',xmomcmcl(1,:,:,:)

!      open(12,file="./zz-tau_ext_ssa_asy.txt",form="formatted")
!        write(12,*) 'ice effective diameter (micons)'
!        write(12,*) dei
!        write(12,*) 'in-cloud ice water path (kg/m2)'
!        write(12,*) iciwp
!        write(12,*) ''
!        write(12,*) '-----------------'
!        write(12,*) 'band 1 to 16'
!        write(12,*) '-----------------'
!        write(12,*) ''
!        write(12,*) 'cloud optical depth'
!        write(12,'(16(F7.4,2X))') taucmcl
!        write(12,*) ''
!        write(12,*) 'extinction coefficient (m2/g)'
!        write(12,'(16(F7.4,2X))') extcmcl
!        write(12,*) ''
!        write(12,*) 'single scattering albedo'
!        write(12,'(16(F7.4,2X))') ssacmcl
!        write(12,*) ''
!        write(12,*) 'Asymmetric factor'
!        write(12,'(16(F7.4,2X))') xmomcmcl(1,:,:,:)
!      close(12)

9999  continue
      stop
      end

      !subroutine mc6_ice_get_rad_props_lw_OFFLINE(pcols, pver, nsubcol, dei, iciwp, ext_od, abs_od, ssa_od, xmomc_od)
      subroutine mc6_ice_get_rad_props_lw_OFFLINE(pcols, pver, nsubcol, dei, iciwp, ext_od0, abs_od0, ssa_od0, xmomc_od0)
       
       !-----------------------------------------------------
       ! yihsuan@UMich
       !
       ! Description:
       !   MC6 ice cloud optics scheme
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
          real(r8), intent(in) :: dei(pcols,pver)     ! ice effective diameter (microns)
          real(r8), intent(in) :: iciwp(pcols,pver) ! in-cloud ice water path (kg/m2)

       !--- Output arguments ---!
          real(r8) :: ext_od  (nsubcol,pcols,pver)        ! cloud ice extinction optical depth, i.e. including absorption and scattering
          real(r8) :: abs_od  (nsubcol,pcols,pver)        ! cloud ice absorption optical depth
          real(r8) :: ssa_od  (nsubcol,pcols,pver)        ! single scattering albedo
          real(r8) :: xmomc_od(0:16,nsubcol,pcols,pver)   ! phase function
       
       !--- Output arguments ---!
          real(r8), intent(out) :: ext_od0  (nbndlw,pcols,pver)        ! cloud ice extinction optical depth, i.e. including absorption and scattering
          real(r8), intent(out) :: abs_od0  (nbndlw,pcols,pver)        ! cloud ice absorption optical depth
          real(r8), intent(out) :: ssa_od0  (nbndlw,pcols,pver)        ! single scattering albedo
          real(r8), intent(out) :: xmomc_od0(0:16,nbndlw,pcols,pver)   ! phase function
       
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
       
          !integer :: dei_idx, iciwp_idx, cld_idx, err, itim
          integer :: dei_idx, iciwp_idx, err, itim
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
               diaice = dei(i,k)  ! effective diameter is defined by that particle-volume divided by particle-projected-area (i.e. effective diameter) times 1.5.
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
!print*,'iciwp=',iciwp(i,k),', dei=',dei(i,k),', diaice=',diaice

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
       end subroutine mc6_ice_get_rad_props_lw_OFFLINE
       
