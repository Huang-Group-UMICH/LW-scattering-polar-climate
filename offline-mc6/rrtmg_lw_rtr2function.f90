!+
! rrtmg_lw_rtr2function.f90
! This program prepares input parameters for two-stream source function
! technique.
!-
! Reference
! Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam (1989),
! Rapid calculation of radiative heating rates and photodissociation
! rates in inhomogeneous multiple scattering atmospheres, J. Geophys.
! Res., 94(D13), 16287–16301, doi:10.1029/JD094iD13p16287.
!+
! History
! Oct. 17, 2016 Develop the code (Chia-Pang Kuo)
!-
! Contact Information
!  Chia-Pang Kuo @ Texas A&M University
!  email: intelb52158@tamu.edu

      module rrtmg_lw_rtr2function

      use shr_kind_mod, only: rb => shr_kind_r8
      !use parkind, only : im => kind_im, rb => kind_rb
      use rrlw_con, only: heatfac
      use rrlw_wvn, only: ngs,delwave

      implicit none

      contains

!*************************************************************************
      subroutine rtr2function(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                              !cldfrac, taucloud, ssacloud, xmomcloud, &
                              planklay, planklev, plankbnd, &
                              pwvcm, fracs, taut, &
                              totuflux, totdflux, fnet, htr, &
                              totuclfl, totdclfl, fnetc, htrc, totufluxs, totdfluxs)
!*************************************************************************

! ----- Input -----
      integer, intent(in) :: nlayers         ! total number of layers
      integer, intent(in) :: istart          ! beginning band of calculation
      integer, intent(in) :: iend            ! ending band of calculation
      integer, intent(in) :: iout            ! output option flag

      integer, parameter :: nbndlw = 16
      integer :: ncbands

! Atmosphere
      real(kind=rb), intent(in) :: pz(0:)             ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(in) :: pwvcm              ! precipitable water vapor (cm)
      real(kind=rb), intent(in) :: semiss(:)          ! lw surface emissivity
                                                      !    Dimensions: (nbndlw)
      real(kind=rb), intent(in) :: fracs(:,:)         ! 
                                                      !    Dimensions: (nlayers,ngptlw)
      real(kind=rb), intent(in) :: taut(:,:)          ! gaseous + aerosol optical depths
                                                      !    Dimensions: (nlayers,ngptlw)
      real(kind=rb), intent(in) :: planklay(:,:)      ! 
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=rb), intent(in) :: planklev(0:,:)     ! 
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(in) :: plankbnd(:)        ! 
                                                      !    Dimensions: (nbndlw)

! Clouds
!      integer, intent(in) :: ncbands         ! number of cloud spectral bands
!                                                      ! Planck derivative [0=off, 1=on]
!      real(kind=rb), intent(in) :: cldfrac(:)         ! layer cloud fraction
!                                                      !    Dimensions: (nlayers)
!      real(kind=rb), intent(in) :: taucloud(:,:)      ! layer cloud optical depth
!                                                      !    Dimensions: (nlayers,nbndlw)
!      real(kind=rb), intent(in) :: ssacloud(:,:)      ! layer cloud single-scattering albedo
!                                                      !    Dimensions: (nlayers,nbndlw)
!      real(kind=rb), intent(in) :: xmomcloud(0:,:,:)  ! layer cloud expansion coefficients of phase function
                                                      !    Dimensions: (0:16,nlayers,nbndlw)

!flag11
      real(kind=rb) :: cldfrac(nlayers)         ! layer cloud fraction
                                                      !    Dimensions:
                                                      !    (nlayers)
      real(kind=rb) :: taucloud(nlayers,nbndlw)      ! layer cloud optical depth
                                                      !    Dimensions:
                                                      !    (nlayers,nbndlw)
      real(kind=rb) :: ssacloud(nlayers,nbndlw)      ! layer cloud single-scattering albedo
                                                      !    Dimensions:
                                                      !    (nlayers,nbndlw)
      real(kind=rb) :: xmomcloud(17,nlayers,nbndlw)  ! layer cloud expansion coefficients of phase function
                                                      !    Dimensions:
                                                      !    (0:16,nlayers,nbndlw)
!flag11




! ----- Output -----
      real(kind=rb), intent(out) :: totuflux(0:)      ! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdflux(0:)      ! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: fnet(0:)          ! net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: htr(0:)           ! longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuclfl(0:)      ! clear sky upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdclfl(0:)      ! clear sky downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: fnetc(0:)         ! clear sky net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: htrc(0:)          ! clear sky longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      !>>> yihsuan 2017-06-06, add spectral flux >>>
      real(kind=rb), intent(out) :: totufluxs(:,0:)  ! upward longwave spectral flux (w/m2)
                                                          !    Dimensions: (nbndlw,0:nlayers)

      real(kind=rb), intent(out) :: totdfluxs(:,0:)  ! downward longwave spectral flux (W/m2)
                                                          !    Dimensions: (nbndlw,0:nlayers)
      !<<< yihsuan 2017-06-06, add spectral flux <<<
! ----- Local -----
! Declarations for radiative transfer
      integer :: iband, lay, lev, ig   ! loop indices
      integer :: ibcnd

      real(kind=rb) :: wavenumlo, wavenumhi
      real(kind=rb) :: plkavg
      real(kind=rb) :: diffus
      real(kind=rb) :: albedosuf
      real(kind=rb) :: dnftoa
      real(kind=rb) :: scatcld
      real(kind=rb) :: planksuf
      real(kind=rb) :: gama1(nlayers)
      real(kind=rb) :: gama2(nlayers)
      real(kind=rb) :: fluxupcld(0:nlayers) ! Upward flux under cloudy sky
      real(kind=rb) :: fluxupclr(0:nlayers) ! Upward flux under clear sky
      real(kind=rb) :: fluxdncld(0:nlayers) ! Downward flux under cloudy sky
      real(kind=rb) :: fluxdnclr(0:nlayers) ! Downward flux under clear sky
      real(kind=rb) :: taurevcld(nlayers)
      real(kind=rb) :: taurevclr(nlayers)
      real(kind=rb) :: ssarevcld(nlayers)
      real(kind=rb) :: ssarevclr(nlayers)
      real(kind=rb) :: asyrevcld(nlayers)
      real(kind=rb) :: asyrevclr(nlayers)
      real(kind=rb) :: plankrev(0:nlayers)
      real(kind=rb) :: fracsrev(nlayers)

      real(kind=rb), parameter :: pi = 3.1415926535897932_rb


!flag11
      cldfrac = 0._rb
      taucloud = 0._rb
      ssacloud = 0._rb
      xmomcloud = 0._rb

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    ssacloud                     ! layer cloud single-scattering albedo
!    xmomcloud                    ! layer cloud expansion coefficients of phase function
!    semiss                       ! surface emissivities for each band
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    odclr                        ! clear sky (gaseous) optical depth
!    bbdgas                       ! gas-only planck function for downward rt
!    d_urad_dt                    ! upward radiance by layer
!    d_clrurad_dt                 ! clear sky upward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)

      !>>> yihsuan 2017-06-06, add spectral flux >>>
      totufluxs = 0.0_rb 
      totdfluxs = 0.0_rb 
      !<<< yihsuan 2017-06-06, add spectral flux <<<

      totuclfl = 0.0_rb
      totdclfl = 0.0_rb
      totuflux = 0.0_rb
      totdflux = 0.0_rb
      fnetc = 0.0_rb
      htrc = 0.0_rb
      fnet = 0.0_rb
      htr = 0.0_rb
      taurevcld = 0.0_rb
      taurevclr = 0.0_rb
      ssarevcld = 0.0_rb
      ssarevclr = 0.0_rb
      asyrevcld = 0.0_rb
      asyrevclr = 0.0_rb
      plankrev = 0.0_rb
      planksuf = 0.0_rb
      albedosuf = 0.0_rb
      fracsrev = 0.0_rb
      scatcld = 0.0_rb
      gama1 = 0.0_rb
      gama2 = 0.0_rb
      dnftoa = 0.0_rb
      
      ig = 1
! *** loop over frequency bands.
      do iband = istart, iend

! ***    planck function for each level
         do lay = 0, nlayers
            plankrev(nlayers-lay) = planklev(lay,iband) * 1.0e4_rb * delwave(iband)
         enddo

! ***    planck function at the surface
         planksuf = plankbnd(iband) * 1.0e4_rb * delwave(iband)

! ***    set surface albedo for this band
         albedosuf = 1.0_rb - semiss(iband)

! ***    loop over g-channels.
         if (iout.gt.0.and.iband.ge.2) ig = ngs(iband-1)+1

 1000    continue

! ***    downward radiative transfer.

         do lay = nlayers, 1, -1

! ***       fraction of planck function
            fracsrev(nlayers-lay+1) = fracs(lay,ig)

! ***       clear sky
            if (taut(lay,ig) .lt. 0.0_rb) then
               taurevclr(nlayers-lay+1) = 0.0_rb
            else
               taurevclr(nlayers-lay+1) = taut(lay,ig)
            endif
            ssarevclr(nlayers-lay+1) = 0.0_rb
            asyrevclr(nlayers-lay+1) = 0.0_rb
            
! ***       mix optical properties of cloud and gas
            if (taut(lay,ig) .lt. 0.0_rb) then
               taurevcld(nlayers-lay+1) = taucloud(lay,iband)
            else
               taurevcld(nlayers-lay+1) = taut(lay,ig) + taucloud(lay,iband)
            endif
            scatcld = ssacloud(lay,iband) * taucloud(lay,iband)
            asyrevcld(nlayers-lay+1) = xmomcloud(1,lay,iband)

            if (taurevcld(nlayers-lay+1) .ne. 0.0) &
                 ssarevcld(nlayers-lay+1) = scatcld / taurevcld(nlayers-lay+1)
            
            if (ssarevcld(nlayers-lay+1) .gt. 1.0) then
               write(*,'(A,F10.5,A,I10,A,I10)') '!! warning ssarevcld: ', &
                            ssarevcld(nlayers-lay+1),' > 1.0 at layer ', lay, &
                            'at g point', ig
            endif
            if (ssarevcld(nlayers-lay+1) .lt. 0.0) then
               write(*,'(A,F10.5,A,I10,A,I10)') '!! warning ssarevcld: ', &
                            ssarevcld(nlayers-lay+1),' < 0.0 at layer ', lay, &
                            'at g point', ig
            endif

         enddo

! ***    Cloudy sky

! ***    Hemispheric mean
         diffus = 1.66_rb
         gama1 = diffus * 0.5_rb * (2.0_rb - ssarevcld*(1.0_rb+asyrevcld))
         gama2 = diffus * 0.5_rb * ssarevcld * (1.0_rb-asyrevcld)

! ***    boundary conditions
         dnftoa  = 0.0_rb
         
         call TwsFunctionIR(nlayers,diffus,taurevcld,ssarevcld,fracsrev,plankrev,planksuf,&
                            gama1,gama2,albedosuf,dnftoa,fluxupcld,fluxdncld)

         ! Set unphysical values to 0.0
         do lev = 0, nlayers
            if (fluxdncld(lev) < 0.0) fluxdncld(lev) = 0.0_rb
            if (fluxupcld(lev) < 0.0) fluxupcld(lev) = 0.0_rb
         enddo

         do lev = nlayers, 0, -1
            totuflux(lev) = totuflux(lev) + fluxupcld(nlayers-lev)
            totdflux(lev) = totdflux(lev) + fluxdncld(nlayers-lev) 

            !>>> yihsuan 2017-06-06, add spectral flux >>>
            totufluxs(iband,lev)= totufluxs(iband,lev) + fluxupcld(nlayers-lev)  ! upward spectral flux
            totdfluxs(iband,lev)= totdfluxs(iband,lev) + fluxdncld(nlayers-lev)  ! downward spectral flux
            !<<< yihsuan 2017-06-06, add spectral flux <<<
         enddo
         
         if (fluxdncld(0) .gt. 1.e-5) &
              write(*,9000) fluxdncld(0), iband, ig
         do lev = 0, nlayers
            if (fluxdncld(lev) < 0.0) &
                 write(*,9001) fluxdncld(lev), nlayers-lev, iband, ig
         enddo

! *** Clear sky

! ***    Hemispheric mean
         diffus = 1.66_rb
         gama1 = diffus * 0.5_rb * (2.0_rb - ssarevclr*(1.0_rb+asyrevclr))
         gama2 = diffus * 0.5_rb * ssarevclr * (1.0_rb-asyrevclr)

! ***    boundary conditions
         dnftoa  = 0.0_rb
         
         call TwsFunctionIR(nlayers,diffus,taurevclr,ssarevclr,fracsrev,plankrev,planksuf,&
                            gama1,gama2,albedosuf,dnftoa,fluxupclr,fluxdnclr)

         ! Set unphysical values to 0.0
         do lev = 0, nlayers
            if (fluxdnclr(lev) < 0.0) fluxdnclr(lev) = 0.0_rb
            if (fluxupclr(lev) < 0.0) fluxupclr(lev) = 0.0_rb
         enddo

         do lev = nlayers, 0, -1
            totuclfl (lev) = totuclfl(lev) + fluxupclr(nlayers-lev)
            totdclfl (lev) = totdclfl(lev) + fluxdnclr(nlayers-lev) 
            !>>> yihsuan 2017-06-06, add spectral flux >>>
            !totufluxsclr(iband,lev)= fluxupclr(nlayers-lev)  ! upward spectral flux
            !totdfluxsclr(iband,lev)= fluxdnclr(nlayers-lev)  ! downward spectral flux
            !<<< yihsuan 2017-06-06, add spectral flux <<<
         enddo
         
         if (fluxdnclr(0) .gt. 1.e-5) &
              write(*,9000) fluxdnclr(0), iband, ig
         do lev = 0, nlayers
            if (fluxdnclr(lev) < 0.0) &
                 write(*,9001) fluxdnclr(lev), nlayers-lev, iband, ig
         enddo

         ig = ig + 1
         if (ig .le. ngs(iband)) go to 1000
            
      enddo
      
      fnet(nlayers)  = totuflux(nlayers) - totdflux(nlayers)
      fnetc(nlayers) = totuclfl(nlayers) - totdclfl(nlayers)
      htr(nlayers)  = 0.0_rb
      htrc(nlayers) = 0.0_rb
      do lev = nlayers-1, 0, -1
         fnet(lev)  = totuflux(lev) - totdflux(lev)
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
         htr(lev)  = heatfac * (fnet(lev) -fnet(lev+1)) / (pz(lev) - pz(lev+1))
         htrc(lev) = heatfac * (fnetc(lev) -fnetc(lev+1)) / (pz(lev) - pz(lev+1))
      enddo

 9000 format('DOWNWARD FLUX ',ES15.7,' AT TOA GTR THAN 0. IN BAND ',I3, &
      ' AT IG =',I3,'. POSSIBLE', &
      ' INSTABILITY IN TwoStream.')
 9001 format('DOWNWARD FLUX ',ES15.7,' AT LEVEL ',I3,&
      ' GTR THAN 0. IN BAND ',I2, &
      ' AT IG =',I3,'. POSSIBLE',/, &
      ' INSTABILITY IN TwoStream.')

      end subroutine rtr2function
      
!*************************************************************************
      subroutine TwsFunctionIR(nlayers,difactor,taulay,ssalay,planckfracs,planckflev,planckfsuf,&
                               gama1,gama2,sufalb,dnf0,upf,dnf)
!*************************************************************************

!+
! TwsFunction
! This program is two-stream source function technique.
!-
! Reference
! Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam (1989),
! Rapid calculation of radiative heating rates and photodissociation
! rates in inhomogeneous multiple scattering atmospheres, J. Geophys.
! Res., 94(D13), 16287–16301, doi:10.1029/JD094iD13p16287.
!+
! History
! Nov. 11, 2016 Develop the code (Chia-Pang Kuo)
! Mar. 29, 2017 Fix a bug in paramH and paramG (Chia-Pang Kuo)
!-
! Contact Information
!  Chia-Pang Kuo @ Texas A&M University
!  email: intelb52158@tamu.edu

!  Structure of atmosphere
!
!  Top of the atmosphere
!  -------------------- level 0
!         layer 1
!  -------------------- level 1
!         layer 2
!  -------------------- level 2
!            .
!            .
!            .
!  -------------------- level N-2
!         layer N-1
!  -------------------- level N-1
!         layer N
!  -------------------- level N
!         Surface

! ----- Input -----
      integer, intent(in) :: nlayers ! number of atmopheric layers

      real(kind=rb), intent(in) :: difactor        ! diffusivity factor
      real(kind=rb), intent(in) :: dnf0            ! downward flux at the top of the atmosphere (boundary condition)
      real(kind=rb), intent(in) :: sufalb          ! surface albedo
      real(kind=rb), intent(in) :: taulay(:)       ! optical thickness in each layer
                                                   ! dimension(layer)
      real(kind=rb), intent(in) :: ssalay(:)       ! single scattering albedo in each layer
                                                   ! dimension(layer)
      real(kind=rb), intent(in) :: gama1(:)        ! gamma 1
                                                   ! dimension(layer)
      real(kind=rb), intent(in) :: gama2(:)        ! gamma 2
                                                   ! dimension(layer)
      real(kind=rb), intent(in) :: planckflev(0:)  ! plank function in each level
                                                   ! dimension(0:layer)
      real(kind=rb), intent(in) :: planckfracs(:)  ! 
                                                   ! dimensions: (layer)
      real(kind=rb), intent(in) :: planckfsuf      ! plank function on the ground
! ----- Output -----
      real(kind=rb), intent(out) :: upf(0:) ! upward flux at the each level
                                            ! dimension(0:layer)
      real(kind=rb), intent(out) :: dnf(0:) ! downward flux at the each level
                                            ! dimension(0:layer)
! ----- Local -----
      integer, parameter :: nangle = 2 ! number of double gauss quadrature
      integer :: iangle
      integer :: lay      ! atmopheric layer
      integer :: lev      ! atmopheric level
      real(kind=rb) :: difactorev  ! inverse diffusivity factor
      real(kind=rb) :: upfn        ! upward flux at the surface (boundary condition)
      real(kind=rb) :: planktop
      real(kind=rb) :: plankbas
      real(kind=rb) :: lamda(nlayers) 
      real(kind=rb) :: gama(nlayers) 
      real(kind=rb) :: exp1(nlayers)
      real(kind=rb) :: exp2(nlayers)
      real(kind=rb) :: e1n(nlayers)
      real(kind=rb) :: e2n(nlayers)
      real(kind=rb) :: e3n(nlayers)
      real(kind=rb) :: e4n(nlayers)
      real(kind=rb) :: mu1
      real(kind=rb) :: b0(nlayers)
      real(kind=rb) :: b1(nlayers)
      real(kind=rb) :: tempc(nlayers)
      real(kind=rb) :: upc0(nlayers)
      real(kind=rb) :: upcn(nlayers)
      real(kind=rb) :: dnc0(nlayers)
      real(kind=rb) :: dncn(nlayers)
      real(kind=rb) :: upc_tmp(nlayers)
      real(kind=rb) :: dnc_tmp(nlayers)
      real(kind=rb) :: aa(2*nlayers)
      real(kind=rb) :: bb(2*nlayers)
      real(kind=rb) :: cc(2*nlayers)
      real(kind=rb) :: dd(2*nlayers)
      real(kind=rb) :: ee(2*nlayers)
      real(kind=rb) :: xx(2*nlayers)
      real(kind=rb) :: k1n(nlayers)
      real(kind=rb) :: k2n(nlayers)
      real(kind=rb) :: dgausa(nangle)
      real(kind=rb) :: dgausw(nangle)
      real(kind=rb) :: upi(0:nlayers)
      real(kind=rb) :: dni(0:nlayers)
      real(kind=rb) :: cosangle
      real(kind=rb) :: cosanglerev
      real(kind=rb) :: param1(nlayers)
      real(kind=rb) :: param2(nlayers)
      real(kind=rb) :: paramG(nlayers)
      real(kind=rb) :: paramK(nlayers)
      real(kind=rb) :: paramJ(nlayers)
      real(kind=rb) :: paramH(nlayers)
      real(kind=rb) :: alpha1(nlayers)
      real(kind=rb) :: alpha2(nlayers)
      real(kind=rb) :: sigma1(nlayers)
      real(kind=rb) :: sigma2(nlayers)
      real(kind=rb) :: exptmp1(nlayers)
      real(kind=rb) :: exptmp2(nlayers)
      real(kind=rb) :: revpar1(nlayers)
      real(kind=rb) :: revpar2(nlayers)
      
      real(kind=rb), parameter :: pi = 3.1415926535897932_rb

      upf     = 0.0_rb
      dnf     = 0.0_rb
      lamda   = 0.0_rb
      gama    = 0.0_rb
      exp1    = 0.0_rb
      exp2    = 0.0_rb
      e1n     = 0.0_rb
      e2n     = 0.0_rb
      e3n     = 0.0_rb
      e4n     = 0.0_rb
      mu1     = 0.0_rb
      b0      = 0.0_rb
      b1      = 0.0_rb
      tempc   = 0.0_rb
      upc0    = 0.0_rb
      upcn    = 0.0_rb
      dnc0    = 0.0_rb
      dncn    = 0.0_rb
      upc_tmp = 0.0_rb
      dnc_tmp = 0.0_rb
      aa      = 0.0_rb
      bb      = 0.0_rb
      cc      = 0.0_rb
      dd      = 0.0_rb
      ee      = 0.0_rb
      xx      = 0.0_rb
      k1n     = 0.0_rb
      k2n     = 0.0_rb
      upi     = 0.0_rb
      dni     = 0.0_rb
      cosangle = 0.0_rb
      cosanglerev = 0.0_rb
      param1 = 0.0_rb
      param2 = 0.0_rb
      paramG = 0.0_rb
      paramK = 0.0_rb
      paramJ = 0.0_rb
      paramH = 0.0_rb
      alpha1 = 0.0_rb
      alpha2 = 0.0_rb
      sigma1 = 0.0_rb
      sigma2 = 0.0_rb
      exptmp1 = 0.0_rb
      exptmp2 = 0.0_rb
      revpar1 = 0.0_rb
      revpar2 = 0.0_rb
      dgausa  = (/0.2113248_rb, 0.7886752_rb/)
      dgausw  = (/0.5_rb, 0.5_rb/)

      ! lamda in eq. 21
      lamda = sqrt(gama1*gama1 - gama2*gama2)

      ! gamma in eq. 22
      gama = gama2 / (gama1 + lamda)
      
      exp1 = exp(lamda*taulay)
      exp2 = exp(-lamda*taulay)
    
      ! exponential functions in eq. 44
      e1n = 1.0_rb + gama*exp2
      e2n = 1.0_rb - gama*exp2
      e3n = gama + exp2
      e4n = gama - exp2

      ! constant in eq. 18 for Hemispheric mean 
      ! (mu = diffusivity factor * inv(diffusivity factor) = 1.0)
      mu1 = 1.0_rb

      ! approximate Planck function by first two terms of Taylor expansion in eq. 25
      ! bn(tau) = b0n + b1n*tau
      do lay = 1, nlayers 
         ! weighted planck function for a g interval
         planktop = planckflev(lay-1) * planckfracs(lay)
         plankbas = planckflev(lay) * planckfracs(lay)
         if (taulay(lay) .lt. 1.0e-6_rb) then
            b1(lay) = 0.0_rb
            b0(lay) = planktop
         else
            b1(lay) = (plankbas - planktop) / taulay(lay)
            b0(lay) = planktop
         endif
      enddo

      ! particular solution in eq. 27
      tempc = 1.0_rb / (gama1 + gama2)
      upc0 = b0 + b1*(0.0_rb+tempc) ! upward direction at the top of the layer
      upcn = b0 + b1*(taulay+tempc) ! upward direction at the bottom of the layer 
      dnc0 = b0 + b1*(0.0_rb-tempc) ! downward direction at the top of the layer
      dncn = b0 + b1*(taulay-tempc) ! downward direction at the bottom of the layer 

      ! upCn+1(0) - upCn(tau) in eq. 41 and 42
      upc_tmp(1:nlayers-1)  = upc0(2:nlayers) - upcn(1:nlayers-1)
      ! dnCn+1(0) - dnCn(tau) in eq. 41 and 42
      dnc_tmp(1:nlayers-1)  = dnc0(2:nlayers) - dncn(1:nlayers-1)

      ! upward flux at the surface (boundary condition)
      upfn = planckfsuf * planckfracs(nlayers)

      ! develop tridiagonal matrix (aa)Yn-1 + (bb)Yn + (dd)Yn+1 = ee in eq. 39
      aa(1) = 0.0_rb
      bb(1) = e1n(1)
      dd(1) = -1.0_rb*e2n(1)
      ee(1) = dnf0/pi - dnc0(1)

      aa(2*nlayers) = e1n(nlayers) - sufalb*e3n(nlayers)
      bb(2*nlayers) = e2n(nlayers) - sufalb*e4n(nlayers)
      dd(2*nlayers) = 0.0_rb
      ee(2*nlayers) = upfn - upcn(nlayers) + sufalb*dncn(nlayers)

      ! even element in a tridiagonal matrix
      aa(2:2*nlayers-2:2) = e1n(1:nlayers-1)*e2n(2:nlayers) - &
                            e3n(1:nlayers-1)*e4n(2:nlayers)
      bb(2:2*nlayers-2:2) = e2n(1:nlayers-1)*e2n(2:nlayers) - &
                            e4n(1:nlayers-1)*e4n(2:nlayers)
      dd(2:2*nlayers-2:2) = e1n(2:nlayers)*e4n(2:nlayers) - &
                            e2n(2:nlayers)*e3n(2:nlayers)
      ee(2:2*nlayers-2:2) = upc_tmp(1:nlayers-1)*e2n(2:nlayers) - &
                            dnc_tmp(1:nlayers-1)*e4n(2:nlayers)
      ! odd element in a tridiagonal matrix
      aa(3:2*nlayers-1:2) = e2n(1:nlayers-1)*e3n(1:nlayers-1) - &
                            e4n(1:nlayers-1)*e1n(1:nlayers-1)
      bb(3:2*nlayers-1:2) = e1n(1:nlayers-1)*e1n(2:nlayers) - &
                            e3n(1:nlayers-1)*e3n(2:nlayers)
      dd(3:2*nlayers-1:2) = e3n(1:nlayers-1)*e4n(2:nlayers) - &
                            e1n(1:nlayers-1)*e2n(2:nlayers)
      ee(3:2*nlayers-1:2) = upc_tmp(1:nlayers-1)*e3n(1:nlayers-1) - &
                            dnc_tmp(1:nlayers-1)*e1n(1:nlayers-1)

      ! use Thomas algorithm to solve tridiagonal matrix
      call TDMA(2*nlayers,aa,bb,dd,ee,xx)

      ! use eq. 29 and 30 to get coefficients
      k1n = xx(1:2*nlayers-1:2) + xx(2:2*nlayers:2)  
      k2n = xx(1:2*nlayers-1:2) - xx(2:2*nlayers:2)

      ! parameter in Table 3
      difactorev = 1.0_rb / difactor
      param1 = (difactor-lamda) * difactorev
      param2 = gama * (lamda+difactor) * difactorev
      paramH = k2n * param2
      paramG = k1n * param1
      paramK = k2n * param1
      paramJ = k1n * param2
      sigma1 = b0 - b1*(tempc-difactorev)
      sigma2 = b1
      alpha1 = b0 + b1*(tempc-difactorev)
      alpha2 = b1

      ! downward intensity for a given direction
      do iangle = 1, nangle ! loop over double gauss quadrature

         ! double gauss quadrature
         cosangle  = dgausa(iangle)
         cosanglerev = 1.0_rb / cosangle

         exptmp1 = exp(-taulay*cosanglerev)
         exptmp2 = exp(-taulay*(lamda+cosanglerev))
         revpar1 = 1.0_rb / (lamda*cosangle + 1.0_rb)
         revpar2 = 1.0_rb / (lamda*cosangle - 1.0_rb)
        
         ! incident downward intensity 
         dni(0) = dnf0 / pi
         ! isotropic source
         dnf(0) = dnf0

         ! downward radiative transfer in eq. 56
         do lev = 1, nlayers
            dni(lev) = dni(lev-1)*exptmp1(lev) + &
                       paramJ(lev)*revpar1(lev)*(1.0_rb-exptmp2(lev)) + &
                       paramK(lev)*revpar2(lev)*(exptmp1(lev)-exp2(lev)) + &
                       sigma1(lev)*(1.0_rb-exptmp1(lev)) + &
                       sigma2(lev)*(cosangle*exptmp1(lev)+taulay(lev)-cosangle)
            ! downward flux at the bottom of the layer
            dnf(lev) = dnf(lev) + 2.0_rb*pi*dni(lev)*cosangle*dgausw(iangle)
         enddo

         upi(nlayers) = upi(nlayers) + cosangle*dgausw(iangle)*dni(nlayers)

      enddo
      
      ! reflected upward intensity for a given direction at the
      ! lambertian surface
      upi(nlayers) = 2.0_rb*sufalb*upi(nlayers) + planckfsuf*planckfracs(nlayers)

      ! upward intensity for a given direction
      do iangle = 1, nangle ! loop over double gauss quadrature

         ! double gauss quadrature
         cosangle  = dgausa(iangle)
         cosanglerev = 1.0_rb / cosangle

         exptmp1 = exp(-taulay*cosanglerev)
         exptmp2 = exp(-taulay*(lamda+cosanglerev))
         revpar1 = 1.0_rb / (lamda*cosangle + 1.0_rb)
         revpar2 = 1.0_rb / (lamda*cosangle - 1.0_rb)

         ! upward radiative transfer in eq. 55
         do lev = nlayers, 1, -1
            upi(lev-1) = upi(lev)*exptmp1(lev) + &
                         paramG(lev)*revpar2(lev)*(exptmp1(lev)-exp2(lev)) + &
                         paramH(lev)*revpar1(lev)*(1.0_rb-exptmp2(lev)) + &
                         alpha1(lev)*(1.0_rb-exptmp1(lev)) + &
                         alpha2(lev)*(cosangle-(cosangle+taulay(lev))*exptmp1(lev))
            ! upward flux at the bottom of the layer
            upf(lev-1) = upf(lev-1) + 2.0_rb*pi*upi(lev-1)*cosangle*dgausw(iangle)
         enddo

         upf(nlayers) = upf(nlayers) + 2.0_rb*pi*upi(nlayers)*cosangle*dgausw(iangle)

      enddo

      end subroutine TwsFunctionIR

!*************************************************************************
      subroutine TDMA(N,A,B,C,D,X)
!*************************************************************************
!+
! TDMA.f90
! This program is Tridiagonal matrix algorithm (TDMA), also called Thomas 
! algorithm.
! A tridiagonal matrix for N unknowns can be written as
! A(i)X(i-1) + B(i)X(i) + C(i)X(i+1) = D(i), where A(1) = 0, and C(N) =
! 0.
! --                           -- --  --   --  --
! |B1   C1                     0| | X1 |   | D1 |
! |A2   B2   C2                 | | X2 |   | D2 |
! |     A3   B3   C3            | | X3 |   | D3 |
! |      .    .    .            | | :  | = | :  |
! |           .    .    .       | | :  |   | :  |
! |                .    .   CN-1| | :  |   | :  |
! |0                   AN     BN| | XN |   | DN |
! --                           -- --  --   --  --
!-
! Reference
! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
! http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_%28Thomas_algorithm%29
!+
! History
! Oct. 17, 2016 Develop the code (Chia-Pang Kuo)
!-
! Contact Information
!  Chia-Pang Kuo @ Texas A&M University
!  email: intelb52158@tamu.edu

! ----- Input -----
      integer, intent(in) :: N ! N*N of the matrix
      real(kind=rb), intent(in) :: A(N)
      real(kind=rb), intent(in) :: B(N)
      real(kind=rb), intent(in) :: C(N)
      real(kind=rb), intent(in) :: D(N)
! ----- Output -----
      real(kind=rb), intent(out) :: X(N)
! ----- Local -----
      integer :: i 
      real(kind=rb) :: P(0:N)
      real(kind=rb) :: Q(0:N)
      real(kind=rb) :: denominator

      X = 0.0_rb
      P = 0.0_rb
      Q = 0.0_rb
      
      ! forward elimination
      do i = 1, N
         denominator = B(i) + A(i) * P(i-1)
         P(i) = -C(i) / denominator
         Q(i) = (D(i)-A(i)*Q(i-1)) / denominator
      enddo

      ! back substitution
      do i = N, 1, -1
         X(i) = P(i) * X(i+1) + Q(i)
      enddo
      end subroutine TDMA

      end module rrtmg_lw_rtr2function
