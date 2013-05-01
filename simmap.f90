!Program to simulate a very simple cmb map, with constant rms for both
!temperature and polarization.

program simmap
    use healpix_types
    use alm_tools
    use rngmod
    use quiet_fileutils
    use quiet_utils
    use math_tools
    implicit none

    real(dp),   allocatable,    dimension(:, :, :)      :: cls
    real(dp),   allocatable,    dimension(:, :)         :: beam
    real(dp), allocatable,      dimension(:)            :: eta
    real(dp),   allocatable,    dimension(:, :, :)      :: chol
    complex(dpc),   allocatable,    dimension(:, :, :)      :: alms
    integer(i4b)                                        :: lmax
    integer(i4b)                                        :: nside, npix
    integer(i4b)                                        :: nmaps
    real(dp),   allocatable,    dimension(:, :)         :: x
    real(dp),   allocatable,    dimension(:, :)         :: rms

    !POLARIZATION
    real(dp),   allocatable,    dimension(:, :)         :: noise_pol
    real(dp),   allocatable,    dimension(:, :)         :: map_pol
    real(dp),   allocatable,    dimension(:, :, :)         :: map_pol_mult
   !COVARIANCE MATRIX
    real(dp),   allocatable,    dimension(:, :)         :: tempmap
    integer(i4b)                                        :: ordering
    integer(i4b)                                        :: polarization
    integer(i8b)                                        :: n
    !DIAG

!   TEMPERATURE ONLY(DIAG)
    real(dp),   allocatable,    dimension(:)         :: noise
    real(dp),   allocatable,    dimension(:)         :: map


    real(dp)                                            :: Cl_TT, Cl_TE, Cl_EE
    real(dp)                                            :: Cl_BB, chisq
    integer(i4b)                                        :: i, j, m, l, unit
    character(len=512)                                  :: filename
    type(planck_rng)                                    :: rng_handle   
    real(dp),   pointer,        dimension(:, :)         :: pixwin
    real(dp), allocatable, dimension(:, :)              :: sqrt_covmatrix
    logical(lgt)                                        :: inv

    character(len=2)                                    :: simulstring
    integer(i4b)                                        :: k
    integer(i4b)                                        :: ipix, seed
    character(len=128)                                  :: operation


    nside = 64
    lmax  = 128
    npix = 12*nside**2
    seed = 213409
    unit = 39
    allocate(cls(2:lmax, 3, 3))
    allocate(beam(0:lmax, 3))

    operation = 'temp_diag'

    select case (trim(operation))
      case('pol_cov')
         allocate(map_pol(0:npix-1, 3))
!         allocate(map_pol_mult(0:npix-1, 3, 10))
         allocate(noise(0:3*npix-1))
         nmaps = 3
      case('pol_diag')
         allocate(map_pol(0:npix-1, 3))
         allocate(noise_pol(0:npix-1, 3))
         allocate(rms(0:npix-1, 3))
         nmaps = 3
      case('temp_diag')
         allocate(map(0:npix-1))
!         allocate(rms(0:npix-1, 1))
         allocate(noise(0:npix-1))
         nmaps = 1
      case default
         print *, 'unknown operation. Quitting'
         stop
   end select
   allocate(x(nmaps, 2))
   allocate(alms(nmaps, 0:lmax, 0:lmax))
   allocate(chol(nmaps, nmaps, 2:lmax))
   if (trim(operation) == 'pol_cov') then
!      filename = '/mn/svati/d1/eirikgje/data/vault/ctp3_simulation_data/covmats/ctp3_madam_covmat_nside16_143GHz_inv_bin.dat.inverse_440_muK_CMB_regularised_sqrt.unf'
!      filename = '/mn/svati/d1/eirikgje/data/vault/dx7_data/covmats/covmat_nside16_70GHz_DX7_sqrt.unf'
      filename = '/mn/svati/d1/eirikgje/data/vault/ffp4_data/covmats/covmat_nside16_70GHz_ffp4_sqrt.unf'
      call read_covmatrix(unit, filename, ordering, polarization, sqrt_covmatrix, inv, n)
   end if

    call rand_init(rng_handle, seed)
!    filename = '/mn/svati/d1/eirikgje/data/vault/gaussian_beams/beam_14arcmin_lmax2000.fits'
!    filename = '/mn/svati/d1/eirikgje/data/vault/gaussian_beams/beam_440arcmin_lmax47.fits'
    filename = '/mn/svati/d1/eirikgje/data/vault/gaussian_beams/beam_900arcmin_1percenttruncated_lmax47.fits'
!    filename = '/mn/svati/d1/eirikgje/data/vault/cobe_dmr_data/beams/cobe_dmr_beam.fits'

    call read_beam(filename, beam)
    call read_pixwin(nside, 3, pixwin)
    beam = beam * pixwin(0:lmax, 1:3)
    deallocate(pixwin)

!      filename = '/mn/svati/d1/eirikgje/data/vault/tau_likelihood_data/tau_powerspec/taucls_tau' // taustring // '_scalCls.dat'
   filename = '/mn/svati/d1/eirikgje/data/vault/general/WMAP_bestfit_scalCls_lmax2000.dat'
   open(unit, file=trim(filename))
   cls = 0.d0
   do l = 2, lmax
!        read(unit, *) i, cls(l, 1, 1), cls(l, 2, 2), cls(l, 3, 3), cls(l, 1, 2)
      read(unit, *) i, cls(l, 1, 1), cls(l, 2, 2), cls(l, 1, 2)
      cls(l, 2, 1) = cls(l, 1, 2)
      cls(l, :, :) = cls(l, :, :)*2.d0*pi/(real(l, dp)*(real(l, dp) + 1.d0))
   end do
   close(unit)

   chol = 0.d0
   chol(1, 1, :) = sqrt(cls(:, 1, 1))
   if (trim(operation) == 'pol_cov' .or. trim(operation) == 'pol_diag') then
      chol(2, 1, :) = cls(:, 2, 1)/chol(1, 1, :)
      chol(2, 2, :) = sqrt(cls(:, 2, 2)-chol(2, 1, :)*chol(2, 1, :))
      chol(3, 1, :) = cls(:, 3, 1)/chol(1, 1, :)
      chol(3, 2, :) = (cls(:, 3, 2) -chol(3, 1, :)*chol(2, 1, :))/chol(2, 2, :)
      chol(3, 3, :) = sqrt(cls(:, 3, 3)-chol(3, 1, :)*chol(3, 1, :)-chol(3, 2, :)*chol(3, 2, :))
   end if

!do k = 1, 10
   allocate(eta(nmaps))
    alms = 0.d0
    x = 0.d0
    eta = 0.d0
   do l = 2, lmax
        do m = 0, l
            do i = 1, 2
                do j = 1, nmaps
                    eta(j) = rand_gauss(rng_handle)
                end do
                x(:, i) = matmul(chol(:, :, l), eta)
            end do
            if (m == 0) then
                alms(:, l, m) = cmplx(x(:, 1), 0.d0)
            else
                alms(:, l, m) = 1.d0/sqrt(2.d0)*cmplx(x(:, 1), x(:, 2))
            end if
        end do
        do j = 1, nmaps
            alms(j, l, :) = alms(j, l, :)*beam(l, j)
        end do
   end do

   deallocate(eta)
   if (trim(operation) == 'pol_cov' .or. trim(operation) == 'pol_diag') then
      call alm2map(nside, lmax, lmax, alms, map_pol)
      !call alm2map(nside, lmax, lmax, alms, map_pol_mult(:, :, 1))
      !do k = 2, 10
      !   map_pol_mult(:, :, k) = map_pol_mult(:, :, 1)
      !end do
   else if (trim(operation) == 'temp_diag') then
      call alm2map(nside, lmax, lmax, alms, map)
   end if
do k = 1, 10
!   map = 0
   deallocate(map)
   call read_map(map, ordering, '/mn/svati/d1/eirikgje/data/vault/dm_data/nonoise/map_DM_64_20GeV_3e-26_33GHz_nonoise_muK.fits')

   if (trim(operation) == 'pol_cov') then
    !Convert from ring to nest
    allocate(tempmap(0:npix - 1, 3))
!    do k = 1, 10
      tempmap = 0
      do i = 0, npix - 1
         call ring2nest(nside, i, ipix)
!         tempmap(ipix, :) = map_pol_mult(i, :, k)
         tempmap(ipix, :) = map_pol(i, :)
      end do
!      map_pol_mult(:, :, k) = tempmap
      map_pol = tempmap
   
      !ADD NOISE (COVMATRIX)
      noise = 0
      allocate(eta(0:3*npix-1))
      do j = 0, 3*npix - 1
         eta(j) = rand_gauss(rng_handle)
      end do
      noise = matmul(sqrt_covmatrix, eta)
      do i = 1, 3
         do j = 0, npix - 1
             !map_pol(j, i) = map_pol(j, i) + noise((i - 1) * npix + j)
             map_pol(j, i) = map_pol(j, i) + 2 * noise((i - 1) * npix + j)
!             map_pol_mult(j, i, k) = map_pol_mult(j, i, k) + 2 * noise((i - 1) * npix + j)
         end do
      end do
      deallocate(eta)

      !Convert back to ring
      tempmap = 0
      do i = 0, npix - 1
         call nest2ring(nside, i, ipix)
         tempmap(ipix, :) = map_pol(i, :)
!         tempmap(ipix, :) = map_pol_mult(i, :, k)
      end do
!      map_pol_mult(:, :, k) = tempmap
      map_pol = tempmap
!   end do
   deallocate(tempmap)
   else if (trim(operation) == 'pol_diag' .or. trim(operation) == 'temp_diag') then

      !ADD NOISE (DIAG)
      if (k == 1) then
         call read_map(rms, ordering, '/mn/svati/d1/eirikgje/data/vault/wmap_data/rms/wmap_9yr_rms_Ka_I_muK_ns64.fits')
      end if
      !rms = rms * 0.005
!      rms(:, 1) = 0.30
!      rms(:, 1) = 65.0
!      rms(:, 2:3) = sqrt(0.35d0)

      do i = 1, nmaps
         do j = 0, npix-1
            if (trim(operation) == 'temp_diag') then
               call nest2ring(nside, j, ipix)
!               noise(j) = rms(j, i)*rand_gauss(rng_handle)
               noise(ipix) = rms(j, i)*rand_gauss(rng_handle)
            else if (trim(operation) == 'pol_diag') then
               noise_pol(j, i) = rms(j, i)*rand_gauss(rng_handle)
            end if
        end do
      end do

      if (trim(operation) == 'temp_diag') then
         map = map + noise
      else if (trim(operation) == 'pol_diag') then
         map_pol = map_pol + noise_pol
      end if
   end if

!   do k = 1, 10
   call int2string(k, simulstring)
!   filename = 'ffp4_ownsimul_4timescovmat_70GHz_nside16_440arcmin_WMAP_bestfit_lmax47_cmbsimul' // simulstring // '.fits'
   !filename = 'ffp4_ownsimul_4timescovmat_70GHz_nside16_600arcmin_1percenttruncated_WMAP_bestfit_lmax47_cmbsimul' // simulstring // '.fits'
!   filename = 'ffp4_ownsimul_4timescovmat_70GHz_nside16_900arcmin_1percenttruncated_WMAP_bestfit_lmax47.fits'
!   filename = 'WMAP_bestfit_simul_n32_cobe_beam_diag_0.3muK_I_lmax95.fits'
!   filename = 'dx7_ownsimul_4timesnonsmoothcovmat_70GHz_nside16_440_arcmin_WMAP_bestfit_lmax47_2.fits'
   filename = 'map_DM_64_20GeV_3e-26_33GHz_wmap_Ka_noise_simul' // simulstring // '.fits'
   if (trim(operation) == 'pol_cov' .or. trim(operation) == 'pol_diag') then
!      call write_map(map_pol_mult(:, :, k), ring, filename)
      call write_map(map_pol, ring, filename)
!      map_pol = 1.d0
!      filename = 'mask_fullsky_n16.fits'
!      call write_map(map_pol, ring, filename)
   else if (trim(operation) == 'temp_diag') then
      call write_map(map, ring, filename)
      map = 1.d0
      filename = 'mask_fullsky_n64.fits'
      call write_map(map, ring, filename)
   end if
   if (trim(operation) == 'temp_diag' .or. trim(operation) == 'pol_diag') then
      filename = 'rms_0.3muK_I_nside32.fits'
      call write_map(rms, ring, filename)
   end if
end do
!   currtau = currtau + deltatau


   
   do l = 2, lmax
      do j = 1, nmaps
         alms(j, l, :) = alms(j, l, :) / beam(l, j)
      end do
   end do

!   open(unit, file='clout_dx7simul_tau' // taustring // '.dat', recl=512)
   open(unit, file='clout_WMAP_bestfit_lmax47.dat', recl=512)
   do l = 2, lmax
      Cl_TT = 0
      Cl_TE = 0
      Cl_EE = 0
      Cl_BB = 0

      do m = 0, l
         if (m == 0) then
            Cl_TT = Cl_TT + alms(1, l, m) ** 2
            if (trim(operation) == 'pol_diag' .or. trim(operation) == 'pol_cov')then
               Cl_TE = Cl_TE + alms(1, l, m) * conjg(alms(2, l, m))
               Cl_EE = Cl_EE + alms(2, l, m) **2
               Cl_BB = Cl_BB + alms(3, l, m) **2
            end if
         else
            Cl_TT = Cl_TT + 2.d0 * real(alms(1, l, m) * conjg(alms(1, l, m)))
            if (trim(operation) == 'pol_diag' .or. trim(operation) == 'pol_cov')then
               Cl_TE = Cl_TE + real(alms(1, l, m) * conjg(alms(2, l, m))) &
                  & + real(alms(2, l, m) * conjg(alms(1, l, m)))
               Cl_EE = Cl_EE + 2.d0 * real(alms(2, l, m) * conjg(alms(2, l, m)))
               Cl_BB = Cl_BB + 2.d0 * real(alms(3, l, m) * conjg(alms(3, l, m)))
            end if
         end if
      end do
      Cl_TT = Cl_TT * float(l * (l + 1)) / (2.d0 * pi * float(2 * l + 1))
      Cl_TE = Cl_TE * float(l * (l + 1)) / (2.d0 * pi * float(2 * l + 1))
      Cl_EE = Cl_EE * float(l * (l + 1)) / (2.d0 * pi * float(2 * l + 1))
      Cl_BB = Cl_BB * float(l * (l + 1)) / (2.d0 * pi * float(2 * l + 1))

      if (trim(operation) == 'temp_diag') then
         write(unit, *) l, Cl_TT
      else if (trim(operation) == 'pol_diag' .or. trim(operation) == 'pol_cov') then
         write(unit, *) l, Cl_TT, Cl_EE, Cl_BB, Cl_TE
      end if
   end do
   close(unit)
!end do

end program simmap
