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
    real(dp)                                            :: temp_unif_noise
    real(dp)                                            :: qu_unif_noise
    integer(i4b)                                        :: i, j, m, l, unit
    character(len=512)                                  :: filename, paramfile
    character(len=512)                                  :: cloutfilename
    character(len=512)                                  :: outfilename
    character(len=512)                                  :: rmsoutfilename
    type(planck_rng)                                    :: rng_handle   
    real(dp),   pointer,        dimension(:, :)         :: pixwin
    real(dp), allocatable, dimension(:, :)              :: sqrt_covmatrix
    logical(lgt)                                        :: inv

    character(len=2)                                    :: simulstring
    integer(i4b)                                        :: k
    integer(i4b)                                        :: ipix, seed
    character(len=128)                                  :: operation

    unit = getlun()

    call getarg(1, paramfile)

    call get_parameter(unit, paramfile, 'NSIDE', par_int=nside)
    call get_parameter(unit, paramfile, 'LMAX', par_int=lmax)
    call get_parameter(unit, paramfile, 'SEED', par_int=seed)
    call get_parameter(unit, paramfile, 'OPERATION', par_string=operation)
    call get_parameter(unit, paramfile, 'OUTFILENAME', par_string=outfilename)

    call rand_init(rng_handle, seed)
    npix = 12*nside**2
    allocate(cls(2:lmax, 3, 3))
    allocate(beam(0:lmax, 3))
    call get_parameter(unit, paramfile, 'BEAMFILE', par_string=filename)
    call read_beam(filename, beam)
    call read_pixwin(nside, 3, pixwin)
    beam = beam * pixwin(0:lmax, 1:3)
    deallocate(pixwin)

    select case (trim(operation))
      case('pol_cov')
         allocate(map_pol(0:npix-1, 3))
         allocate(noise(0:3*npix-1))
         nmaps = 3
      case('pol_diag')
         allocate(map_pol(0:npix-1, 3))
         allocate(noise_pol(0:npix-1, 3))
         allocate(rms(0:npix-1, 3))
         nmaps = 3
      case('temp_uniform_diag')
         allocate(map(0:npix-1))
         allocate(noise(0:npix-1))
         allocate(rms(0:npix-1, 1))
         nmaps = 1
      case('temp_rmsmap_diag')
         allocate(map(0:npix-1))
         allocate(noise(0:npix-1))
         allocate(rms(0:npix-1, 1))
      case default
         print *, 'unknown operation. Quitting'
         stop
   end select
   allocate(x(nmaps, 2))
   allocate(alms(nmaps, 0:lmax, 0:lmax))
   allocate(chol(nmaps, nmaps, 2:lmax))
   if (trim(operation) == 'pol_cov') then
      call get_parameter(unit, paramfile, 'COVMATFILE', par_string=filename)
      call read_covmatrix(unit, filename, ordering, polarization, sqrt_covmatrix, inv, n)
   else if (trim(operation) == 'temp_rmsmap_diag') then
      call get_parameter(unit, paramfile, 'RMSFILE', par_string=filename)
      call read_map(rms, ordering, trim(filename))
   else if (trim(operation) == 'temp_uniform_diag' .or. trim(operation) == 'pol_diag') then
      call get_parameter(unit, paramfile, 'RMS_OUTFILE', par_string=rmsoutfilename)
      call get_parameter(unit, paramfile, 'TEMP_NOISE', par_dp=temp_unif_noise)
      rms(:, 1) = temp_unif_noise
      if (trim(operation) == 'pol_diag') then
         call get_parameter(unit, paramfile, 'QU_NOISE', par_dp=qu_unif_noise)
         rms(:, 2:3) = qu_unif_noise
      end if 
   end if

   call get_parameter(unit, paramfile, 'CLFILE', par_string=filename)
   call get_parameter(unit, paramfile, 'CL_OUTFILE', par_string=cloutfilename)
   open(unit, file=trim(filename))
   cls = 0.d0
   do l = 2, lmax
      !There is an assumed ordering in the cls, no matter whether we want
      !polarization or not. Can be fixed later.
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
   else
      call alm2map(nside, lmax, lmax, alms, map)
   end if

   if (trim(operation) == 'pol_cov') then
      if (ordering == nest) then
         allocate(tempmap(0:npix - 1, 3))
         tempmap = 0
         do i = 0, npix - 1
            call ring2nest(nside, i, ipix)
            tempmap(ipix, :) = map_pol(i, :)
         end do
         map_pol = tempmap
      end if
      !ADD NOISE (COVMATRIX)
      noise = 0
      allocate(eta(0:3*npix-1))
      do j = 0, 3*npix - 1
         eta(j) = rand_gauss(rng_handle)
      end do
      noise = matmul(sqrt_covmatrix, eta)
      do i = 1, 3
         do j = 0, npix - 1
               map_pol(j, i) = map_pol(j, i) + noise((i - 1) * npix + j)
         end do
      end do
      deallocate(eta)

      if (ordering == nest) then
         !Convert back to ring
         tempmap = 0
         do i = 0, npix - 1
            call nest2ring(nside, i, ipix)
            tempmap(ipix, :) = map_pol(i, :)
         end do
         map_pol = tempmap
         deallocate(tempmap)
      end if
   else
      do i = 1, nmaps
         do j = 0, npix-1
            if (trim(operation) == 'temp_rmsmap_diag' .and. ordering == nest) then
               call nest2ring(nside, j, ipix)
               noise(ipix) = rms(j, i)*rand_gauss(rng_handle)
            else if (trim(operation) == 'pol_diag') then
               noise_pol(j, i) = rms(j, i)*rand_gauss(rng_handle)
            else 
               noise(j) = rms(j, i)*rand_gauss(rng_handle)
            end if
        end do
      end do

      if (trim(operation) == 'pol_diag') then
         map_pol = map_pol + noise_pol
      else
         map = map + noise
      end if
   end if

   if (trim(operation) == 'pol_cov' .or. trim(operation) == 'pol_diag') then
      call write_map(map_pol, ring, outfilename)
   else
      call write_map(map, ring, outfilename)
   end if
   if (trim(operation) == 'temp_uniform_diag' .or. trim(operation) == 'pol_diag') then
      call write_map(rms, ring, rmsoutfilename)
   end if
   
   do l = 2, lmax
      do j = 1, nmaps
         alms(j, l, :) = alms(j, l, :) / beam(l, j)
      end do
   end do

   open(unit, file=cloutfilename, recl=512)
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

      if (trim(operation) == 'temp_uniform_diag' .or. trim(operation) == 'temp_rmsmap_diag') then
         write(unit, *) l, Cl_TT
      else if (trim(operation) == 'pol_diag' .or. trim(operation) == 'pol_cov') then
         write(unit, *) l, Cl_TT, Cl_EE, Cl_BB, Cl_TE
      end if
   end do
   close(unit)

end program simmap
