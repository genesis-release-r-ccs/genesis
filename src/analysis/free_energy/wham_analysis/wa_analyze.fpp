!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   wa_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module wa_analyze_mod

  use wa_option_str_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type s_data_k
    real(wp),      allocatable :: v(:,:)       ! (nstep, ndim)
  end type s_data_k

  type s_pmf
    real(wp),      allocatable :: v(:)         ! (nbin)
  end type s_pmf

  ! constants
  real(wp),        parameter   :: KB   = 0.00198719168260038_wp

  ! subroutines
  public  :: analyze

  private :: build_data_k_cv
  private :: build_data_k_dcd
  private :: build_bias_km
  private :: solve_wham
  private :: solve_wham_block
  private :: output_wham

  private :: compute_grid_center
  private :: exec_fhandle_1d
  private :: exec_fhandle_2d
  private :: get_dcd_cv
  private :: get_com_dist
  private :: get_com_angl
  private :: get_com_dihe
  private :: check_cvfile
  private :: check_dcdfile
  private :: get_replicate_name1
  private :: get_replicate_name2
  private :: periodic
  private :: logsumexp
  private :: std

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    input    : input information
  !! @param[in]    output   : output information
  !! @param[in]    option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, input, output, option)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option

    ! local variables
    type(s_data_k), allocatable :: data_k(:)     ! (nbrella)
    real(wp),       allocatable :: bias_km(:,:)  ! (nbin,nbrella)
    integer,        allocatable :: h_km(:,:)     ! (nbin,nbrella)
    type(s_pmf),    allocatable :: pmf_m(:)      ! (nblocks)


    ! check check only
    !
    if (option%check_only) &
      return


    ! build data_k
    !
    if (input%cvfile /= '') then

      call build_data_k_cv (input%cvfile, option, data_k)

    else if (input%dcdfile /= '') then

      call build_data_k_dcd(input%dcdfile, molecule, option, data_k)

    end if


    ! build bias_km
    !
    call build_bias_km(option, data_k, bias_km)


    ! solve WHAM
    !
    call solve_wham(option, data_k, bias_km, pmf_m)


    ! output f_k and pmf
    !
    call output_wham(option, output, pmf_m)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_data_k_cv(cvfile, option, data_k)

    ! formal arguments
    character(*),            intent(in)    :: cvfile
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)

    ! local variables
    integer                  :: file, tim, i, j, k, istep
    integer                  :: nbrella, ndim, nstep, nfunc
    character(MaxFilename)   :: filename


    write(MsgOut,'(a)') 'Build_Data_K_Cv> '

    ndim    = option%dimension
    nfunc   = size(option%rest_func_no)
    nbrella = 1
    do i = 1, ndim
      nbrella = nbrella * option%rest_nreplica(option%rest_func_no(i))
    end do

    ! allocate data_k
    !
    allocate(data_k(nbrella))


    ! check nstep
    !
    call check_cvfile(cvfile, nstep)


    ! read cv file and setup data_k
    !
    nbrella = 0

    if (ndim == 1) then

      do i = 1, option%rest_nreplica(option%rest_func_no(1))

        nbrella = nbrella + 1
        allocate(data_k(nbrella)%v(nstep,nfunc))

        filename = get_replicate_name1(cvfile, nbrella)
        write(MsgOut,'(a,a)') '  read cv file: ',trim(filename)

        call open_file(file, &
                       filename, &
                       IOFileInput)

        do istep = 1, nstep
          read(file,*) tim, (data_k(nbrella)%v(istep,k),k=1,nfunc)
        end do

        call close_file(file)

      end do

    else

      if (nfunc == 1) &
        call error_msg('Build_Data_K_Cv> # of rest_func must be 2 on dim=2')

      do i = 1, option%rest_nreplica(option%rest_func_no(1))
        do j = 1, option%rest_nreplica(option%rest_func_no(2))

          nbrella = nbrella + 1
          allocate(data_k(nbrella)%v(nstep,nfunc))

          filename = get_replicate_name1(cvfile, nbrella)
          write(MsgOut,'(a,a)') '  read cv file: ',trim(filename)

          call open_file(file, &
                         filename, &
                         IOFileInput)

          do istep = 1, nstep
            read(file,*) tim, (data_k(nbrella)%v(istep,k),k=1,nfunc)
          end do

          call close_file(file)

        end do
      end do

    end if

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

  end subroutine build_data_k_cv

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_data_k_dcd(dcdfile, molecule, option, data_k)

    ! formal arguments
    character(*),            intent(in)    :: dcdfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)

    ! local variables
    type(s_trajectory)       :: trajectory
    type(s_trj_file)         :: file
    integer                  :: i, j, k, istep
    integer                  :: nbrella, ndim, nstep, natom, nfunc
    character(MaxFilename)   :: filename


    write(MsgOut,'(a)') 'Build_Data_K_Dcd> '

    ndim    = option%dimension
    nfunc   = size(option%rest_func_no)
    nbrella = 1
    do i = 1, ndim
      nbrella = nbrella * option%rest_nreplica(option%rest_func_no(i))
    end do

    ! allocate data_k
    !
    allocate(data_k(nbrella))


    ! check nstep and allocate trajectory
    !
    call check_dcdfile(dcdfile, nstep, natom)

    if (natom /= option%num_atoms) &
      call error_msg( &
      'Build_Data_K_Dcd> Dcd atom count is different from PSF/PRMTOP.')

    call alloc_trajectory(trajectory, natom)


    ! read trajectory and setup data_k
    !
    nbrella = 0

    if (ndim == 1) then

      do i = 1, option%rest_nreplica(option%rest_func_no(1))

        nbrella = nbrella + 1
        allocate(data_k(nbrella)%v(nstep,nfunc))

        filename = get_replicate_name1(dcdfile, nbrella)
        write(MsgOut,'(a,a)') '  read and analyze trajectory: ',trim(filename)

        call open_trj(file, &
                    filename,&
                    TrjFormatDCD,   &
                    TrjTypeCoorBox, &
                    IOFileInput)

        do istep = 1, nstep

          call read_trj(file, trajectory)

          do k = 1, nfunc
            data_k(nbrella)%v(istep,k) = &
                 get_dcd_cv(molecule, option, trajectory, k)

          end do

        end do

        call close_trj(file)

      end do

    else

      if (nfunc == 1) &
        call error_msg('Build_Data_K_Cv> # of rest_func must be 2 on dim=2')

      do i = 1, option%rest_nreplica(option%rest_func_no(1))
        do j = 1, option%rest_nreplica(option%rest_func_no(2))

          nbrella = nbrella + 1
          allocate(data_k(nbrella)%v(nstep,nfunc))

          filename = get_replicate_name1(dcdfile, nbrella)
          write(MsgOut,'(a,a)') '  read and analyze trajectory: ',trim(filename)

          call open_trj(file, &
                      filename,&
                      TrjFormatDCD,   &
                      TrjTypeCoorBox, &
                      IOFileInput)

          do istep = 1, nstep

            call read_trj(file, trajectory)

            do k = 1, nfunc
              data_k(nbrella)%v(istep,k) = &
                   get_dcd_cv(molecule, option, trajectory, k)
            end do

          end do

          call close_trj(file)

        end do
      end do

    end if

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine build_data_k_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_bias_km(option, data_k, bias_km)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)
    real(wp),                allocatable   :: bias_km(:,:)

    ! local variables
    integer                  :: nbrella, ndim, nbin, nbin_x, nbin_y, nfunc
    integer                  :: i, j, mx, my, m

    real(wp),    allocatable :: umbrella_center(:,:), k_constant(:,:)
    real(wp),    allocatable :: grid_x(:), center_x(:), grid_y(:), center_y(:)


    write(MsgOut,'(a)') 'Build_Bias_Km> '

    ndim    = option%dimension
    !nfunc   = size(option%rest_func_no)
    nbrella = size(data_k)

    if (ndim == 1) then

      nbin = option%num_grids(1) - 1

      allocate(bias_km(nbin, nbrella))
      allocate(umbrella_center(nbrella, 1), &
               k_constant     (nbrella, 1))

      umbrella_center(:,1) = option%rest_references(:,option%rest_func_no(1))
      k_constant(:,1)      = option%rest_constants (:,option%rest_func_no(1))

      call compute_grid_center(option, 1, grid_x, center_x)
      
      write(MsgOut, '(a,i3)') 'Build_Bias_Km> centers of grids in dimension ',1
      do j = 1, size(grid_x) - 1
        write(MsgOut,'(F10.4,$)') center_x(j)
      end do
      write(MsgOut, *)
      write(MsgOut, *)

      do i = 1, nbrella
        do j = 1, nbin
          call exec_fhandle_1d(option, umbrella_center(i,:), k_constant(i,:), &
                               center_x(j), bias_km(j, i))
        end do
      end do

      deallocate(grid_x, center_x)

    else

      nbin_x = (option%num_grids(1) - 1)
      nbin_y = (option%num_grids(2) - 1)
      nbin = nbin_x * nbin_y

      allocate(bias_km(nbin, nbrella))
      allocate(umbrella_center(nbrella, 2), &
               k_constant     (nbrella, 2))

      nbrella = 0
      do i = 1, option%rest_nreplica(option%rest_func_no(1))
        do j = 1, option%rest_nreplica(option%rest_func_no(2))
          nbrella = nbrella + 1
          umbrella_center(nbrella,1:2) = &
               (/option%rest_references(i, option%rest_func_no(1)), &
                 option%rest_references(j, option%rest_func_no(2))/)
          k_constant(nbrella,1:2) = &
               (/option%rest_constants(i, option%rest_func_no(1)), &
                 option%rest_constants(j, option%rest_func_no(2))/)
        end do
      end do

      call compute_grid_center(option, 1, grid_x, center_x)
      call compute_grid_center(option, 2, grid_y, center_y)

      write(MsgOut, '(a,i3)') 'Build_Bias_Km> centers of grids in dimension ',1
      do j = 1, size(grid_x) - 1
        write(MsgOut,'(F10.4,$)') center_x(j)
      end do
      write(MsgOut, *)
      write(MsgOut, *)
      write(MsgOut, '(a,i3)') 'Build_Bias_Km> centers of grids in dimension ',2
      do j = 1, size(grid_y) - 1
        write(MsgOut,'(F10.4,$)') center_y(j)
      end do
      write(MsgOut, *)
      write(MsgOut, *)

      do i = 1, nbrella
        do mx = 1, nbin_x
          do my = 1, nbin_y
            m = (my-1)*nbin_x + mx
            call exec_fhandle_2d(option, umbrella_center(i,:), k_constant(i,:),&
                                 center_x(mx), center_y(my), bias_km(m, i))
          end do
        end do
      end do

      deallocate(grid_y, center_y)
      deallocate(grid_x, center_x)

    end if

    bias_km(1:nbin,1:nbrella) = bias_km(1:nbin, 1:nbrella) / &
                                (KB * option%temperature)

    deallocate(umbrella_center, k_constant)


    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine build_bias_km

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine solve_wham(option, data_k, bias_km, pmf_m)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)
    real(wp),                intent(in)    :: bias_km(:,:)
    type(s_pmf),             allocatable   :: pmf_m(:)

    ! local variables
    real(wp)                 :: rstep, pmf_min
    integer                  :: nblock, iblock, nstep, istep0, istep1
    integer                  :: nbin, ibin, pmf_min_idx

    type(s_pmf), allocatable :: pmf0(:)
    real(wp),    allocatable :: v(:)


    write(MsgOut,'(a)') 'Solve_WHAM> '
    write(MsgOut,'(a)') ''

    nblock = option%nblocks
    nstep  = size(data_k(1)%v(:,1))

    allocate(pmf0(nblock))

    istep0 = 1
    rstep  = 0.0_wp

    do iblock = 1, nblock

      rstep = rstep + (real(nstep, wp) / real(nblock, wp))

      istep1 = nint(rstep)

      call solve_wham_block(option, data_k, bias_km, &
                            iblock, istep0-1, istep1-istep0+1, &
                            pmf0(iblock)%v)

      istep0 = istep1 + 1

    end do

    nbin = size(pmf0(1)%v)

    if (nblock > 1) then

      allocate(pmf_m(2))
      allocate(pmf_m(1)%v(nbin))
      allocate(pmf_m(2)%v(nbin))
      allocate(v(nblock))

      pmf_min = 1000000.0_wp
      do ibin = 1, nbin
        if (pmf_min > pmf0(1)%v(ibin)) then
          pmf_min_idx = ibin
          pmf_min = pmf0(1)%v(ibin)
        end if
      end do
      do iblock = 1, nblock
        v(iblock) = pmf0(iblock)%v(pmf_min_idx)
      end do
      do ibin=1,nbin
        do iblock=1,nblock
          pmf0(iblock)%v(ibin) = pmf0(iblock)%v(ibin) - v(iblock)
        end do
      end do

      do ibin = 1, nbin
        do iblock = 1, nblock
          v(iblock) = pmf0(iblock)%v(ibin)
        end do
        pmf_m(1)%v(ibin) = sum(v(1:nblock))/real(nblock,wp)
        pmf_m(2)%v(ibin) = 2.0_wp * std(v)
      end do

      deallocate(v)

    else

      allocate(pmf_m(1))
      allocate(pmf_m(1)%v(nbin))

      pmf_min = minval(pmf0(1)%v)
      pmf_m(1)%v(1:nbin) = pmf0(1)%v(1:nbin) - pmf_min

    end if

    do iblock = 1, nblock
      deallocate(pmf0(iblock)%v)
    end do
    deallocate(pmf0)

    write(MsgOut,'(a)') ''

    return

  end subroutine solve_wham

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine solve_wham_block(option, data_k, bias_km, &
                              iblock, step0, nstep, pmf_m)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)
    real(wp),                intent(in)    :: bias_km(:,:)
    integer,                 intent(in)    :: iblock
    integer,                 intent(in)    :: step0
    integer,                 intent(in)    :: nstep
    real(wp),                allocatable   :: pmf_m(:)

    ! local variables
    real(wp)                 :: check_convergence, fk, k, k1, k2
    integer                  :: iter
    integer                  :: nbrella, nbin, nbin_x, nbin_y
    integer                  :: ibrella, ibin, ibin_x, ibin_y, istep

    real(wp),    allocatable :: g_km(:,:), heff_m(:), log_heff_m(:)
    real(wp),    allocatable :: Neff_km(:,:), log_Neff_km(:,:)
    real(wp),    allocatable :: log_den_km(:,:), log_den_m(:), log_num_m(:)
    real(wp),    allocatable :: f_k(:), f_k_new(:), log_prob_m(:), tmp_km(:,:)
    integer,     allocatable :: h_km(:,:), N_k(:)
    real(wp),    allocatable :: grid_x(:), center_x(:), grid_y(:), center_y(:)


    ! h_km
    !
    nbrella = size(data_k)

    if (option%dimension == 1) then

      nbin = option%num_grids(1) - 1

      allocate(h_km(nbin, nbrella))

      call compute_grid_center(option, 1, grid_x, center_x)

      h_km(1:nbin,1:nbrella) = 0

      do ibrella = 1, nbrella
        do istep = step0+1, step0+nstep
          k = data_k(ibrella)%v(istep,1)
          do ibin = 1, nbin
            if (grid_x(ibin) <= k .and. k < grid_x(ibin+1)) then
              h_km(ibin,ibrella) = h_km(ibin,ibrella) + 1
              exit
            end if
          end do
        end do
      end do

      deallocate(grid_x, center_x)

    else

      nbin_x = (option%num_grids(1) - 1)
      nbin_y = (option%num_grids(2) - 1)
      nbin = nbin_x * nbin_y

      allocate(h_km(nbin, nbrella))

      call compute_grid_center(option, 1, grid_x, center_x)
      call compute_grid_center(option, 2, grid_y, center_y)

      do ibrella = 1, nbrella
        do istep = step0+1, step0+nstep
          k1 = data_k(ibrella)%v(istep,1)
          k2 = data_k(ibrella)%v(istep,2)
          do ibin_x = 1, nbin_x
            do ibin_y = 1, nbin_y
              if (grid_x(ibin_x) <= k1 .and. k1 < grid_x(ibin_x+1) .and. &
                  grid_y(ibin_y) <= k2 .and. k2 < grid_y(ibin_y+1)) then
                ibin = (ibin_x-1)+(ibin_y-1)*nbin_x+1
                h_km(ibin,ibrella) = h_km(ibin,ibrella) + 1
              end if
            end do
          end do
        end do
      end do

      deallocate(grid_y, center_y)
      deallocate(grid_x, center_x)

    end if


    ! N_k
    !
    allocate(N_k(nbrella))

    do ibrella = 1, nbrella
      N_k(ibrella) = sum(h_km(1:nbin,ibrella))
    end do


    ! calculate statistical inefficiency (g_km)
    !
    allocate(g_km(nbin,nbrella))

    g_km(1:nbin,1:nbrella) = 1.0_wp


    ! calculate effective histogram (heff_m)
    !
    allocate(heff_m(nbin), log_heff_m(nbin))

    do ibin = 1, nbin

      ! heff_m
      heff_m(ibin) = 0.0_wp
      do ibrella = 1, nbrella
        heff_m(ibin) = heff_m(ibin) + h_km(ibin,ibrella) / g_km(ibin,ibrella)
      end do

      ! log_heff_m
      log_heff_m(ibin) = real(log(1.1921e-07),wp)
      if (heff_m(ibin) > 0.0_wp) &
        log_heff_m(ibin) = log(heff_m(ibin))

    end do


    ! calculate effective counts (Neff_km)
    !
    allocate(Neff_km(nbin,nbrella),log_Neff_km(nbin,nbrella))

    do ibin = 1, nbin

      do ibrella = 1, nbrella

        ! Neff_km
        Neff_km(ibin,ibrella) = N_k(ibrella) / g_km(ibin,ibrella)

        ! log_Neff_km
        log_Neff_km(ibin,ibrella) = real(log(1.1921e-07),wp)
        if (Neff_km(ibin,ibrella) > 0.0_wp) &
          log_Neff_km(ibin,ibrella) = log(Neff_km(ibin,ibrella))

      end do

    end do


    ! solve
    !
    allocate(pmf_m(nbin))

    allocate(log_den_km(nbin,nbrella), &
             log_den_m (nbin),         &
             log_num_m (nbin),         &
             f_k       (nbrella),      &
             f_k_new   (nbrella),      &
             log_prob_m(nbin),         &
             tmp_km    (nbrella,nbin))

    f_k(1:nbrella) = 0.0_wp

    check_convergence = 1.0e+8_wp

    iter = 0

    do while(check_convergence > option%tolerance)

      ! 1st equation
      do ibin = 1, nbin
        do ibrella = 1, nbrella
          log_den_km(ibin,ibrella) = f_k(ibrella) - bias_km(ibin,ibrella)
          log_den_km(ibin,ibrella) = log_den_km (ibin,ibrella) + &
                                     log_Neff_km(ibin,ibrella)
        end do
      end do

      call logsumexp(log_den_km, log_den_m)

      do ibin = 1, nbin
        log_num_m (ibin) = log_heff_m(ibin)
        log_prob_m(ibin) = log_num_m (ibin) - log_den_m(ibin)
      end do


      ! 2nd equation
      do ibin = 1, nbin
        do ibrella = 1, nbrella
          tmp_km(ibrella,ibin) = log_prob_m(ibin) - bias_km(ibin,ibrella)
        end do
      end do
      call logsumexp(tmp_km, f_k_new)
      f_k_new(1:nbrella) = -f_k_new(1:nbrella)

      fk = f_k_new(1)
      f_k_new(1:nbrella) = f_k_new(1:nbrella) - fk


      ! check convergence
      check_convergence = maxval(abs(f_k_new(1:nbrella) - f_k(1:nbrella))) / &
           std(f_k_new)

      f_k(1:nbrella) = f_k_new(1:nbrella)

      iter = iter + 1
      if (mod(iter, 100) == 0) then
        write(MsgOut,'(a,i5,a,i5,a,e17.9,a,e17.9)') &
           '[block ', iblock, '] ', &
           iter,   'th iteration  delta = ', &
           check_convergence, '  tolerance = ', option%tolerance
      end if

    end do

    if (mod(iter, 100) /= 0) &
      write(MsgOut,'(a,i5,a,i5,a)') &
         '[block ', iblock, '] ', &
         iter,   'th iteration (finish)'

    pmf_m(1:nbin) = - KB * option%temperature * log_prob_m(1:nbin)


    ! deallocate memory
    !

    deallocate(log_den_km, log_den_m)
    deallocate(log_num_m, f_k_new, f_k)
    deallocate(log_prob_m, tmp_km)

    deallocate(Neff_km, log_Neff_km)
    deallocate(heff_m, log_heff_m, g_km)
    deallocate(N_k, h_km)

    return

  end subroutine solve_wham_block

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_wham(option, output, pmf_m)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_output),          intent(in)    :: output
    type(s_pmf),             intent(in)    :: pmf_m(:)

    ! local variables
    integer                  :: file, ncol, j
    integer                  :: nbin, nbin_x, nbin_y
    integer                  :: ibin, ibin_x, ibin_y
    character(64)            :: fmt

    real(wp),    allocatable :: grid(:), grid2(:), center(:), center2(:)


    ! output pmf file
    !
    if (output%pmffile /= '') then

      call open_file(file, output%pmffile, IOFileOutputNew)

      if (option%dimension == 1) then

        call compute_grid_center(option, 1, grid, center)

        nbin = size(center)

        if (option%nblocks > 1) then
          ncol = 2
        else
          ncol = 1
        end if

        write(fmt,'(a,i0,a)') '(', ncol+1,'es25.16e3)'

        do ibin = 1, nbin
          write(file,fmt=fmt) center(ibin), (pmf_m(j)%v(ibin),j=1,ncol)
        end do

        deallocate(grid, center)

      else

        ! call compute_grid_center(option, 1, grid, center)
        ! call compute_grid_center(option, 2, grid2, center2)

        ! nbin_x = option%num_grids(1)-1
        ! nbin_y = option%num_grids(2)-1

        ! do ibin_x = 1, nbin_x
        !   do ibin_y = 1, nbin_y
        !     write(file,*) center(ibin_x), center2(ibin_y), &
        !       pmf_m(1)%v((ibin_x-1)+(ibin_y-1)*nbin_x+1)
        !   end do
        !   write(file,*)
        ! end do

        if (option%nblocks > 1) &
          call error_msg('Output_Wham> n-block is not supported in 2D')

        nbin_x = option%num_grids(1)-1
        nbin_y = option%num_grids(2)-1
        
        write(fmt,'(a,i0,a)') '(',nbin_x, 'es25.16e3)'
        
        do ibin_y = 1, nbin_y
          write(file,fmt=fmt) &
               (pmf_m(1)%v((ibin_x-1)+(ibin_y-1)*nbin_x+1),ibin_x=1,nbin_x)
        end do

      end if

      call close_file(file)

    end if


    ! Output summary
    !
    if (output%pmffile /= '') then
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '  [pmffile] ' // trim(output%pmffile)
      if (option%dimension == 1) then
        write(MsgOut,'(A)') '    Column 1: coordinates of bin centers'
        write(MsgOut,'(A)') '    Column 2: Free energy profile at the corresponding bin estimated using 1st block of input data'
        write(MsgOut,'(A)') '    Column X: Free energy profile at the corresponding bin estimated using (X-1)th block of input data'
        write(MsgOut,'(A)') ''
      else
        write(MsgOut,'(A)') '    Row X and Column Y: free energy profile at a bin center located at center(X,Y)'
        write(MsgOut,'(A)') ''
      end if
    end if

    return

  end subroutine output_wham

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_grid_center(option, idim, grid, center)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    integer,                 intent(in)    :: idim
    real(wp),                allocatable   :: grid(:)
    real(wp),                allocatable   :: center(:)

    ! local variables
    real(wp)                 :: grid_min, grid_max, grid_range
    integer                  :: num_bins, num_grid
    integer                  :: i


    grid_min = option%grid_min(idim)
    grid_max = option%grid_max(idim)
    num_bins = option%num_grids(idim) -1
    num_grid = option%num_grids(idim)
    grid_range = grid_max - grid_min + num_grid * EPS

    if (num_grid <= 1) &
      call error_msg('Assign_Bin> ERROR: # of grid must be > 1.')

    allocate(grid(num_grid), center(num_bins))

    ! grid
    do i = 1, num_grid
      grid(i) = grid_min + grid_range * ((i-1) / real(num_bins,wp))
    end do

    ! center
    do i = 1, num_bins
      center(i) = (grid(i) + grid(i+1)) * 0.5_wp
    end do

    return

  end subroutine compute_grid_center

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine exec_fhandle_1d(option, umbrella_center, k_constant, &
                             center, bias_km)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: umbrella_center(:)
    real(wp),                intent(in)    :: k_constant(:)
    real(wp),                intent(in)    :: center
    real(wp),                intent(inout) :: bias_km

    ! local variables
    real(wp)                 :: K
    integer                  :: idim, ndim


    ndim = size(umbrella_center)

    bias_km = 0.0_wp
    do idim = 1, ndim
      K = k_constant(idim)
      if (option%is_periodic(1)) then
        bias_km = bias_km + K*(periodic(center, umbrella_center(idim), option%box_size(1))**2.0_wp)
      else
        bias_km = bias_km + K*((center - umbrella_center(idim))**2.0_wp)
      end if
    end do

    return

  end subroutine exec_fhandle_1d

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine exec_fhandle_2d(option, umbrella_center, k_constant, &
                             center_x, center_y, bias_km)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: umbrella_center(:)
    real(wp),                intent(in)    :: k_constant(:)
    real(wp),                intent(in)    :: center_x
    real(wp),                intent(in)    :: center_y
    real(wp),                intent(inout) :: bias_km

    ! local variables
    real(wp)                 :: K1, K2
    real(wp)                 :: d1, d2


    K1 = k_constant(1)
    K2 = k_constant(2)

    if (option%is_periodic(1)) then
      d1 = periodic(center_x, umbrella_center(1), option%box_size(1))
    else
      d1 = center_x - umbrella_center(1)
    end if

    if (option%is_periodic(2)) then
      d2 = periodic(center_y, umbrella_center(2), option%box_size(2))
    else
      d2 = center_y - umbrella_center(2)
    end if

    bias_km = K1*(d1**2.0_wp) + K2*(d2**2.0_wp)

    return

  end subroutine exec_fhandle_2d

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_dcd_cv(molecule, option, trajectory, func_idx)

    ! function
    real(wp)                 :: get_dcd_cv

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(in)    :: option
    type(s_trajectory),      intent(in)    :: trajectory
    integer,                 intent(in)    :: func_idx

    ! local variables
    integer                  :: func_no, func, func_sel(4)


    ! select restraint function
    !
    func_no = option%rest_func_no(func_idx)
    if (func_no > size(option%rest_funcs(:))) &
      call error_msg('Get_Dcd_Cv> bad restraint function No.')

    func          = option%rest_funcs(func_no)
    func_sel(1:4) = option%rest_sel_index(1:4, func_no)


    ! get cv
    !
    select case (func)

    case (RestraintsFuncPOSI)
      call error_msg( &
      'Get_Dcd_Cv> ERROR : RestraintsFuncPOSI : not supprted.')

    case (RestraintsFuncDIST, RestraintsFuncDISTCOM)
      get_dcd_cv = get_com_dist(molecule, trajectory, option%selatoms, func_sel)

    case (RestraintsFuncRMSD, RestraintsFuncRMSDCOM)
      call error_msg( &
      'Get_Dcd_Cv> ERROR : RestraintsFuncRMSD/RMSDCOM : not supprted.')

    case (RestraintsFuncANGLE, RestraintsFuncANGLECOM)
      get_dcd_cv = get_com_angl(molecule, trajectory, option%selatoms, func_sel)

    case (RestraintsFuncDIHED, RestraintsFuncDIHEDCOM)
      get_dcd_cv = get_com_dihe(molecule, trajectory, option%selatoms, func_sel)

    end select

    return

  end function get_dcd_cv

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_com_dist(molecule, trajectory, selatoms, sel_index)

    ! return value
    real(wp)                 :: get_com_dist

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_selatoms),        intent(in)    :: selatoms(:)
    integer,                 intent(in)    :: sel_index(1:4)

    ! local variables
    real(wp)                 :: c1(3), c2(3)


    ! atom group1
    c1 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(1))%idx)

    ! atom group2
    c2 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(2))%idx)

    ! compute distance
    get_com_dist = compute_dis(c1, c2)

    return

  end function get_com_dist

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_com_angl(molecule, trajectory, selatoms, sel_index)

    ! return value
    real(wp)                 :: get_com_angl

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_selatoms),        intent(in)    :: selatoms(:)
    integer,                 intent(in)    :: sel_index(1:4)

    ! local variables
    real(wp)                 :: c1(3), c2(3), c3(3)


    ! atom group1
    c1 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(1))%idx)

    ! atom group2
    c2 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(2))%idx)

    ! atom group3
    c3 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(3))%idx)

    ! compute angle
    get_com_angl = compute_ang(c1, c2, c3)

!TODO
!    if (get_com_angl < 0.0_wp) &
!      get_com_angl = get_com_angl + 360.0_wp

    return

  end function get_com_angl

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_com_dihe(molecule, trajectory, selatoms, sel_index)

    ! return value
    real(wp)                 :: get_com_dihe

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_selatoms),        intent(in)    :: selatoms(:)
    integer,                 intent(in)    :: sel_index(1:4)

    ! local variables
    real(wp)                 :: c1(3), c2(3), c3(3), c4(3)


    ! atom group1
    c1 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(1))%idx)

    ! atom group2
    c2 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(2))%idx)

    ! atom group3
    c3 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(3))%idx)

    ! atom group4
    c4 = compute_com(trajectory%coord, molecule%mass, &
                     selatoms(sel_index(4))%idx)

    ! compute dihedral angle
    get_com_dihe = compute_dih(c1, c2, c3, c4)

!TODO
!    if (get_com_dihe < 0.0_wp) &
!      get_com_dihe = get_com_dihe + 360.0_wp

    return

  end function get_com_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_cvfile(cvfile, nsteps)

    ! formal arguments
    character(*),            intent(in)    :: cvfile
    integer,                 intent(inout) :: nsteps

    ! local variables
    integer                  :: file
    character(MaxFilename)   :: filename
    character(MaxLine)       :: line


    filename = get_replicate_name1(cvfile, 1)

    call open_file(file, filename, IOFileInput)

    nsteps = 0
    do while(.true.)
      read(file,'(a)',end=10) line
      nsteps = nsteps + 1
    end do

10  call close_file(file)

    write(MsgOut,'(a,i12)') '  number of time steps       : ', nsteps
    write(MsgOut,'(a)')   ''

    return

  end subroutine check_cvfile

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_dcdfile(dcdfile, nsteps, natom)

    ! formal arguments
    character(*),            intent(in)    :: dcdfile
    integer,                 intent(inout) :: nsteps
    integer,                 intent(inout) :: natom

    ! local variables
    type(s_trj_file)         :: file
    integer                  :: i, file_size, hdr_size, step_size
    integer(4)               :: icntrl(20), ntitle
    character(MaxFilename)   :: filename
    character(80)            :: title(10)
    character(4)             :: hdr
    logical                  :: exist

    integer(8)               :: ftell


    filename = get_replicate_name1(dcdfile, 1)

    call open_trj(file, filename, TrjFormatDCD, TrjTypeCoorBox, IOFileInput)

    read(file%unit_no) hdr, icntrl(1:20)
    read(file%unit_no) ntitle,(title(i),i=1,ntitle)
    read(file%unit_no) natom

    ! check header size
    hdr_size = &
         4 + 4 + 20*4      + 4 + &  ! => read() hdr, icntrl
         4 + 4 + 80*ntitle + 4 + &  ! => read() ntitle, title(:)
         4 + 4             + 4      ! => read() natom

    ! check trajectory step size
    step_size = (4 + 4 * natom + 4) * 3

    ! check file size
#ifdef KCOMP
    inquire(file%unit_no,flen=file_size)
#else
    inquire(file%unit_no,size=file_size)
#endif

    if (mod(file_size - hdr_size, step_size) /= 0) &
      step_size = step_size + 4 + 8 * 6 + 4   ! box

    nsteps = (file_size - hdr_size) / step_size

    write(MsgOut,'(a,i12)') '  number of trajectory steps       : ', nsteps
    write(MsgOut,'(a)')   ''

    call close_trj(file)

    return

  end subroutine check_dcdfile

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name1(filename, no)

    ! return
    character(Maxfilename)   :: get_replicate_name1

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Replicate_Name1> Syntax error.')

    write(get_replicate_name1, '(a,i0,a)') &
          filename(:bl-1),no,filename(br+1:len_trim(filename))

    return

  end function get_replicate_name1

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name2(filename, no1, no2)

    ! return
    character(Maxfilename)   :: get_replicate_name2

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no1
    integer,                 intent(in)    :: no2

    ! local variables
    integer                  :: bl1, br1, bl2, br2


    bl1 = index(filename, '{')
    br1 = index(filename, '}')
    bl2 = index(filename, '{', back=.true.)
    br2 = index(filename, '}', back=.true.)

    if (bl1 == 0 .or. br1 == 0 .or. bl1 == bl2 .or. br1 == br2 .or. bl1 > br1) &
      call error_msg('Get_Replicate_Name2> Syntax error.')

    write(get_replicate_name2, '(a,i0,a,i0,a)') &
         filename(     :bl1-1),no1, &
         filename(br1+1:bl2-1),no2, &
         filename(br2+1:     )

    return

  end function get_replicate_name2

  !======1=========2=========3=========4=========5=========6=========7=========8

  function periodic(x, center, box_size)

    ! return value
    real(wp)                 :: periodic

    ! formal arguments
    real(wp)                 :: x
    real(wp)                 :: center
    real(wp)                 :: box_size


    periodic = x - center
    periodic = periodic - nint(periodic/box_size)*box_size

    return

  end function periodic

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine logsumexp(x, s)

    ! formal arguments
    real(wp),                intent(in)    :: x(:,:)
    real(wp),                intent(inout) :: s(:)

    ! local variables
    integer                  :: irow, icol, nrow, ncol
    real(wp)                 :: max_x, exp_total


    nrow = size(x(1,:))
    ncol = size(x(:,1))

    !$omp parallel do private(max_x, irow, exp_total)
    !

    do icol = 1, ncol
      max_x = maxval(x(icol,:))

      exp_total = 0.0_wp
      do irow = 1, nrow
        exp_total = exp_total + exp(x(icol,irow)-max_x)
      end do

      s(icol) = log(exp_total) + max_x
    end do

    !
    !$omp end parallel do

    return

  end subroutine logsumexp

  !======1=========2=========3=========4=========5=========6=========7=========8

  function std(x)

    ! return value
    real(wp)                 :: std

    ! formal arguments
    real(wp),                intent(in)    :: x(:)

    ! local variables
    integer                  :: i, n
    real(wp)                 :: xa, xc


    n = size(x(:))

    xa = 0.0_wp
    do i = 1, n
      xa = xa + x(i)
    end do
    xa = xa / real(n,wp)

    xc = 0.0_wp
    do i = 1, n
      xc = xc + (x(i) - xa)**2
    end do
    xc = xc / real(n-1,wp)

    std = sqrt(xc)

    return

  end function std

end module wa_analyze_mod
