!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ma_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ma_analyze_mod

  use ma_option_str_mod
  use ma_matrix_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_pdb_mod
  use string_mod
  use messages_mod
  use constants_mod
#ifdef OMP
  use omp_lib
#endif

  implicit none
  private

  ! structures
  type s_data_k
    real(wp),      allocatable :: v(:,:)       ! (nstep, ndim)
    real(wp),      allocatable :: vcrd(:,:)    ! (natom*3, nstep)
    real(wp),      allocatable :: vtarget(:)   ! (natom*3)
    real(wp),      allocatable :: vrefene(:)   ! (natom*3)
  end type s_data_k

  type s_u_kl
    real(wp),      allocatable :: v(:)         ! (nstep)
  end type s_u_kl

  type s_f_k
    real(wp),      allocatable :: v(:)         ! (nbrella)
  end type s_f_k

  type s_bin_k
    integer,       allocatable :: v(:)         ! (nstep)
    real(wp),      allocatable :: center_x(:)  ! (nbin_x,ndim)
    real(wp),      allocatable :: center_y(:)  ! (nbin_y,ndim)
  end type s_bin_k

  type s_pmf
    real(wp),      allocatable :: v(:)         ! (maximum nbin)
  end type s_pmf

  ! constants
  real(wp),        parameter   :: KB   = 0.00198719168260038_wp

  ! module variables
  real(wp),    allocatable :: g_rep_N_k  (:,:) ! (nstp, nrep)
  real(wp),    allocatable :: g_rep_f_k  (:,:) ! (nstp, nrep)
  real(wp),    allocatable :: g_sqz_u_kln(:,:) ! (nstp, nrep)
  real(wp),    allocatable :: g_rep_u_kn (:,:) ! (nstp, nrep)
  real(wp),    allocatable :: g_temp     (:,:) ! (nstp, nrep)

  !$omp threadprivate (g_rep_N_k)
  !$omp threadprivate (g_rep_f_k)
  !$omp threadprivate (g_sqz_u_kln)
  !$omp threadprivate (g_rep_u_kn)
  !$omp threadprivate (g_temp)

  ! subroutines
  public  :: analyze

  private :: build_data_k_from_cv
  private :: build_data_k_from_ene
  private :: add_data_k_from_target
  private :: build_data_k_from_dcd
  private :: build_data_k_from_dcd_posi_readref 
  private :: build_u_kl_from_cv
  private :: build_u_kl_from_ene
  private :: build_u_kl_from_posi_readref
  private :: solve_mbar
  private :: solve_mbar_block
  private :: assign_bin
  private :: compute_pmf
  private :: compute_weight
  private :: compute_pmf_block
  private :: compute_weight_block
  private :: output_mbar

  private :: exec_fhandle
  private :: get_dcd_cv
  private :: get_com_dist
  private :: get_com_angl
  private :: get_com_dihe
  private :: check_cvfile
  private :: check_dcdfile
  private :: get_replicate_name1
  private :: get_replicate_name2
  private :: periodic
  private :: mbar_log_wi_jn
  private :: logsumexp2d
  private :: logsumexp1dn
  private :: logsumexp1d2
  private :: std
  private :: minval_from_not_nan
  private :: assign_bin_1d

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
    type(s_option),          intent(inout) :: option

    ! local variables
    type(s_data_k),   allocatable :: data_k(:)     ! (nbrella)
    type(s_u_kl),     allocatable :: u_kl(:,:)     ! (nbrella, nbrella)
    type(s_u_kl),     allocatable :: u_k(:)        ! (nbrella)
    type(s_f_k),      allocatable :: f_k(:), f_k_tmp(:) ! (nblocks)
    type(s_bin_k),    allocatable :: bin_k(:)      ! (nbrella)
    type(s_pmf),      allocatable :: pmf(:)        ! (nblocks)
    real(wp),         allocatable :: weight_k(:,:) ! (nstep, nbrella)
    real(wp),         allocatable :: time_k(:,:)   ! (nstep, nbrella)
    integer                       :: i, j, iblock, nstep
    real(wp)                      :: f_sum

    real(8):: time_start, time_end


    ! check only
    !
    if (option%check_only) &
      return

    ! read reference files if Cartesian CV
    !
    if (option%read_ref_path) then
      if (input%pathfile /= '') then
        call readref_path(input%pathfile, option)
      else
        call error_msg('Analyze> Please set pathfile in [INPUT]')
      endif
    else if (option%read_ref_pdb) then
      if (input%reffile /= '') then
        call readref_pdb(input%pdbfile, molecule, option)
      else
        call error_msg('Analyze> Please set reffile in [INPUT]')
      endif
    endif

    ! build data_k
    !
    if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then

      if (input%cvfile /= '') then

        call build_data_k_from_cv(input%cvfile, option, data_k, time_k)
    
      else if (input%dcdfile /= '') then
    
        if (option%read_ref_path .or. option%read_ref_pdb) then
          call build_data_k_from_dcd_posi_readref(input%dcdfile, molecule, option, data_k, &
                                           time_k)
        else
          call build_data_k_from_dcd(input%dcdfile, molecule, option, data_k, time_k)
        endif
    
      end if

        if (input%refenefile /= '') &
         call build_data_k_from_refene(input%refenefile, option, data_k, time_k)

    else

      call build_data_k_from_ene(input%cvfile, option, data_k, time_k)
      
    end if


    ! build target data data_k()%vtarget() for weight calculation
    !
    do i = 1, option%num_replicas
      if (input%dcdfile /= '' .and. (option%read_ref_path .or. option%read_ref_pdb)) then
        nstep = size(data_k(i)%vcrd(1,:))
      else
        nstep = size(data_k(i)%v(:,1))
      end if
      allocate(data_k(i)%vtarget(nstep))
    end do

    if (input%targetfile /= '') then    
      call add_data_k_from_target(input%targetfile, option, data_k)
    else
      do i = 1, option%num_replicas
        nstep = size(data_k(i)%vtarget(:))

        if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then        
          data_k(i)%vtarget(:) = 0.0_wp

        else if (option%input_type == InputTypeEneSingle .or. option%input_type == InputTypeREMD) then
          do j = 1, nstep
            data_k(i)%vtarget(j) = data_k(i)%v(j, 1)
          end do

        else if (option%input_type == InputTypeEnePair .or. option%input_type == InputTypeFEP) then
          data_k(i)%vtarget(:) = 0.0_wp

        else if (option%input_type == InputTypeEneAll .or. option%input_type == InputTypeREST .or. option%input_type == InputTypeMBGO) then
          do j = 1, nstep
            data_k(i)%vtarget(j) = data_k(i)%v(j, i)
          end do
        
        end if

      end do
    end if

    ! build u_kl
    !
    if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then

      if (option%read_ref_path .or. option%read_ref_pdb) then

        call build_u_kl_from_posi_readref(option, data_k, u_kl)

      else

        if (input%refenefile /= '') then
          call build_u_kl_from_cv_ene(option, data_k, u_kl)        
        else
          call build_u_kl_from_cv(option, data_k, u_kl)        
        end if

      end if

    else

        call build_u_kl_from_ene(option, data_k, u_kl)

    end if

    ! build u_k
    !
    allocate(u_k(option%num_replicas))
    do i = 1, option%num_replicas
      nstep = size(data_k(i)%vtarget(:))
      allocate(u_k(i)%v(nstep))
      do j = 1, nstep
        u_k(i)%v(j) = data_k(i)%vtarget(j) / (KB * option%target_temperature)
      end do
    end do

    ! solve mbar
    !
    if (option%input_type == InputTypeEnePair .or. option%input_type == InputTypeFEP) then
      allocate(f_k(option%nblocks))
      do iblock = 1, option%nblocks
        allocate(f_k(iblock)%v(option%num_replicas))
        do i = 1, option%num_replicas
          f_k(iblock)%v(i) = 0.0_wp
        end do
      end do
      !f_sum = 0.0_wp
      do i = 1, option%num_replicas-1
        call solve_mbar(option, u_kl(i:i+1, i:i+1), f_k_tmp)
        !f_sum = f_sum + f_k(1)%v(2)
        do iblock = 1, option%nblocks
          f_k(iblock)%v(i+1) = f_k(iblock)%v(i) + f_k_tmp(iblock)%v(2)
        end do
        write(*, *)'f_k_tmp = ', f_k_tmp(1)%v
        write(*, *)'f_k = ', f_k(1)%v
        deallocate(f_k_tmp)
      end do
    else
!      time_start = omp_get_wtime()
      call solve_mbar(option, u_kl, f_k)
!      time_end = omp_get_wtime()
      ! write(*,*) 'Elapsed time in sec = ', time_end - time_start
    end if


    ! solve mbar weight
    !
    if (option%input_type /= InputTypeEnePair .and. option%input_type /= InputTypeFEP) then
      call compute_weight(option, u_kl, u_k, bin_k, f_k, weight_k)
    end if


    ! assign bin
    !
    if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then
      call assign_bin(option, data_k, bin_k)
    end if


    ! solve mbar pmf
    !
    if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then
      call compute_pmf(option, u_kl, u_k, bin_k, f_k, pmf)
    end if


    ! output f_k and pmf
    !
    call output_mbar(option, output, bin_k(1), f_k, pmf, weight_k, time_k)


    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine readref_path(pathfile, option)
    ! formal arguments
    character(*),            intent(in)    :: pathfile
    type(s_option),          intent(inout) :: option

    ! local variables
    real(wp)                               :: dummy

    integer                  :: nbrella, ibrella, ndim, nfunc, ifunc, igroup
    integer                  :: natoms3, file
    integer                  :: i, j, k, l


    write(MsgOut,'(a)') 'Read_Ref_Rpath> '

    ndim    = option%dimension
    nfunc   = 1
    igroup  = option%rest_sel_index(1,nfunc)
    natoms3 = size(option%selatoms(igroup)%idx)*3
    nbrella = option%rest_nreplica(nfunc)

    allocate(option%rest_ref_coord(1:natoms3,1:nbrella))

    option%rest_ref_coord(:,:)=0.0_wp

    write(MsgOut,'(a,a)') '  read pathfile: ',trim(pathfile)

    call open_file(file, &
                   pathfile, &
                   IOFileInput)

    do i = 1, nbrella
      read(file,*) dummy, (option%rest_ref_coord(j,i),j=1,natoms3)
    end do

    call close_file(file)

    return

  end subroutine readref_path

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine readref_pdb(pdbfile, molecule, option)
    ! formal arguments
    character(*),            intent(in)    :: pdbfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(inout) :: option

    ! local variables
    type(s_pdb)                            :: pdb
    character(MaxFilename)                 :: filename
    character(MaxLineLong_CV)              :: line
    real(wp)                               :: dummy

    integer                  :: nbrella, ibrella, ndim, nfunc, ifunc, igroup
    integer                  :: natoms
    integer                  :: i, j, jj, iatom

    write(MsgOut,'(a)') 'Read_Ref_Pdb> '

    ndim    = option%dimension
    nfunc   = 1
    igroup  = option%rest_sel_index(1,nfunc)
    natoms  = size(option%selatoms(igroup)%idx(:))
    nbrella = option%rest_nreplica(nfunc)

    allocate(option%rest_ref_coord(1:natoms*3,1:nbrella))

    option%rest_ref_coord(:,:)=0.0_wp


    do i = 1, nbrella

      filename = get_replicate_name1(pdbfile, i)
      write(MsgOut,'(a,a)') '  read pdbfile: ',trim(filename)

      call input_pdb(filename, pdb)
      if (molecule%num_atoms /= pdb%num_atoms) &
        call error_msg('Readref_PDB> # of atoms is mismatch ref/system')

      jj = 0
      do j = 1, natoms
        
        iatom = option%selatoms(igroup)%idx(j)
        jj = jj + 1 
        option%rest_ref_coord(jj, i) = pdb%atom_coord(1,iatom)
        jj = jj + 1 
        option%rest_ref_coord(jj, i) = pdb%atom_coord(2,iatom)
        jj = jj + 1 
        option%rest_ref_coord(jj, i) = pdb%atom_coord(3,iatom)

      end do

    end do

   return

  end subroutine readref_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine build_data_k_from_cv(cvfile, option, data_k, time_k)

    ! formal arguments
    character(*),            intent(in)    :: cvfile
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)
    real(wp),                allocatable   :: time_k(:,:)

    ! local variables
    integer                  :: file, tim, i, j, k, istep
    integer                  :: nbrella, ibrella, ndim, nstep, nfunc
    character(MaxFilename)   :: filename
    character(MaxLineLong_CV) :: line


    write(MsgOut,'(a)') 'Build_Data_K_From_Cv> '

    ndim    = option%dimension

    if (ndim == 1) then
      nbrella = option%rest_nreplica(option%rest_func_no(1, 1))
    else if (ndim == 2) then
      nbrella = 1
      do i = 1, 2
        nbrella = nbrella * option%rest_nreplica(option%rest_func_no(i, 1))
      end do
    else
      call error_msg('Build_Data_K_From_Cv> dimension > 2 not supported')
    end if

    
    ! check nstep
    !
    call check_cvfile(cvfile, nstep)


    ! allocate data_k and time_k
    !
    allocate(data_k(nbrella))
    allocate(time_k(nstep, nbrella))


    ! read cv file and setup data_k
    !
    ibrella = 0

    if (ndim == 1) then

      nfunc = option%rest_func_num(1)

      do i = 1, nbrella

        ibrella = ibrella + 1
        allocate(data_k(ibrella)%v(nstep,nfunc))

        filename = get_replicate_name1(cvfile, ibrella)
        write(MsgOut,'(a,a)') '  read cv file: ',trim(filename)

        call open_file(file, &
                       filename, &
                       IOFileInput)

        ! do istep = 1, nstep          
        !   read(file,*) time_k(istep, ibrella), (data_k(ibrella)%v(istep,k),k=1,nfunc)
        ! end do
        istep = 0
        do while(.true.)
          read(file,'(a)',end=12) line
          if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
            istep = istep + 1
            read(line,*) time_k(istep, ibrella), (data_k(ibrella)%v(istep,k),k=1,nfunc)
          end if
        end do

12      call close_file(file)

      end do

    else if (ndim == 2) then

      nfunc = option%rest_func_num(1) + option%rest_func_num(2)
      ! if (nfunc == 1) &
      !   call error_msg('Build_Data_K_From_Cv> # of rest_func must be 2 on dim=2')

      do i = 1, option%rest_nreplica(option%rest_func_no(1, 1))
        do j = 1, option%rest_nreplica(option%rest_func_no(2, 1))

          ibrella = ibrella + 1
          allocate(data_k(ibrella)%v(nstep,nfunc))

          filename = get_replicate_name1(cvfile, ibrella)
          write(MsgOut,'(a,a)') '  read cv file: ',trim(filename)

          call open_file(file, &
                         filename, &
                         IOFileInput)

          do istep = 1, nstep
            read(file,*) time_k(istep, ibrella), (data_k(ibrella)%v(istep,k),k=1,nfunc)
          end do

          call close_file(file)

        end do
      end do

    end if

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

  end subroutine build_data_k_from_cv

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine build_data_k_from_ene(cvfile, option, data_k, time_k)

    ! formal arguments
    character(*),            intent(in)    :: cvfile
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)
    real(wp),                allocatable   :: time_k(:,:)

    ! local variables
    integer                  :: file, tim, i, j, k, istep
    integer                  :: nbrella, ibrella, nstep, nfunc
    character(MaxFilename)   :: filename
    character(MaxLineLong_CV) :: line


    write(MsgOut,'(a)') 'Build_Data_K_From_Ene> '

    
    ! check nstep
    !
    call check_cvfile(cvfile, nstep)


    ! read cv file and setup data_k
    !
    nbrella = option%num_replicas


    ! allocate data_k and time_k
    !
    allocate(data_k(nbrella))
    allocate(time_k(nstep, nbrella))


    if (option%input_type == InputTypeEneSingle .or. option%input_type == InputTypeREMD) then
      nfunc = 1
    else if (option%input_type == InputTypeEnePair .or. option%input_type == InputTypeFEP) then
      nfunc = 3
    else if (option%input_type == InputTypeEneAll .or. option%input_type == InputTypeREST .or. option%input_type == InputTypeMBGO) then
      nfunc = nbrella
    else
      call error_msg('Build_Data_K_From_Ene> unsupported input type')
    end if


    ibrella = 0
    do i = 1, nbrella

      ibrella = ibrella + 1
      allocate(data_k(ibrella)%v(nstep,nfunc))

      filename = get_replicate_name1(cvfile, ibrella)
      write(MsgOut,'(a,a)') '  read cv file: ',trim(filename)

      call open_file(file, &
                     filename, &
                     IOFileInput)

      istep = 0
      !do istep = 1, nstep
      do while(.true.)
        read(file,'(a)',end=11) line
        if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
          istep = istep + 1
          read(line,*) time_k(istep, ibrella), (data_k(ibrella)%v(istep,k),k=1,nfunc)
        end if
        ! read(file,*) time_k(istep, ibrella), (data_k(ibrella)%v(istep,k),k=1,nfunc)
      end do

11    call close_file(file)
!      call close_file(file)

    end do

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

  end subroutine build_data_k_from_ene

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine build_data_k_from_refene(refenefile, option, data_k, time_k)

    ! formal arguments
    character(*),            intent(in)    :: refenefile
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)
    real(wp),                allocatable   :: time_k(:,:)


    ! local variables
    integer                   :: file, tim, i, j, k, istep
    integer                   :: nbrella, ibrella, nstep, nfunc
    character(MaxLine)        :: filename
    character(MaxLineLong_CV) :: line


    write(MsgOut,'(a)') 'Build_Data_K_From_Ref_Ene> '

    
    ! check nstep
    !
    call check_cvfile(refenefile, nstep)


    ! read cv file and setup data_k
    !
    nbrella = option%num_replicas

    ibrella = 0
    do i = 1, nbrella

      ibrella = ibrella + 1
      allocate(data_k(ibrella)%vrefene(nstep))

      filename = get_replicate_name1(refenefile, ibrella)
      write(MsgOut,'(a,a)') '  read refene file: ',trim(filename)

      call open_file(file, &
                     filename, &
                     IOFileInput)

      istep = 0
      do while(.true.)
        read(file,'(a)',end=11) line
        if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
          istep = istep + 1
          read(line,*) time_k(istep, ibrella), data_k(ibrella)%vrefene(istep)
          data_k(ibrella)%vrefene(istep) = data_k(ibrella)%vrefene(istep) &
            /(KB * option%temperature(1))
        end if
      end do

11    call close_file(file)
!      call close_file(file)

    end do

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

  end subroutine build_data_k_from_refene

  !======1=========2=========3=========4=========5=========6=========7=========8  
  subroutine add_data_k_from_target(targetfile, option, data_k)

    ! formal arguments
    character(*),            intent(in)    :: targetfile
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)

    ! local variables
    integer                  :: file, tim, i, j, k, istep
    integer                  :: nbrella, ibrella, nstep, nfunc
    character(MaxFilename)   :: filename
    character(MaxLineLong_CV) :: line
    real(wp)                 :: dummy


    write(MsgOut,'(a)') 'Add_Data_K_From_Target> '

    
    ! check nstep
    !
    call check_cvfile(targetfile, nstep)


    ! read cv file and setup data_k
    !
    nbrella = option%num_replicas

    nfunc = 1

    ibrella = 0
    do i = 1, nbrella

      ibrella = ibrella + 1
      !allocate(data_k(ibrella)%vtarget(nstep))

      filename = get_replicate_name1(targetfile, ibrella)
      write(MsgOut,'(a,a)') '  read target energy file: ',trim(filename)

      call open_file(file, &
                     filename, &
                     IOFileInput)

      istep = 0
      !do istep = 1, nstep
      do while(.true.)
        read(file,'(a)',end=11) line
        if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
          istep = istep + 1
          read(line,*) dummy, data_k(ibrella)%vtarget(istep)
        end if
      end do

11    call close_file(file)

    end do

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

  end subroutine add_data_k_from_target

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_data_k_from_dcd(dcdfile, molecule, option, data_k, time_k)

    ! formal arguments
    character(*),            intent(in)    :: dcdfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)
    real(wp),                allocatable   :: time_k(:,:)

    ! local variables
    type(s_trajectory)       :: trajectory
    type(s_trj_file)         :: file
    integer                  :: i, j, k, istep
    integer                  :: nbrella, ibrella, ndim, nstep, natom, nfunc
    character(MaxFilename)   :: filename


    write(MsgOut,'(a)') 'Build_Data_K_From_Dcd> '

    ndim    = option%dimension

    if (ndim == 1) then
      nbrella = option%rest_nreplica(option%rest_func_no(1, 1))
    else if (ndim == 2) then
      nbrella = 1
      do i = 1, 2
        !!! todo: support for temperature
        nbrella = nbrella * option%rest_nreplica(option%rest_func_no(i, 1))
      end do
    else
      call error_msg('Build_Data_K_From_Cv> dimension > 2 not supported')
    end if


    ! check nstep and allocate trajectory
    !
    call check_dcdfile(dcdfile, nstep, natom)


    ! allocate data_k
    !
    allocate(data_k(nbrella))
    allocate(time_k(nstep, nbrella))


    if (natom /= option%num_atoms) &
      call error_msg( &
      'Build_Data_K_From_Dcd> Dcd atom count is different from PSF/PRMTOP.')

    call alloc_trajectory(trajectory, natom)


    ! read trajectory and setup data_k
    !
    ibrella = 0

    if (ndim == 1) then

      nfunc = option%rest_func_num(1)

      do i = 1, nbrella

        ibrella = ibrella + 1
        allocate(data_k(ibrella)%v(nstep,nfunc))

        filename = get_replicate_name1(dcdfile, ibrella)
        write(MsgOut,'(a,a)') '  read and analyze trajectory: ',trim(filename)

        call open_trj(file, &
                    filename,&
                    TrjFormatDCD,   &
                    TrjTypeCoor, &
                    IOFileInput)

!        call open_trj(file, &
!                    filename,&
!                    TrjFormatDCD,   &
!                    TrjTypeCoorBox, &
!                    IOFileInput)

        do istep = 1, nstep

          call read_trj(file, trajectory)

          do k = 1, nfunc
            data_k(ibrella)%v(istep,k) = &
                 get_dcd_cv(molecule, option, trajectory, k)
          end do

          time_k(istep, ibrella) = istep

        end do

        call close_trj(file)

      end do

    else if (ndim == 2) then

      nfunc = option%rest_func_num(1) + option%rest_func_num(2)
      ! if (nfunc == 1) &
      !   call error_msg('Build_Data_K_From_Cv> # of rest_func must be 2 on dim=2')

      do i = 1, option%rest_nreplica(option%rest_func_no(1, 1))
        do j = 1, option%rest_nreplica(option%rest_func_no(2, 1))

          ibrella = ibrella + 1
          allocate(data_k(ibrella)%v(nstep,nfunc))

          filename = get_replicate_name1(dcdfile, ibrella)
          write(MsgOut,'(a,a)') '  read and analyze trajectory: ',trim(filename)

          call open_trj(file, &
                      filename,&
                      TrjFormatDCD,   &
                      TrjTypeCoor, &
                      IOFileInput)

!          call open_trj(file, &
!                      filename,&
!                      TrjFormatDCD,   &
!                      TrjTypeCoorBox, &
!                      IOFileInput)

          do istep = 1, nstep

            call read_trj(file, trajectory)

            do k = 1, nfunc
              data_k(ibrella)%v(istep,k) = &
                   get_dcd_cv(molecule, option, trajectory, k)
            end do

            time_k(istep, ibrella) = istep

          end do

          call close_trj(file)

        end do
      end do

    end if

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine build_data_k_from_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine build_data_k_from_dcd_posi_readref &
                (dcdfile, molecule, option, data_k, time_k)

    ! formal arguments
    character(*),            intent(in)    :: dcdfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(in)    :: option
    type(s_data_k),          allocatable   :: data_k(:)
    real(wp),                allocatable   :: time_k(:,:)

    ! local variables
    type(s_trajectory)       :: trajectory
    type(s_trj_file)         :: file
    integer                  :: i, j, k, istep, jj, iatom
    integer                  :: nbrella, ibrella, ndim, nstep, natom, nfunc, igroup
    integer                  :: nselatom
    character(MaxFilename)   :: filename


    write(MsgOut,'(a)') 'Build_Data_K_From_Dcd_Posi_Readref> '

    ndim    = option%dimension
    nfunc   = 1
    igroup  = option%rest_sel_index(1,nfunc)
    nbrella = option%rest_nreplica(nfunc)

    nselatom  = size(option%selatoms(igroup)%idx(:))

    ! check nstep and allocate trajectory
    !
    call check_dcdfile(dcdfile, nstep, natom)

    ! allocate data_k
    !
    allocate(data_k(nbrella))
    allocate(time_k(nstep, nbrella))

    if (natom /= option%num_atoms) &
      call error_msg( &
      'Build_Data_K_From_Dcd> Dcd atom count is different from PSF/PRMTOP.')

    call alloc_trajectory(trajectory, natom)


    ! read trajectory and setup data_k
    !

    nfunc = option%rest_func_num(1)

    do i = 1, nbrella

      allocate(data_k(i)%vcrd(1:nselatom*3, 1:nstep))

      filename = get_replicate_name1(dcdfile, i)
      write(MsgOut,'(a,a)') '  read and analyze trajectory: ',trim(filename)

      call open_trj(file, &
                  filename,&
                  TrjFormatDCD,   &
                  TrjTypeCoor, &
                  IOFileInput)

!      call open_trj(file, &
!                  filename,&
!                  TrjFormatDCD,   &
!                  TrjTypeCoorBox, &
!                  IOFileInput)

!      do istep = 1, nstep
      do istep = 1, nstep

        call read_trj(file, trajectory)

        jj = 0
        do j = 1, nselatom

          iatom = option%selatoms(igroup)%idx(j)

          jj = jj + 1 
          data_k(i)%vcrd(jj, istep) = trajectory%coord(1,iatom)
          jj = jj + 1 
          data_k(i)%vcrd(jj, istep) = trajectory%coord(2,iatom)
          jj = jj + 1 
          data_k(i)%vcrd(jj, istep) = trajectory%coord(3,iatom)

        end do

        time_k(istep, i) = istep

      end do

      call close_trj(file)

    end do

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine build_data_k_from_dcd_posi_readref

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_u_kl_from_cv(option, data_k, u_kl)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)
    type(s_u_kl),            allocatable   :: u_kl(:,:)
 
    ! local variables
    real(wp),    allocatable :: umbrella_center(:,:), k_constant(:,:)
    logical,     allocatable :: is_periodic(:)
    real(wp),    allocatable :: box_size(:)
    integer                  :: nbrella, ibrella, ndim, nfunc, ifunc
    integer                  :: i, j, k, l


    write(MsgOut,'(a)') 'Build_U_Kl_From_CV> '

    ndim    = option%dimension
    !nfunc   = size(option%rest_func_no)

    if (ndim == 1) then
      nbrella = option%rest_nreplica(option%rest_func_no(1, 1))
    else if (ndim == 2) then
      nbrella = 1
      do i = 1, 2
        !!! todo: support for temperature
        nbrella = nbrella * option%rest_nreplica(option%rest_func_no(i, 1))
      end do
    else
      call error_msg('Build_U_Kl_From_CV> dimension > 2 not supported')
    end if

    allocate(u_kl(nbrella, nbrella))

    if (ndim == 1) then

      nfunc = option%rest_func_num(1)

      allocate(umbrella_center(nbrella, nfunc), &
               k_constant(nbrella, nfunc), &
               is_periodic(nfunc), &
               box_size(nfunc))

      do i = 1, nfunc
        umbrella_center(:,i) = option%rest_references(:,option%rest_func_no(1, i))
        k_constant(:,i)      = option%rest_constants (:,option%rest_func_no(1, i))
        is_periodic(i)       = option%is_periodic(option%rest_func_no(1, i))
        box_size(i)          = option%box_size(option%rest_func_no(1, i))
      end do

    else

      nfunc = option%rest_func_num(1) + option%rest_func_num(2)

      allocate(umbrella_center(nbrella, nfunc), &
               k_constant(nbrella, nfunc), &
               is_periodic(nfunc), &
               box_size(nfunc))

      ibrella = 0
      do i = 1, option%rest_nreplica(option%rest_func_no(1, 1))
        do j = 1, option%rest_nreplica(option%rest_func_no(2, 1))
          ibrella = ibrella + 1
          umbrella_center(ibrella,1:2) = &
               (/option%rest_references(i, option%rest_func_no(1, 1)), &
                 option%rest_references(j, option%rest_func_no(2, 1))/)
          k_constant(ibrella,1:2) = &
               (/option%rest_constants(i, option%rest_func_no(1, 1)), &
                 option%rest_constants(j, option%rest_func_no(2, 1))/)
        end do
      end do

      ifunc = 0
      do i = 1, option%rest_func_num(1)
        ifunc = ifunc + 1
        is_periodic(ifunc) = option%is_periodic(option%rest_func_no(1, i))
        box_size(ifunc)    = option%box_size(option%rest_func_no(1, i))
      end do
      do i = 1, option%rest_func_num(2)
        ifunc = ifunc + 1
        is_periodic(ifunc) = option%is_periodic(option%rest_func_no(2, i))
        box_size(ifunc)    = option%box_size(option%rest_func_no(2, i))
      end do

    end if


    do k = 1, nbrella
      do l = 1, nbrella
        call exec_fhandle(option, umbrella_center(l,:), k_constant(l,:), &
                          data_k(k), u_kl(k,l), is_periodic, box_size)
      end do
    end do

    deallocate(umbrella_center, k_constant, is_periodic, box_size)

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine build_u_kl_from_cv

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine build_u_kl_from_cv_ene(option, data_k, u_kl)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)
    type(s_u_kl),            allocatable   :: u_kl(:,:)
 
    ! local variables
    real(wp),    allocatable :: umbrella_center(:,:), k_constant(:,:)
    logical,     allocatable :: is_periodic(:)
    real(wp),    allocatable :: box_size(:)
    integer                  :: nbrella, ibrella, ndim, nfunc, ifunc
    integer                  :: i, j, k, l


    write(MsgOut,'(a)') 'Build_U_Kl_From_CV_with_Ene> '

    ndim    = option%dimension

    if (ndim == 1) then
      nbrella = option%rest_nreplica(option%rest_func_no(1, 1))
    else if (ndim == 2) then
      nbrella = 1
      do i = 1, 2
        !!! todo: support for temperature
        nbrella = nbrella * option%rest_nreplica(option%rest_func_no(i, 1))
      end do
    else
      call error_msg('Build_U_Kl_From_CV> dimension > 2 not supported')
    end if

    allocate(u_kl(nbrella, nbrella))

    if (ndim == 1) then

      nfunc = option%rest_func_num(1)

      allocate(umbrella_center(nbrella, nfunc), &
               k_constant(nbrella, nfunc), &
               is_periodic(nfunc), &
               box_size(nfunc))

      do i = 1, nfunc
        umbrella_center(:,i) = option%rest_references(:,option%rest_func_no(1, i))
        k_constant(:,i)      = option%rest_constants (:,option%rest_func_no(1, i))
        is_periodic(i)       = option%is_periodic(option%rest_func_no(1, i))
        box_size(i)          = option%box_size(option%rest_func_no(1, i))
      end do

    else

      nfunc = option%rest_func_num(1) + option%rest_func_num(2)

      allocate(umbrella_center(nbrella, nfunc), &
               k_constant(nbrella, nfunc), &
               is_periodic(nfunc), &
               box_size(nfunc))

      ibrella = 0
      do i = 1, option%rest_nreplica(option%rest_func_no(1, 1))
        do j = 1, option%rest_nreplica(option%rest_func_no(2, 1))
          ibrella = ibrella + 1
          umbrella_center(ibrella,1:2) = &
               (/option%rest_references(i, option%rest_func_no(1, 1)), &
                 option%rest_references(j, option%rest_func_no(2, 1))/)
          k_constant(ibrella,1:2) = &
               (/option%rest_constants(i, option%rest_func_no(1, 1)), &
                 option%rest_constants(j, option%rest_func_no(2, 1))/)
        end do
      end do

      ifunc = 0
      do i = 1, option%rest_func_num(1)
        ifunc = ifunc + 1
        is_periodic(ifunc) = option%is_periodic(option%rest_func_no(1, i))
        box_size(ifunc)    = option%box_size(option%rest_func_no(1, i))
      end do
      do i = 1, option%rest_func_num(2)
        ifunc = ifunc + 1
        is_periodic(ifunc) = option%is_periodic(option%rest_func_no(2, i))
        box_size(ifunc)    = option%box_size(option%rest_func_no(2, i))
      end do

    end if

    do k = 1, nbrella
      do l = 1, nbrella
        call exec_fhandle(option, umbrella_center(l,:), k_constant(l,:), &
                          data_k(k), u_kl(k,l), is_periodic, box_size)
        u_kl(k,l)%v(:) = (u_kl(k,l)%v(:) + data_k(k)%vrefene(:))
      end do
    end do

    deallocate(umbrella_center, k_constant, is_periodic, box_size)

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''         

    return

  end subroutine build_u_kl_from_cv_ene

  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine build_u_kl_from_ene(option, data_k, u_kl)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)
    type(s_u_kl),            allocatable   :: u_kl(:,:)
 
    ! local variables
    integer                  :: nbrella, ibrella, ndim, nfunc, ifunc, nstep, istep
    integer                  :: i, j, k, l


    write(MsgOut,'(a)') 'Build_U_Kl_From_Ene> '

    nbrella = option%num_replicas
    allocate(u_kl(nbrella, nbrella))

    if (option%input_type == InputTypeEneSingle .or. option%input_type == InputTypeREMD) then
      nfunc = 1
    else if (option%input_type == InputTypeEnePair .or. option%input_type == InputTypeFEP) then
      nfunc = 2
    else if (option%input_type == InputTypeEneAll .or. option%input_type == InputTypeREST .or. option%input_type == InputTypeMBGO) then
      nfunc = nbrella
    else
      call error_msg('Build_U_KL_From_Ene> unsupported input type')
    end if

    nstep = size(data_k(1)%v(:,1))

    if (option%input_type == InputTypeEneSingle .or. option%input_type == InputTypeREMD) then
      do k = 1, nbrella
        do l = 1, nbrella
          allocate(u_kl(k,l)%v(nstep))
          do istep = 1, nstep
            u_kl(k,l)%v(istep) = data_k(k)%v(istep,1) / (KB * option%temperature(l))
          end do
        end do
      end do
    else if (option%input_type == InputTypeEnePair .or. option%input_type == InputTypeFEP) then
      do k = 1, nbrella
        if (k /= 1) then
          allocate(u_kl(k,k-1)%v(nstep))
          do istep = 1, nstep
            u_kl(k,k-1)%v(istep) = (data_k(k)%v(istep,1) + data_k(k)%v(istep,2)) / (KB * option%temperature(k-1))
          end do
        end if
        allocate(u_kl(k,k)%v(nstep))
        do istep = 1, nstep
          u_kl(k,k)%v(istep) = data_k(k)%v(istep,1) / (KB * option%temperature(k))
        end do
        if (k /= nbrella) then
          allocate(u_kl(k,k+1)%v(nstep))
          do istep = 1, nstep
            u_kl(k,k+1)%v(istep) = (data_k(k)%v(istep,1) + data_k(k)%v(istep,3)) / (KB * option%temperature(k+1))
          end do
        end if
      end do
    else if (option%input_type == InputTypeEneAll .or. option%input_type == InputTypeREST .or. option%input_type == InputTypeMBGO) then
      do k = 1, nbrella
        do l = 1, nbrella
          allocate(u_kl(k,l)%v(nstep))
          do istep = 1, nstep
            u_kl(k,l)%v(istep) = data_k(k)%v(istep,l) / (KB * option%temperature(1))
          end do
        end do
      end do
    end if

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine build_u_kl_from_ene

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine build_u_kl_from_posi_readref(option, data_k, u_kl)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)
    type(s_u_kl),            allocatable   :: u_kl(:,:)
 
    ! local variables
    integer                  :: nbrella, ibrella, ndim, nstep, natom, nfunc, igroup
    integer                  :: nselatom3
    integer                  :: i, j, k, l
    real(wp) , allocatable   :: ref(:)
    real(wp)                 :: k_constant, v, dtmp
    real(wp)                 :: KBT, BETA


    write(MsgOut,'(a)') 'Build_U_Kl_From_Posi_Readref> '

    ndim    = option%dimension
    nfunc   = 1
    igroup  = option%rest_sel_index(1,nfunc)
    nbrella = option%rest_nreplica(nfunc)

    !KBT  = KB * option%temperature
    !BETA = 1.0_wp / KBT

    nselatom3  = size(option%selatoms(igroup)%idx(:))*3

    nstep = size(data_k(1)%vcrd(1,:))

    allocate(u_kl(1:nbrella, 1:nbrella), &
             ref(1:nselatom3))

    do l = 1, nbrella

      !k_constant       = option%rest_constants (l, nfunc) * BETA
      k_constant       = option%rest_constants(l, nfunc) / (KB * option%temperature(l))
      ref(1:nselatom3) = option%rest_ref_coord(1:nselatom3, l)

      do k = 1, nbrella
        allocate(u_kl(k, l)%v(1:nstep))

        do i = 1, nstep
          v = 0.0_wp
          do j = 1, nselatom3
            dtmp = data_k(k)%vcrd(j,i)-ref(j)
            v = v + k_constant * dtmp * dtmp
          end do
          u_kl(k, l)%v(i) = v
        end do
      end do
    end do

    deallocate(ref)

    write(MsgOut,'(a)') '  ..done'
    write(MsgOut,'(a)') ''

    return

  end subroutine build_u_kl_from_posi_readref

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine solve_mbar(option, u_kl, f_k)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_u_kl),            intent(in)    :: u_kl(:,:)
    type(s_f_k),             allocatable   :: f_k(:)

    ! local variables
    real(wp)                 :: rstep
    integer                  :: iblock, nblock, nstep, istep0, istep1


    write(MsgOut,'(a)') 'Solve_Mbar> '
    write(MsgOut,'(a)') ''

    nblock = option%nblocks
    nstep  = size(u_kl(1,1)%v)

    allocate(f_k(nblock))

    istep0 = 1
    rstep  = 0.0_wp

    do iblock = 1, nblock

      rstep = rstep + (real(nstep, wp) / real(nblock, wp))

      istep1 = nint(rstep)

      call solve_mbar_block(option, u_kl, iblock, &
                            istep0-1, istep1-istep0+1, f_k(iblock)%v)

      istep0 = istep1 + 1

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine solve_mbar

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine solve_mbar_block(option, u_kl, iblock, step0, nstep, f_k)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_u_kl),            intent(in)    :: u_kl(:,:)
    integer,                 intent(in)    :: iblock
    integer,                 intent(in)    :: step0
    integer,                 intent(in)    :: nstep
    real(wp),                allocatable   :: f_k(:)

    ! local variables
    real(wp)                 :: check_convergence, N_max, f_kc
    integer                  :: i, j, ibrella, count_iteration
    integer                  :: nbrella

    real(wp),  allocatable   :: u_kln(:,:,:)
    integer,   allocatable   :: N_k(:)
    real(wp),  allocatable   :: f_k_new(:)
    real(wp),  allocatable   :: log_wi_jn(:,:)
    real(wp),  allocatable   :: W_nk(:,:)
    real(wp),  allocatable   :: g(:), ginv(:)
    real(wp),  allocatable   :: H(:,:), Hinv(:,:)


    nbrella = size(u_kl(:,1))

    ! convert data
    !
    allocate(u_kln(nstep,nbrella,nbrella))
    do i = 1, nbrella
      do j = 1, nbrella
        u_kln(1:nstep,j,i) = u_kl(i,j)%v(step0+1:step0+nstep)
      end do
    end do


    ! N_k: number of data in k-th umbrella window
    !
    allocate(N_k(nbrella))

    do ibrella = 1, nbrella
      N_k(ibrella) = size(u_kln(:,ibrella,1))
    end do

    N_max = maxval(N_k(1:nbrella))
    if (int(N_max) /= nstep) &
      call error_msg('Solve_Mbar> must be N_max eq nstep')


    ! solve the MBAR equation by self-consistent iteration
    !
    allocate(f_k(nbrella), f_k_new(nbrella), log_wi_jn(nstep,nbrella))

    f_k    (1:nbrella) = 0.0_wp
    f_k_new(1:nbrella) = 0.0_wp

    log_wi_jn(1:nstep,1:nbrella)  = 0.0_wp
    do ibrella = 1, nbrella
      log_wi_jn(1:N_k(ibrella), ibrella) = 1.0_wp
    end do

    check_convergence = 1.0e+8_wp

    do count_iteration = 1, option%self_iteration

      do ibrella = 1, nbrella

        call mbar_log_wi_jn(N_k, f_k, u_kln, u_kln(:,ibrella,:), &
                            nbrella, nstep, log_wi_jn)
        f_k_new(ibrella) = - logsumexp2d(log_wi_jn)

      end do

      f_kc = f_k_new(1)
      f_k_new(1:nbrella) = f_k_new(1:nbrella) - f_kc

      check_convergence = maxval(abs(f_k_new(1:nbrella) - f_k(1:nbrella))) / &
           std(f_k_new)

      f_k(1:nbrella) = f_k_new(1:nbrella)

      write(MsgOut,'(a,i5,a,i5,a,e17.9,a,e17.9)') &
           '[block ', iblock, '] ', &
           count_iteration,   'th iteration  delta = ', &
           check_convergence, '  tolerance = ', option%tolerance

    end do


    ! solve the MBAR equation by the Newton-Raphson method
    !
    check_convergence = 10000000.0_wp

    log_wi_jn(1:nstep,1:nbrella) = 0.0_wp
    do ibrella = 1, nbrella
      log_wi_jn(1:N_k(ibrella), ibrella) = 1.0_wp
    end do

    allocate(W_nk(nbrella*nstep,nbrella))
    allocate(g(nbrella-1), ginv(nbrella-1))
    allocate(H(nbrella-1,nbrella-1), Hinv(nbrella-1,nbrella-1))

    W_nk(1:nbrella*nstep,1:nbrella) = 0.0_wp

    count_iteration = option%self_iteration
    do while (check_convergence > option%tolerance .and.            &  
              (count_iteration <                                    &
              option%newton_iteration + option%self_iteration)) 

      f_k_new(1:nbrella) = f_k(1:nbrella)

      do ibrella = 1, nbrella

        call mbar_log_wi_jn(N_k, f_k, u_kln, u_kln(:,ibrella,:), &
                            nbrella, nstep, log_wi_jn)
        W_nk(:,ibrella) = reshape(exp(log_wi_jn+f_k(ibrella)),(/nstep*nbrella/))

      end do

      g(1:nbrella-1)             = 0.0_wp
      H(1:nbrella-1,1:nbrella-1) = 0.0_wp

      do i = 1, nbrella - 1
        g(i)   = N_k(i+1) - N_k(i+1) * sum(W_nk(:,i+1))
        H(i,i) = - sum(N_k(i+1)*W_nk(:,i+1)*(1.0_wp - N_k(i+1)*W_nk(:,i+1)))
        do j = 1, i-1
          H(i,j) = sum((N_k(i+1)*W_nk(:,i+1))*(N_k(j+1)*W_nk(:,j+1)))
          H(j,i) = H(i,j)
        end do
      end do

      call inv_matrix(H, Hinv)
      call mul_matrix(Hinv, g, ginv)

      do i = 1, nbrella - 1
        if (count_iteration == option%self_iteration) then
          f_k_new(i+1) = f_k_new(i+1) - 0.1_wp * ginv(i)
        else
          f_k_new(i+1) = f_k_new(i+1) - 1.0_wp * ginv(i)
        end if

      end do

      check_convergence = maxval(abs(f_k_new(1:nbrella) - f_k(1:nbrella))) / &
           std(f_k_new)
      f_k(1:nbrella) = f_k_new(1:nbrella)

      count_iteration = count_iteration + 1
      write(MsgOut,'(a,i5,a,i5,a,e17.9,a,e17.9)') &
           '[block ', iblock, '] ', &
           count_iteration,   'th iteration  delta = ', &
           check_convergence, '  tolerance = ', option%tolerance

      if (count_iteration ==  &
            option%newton_iteration + option%self_iteration) then
        if  (check_convergence > option%tolerance) then
          write(MsgOut,'(a)') 'WARNING: Not converged'
        endif
      endif
    end do


    ! recompute all free energies
    !
    do ibrella = 1, nbrella

      call mbar_log_wi_jn(N_k, f_k, u_kln, u_kln(:,ibrella,:), &
                          nbrella, nstep, log_wi_jn)
      f_k(ibrella) = - logsumexp2d(log_wi_jn)

    end do

    f_kc = f_k(1)
    f_k(1:nbrella) = f_k(1:nbrella) - f_kc


    ! deallocate memory
    !
    deallocate(H, Hinv, g, ginv, W_nk)
    deallocate(log_wi_jn, f_k_new, N_k, u_kln)

    return

  end subroutine solve_mbar_block

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_bin(option, data_k, bin_k)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_data_k),          intent(in)    :: data_k(:)  ! nbrella
    type(s_bin_k),           allocatable   :: bin_k(:)   ! nbrella

    ! local variables
    integer                  :: nbrella, nstep, i, j, nbin_y

    integer,     allocatable :: bin_x(:), bin_y(:)


    write(MsgOut,'(a)') 'Assign_Bin>'

    if (option%num_grids(1) == 0) then

      write(MsgOut,'(a)') 'Grids are not specified. Constructing histogram is skipped.'
      write(MsgOut,'(a)') ''

      return

    else

      write(MsgOut,'(a)') ''
    
    end if

    nbrella = size(data_k)
    allocate(bin_k(nbrella))

    do i = 1, nbrella

      nstep = size(data_k(i)%v(:,1))

      if (option%dimension == 1) then

        call assign_bin_1d(option, 1, data_k(i)%v(:,1), &
                           bin_k(i)%v, bin_k(i)%center_x)

        if ((i == 1) .and. (allocated(bin_k(i)%center_x))) then
          write(MsgOut, '(a,i3)') 'Assign_Bin> centers of grids in dimension ',1
          do j = 1, size(bin_k(1)%center_x)
            write(MsgOut,'(F10.4,$)') bin_k(1)%center_x(j)
          end do
          write(MsgOut, *)
          write(MsgOut, *)
        end if

      else

        call assign_bin_1d(option, 1, data_k(i)%v(:,1), &
                           bin_x, bin_k(i)%center_x)

        call assign_bin_1d(option, 2, data_k(i)%v(:,2), &
                           bin_y, bin_k(i)%center_y)

        if ((i == 1) .and. (allocated(bin_k(i)%center_x))) then
          write(MsgOut, '(a,i3)') 'Assign_Bin> centers of grids in dimension ',1
          do j = 1, size(bin_k(1)%center_x)
            write(MsgOut,'(F10.4,$)') bin_k(1)%center_x(j)
          end do
          write(MsgOut, *)
          write(MsgOut, *)
        end if

        if ((i == 1) .and. (allocated(bin_k(i)%center_y))) then
          write(MsgOut, '(a,i3)') 'Assign_Bin> centers of grids in dimension ',2
          do j = 1, size(bin_k(1)%center_y)
            write(MsgOut,'(F10.4,$)') bin_k(1)%center_y(j)
          end do
          write(MsgOut, *)
          write(MsgOut, *)
        end if

        allocate(bin_k(i)%v(nstep))

        nbin_y = option%num_grids(2)-1

        do j = 1, nstep
          bin_k(i)%v(j) = nbin_y * (bin_x(j)-1) + bin_y(j)
        end do

        deallocate(bin_x, bin_y)

      end if

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine assign_bin

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_pmf(option, u_kl, u_k, bin_k, f_k, pmf_i)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_u_kl),            intent(in)    :: u_kl(:,:)
    type(s_u_kl),            intent(in)    :: u_k(:)
    type(s_bin_k),           intent(in)    :: bin_k(:)
    type(s_f_k),             intent(in)    :: f_k(:)
    type(s_pmf),             allocatable   :: pmf_i(:)

    ! local variables
    real(wp)                 :: rstep, pmf_min
    integer                  :: nblock, iblock, nstep, istep0, istep1
    integer                  :: nbin, ibin, pmf_min_idx

    type(s_pmf), allocatable :: pmf0(:)
    real(wp),    allocatable :: v(:)

    write(MsgOut,'(a)') 'Compute_Pmf>'

    if (option%num_grids(1) == 0) then

      write(MsgOut,'(a)') 'Grids are not specified. PMF evaluation is skipped.'
      write(MsgOut,'(a)') ''

      return

    else

      write(MsgOut,'(a)') ''
    
    end if

    nblock = option%nblocks
    nstep  = size(u_kl(1,1)%v)

    allocate(pmf0(nblock))

    istep0 = 1
    rstep  = 0.0_wp

    do iblock = 1, nblock

      rstep = rstep + (real(nstep, wp) / real(nblock, wp))
      istep1 = nint(rstep)

      call compute_pmf_block(option, u_kl, u_k, bin_k, f_k(iblock)%v, &
                               iblock, istep0-1, istep1-istep0+1, &
                               pmf0(iblock)%v)

      istep0 = istep1 + 1

    end do

    nbin = size(pmf0(1)%v)

    if (nblock > 1) then

      allocate(pmf_i(2))
      allocate(pmf_i(1)%v(nbin))
      allocate(pmf_i(2)%v(nbin))
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
        pmf_i(1)%v(ibin) = sum(v(1:nblock))/real(nblock,wp)
        pmf_i(2)%v(ibin) = 2.0_wp * std(v)
      end do

      deallocate(v)

    else

      allocate(pmf_i(1))
      allocate(pmf_i(1)%v(nbin))

      !pmf_min = minval_from_not_nan(pmf0(1)%v)
      !pmf_i(1)%v(1:nbin) = pmf0(1)%v(1:nbin) - pmf_min
      pmf_i(1)%v(1:nbin) = pmf0(1)%v(1:nbin)

    end if

    do iblock = 1, nblock
      deallocate(pmf0(iblock)%v)
    end do
    deallocate(pmf0)

    write(MsgOut,'(a)') ''

    return

  end subroutine compute_pmf

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_weight(option, u_kl, u_k, bin_k, f_k, weight_k)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_u_kl),            intent(in)    :: u_kl(:,:)
    type(s_u_kl),            intent(in)    :: u_k(:)
    type(s_bin_k),           intent(in)    :: bin_k(:)
    type(s_f_k),             intent(in)    :: f_k(:)
    real(wp),                allocatable   :: weight_k(:,:)

    ! local variables
    real(wp)                 :: rstep
    integer                  :: nblock, iblock, nstep, istep0, istep1
    integer                  :: nbin, ibin


    write(MsgOut,'(a)') 'Compute_Weight>'
    write(MsgOut,'(a)') ''

    nblock = option%nblocks
    nstep  = size(u_kl(1,1)%v)

    istep0 = 1
    rstep  = 0.0_wp

    do iblock = 1, nblock

      rstep = rstep + (real(nstep, wp) / real(nblock, wp))
      istep1 = nint(rstep)

      if (iblock == 1) then
        if (allocated(weight_k)) &
          deallocate(weight_k)
        call compute_weight_block(option, u_kl, u_k, bin_k, f_k(iblock)%v, &
                                 iblock, istep0-1, istep1-istep0+1, &
                                 weight_k)
      end if

      istep0 = istep1 + 1

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine compute_weight

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_pmf_block(option, u_kl, u_k, bin_k, f_k, &
                                 iblock, step0, nstep, pmf_i)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_u_kl),            intent(in)    :: u_kl(:,:)
    type(s_u_kl),            intent(in)    :: u_k(:)
    type(s_bin_k),           intent(in)    :: bin_k(:)
    real(wp),                intent(in)    :: f_k(:)
    integer,                 intent(in)    :: iblock
    integer,                 intent(in)    :: step0
    integer,                 intent(in)    :: nstep
    real(wp),                allocatable   :: pmf_i(:)

    ! local variables
    real(wp)                 :: N_max, s
    integer                  :: i, j, k
    integer                  :: ibrella, istep, ibin
    integer                  :: nbrella, nbin, nsum

    real(wp)                 :: ZERO = 0.0_wp
    real(wp)                 :: NaN

    real(wp),  allocatable   :: u_kln(:,:,:)
    real(wp),  allocatable   :: u_kn(:,:)
    integer,   allocatable   :: N_k(:)
    real(wp),  allocatable   :: log_w_kn(:,:)
    real(wp),  allocatable   :: log_w_n(:)
    integer,   allocatable   :: bin_kn(:,:)


    NaN  = 0.0_wp/ZERO

    nbrella = size(u_kl(:,1))

    ! convert data
    !
    allocate(u_kln(nstep,nbrella,nbrella))
    do i = 1, nbrella
      do j = 1, nbrella
        u_kln(1:nstep,j,i) = u_kl(i,j)%v(step0+1:step0+nstep)
      end do
    end do

    allocate(bin_kn(nstep,nbrella))
    do i = 1, nbrella
      bin_kn(1:nstep,i) = bin_k(i)%v(step0+1:step0+nstep)
    end do

    nbin = maxval(bin_kn)


    ! N_k: number of data in k-th umbrella window
    !
    allocate(N_k(nbrella))

    do ibrella = 1, nbrella
      N_k(ibrella) = size(u_kln(:,ibrella,1))
    end do

    N_max = maxval(N_k(1:nbrella))
    if (int(N_max) /= nstep) &
      call error_msg('Solve_Mbar> must be N_max eq n-step')


    ! !
    !
    allocate(u_kn(nstep,nbrella))
    do k = 1, nbrella
      u_kn(1:nstep, k) = u_k(k)%v(step0+1:step0+nstep)
    end do
    
    !if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then
    !  u_kn(1:nstep,1:nbrella) = 0.0_wp
    !else
    !  do k = 1, nbrella
    !    do i = 1, nstep
    !      u_kn(i, k) = u_kln(i, k, k) * option%temperature(k) / option%target_temperature
    !    end do
    !  end do
    !end if

    allocate(log_w_kn(nstep,nbrella))

    log_w_kn(1:nstep,1:nbrella)  = 0.0_wp
    do ibrella = 1, nbrella
      log_w_kn(1:N_k(ibrella), ibrella) = 1.0_wp
    end do


    call mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, &
                        nbrella, nstep, log_w_kn)

    ! calc PMF
    !
    allocate(pmf_i(nbin), log_w_n(nstep*nbrella))

    pmf_i(1:nbin) = 0.0_wp

    do ibin = 1, nbin

      nsum = 0
      do ibrella = 1, nbrella
        do istep = 1, nstep
          if (bin_kn(istep,ibrella) == ibin) then
            nsum = nsum + 1
            log_w_n(nsum) = log_w_kn(istep,ibrella)
          end if
        end do
      end do

      if (nsum > 0) then
        pmf_i(ibin) = - logsumexp1dn(log_w_n, nsum)
      else
        pmf_i(ibin) = NaN
      end if

    end do


    ! deallocate memory
    !
    deallocate(log_w_n, log_w_kn, u_kn, N_k, u_kln, bin_kn)

    return

  end subroutine compute_pmf_block

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_weight_block(option, u_kl, u_k, bin_k, f_k, &
                            iblock, step0, nstep, weight_k)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_u_kl),            intent(in)    :: u_kl(:,:)
    type(s_u_kl),            intent(in)    :: u_k(:)
    type(s_bin_k),           intent(in)    :: bin_k(:)
    real(wp),                intent(in)    :: f_k(:)
    integer,                 intent(in)    :: iblock
    integer,                 intent(in)    :: step0
    integer,                 intent(in)    :: nstep
    real(wp),                allocatable   :: weight_k(:,:)

    ! local variables
    real(wp)                 :: N_max, s
    integer                  :: i, j, k
    integer                  :: ibrella, istep, ibin
    integer                  :: nbrella, nbin, nsum

    real(wp)                 :: ZERO = 0.0_wp
    real(wp)                 :: NaN

    real(wp),  allocatable   :: u_kln(:,:,:)
    real(wp),  allocatable   :: u_kn(:,:)
    integer,   allocatable   :: N_k(:)
    real(wp),  allocatable   :: log_w_kn(:,:)
    real(wp),  allocatable   :: log_w_n(:)


    NaN  = 0.0_wp/ZERO

    nbrella = size(u_kl(:,1))

    ! convert data
    !
    allocate(u_kln(nstep,nbrella,nbrella))
    do i = 1, nbrella
      do j = 1, nbrella
        u_kln(1:nstep,j,i) = u_kl(i,j)%v(step0+1:step0+nstep)
      end do
    end do

    ! N_k: number of data in k-th umbrella window
    !
    allocate(N_k(nbrella))

    do ibrella = 1, nbrella
      N_k(ibrella) = size(u_kln(:,ibrella,1))
    end do

    N_max = maxval(N_k(1:nbrella))
    if (int(N_max) /= nstep) &
      call error_msg('Solve_Mbar> must be N_max eq n-step')


    ! !
    !
    allocate(u_kn(nstep,nbrella))
    do k = 1, nbrella
      u_kn(1:nstep, k) = u_k(k)%v(step0+1:step0+nstep)
    end do
    
    !if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then
    !  u_kn(1:nstep,1:nbrella) = 0.0_wp
    !else
    !  do k = 1, nbrella
    !    do i = 1, nstep
    !      u_kn(i, k) = u_kln(i, k, k) * option%temperature(k) / option%target_temperature
    !    end do
    !  end do
    !end if


    allocate(log_w_kn(nstep,nbrella))

    log_w_kn(1:nstep,1:nbrella)  = 0.0_wp
    do ibrella = 1, nbrella
      log_w_kn(1:N_k(ibrella), ibrella) = 1.0_wp
    end do

    call mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, &
                        nbrella, nstep, log_w_kn)

    allocate(log_w_n(nstep*nbrella))

    nsum = 0
    do ibrella = 1, nbrella
      do istep = 1, nstep
        nsum = nsum + 1
        log_w_n(nsum) = log_w_kn(istep,ibrella)
      end do
    end do


    ! calc Weight
    !

    allocate(weight_k(nstep, nbrella))

    log_w_n = reshape(log_w_kn, (/nstep*nbrella/))

    s = logsumexp1dn(log_w_n, size(log_w_n))

    do ibrella = 1,nbrella
      weight_k(1:nstep, ibrella) = exp(log_w_kn(1:nstep,ibrella)-s)
      !weight_k(1:nstep, ibrella) = log_w_kn(1:nstep,ibrella)-s
    end do


    ! deallocate memory
    !
    deallocate(log_w_n, log_w_kn, u_kn, N_k, u_kln)

    return

  end subroutine compute_weight_block

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_mbar(option, output, bin_k, f_k, pmf, weight_k, time_k)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_output),          intent(in)    :: output
    type(s_bin_k),           intent(in)    :: bin_k
    type(s_f_k),             intent(in)    :: f_k(:)
    type(s_pmf),             intent(in)    :: pmf(:)
    real(wp),                intent(in)    :: weight_k(:,:)
    real(wp),                intent(in)    :: time_k(:,:) 

    ! local variables
    integer                  :: file, i, j, nbrella, nstep, irep_x, irep_y
    integer                  :: nrep_x, nrep_y
    integer                  :: nbin_pmf, nbin_x, nbin_y
    integer                  :: ibin_pmf, ibin_x, ibin_y
    integer                  :: ncol
    character(64)            :: fmt

    real(wp)                 :: ZERO = 0.0_wp
    real(wp)                 :: NaN, KBT
    real(wp)                 :: out_unit

    real(wp),    allocatable :: f_k2d(:,:), pmf1d(:), pmf2d(:,:)


    NaN  = 0.0_wp/ZERO
    KBT  = KB * option%target_temperature
    if (option%out_unit == 'kcal/mol') then
      out_unit = KBT
    else if (option%out_unit == 'NONE') then
      out_unit = 1.0_wp
    end if

    ! output fene file
    !
    if (output%fenefile /= '') then

      call open_file(file, output%fenefile, IOFileOutputNew)

      nbrella  = size(f_k(1)%v)

      if (option%dimension == 1) then

        write(fmt,'(a,i0,a)') '(', option%nblocks,'f25.16)'

        do i = 1, nbrella
          write(file,fmt=fmt) (f_k(j)%v(i)*out_unit,j=1,option%nblocks)
        end do

      else

        if (option%nblocks > 1) &
          call error_msg('Output_MBar> n-block > 1 is not supported in 2D')

        nrep_x = option%rest_nreplica(option%rest_func_no(1, 1))
        nrep_y = option%rest_nreplica(option%rest_func_no(2, 1))

        allocate(f_k2d(nrep_y, nrep_x))

        f_k2d = reshape(f_k(1)%v, (/nrep_y, nrep_x/))

        write(fmt,'(a,i0,a)') '(',nrep_x, 'f25.16)'

        do j = 1, nrep_y
          write(file,fmt=fmt) (f_k2d(j,i)*out_unit,i=1,nrep_x)
        end do

        deallocate(f_k2d)

      end if

      call close_file(file)

    end if


    ! output pmf file
    !
    if ((output%pmffile /= '') .and. (option%num_grids(1) > 0)) then

      call open_file(file, output%pmffile, IOFileOutputNew)

      nbin_pmf = size(pmf(1)%v)

      if (option%dimension == 1) then

        if (option%nblocks > 1) then
          ncol = 2
        else
          ncol = 1
        end if

        write(fmt,'(a,i0,a)') '(', ncol+1,'f25.16)'

        do i = 1, nbin_pmf
          write(file,fmt=fmt) bin_k%center_x(i), (pmf(j)%v(i)*out_unit,j=1,ncol)
        end do

      else

        if (option%nblocks > 1) &
          call error_msg('Output_MBar> n-block is not supported in 2D')

        nbin_x = option%num_grids(1)-1
        nbin_y = option%num_grids(2)-1

        allocate(pmf1d(nbin_y*nbin_x), pmf2d(nbin_y, nbin_x))

        pmf1d(:) = NaN
        pmf1d(1:nbin_pmf) = pmf(1)%v(1:nbin_pmf)
        pmf2d = reshape(pmf1d, (/nbin_y, nbin_x/))

        write(fmt,'(a,i0,a)') '(',nbin_x, 'f25.16)'

        do j = 1, nbin_y
          write(file,fmt=fmt) (pmf2d(j,i)*out_unit,i=1,nbin_x)
        end do

        ! do ibin_x = 1, nbin_x
        !   do ibin_y = 1, nbin_y
        !     write(file,*) bin_k%center_x(ibin_x), bin_k%center_y(ibin_y), &
        !                   pmf2d(ibin_y, ibin_x)
        !   end do
        !   write(file,*)
        ! end do

        deallocate(pmf2d, pmf1d)

      end if

      call close_file(file)

    end if


    ! output weight file
    !
    if (output%weightfile /= '') then

      if (option%nblocks > 1) &
     call error_msg('Output_MBar> # of block must be 1 when given weight file.')

      nstep = size(weight_k(:,1))

      if (option%dimension == 1) then

        if (option%input_type == InputTypeCV .or. option%input_type == InputTypeUS) then
          nrep_x = option%rest_nreplica(option%rest_func_no(1, 1))
        else
          nrep_x = option%num_replicas
        end if

        do irep_x = 1, nrep_x

          call open_file(file, &
               get_replicate_name1(output%weightfile, irep_x), IOFileOutputNew)

          do i = 1, nstep
            write(file,*) int(time_k(i,irep_x)), weight_k(i,irep_x)
          end do

          call close_file(file)

        end do

      else

        nrep_x = option%rest_nreplica(option%rest_func_no(1, 1))
        nrep_y = option%rest_nreplica(option%rest_func_no(2, 1))

        j = 0
        do irep_x = 1, nrep_x
          do irep_y = 1, nrep_y

            j = j + 1

            call open_file(file, &
               get_replicate_name1(output%weightfile, j), &
               IOFileOutputNew)

            do i = 1, nstep
              write(file,*) &
                time_k(i,(irep_y-1)+(irep_x-1)*nrep_x+1), &
                weight_k(i,(irep_y-1)+(irep_x-1)*nrep_x+1)
            end do

            call close_file(file)

          end do
        end do

      end if

    end if

    return

  end subroutine output_mbar

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine exec_fhandle(option, umbrella_center, k_constant, data_k, u_kl, is_periodic, box_size)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: umbrella_center(:)
    real(wp),                intent(in)    :: k_constant(:)
    type(s_data_k),          intent(in)    :: data_k
    logical,                 intent(in)    :: is_periodic(:)
    real(wp),                intent(in)    :: box_size(:)
    type(s_u_kl),            intent(inout) :: u_kl

    ! local variables
    real(wp)                 :: K, v
    real(wp)                 :: KBT, BETA
    integer                  :: nstep, nfunc, ifunc, istep


    !KBT  = KB * option%temperature
    !BETA = 1.0_wp / KBT

    nstep = size(data_k%v(:,1))
    nfunc = size(data_k%v(1,:))

    allocate(u_kl%v(nstep))

    do istep = 1, nstep
      v = 0
      do ifunc = 1, nfunc
        K = k_constant(ifunc)
        KBT  = KB * option%temperature(1)
        BETA = 1.0_wp / KBT
        if (is_periodic(ifunc)) then
          v = v + BETA * K * &
            (periodic(data_k%v(istep,ifunc), umbrella_center(ifunc), box_size(ifunc)) ** 2.0_wp)
        else
          v = v + BETA * K * &
            ((data_k%v(istep,ifunc) - umbrella_center(ifunc))**2.0_wp)
        end if
      end do
      u_kl%v(istep) = v
    end do

    return

  end subroutine exec_fhandle

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
    if (func_idx <= option%rest_func_num(1)) then
      func_no = option%rest_func_no(1, func_idx)
    else
      func_no = option%rest_func_no(2, func_idx-option%rest_func_num(1))
    end if
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
    character(MaxLineLong_CV) :: line


    filename = get_replicate_name1(cvfile, 1)

    call open_file(file, filename, IOFileInput)

    nsteps = 0
    do while(.true.)
      read(file,'(a)',end=10) line
      if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
        nsteps = nsteps + 1
      end if
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

!    call open_trj(file, filename, TrjFormatDCD, TrjTypeCoorBox, IOFileInput)
    call open_trj(file, filename, TrjFormatDCD, TrjTypeCoor, IOFileInput)

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

  subroutine mbar_log_wi_jn(N_k, f_k, u_kln, u_kn, nrep, nstp, log_wi_jn)

    ! formal arguments
    integer,                 intent(in)    :: N_k(:)
    real(wp),                intent(in)    :: f_k(:)
    real(wp),                intent(in)    :: u_kln(:,:,:)
    real(wp),                intent(in)    :: u_kn(:,:)
    integer,                 intent(in)    :: nrep
    integer,                 intent(in)    :: nstp
    real(wp),                intent(inout) :: log_wi_jn(:,:)

    ! local variables
    integer                  :: irep, istp, irep2


    log_wi_jn(1:nstp,1:nrep) = 0.0_wp

    !$omp parallel private(istp, irep2)
    !
    if (.not. allocated(g_temp)) then
      allocate(g_temp (nstp,nrep))
    end if

    !$omp do
    !
    ! for K = 1:K
    do irep = 1, nrep

      if (nstp /= N_k(irep)) &
        call error_msg('MBar_Log_Wi_Jn> ERROR ')

      do irep2 = 1, nrep
        do istp = 1, nstp
          g_temp(istp,irep2) = log(real(N_k(irep2), wp)) + f_k(irep2) &
            - u_kln(istp,irep2,irep) + u_kn(istp, irep)
        end do
      end do

      call logsumexp1d2(g_temp(1:nstp,1:nrep), log_wi_jn(:,irep))
      log_wi_jn(1:nstp,irep) = -log_wi_jn(1:nstp,irep)

    end do
    !
    !$omp end do
    !
    !$omp end parallel


!    !$omp parallel private(irep, istp, irep2)
!    !
!
!    if (.not. allocated(g_rep_N_k)) then
!      allocate(g_rep_N_k  (nstp,nrep))
!      allocate(g_rep_f_k  (nstp,nrep))
!      allocate(g_sqz_u_kln(nstp,nrep))
!      allocate(g_rep_u_kn (nstp,nrep))
!    end if
!
!    !$omp do
!    !
!
!    ! for K = 1:K
!    do irep = 1, nrep
!
!      if (nstp /= N_k(irep)) &
!        call error_msg('MBar_Log_Wi_Jn> ERROR ')
!
!      ! repeat(log(N_k), 1, N_k(k))
!      do istp = 1, nstp
!        g_rep_N_k(istp,1:nrep) = log(real(N_k(1:nrep), wp))
!      end do
!
!      ! repeat(f_k, 1, N_k(k))
!      do istp = 1, nstp
!        g_rep_f_k(istp,1:nrep) = f_k(1:nrep)
!      end do
!
!      ! squeeze(u_kln(k, :, 1:N_k(k)))
!      g_sqz_u_kln(1:nstp,1:nrep) = u_kln(1:nstp,1:nrep,irep)
!
!      ! repeat(u_kn(k, 1:N_k(k)), K, 1)
!      do irep2 = 1, nrep
!        g_rep_u_kn(1:nstp,irep2) = u_kn(1:nstp, irep)
!      end do
!
!      call logsumexp1d2(g_rep_N_k   (1:nstp,1:nrep) + &
!                        g_rep_f_k   (1:nstp,1:nrep) - &
!                        (g_sqz_u_kln(1:nstp,1:nrep) - &
!                         g_rep_u_kn (1:nstp,1:nrep)), &
!                        log_wi_jn(:,irep))
!      log_wi_jn(1:nstp,irep) = -log_wi_jn(1:nstp,irep)
!
!    end do
!
!    !
!    !$omp end do
!    !
!    !$omp end parallel

    return

  end subroutine mbar_log_wi_jn

  !======1=========2=========3=========4=========5=========6=========7=========8

  function logsumexp2d(x)

    ! return value
    real(wp)                 :: logsumexp2d

    ! formal arguments
    real(wp),                intent(in)    :: x(:,:)

    ! local variables
    real(wp)                 :: max_x, exp_x
    integer                  :: i, j, ni, nj


    ni = size(x(1,:))
    nj = size(x(:,1))

    max_x = maxval(x(1:nj,1:ni))

    exp_x = 0.0_wp

    do i = 1, ni
      do j = 1, nj
        exp_x = exp_x + exp(x(j,i) - max_x)
      end do
    end do

    logsumexp2d = log(exp_x) + max_x

    return

  end function logsumexp2d

  !======1=========2=========3=========4=========5=========6=========7=========8

  function logsumexp1dn(x,n)

    ! return value
    real(wp)                 :: logsumexp1dn

    ! formal arguments
    real(wp),                intent(in)    :: x(:)
    integer,                 intent(in)    :: n

    ! local variables
    real(wp)                 :: max_x, exp_x
    integer                  :: i


    max_x = maxval(x(1:n))

    exp_x = 0.0_wp

    do i = 1, n
      exp_x = exp_x + exp(x(i) - max_x)
    end do

    logsumexp1dn = log(exp_x) + max_x

    return

  end function logsumexp1dn

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine logsumexp1d2(x, s)

    ! formal arguments
    real(wp),                intent(in)    :: x(:,:)
    real(wp),                intent(inout) :: s(:)

    ! local variables
    integer                  :: icol, ncol, nrow


    nrow = size(x(1,:))
    ncol = size(x(:,1))

    do icol = 1, ncol
      s(icol) = logsumexp1dn(x(icol,:), nrow)
    end do

    return

  end subroutine logsumexp1d2

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

  !======1=========2=========3=========4=========5=========6=========7=========8

  function minval_from_not_nan(x)

    ! return value
    real(wp)                 :: minval_from_not_nan

    ! formal arguments
    real(wp),                intent(in)    :: x(:)

    ! local variables
    integer                  :: i, j, n
    real(wp)                 :: ZERO = 0.0_wp
    real(wp)                 :: NaN


    NaN = 0.0_wp/ZERO
    n = size(x(:))

    i = 1
    do while(x(i) == NaN)
      i = i + 1
      if (i > n) then
        minval_from_not_nan = NaN
        return
      end if
    end do
    minval_from_not_nan = x(i)

    do j = i+1, n
      if (x(j) .ne. NaN) then
        if (x(j) < minval_from_not_nan) then
          minval_from_not_nan = x(j)
        end if
      end if
    end do

    return

  end function minval_from_not_nan

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_bin_1d(option, idim, data_k, bin_k, center)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    integer,                 intent(in)    :: idim
    real(wp),                intent(in)    :: data_k(:)
    integer,                 allocatable   :: bin_k(:)
    real(wp),                allocatable   :: center(:)

    ! local variables
    real(wp)                 :: grid_min, grid_max, grid_range, k
    integer                  :: num_bins, num_grid, num_step
    integer                  :: i, istp, ibin

    real(wp),  allocatable   :: grid(:)


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

    ! bin
    num_step = size(data_k)

    allocate(bin_k(num_step))

    bin_k(1:num_step) = 0

    do istp = 1, num_step

      k = data_k(istp)

      do ibin = 1, num_bins
        if (grid(ibin) <= k .and. k < grid(ibin+1)) then
          bin_k(istp) = ibin
          exit
        end if
      end do

    end do

    deallocate(grid)

    return

  end subroutine assign_bin_1d

end module ma_analyze_mod
