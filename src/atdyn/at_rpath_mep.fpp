!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_rpath_mep_mod
!> @brief   Minimum energy path search using replicas
!! @authors Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_rpath_mep_mod

  use at_energy_mod
  use at_md_leapfrog_mod
  use at_md_vverlet_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_dynvars_mod
  use at_ensemble_str_mod
  use at_rpath_str_mod
  use at_constraints_str_mod
  use at_constraints_mod
  use at_boundary_str_mod
  use at_boundary_mod
  use at_pairlist_str_mod
  use at_pairlist_mod
  use at_enefunc_str_mod
  use at_minimize_str_mod
  use at_minimize_mod
  use at_output_str_mod
  use at_output_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_rstmep_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! local variables
  logical,                private :: vervose = .true.

  ! subroutines
  public  :: setup_rpath_mep
  public  :: setup_mepatoms
  public  :: setup_mepatoms_qmmm
  public  :: run_rpath_mep
  private :: run_rpath_string
  private :: path_update_string
  private :: run_rpath_neb_glbfgs
  private :: calc_neb_force
  private :: minimize_micro
  private :: energy_and_force
  private :: check_convergence
  private :: check_convergence_neb
  private :: print_images
  private :: trans_mass_weight_coord
  private :: backtrans_mass_weight_coord
  private :: calc_pathlength

  private :: gradient_correction

  private :: printout_mep_replica        
  private :: energy_and_force_replica    
  private :: calc_pathlength_replica     
  private :: run_rpath_string_replica    
  private :: path_update_string_replica  
  private :: run_rpath_neb_glbfgs_replica  
  private :: trans_mass_weight_coord_replica
  private :: backtrans_mass_weight_coord_replica
  
  public  :: setup_rpath_mepmd
  public  :: run_rpath_mepmd
  private :: run_rpath_string_mepmd
  private :: ene_force_mepmd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath_mep
  !> @brief        setup RPATH
  !! @authors      YA, KY, SI
  !! @param[in]    mepatm_select_index : index of MEP atoms
  !! @param[in]    sel_info   : selector input information
  !! @param[in]    rst        : restart information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    qmmm       : qmmm information
  !! @param[inout] minimize   : minimize information
  !! @param[inout] dynvars    : dynamic variables
  !! @param[inout] rpath      : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath_mep(mepatm_select_index, sel_info, rst, &
                             molecule, qmmm, minimize, dynvars, rpath, rstmep)

    ! formal arguments
    character(256),           intent(in)    :: mepatm_select_index
    type(s_sel_info),         intent(in)    :: sel_info
    type(s_rst),              intent(in)    :: rst
    type(s_molecule),         intent(in)    :: molecule
    type(s_qmmm),             intent(inout) :: qmmm
    type(s_minimize),         intent(inout) :: minimize
    type(s_dynvars),          intent(inout) :: dynvars
    type(s_rpath),            intent(inout) :: rpath
    type(s_rstmep), optional, intent(in)    :: rstmep

    ! local variables
    integer   :: i, ii
    integer   :: replicaid
    integer   :: atomid, iatom
    integer   :: errorcode, ierr
    logical   :: err


    if (nrep_per_proc > 1) replicaid = my_replica_no

    ! define atoms
    !
    call setup_mepatoms(mepatm_select_index, sel_info, molecule, rpath)

    if (qmmm%do_qmmm) then
      call setup_mepatoms_qmmm(molecule, qmmm, rpath)

      if (minimize%macro) then
        ! Since the QM region is included in the MEP region,
        ! opt_micro is set to true when miminize%macro is true.
        rpath%opt_micro = .true.
        qmmm%qm_get_esp = .true.
        minimize%nsteps_micro = minimize%nsteps

      else
        rpath%opt_micro = .false.

      end if

    else
      rpath%opt_micro        = .false.
      rpath%lbfgs_bnd_qmonly = .false.

    end if
    call add_fixatm(rpath%mep_natoms, rpath%mepatom_id, molecule, minimize)

    ! allocate memory
    rpath%dimension = rpath%mep_natoms * 3
    if (qmmm%do_qmmm) then
      call alloc_rpath_mep(rpath, rpath%dimension, rpath%nreplica, &
                           minimize%num_optatoms_micro, nrep_per_proc, &
                           qmmm%qm_natoms)
    else
      call alloc_rpath_mep(rpath, rpath%dimension, rpath%nreplica, &
                           minimize%num_optatoms_micro, nrep_per_proc)
    end if

    ! set initial images
    replicaid = my_country_no + 1

    if (nrep_per_proc > 1) replicaid = my_replica_no 
    if (rpath%first_iter) then
      ii = 0
      do i = 1, rpath%mep_natoms
        iatom = rpath%mepatom_id(i)
        rpath%mep_coord(ii+1:ii+3,replicaid) = dynvars%coord(1:3,iatom)
        ii = ii + 3
      end do
    end if

    !dbg write(MsgOut,'("ID",i4,"> ",100f8.4)') my_country_no, qmmm%qm_charge
    !dbg call mpi_barrier(mpi_comm_world, ierr)
    !dbg call mpi_abort(mpi_comm_world, errorcode, ierr)

    ! write the summary of setup
    !
    if (main_rank .and. vervose) then
      write(MsgOut,'(A)') 'Setup_Rpath_MEP> Rpath information'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A20,I10)') '  dimension       = ', rpath%dimension
      write(MsgOut,'(A)') ''

      ! Punch out atom info.
      !
      write(MsgOut,'(a)') "  Atoms involved in MEP search"
      do i = 1, rpath%mep_natoms
        atomid   = rpath%mepatom_id(i)
        write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
            atomid, &
            molecule%segment_name(atomid), &
            molecule%residue_no(atomid),   &
            molecule%residue_name(atomid), &
            molecule%atom_name(atomid),    &
            molecule%atom_cls_name(atomid)
      end do
      write(MsgOut,'(a,i0)') "  number of atoms in MEP search = ", rpath%mep_natoms
      write(MsgOut, '(a)') ' '
      vervose = .false.
    end if

    ! supress output of minimize
    minimize%eneout_short = .true.
    if (rpath%eneout_period == 0) then
      minimize%eneout_period = 0
    end if

    if (minimize%crdout_period > 0) then
      minimize%crdout_period = 0
      if (main_rank) &
        write(MsgOut,'(A,/)') &
          ' (NOTICE) crdout_period of [MINIMIZE] is reset to zero.'
    end if

    if (minimize%rstout_period > 0) then
      minimize%rstout_period = 0
      if (main_rank) &
        write(MsgOut,'(A,/)') &
          ' (NOTICE) rstout_period of [MINIMIZE] is reset to zero.'
    end if

    rpath%eneout = .false.
    rpath%crdout = .false.
    rpath%rstout = .false.

    return

  end subroutine setup_rpath_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath_mepmd
  !> @brief        setup RPATH
  !! @authors      KY
  !! @param[in]    mepatm_select_index : index of MEP atoms
  !! @param[in]    sel_info    : selector input information
  !! @param[in]    rst         : restart information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    constraints : information of constraints
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath_mepmd(mepatm_select_index, sel_info, rst, dynamics,   &
                               constraints, molecule, enefunc, dynvars, rpath, &
                               rstmep)

    ! formal arguments
    character(256),           intent(in)    :: mepatm_select_index
    type(s_sel_info),         intent(in)    :: sel_info
    type(s_rst),              intent(in)    :: rst
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_constraints),      intent(in)    :: constraints
    type(s_molecule),         intent(in)    :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),          intent(inout) :: dynvars
    type(s_rpath),            intent(inout) :: rpath
    type(s_rstmep), optional, intent(in)    :: rstmep

    ! local variables
    integer   :: i, ii
    integer   :: replicaid
    integer   :: atomid, iatom
    integer   :: errorcode, ierr
    logical   :: err
    type(s_qmmm), pointer  :: qmmm

    ! Pointers
    !
    qmmm => enefunc%qmmm

    if (nrep_per_proc > 1) replicaid = my_replica_no

    ! setup reparameterization period
    !
    if (rpath%rpath_period > 0) then
      if (mod(dynamics%nsteps, rpath%rpath_period) /= 0) then
        call error_msg('Setup_Rpath> mod(nsteps,rpath_period) must be zero')
      else
        rpath%equilibration_only = .false.
        rpath%ncycle = dynamics%nsteps/rpath%rpath_period
      end if
    else if (rpath%rpath_period == 0) then
      rpath%equilibration_only = .true.
      rpath%rpath_period = dynamics%nsteps
      rpath%ncycle = 1
    else
      call error_msg('Setup_Rpath> rpath_period should be non-negative integer')
    end if

    ! define atoms
    !
    call setup_mepatoms(mepatm_select_index, sel_info, molecule, rpath)

    if (qmmm%do_qmmm) then
      call setup_mepatoms_qmmm(molecule, qmmm, rpath)
      qmmm%qm_get_esp = .true.

    else
      rpath%lbfgs_bnd_qmonly = .false.

    end if

    ! allocate memory
    rpath%dimension = rpath%mep_natoms * 3
    if (qmmm%do_qmmm) then
      call alloc_rpath_mepmd(rpath, rpath%dimension, rpath%nreplica, &
                           qmmm%qm_natoms)
    else
      call alloc_rpath_mepmd(rpath, rpath%dimension, rpath%nreplica)
    end if

    ! setup qm force
    allocate(rpath%mepmd_qmforce(3,rpath%mep_natoms))
    rpath%mepmd_qmforce = 0.0_wp

    ! setup average force
    enefunc%rpath_sum_avforce = .false.
    allocate(enefunc%mepmd_avforce(3,rpath%mep_natoms))
    enefunc%mepmd_avforce = 0.0_wp

    ! set initial images
    replicaid = my_country_no + 1

    if (nrep_per_proc > 1) replicaid = my_replica_no 
    if (rpath%first_iter) then
      ii = 0
      do i = 1, rpath%mep_natoms
        iatom = rpath%mepatom_id(i)
        rpath%mep_coord(ii+1:ii+3,replicaid) = dynvars%coord(1:3,iatom)
        ii = ii + 3
      end do
    end if

    !dbg write(MsgOut,'("ID",i4,"> ",100f8.4)') my_country_no, qmmm%qm_charge
    !dbg call mpi_barrier(mpi_comm_world, ierr)
    !dbg call mpi_abort(mpi_comm_world, errorcode, ierr)

    ! write the summary of setup
    !
    if (main_rank .and. vervose) then
      write(MsgOut,'(A)') 'Setup_Rpath_MEP> Rpath information'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A20,I10)') '  dimension       = ', rpath%dimension
      write(MsgOut,'(A)') ''

      ! Punch out atom info.
      !
      write(MsgOut,'(a)') "  Atoms involved in MEP search"
      do i = 1, rpath%mep_natoms
        atomid   = rpath%mepatom_id(i)
        write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
            atomid, &
            molecule%segment_name(atomid), &
            molecule%residue_no(atomid),   &
            molecule%residue_name(atomid), &
            molecule%atom_name(atomid),    &
            molecule%atom_cls_name(atomid)
      end do
      write(MsgOut,'(a,i0)') "  number of atoms in MEP search = ", rpath%mep_natoms
      write(MsgOut, '(a)') ' '
      vervose = .false.
    end if

    ! Check if MEP atoms are fixed in MD
    !
    ierr = 0
    do i = 1, rpath%mep_natoms
      atomid = rpath%mepatom_id(i)
      if (.not. constraints%fixatm(atomid)) then
        ierr = 1
        if (main_rank) then
          write(MsgOut,'(A,i0,A)') 'Setup_rpath_mepmd> Error: MEP AtomID ',atomid,' is not fixed!'
        end if
      end if
    end do
    if (ierr == 1) call error_msg('Setup_rpath_mepmd> Some MEP atoms are not fixed in [CONSTRAINTS]')

    rpath%eneout = .false.
    rpath%crdout = .false.
    rpath%rstout = .false.

    return

  end subroutine setup_rpath_mepmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mepatoms
  !> @brief        define MEP atoms
  !! @authors      YA
  !! @param[in]    mepatm_select_index : index of MEP atoms
  !! @param[in]    sel_info   : selector input information
  !! @param[in]    molecule   : molecule information
  !! @param[out]   rpath      : rpath parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mepatoms(mepatm_select_index, sel_info, molecule, rpath)

    ! formal arguments
    character(256),      intent(in)    :: mepatm_select_index
    type(s_sel_info),    intent(in)    :: sel_info
    type(s_molecule),    intent(in)    :: molecule
    type(s_rpath),       intent(inout) :: rpath

    ! local variables
    integer                :: igroup, ngroup, natom, i, j, offset, temp
    integer, allocatable   :: group_list(:)
    type(s_selatoms), allocatable :: selatoms(:)


    ! Number of atoms in MEP
    !
    ngroup = split_num(trim(mepatm_select_index))

    allocate(group_list(ngroup))
    call split(ngroup, ngroup, mepatm_select_index, group_list)

    allocate(selatoms(ngroup))

    natom = 0
    do i = 1, ngroup
      igroup = group_list(i)
      call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
      natom = natom + size(selatoms(i)%idx)
    end do

    rpath%mep_natoms = natom

    ! List of atoms
    !
    if (allocated(rpath%mepatom_id)) deallocate(rpath%mepatom_id)    
    allocate(rpath%mepatom_id(rpath%mep_natoms))

    offset = 0
    do i = 1, ngroup
      igroup = group_list(i)
      natom = size(selatoms(i)%idx)
      rpath%mepatom_id(offset+1:offset+natom) = selatoms(i)%idx(1:natom)
      offset = offset + natom
    end do

    deallocate(selatoms)
    deallocate(group_list)

    ! sort atom indices in ascending order
    !
    do i = rpath%mep_natoms, 2, -1
      do j = 1, i - 1
        if (rpath%mepatom_id(j) > rpath%mepatom_id(j+1)) then
          temp = rpath%mepatom_id(j)
          rpath%mepatom_id(j)   = rpath%mepatom_id(j+1)
          rpath%mepatom_id(j+1) = temp
        end if
      end do
    end do

    return

  end subroutine setup_mepatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mepatoms_qmmm
  !> @brief        define MEP atoms for QM/MM jobs
  !! @authors      KY
  !! @param[in]    molecule   : molecule information
  !! @param[in]    qmmm       : QMMM information
  !! @param[out]   rpath      : rpath parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mepatoms_qmmm(molecule, qmmm, rpath)

    ! formal arguments
    type(s_molecule),  intent(in)    :: molecule
    type(s_qmmm),      intent(in)    :: qmmm
    type(s_rpath),     intent(inout) :: rpath

    ! local variables
    integer              :: i, j, k, id, nad, ntmp
    logical              :: error
    logical, allocatable :: found(:)
    integer, allocatable :: add(:), tmp(:)


    allocate(found(qmmm%qm_natoms))
    error = .false.
    do i = 1, qmmm%qm_natoms
      id = qmmm%qmatom_id(i)
      found(i) = .false.
      do j = 1, rpath%mep_natoms
        if (id == rpath%mepatom_id(j)) then
          found(i) = .true.
          exit
        end if
      end do
      if (.not. found(i)) error = .true.
    end do

    if (error) then
      if (main_rank) then
        write(MsgOut,'(a)') "Setup_MEP> Fatal error while setting up MEP."
        write(MsgOut,'(a)') "Setup_MEP> These QM atoms are not in MEP:"
        do i = 1, qmmm%qm_natoms
          if (.not. found(i)) then
            id = qmmm%qmatom_id(i)
            write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
                id, &
                molecule%segment_name(id), &
                molecule%residue_no(id),   &
                molecule%residue_name(id), &
                molecule%atom_name(id),    &
                molecule%atom_cls_name(id)
          end if
        end do
        write(MsgOut,'(a)') "Setup_MEP> check mepatm_select_index."
      end if
      call error_msg("Setup_MEP> Abort with error.")
    end if

    deallocate(found)

    if (qmmm%num_qmmmbonds > 0) then

      ! search boundary MM atoms in MEP atoms
      allocate(found(qmmm%num_qmmmbonds), add(qmmm%num_qmmmbonds))
      nad = 0
      do i = 1, qmmm%num_qmmmbonds
        id = qmmm%qmmmbond_list(2,i)
        found(i) = .false.
        do j = 1, rpath%mep_natoms
          if (id == rpath%mepatom_id(j)) then
            found(i) = .true.
            exit
          end if
        end do
        if (.not. found(i)) then
          nad = nad + 1
          add(nad) = id
        end if
      end do

      ! add boundary MM atoms to MEP atoms
      if (nad > 0) then
        ntmp = rpath%mep_natoms
        allocate(tmp(rpath%mep_natoms))
        tmp = rpath%mepatom_id

        deallocate(rpath%mepatom_id)
        rpath%mep_natoms = rpath%mep_natoms + nad
        allocate(rpath%mepatom_id(rpath%mep_natoms))

        j = 1
        k = 1
        do i = 1, rpath%mep_natoms
          if (add(j) < tmp(k)) then
            rpath%mepatom_id(i) = add(j)
            j = j + 1
            if (j > nad) then
              rpath%mepatom_id(i+1:) = tmp(k:)
              exit
            end if

          else if (tmp(k)  < add(j)) then
            rpath%mepatom_id(i) = tmp(k)
            k = k + 1
            if (k > ntmp) then
              rpath%mepatom_id(i+1:) = add(j:)
              exit
            end if 
          end if 
       
        end do

        deallocate(tmp)
      end if

      deallocate(found, add)
    end if

    return

  end subroutine setup_mepatoms_qmmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    printout_mep_replica
  !> @brief        print out the data when nrepica > nproc
  !! @authors      SI
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine printout_mep_replica(output, molecule, enefunc, dynvars, &
                                  dynamics, boundary, rpath, icycle)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    integer,                  intent(in)    :: icycle

    ! local variables
    integer                                 :: save_ene, save_crd


    output%rpathout     = .false.
    rpath%rstout_period = 0
    save_ene            = rpath%eneout_period
    save_crd            = rpath%crdout_period

    if (rpath%method == MEPmethod_NEB) then
      dynvars%step = rpath%neb_cycle - 1
    else 
      dynvars%step = icycle - 1
    end if 

    if (mod(dynvars%step,rpath%eneout_period) == 0) then
      rpath%eneout = .false.
      rpath%eneout_period = 0
    else
      rpath%eneout = .true.
    end if  

    if (mod(dynvars%step,rpath%crdout_period) == 0) then
      rpath%crdout = .false.
      rpath%crdout_period = 0
    else
      rpath%crdout = .true.
    end if  

    ! output variables
    !
    call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                      boundary, rpath)
    rpath%eneout_period = save_ene
    rpath%crdout_period = save_crd  

    return

  end subroutine printout_mep_replica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_mep
  !> @brief        control string method
  !! @authors      TM, YK, YA, KY, SI
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_mep(output, molecule, enefunc, dynvars, minimize, &
                           dynamics, pairlist, boundary, rpath, conv, icycle, &
                           ireplica)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath

    integer, optional,        intent(in)    :: icycle
    integer, optional,        intent(in)    :: ireplica
    logical,                  intent(out)   :: conv   

    ! local variables
    integer                :: niter
    integer                :: replicaid

    real(wp)    , pointer  :: coord(:,:)


    ! Pointers
    !
    coord => dynvars%coord

    ! My replica ID
    !
    replicaid = my_country_no + 1

    if (nrep_per_proc > 1) then
      replicaid = my_replica_no
      niter = icycle
    end if

    ! Open output files
    !
    call open_output(output)
    DynvarsOut = output%logunit
    
    ! Print out mep data if iteration is converged or reach limit
    ! Only work when nreplica > nproc
    !
    if (nrep_per_proc > 1) then
      if (conv .or. ((.not. rpath%method == MEPmethod_NEB) .and. &
      	  icycle == rpath%ncycle+1) .or. &
          (rpath%neb_cycle == rpath%ncycle+1)) then
        call printout_mep_replica(output, molecule, enefunc, dynvars, &
                                  dynamics, boundary, rpath, icycle)
        rpath%mep_exit = .true.
        call close_output(output)
        return
      end if
    end if

    ! update boundary and pairlist
    !
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(enefunc%cg_pairlistdist_ele, &
          enefunc%cg_pairlistdist_126,                     &
          enefunc%cg_pairlistdist_PWMcos,                  &
          enefunc%cg_pairlistdist_DNAbp,                   &
          enefunc%cg_pairlistdist_exv,                     &
          boundary)
    else
      call update_boundary(enefunc%table%table,               &
          enefunc%pairlistdist,                               &
          boundary)
    end if
    if (real_calc) then
      pairlist%allocation = .true.
      call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
          dynvars%coord_pbc, pairlist)

    end if

    ! start minimum energy path search
    !
    if (rpath%method == MEPmethod_String) then
      if (nrep_per_proc == 1) then
        call run_rpath_string(output, molecule, enefunc, dynvars, &
                              minimize, dynamics, pairlist, boundary, &
                              rpath, conv, niter)
      else
        call run_rpath_string_replica(output, molecule, enefunc, dynvars, &
                                      minimize, dynamics, pairlist, boundary, &
                                      rpath, conv, niter, ireplica)
      end if
    else if (rpath%method == MEPmethod_NEB) then
      if (nrep_per_proc == 1) then
        call run_rpath_neb_glbfgs(output, molecule, enefunc, dynvars, &
                                  minimize, dynamics, pairlist, boundary, &
                                  rpath, conv, niter)
      else
        call run_rpath_neb_glbfgs_replica(output, molecule, enefunc, &
                                          dynvars, minimize, dynamics, &
                                          pairlist, boundary, rpath, conv, &
                                          niter, ireplica)
      end if
    end if

    if (conv) then
      if (main_rank) write(MsgOut, '(A,I0,A)') &
             'Convergence achieved in ', niter, ' iterations'
    else 
      if (nrep_per_proc == 1) then
        if (main_rank) write(MsgOut, '(A,I0,A)') &
               'No convergence in ', rpath%ncycle, ' iterations'
      else if (rpath%method == MEPmethod_NEB) then
        if (rpath%neb_cycle == rpath%ncycle+1) then
          if (main_rank) write(MsgOut, '(A,I0,A)') &
           'No convergence in ', rpath%ncycle, ' iterations'
        end if           
      else
        if (icycle == rpath%ncycle) then
          if (main_rank) write(MsgOut, '(A,I0,A)') &
           'No convergence in ', rpath%ncycle, ' iterations'
        end if   
      end if 
    end if 

    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_rpath_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_mepmd
  !> @brief        control string method
  !! @authors      KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_mepmd(output, molecule, enefunc, dynvars, dynamics,   &
                           pairlist, boundary, constraints, ensemble, rpath, &
                           conv, icycle, ireplica)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_rpath),            intent(inout) :: rpath

    integer, optional,        intent(in)    :: icycle
    integer, optional,        intent(in)    :: ireplica
    logical,                  intent(out)   :: conv   

    ! local variables
    integer                :: niter
    integer                :: replicaid

    real(wp)    , pointer  :: coord(:,:)

    ! Pointers
    !
    coord => dynvars%coord

    ! My replica ID
    !
    replicaid = my_country_no + 1

    ! Open output files
    !
    call open_output(output)
    DynvarsOut = output%logunit
    
    ! update boundary and pairlist
    !
    call update_boundary(enefunc%table%table, enefunc%pairlistdist, boundary)
    if (real_calc) then
      pairlist%allocation = .true.
      call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                           dynvars%coord_pbc, pairlist)

    end if

    ! start minimum energy path search
    !
    if (rpath%method == MEPmethod_String) then
      call run_rpath_string_mepmd(output, molecule, enefunc, dynvars,       &
                    dynamics, pairlist, boundary, constraints, ensemble, &
                    rpath, conv)
    !else if (rpath%method == MEPmethod_NEB) then
    !  call run_rpath_neb_glbfgs(output, molecule, enefunc, dynvars, &
    !                            minimize, dynamics, pairlist, boundary, &
    !                            rpath, conv, niter)
    end if

    if (.not. rpath%equilibration_only) then
      if (conv) then
        if (main_rank) write(MsgOut, '(A,I0,A)') &
               'Convergence achieved in ', niter, ' iterations'
      else 
        if (nrep_per_proc == 1) then
          if (main_rank) write(MsgOut, '(A,I0,A)') &
                 'No convergence in ', rpath%ncycle, ' iterations'
        else if (rpath%method == MEPmethod_NEB) then
          if (rpath%neb_cycle == rpath%ncycle+1) then
            if (main_rank) write(MsgOut, '(A,I0,A)') &
             'No convergence in ', rpath%ncycle, ' iterations'
          end if           
        else
          if (icycle == rpath%ncycle) then
            if (main_rank) write(MsgOut, '(A,I0,A)') &
             'No convergence in ', rpath%ncycle, ' iterations'
          end if   
        end if 
      end if 
    end if 

    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_rpath_mepmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_string
  !> @brief        string method
  !! @authors      TM, YK, YA, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_string(output, molecule, enefunc, dynvars, minimize, &
                              dynamics, pairlist, boundary, rpath, conv, niter)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    logical,                  intent(out)   :: conv
    integer,                  intent(out)   :: niter

    ! local variables
    integer                :: n
    integer                :: min_ene_period, save_step


    ! save eneout_period for minimize
    !
    min_ene_period = minimize%eneout_period

    ! initialize
    !
    rpath%energy_prev     = 0.0_wp
    rpath%mep_force       = 0.0_wp
    rpath%mep_length      = 0.0_wp
    rpath%pathlength_prev = 0.0_wp

    ! run rpath
    !
    dynvars%step = 0
    if (main_rank) then
       write(MsgOut, *) "--- Start string iteration ---"
       write(MsgOut, '(/,"Iter. ",i5)') 0
    end if

    ! compute the initial energy and force
    !
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        .false.,           &
                        dynvars%coord,     &
                        dynvars%trans,     &
                        dynvars%coord_pbc, &
                        dynvars%energy,    &
                        dynvars%temporary, &
                        dynvars%force,     &
                        dynvars%force_omp, &
                        dynvars%virial,    &
                        dynvars%virial_extern)
    call energy_and_force(dynvars, rpath)

    ! calc initial path length
    !
    if (rpath%massWeightCoord) then
      call trans_mass_weight_coord(molecule, rpath)
      call calc_pathlength(rpath)
      call backtrans_mass_weight_coord(molecule, rpath)

    else
      call calc_pathlength(rpath)

    end if

    ! output variables
    !
    call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                      boundary, rpath)
    call print_images(rpath)

    rpath%first_iter = .false.

    ! main loop
    !
    do n = 1, rpath%ncycle
      niter        = n
      dynvars%step = n
      ! output settings
      !
      if (main_rank) write(MsgOut, '(/,"Iter. ",i5)') n
      if (rpath%eneout_period > 0) then
        if (mod(n, rpath%eneout_period) == 0) then
          if (replica_main_rank) write(output%logunit, '("Iter. ",i5,/)') n
          minimize%eneout_period = min_ene_period

        else
          minimize%eneout_period = 0

        end if
      end if

      ! optimize surrounding atoms
      !
      if (rpath%mep_partial_opt) then
        save_step = dynvars%step

        if (rpath%opt_micro) then
          call minimize_micro (output, molecule, enefunc, dynvars, &
                             minimize, pairlist, boundary, rpath, n)

        else
          if (minimize%method == MinimizeMethodSD) then
            call steepest_descent (output, molecule, enefunc, dynvars, &
                                   minimize, pairlist, boundary)

          else if (minimize%method == MinimizeMethodLBFGS) then
            call minimize_lbfgs (output, molecule, enefunc, dynvars,   &
                                 minimize, pairlist, boundary)
          end if

        end if
        dynvars%step = save_step

      end if

      ! set the energy and force
      !
      call energy_and_force(dynvars, rpath)

      ! check convergence
      !
      call check_convergence(rpath, conv)
      if (conv)  exit

      ! output variables
      !
      call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                        boundary, rpath)
      if (.not. (n == rpath%ncycle)) call print_images(rpath)

      ! path update
      !
      call path_update_string(molecule, rpath, dynvars)

    end do

    ! output final status
    !
    rpath%eneout = .true.
    rpath%crdout = .true.
    rpath%rstout = .true.
    call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                      boundary, rpath)
    rpath%eneout = .false.
    rpath%crdout = .false.
    rpath%rstout = .false.

    call print_images(rpath)

    return

  end subroutine run_rpath_string

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_string_mepmd
  !> @brief        string method
  !! @authors      KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_string_mepmd(output, molecule, enefunc, dynvars, &
                   dynamics, pairlist, boundary, constraints, ensemble, &
                   rpath, conv)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_rpath),    target, intent(inout) :: rpath
    logical,                    intent(out) :: conv

    ! local variables
    integer                :: n
    integer                :: save_step


    ! initialize
    !
    rpath%energy_prev     = 0.0_wp
    rpath%mep_force       = 0.0_wp
    rpath%mep_length      = 0.0_wp
    rpath%pathlength_prev = 0.0_wp
    rpath%mepmd_qmforce   = 0.0_wp

    if(enefunc%qmmm%do_qmmm) then
      enefunc%qmmm%qm_classical = .true.
    end if
    enefunc%mep_natoms = rpath%mep_natoms
    enefunc%mepatom_id => rpath%mepatom_id

    ! run rpath
    !
    dynvars%step = 0
    if (main_rank) then
      if (.not. rpath%equilibration_only) then
        write(MsgOut, *) "--- Start string iteration ---"
      else
        write(MsgOut, *) "--- Equilibration run ---"
      end if
    end if

    ! calc initial path length
    !
    if (rpath%massWeightCoord) then
      call trans_mass_weight_coord(molecule, rpath)
      call calc_pathlength(rpath)
      call backtrans_mass_weight_coord(molecule, rpath)

    else
      call calc_pathlength(rpath)

    end if

    ! output variables
    !
    call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                      boundary, rpath)

    ! main loop
    !
    do n = 1, rpath%ncycle

      if (main_rank .and. .not. rpath%equilibration_only) &
        write(MsgOut, '(/,"Iter. ",i5)') n

      if(enefunc%qmmm%do_qmmm) &
          call update_qmforce(molecule, enefunc, pairlist, boundary, &
                              dynvars, constraints, rpath)

      ! MD loop
      !
      dynamics%istart_step  = (n-1)*rpath%rpath_period + 1
      dynamics%iend_step    =  n   *rpath%rpath_period
      dynamics%initial_time = dynvars%time

      enefunc%rpath_sum_avforce = .true.
      enefunc%mepmd_avforce     = 0.0_wp
      if (dynamics%integrator == IntegratorLEAP) then
        call leapfrog_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      else if (dynamics%integrator == IntegratorVVER) then
        call vverlet_dynamics (output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      end if
      enefunc%rpath_sum_avforce = .false.
      enefunc%mepmd_avforce = enefunc%mepmd_avforce / &
                              real(rpath%rpath_period,wp)

      ! set the energy and force
      !
      call ene_force_mepmd(dynvars, enefunc, rpath)

      ! output variables
      !
      call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                        boundary, rpath)

      if (main_rank .and. .not. rpath%equilibration_only) &
        call print_images(rpath)

      ! check convergence
      !
      call check_convergence(rpath, conv)
      if (conv)  exit

      ! path update
      !
      if (.not. rpath%equilibration_only) &
        call path_update_string(molecule, rpath, dynvars)

    end do


    return

  end subroutine run_rpath_string_mepmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_qmforce
  !> @brief        update QM force
  !! @authors      KY
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_qmforce(molecule, enefunc, pairlist, boundary, dynvars, &
                            constraints, rpath)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_rpath),    target, intent(inout) :: rpath

    ! local variables
    integer :: i, id

      ! compute the QM/MM energy and force
      !
      enefunc%qmmm%qm_classical = .false.
      call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                          .false.,           &
                          dynvars%coord,     &
                          dynvars%trans,     &
                          dynvars%coord_pbc, &
                          dynvars%energy,    &
                          dynvars%temporary, &
                          dynvars%force,     &
                          dynvars%force_omp, &
                          dynvars%virial,    &
                          dynvars%virial_extern)

      do i = 1, rpath%mep_natoms
        id = rpath%mepatom_id(i)
        rpath%mepmd_qmforce(:,i) = dynvars%force(:,id)
      end do

      !dbg if (main_rank) then
      !dbg   write(MsgOut,*) 'debug: QM/MM force'
      !dbg   do i = 1, rpath%mep_natoms
      !dbg     write(MsgOut,'(i8, 3f12.4)') rpath%mepatom_id(i), rpath%mepmd_qmforce(:,i)
      !dbg   end do
      !dbg end if

      ! compute the ESP/MM energy and force
      !
      enefunc%qmmm%qm_classical = .true.
      call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                          .false.,           &
                          dynvars%coord,     &
                          dynvars%trans,     &
                          dynvars%coord_pbc, &
                          dynvars%energy,    &
                          dynvars%temporary, &
                          dynvars%force,     &
                          dynvars%force_omp, &
                          dynvars%virial,    &
                          dynvars%virial_extern)

      do i = 1, rpath%mep_natoms
        id = rpath%mepatom_id(i)
        rpath%mepmd_qmforce(:,i) = rpath%mepmd_qmforce(:,i) - dynvars%force(:,id)
      end do

      !dbg if (main_rank) then
      !dbg   write(MsgOut,*) 'debug: QM force'
      !dbg   do i = 1, rpath%mep_natoms
      !dbg     write(MsgOut,'(i8, 3f12.4)') rpath%mepatom_id(i), rpath%mepmd_qmforce(:,i)
      !dbg   end do
      !dbg end if

      ! zero clear fixed atom force
      !
      call clear_fixatm_component(constraints, molecule%num_atoms, &
                                  dynvars%force)


  end subroutine update_qmforce

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_string_replica
  !> @brief        string method when nreplica is larger than total MPI proccess
  !! @authors      SI
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !! @param[out]   conv        : convergence flag
  !! @param[in]    niter       : current iteration number
  !! @param[in]    ireplica    : replica number in mpi process
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_string_replica(output, molecule, enefunc, dynvars, &
                                      minimize, dynamics, pairlist, boundary, &
                                      rpath, conv, niter, ireplica)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    logical,                  intent(out)   :: conv
    integer,                  intent(in)    :: niter
    integer,                  intent(in)    :: ireplica

    ! local variables
    integer                :: i, ii, j, repid, iatom
    integer                :: min_ene_period, save_step
    real(wp), pointer      :: coord(:,:)


    repid = my_replica_no

    min_ene_period = minimize%eneout_period

    ! First iteration
    !
    if (niter == 0) then
      if (ireplica == 1) then
        ! save eneout_period for minimize
        !
        min_ene_period = minimize%eneout_period
    
        ! initialize
        !
        rpath%energy_prev     = 0.0_wp
        rpath%mep_force       = 0.0_wp
        rpath%mep_length      = 0.0_wp
        rpath%pathlength_prev = 0.0_wp
    
        ! run rpath
        !
        dynvars%step = 0
        if (main_rank) then
           write(MsgOut, *) "--- Start string iteration ---"
           write(MsgOut, '(/,"Iter. ",i5)') 0
        end if
      end if

      ! compute the initial energy and force
      !
      call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                          .false.,           &
                          dynvars%coord,     &
                          dynvars%trans,     &
                          dynvars%coord_pbc, &
                          dynvars%energy,    &
                          dynvars%temporary, &
                          dynvars%force,     &
                          dynvars%force_omp, &
                          dynvars%virial,    &
                          dynvars%virial_extern)
      call energy_and_force_replica(dynvars, rpath)

      ! Store the old force by micro iteration
      ! 
      do i = 1, minimize%num_optatoms_micro
          iatom = minimize%optatom_micro_id(i)
          rpath%micro_force(1:3, i, ireplica) = dynvars%force(1:3, iatom)
      end do

      if (ireplica == nrep_per_proc) then
        ! calc initial path length
        !
        if (rpath%massWeightCoord) then
          call trans_mass_weight_coord_replica(molecule, rpath)
          call calc_pathlength_replica(rpath)
          call backtrans_mass_weight_coord_replica(molecule, rpath)
    
        else
          call calc_pathlength_replica(rpath)
    
        end if

        ! output collective variables
        !
        call print_images(rpath)
      end if

      ! output variables
      !
      call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                        boundary, rpath)
    end if

    ! If first cycle, then return
    !
    if (rpath%first_iter) return

    ! start string method
    !
    dynvars%step = niter

    ! update the coordinates
    !
    if (niter /= 1) then
      coord     => dynvars%coord 
      ii = 1
      do i = 1, rpath%mep_natoms
        coord(1:3, rpath%mepatom_id(i)) = rpath%mep_coord(ii:ii+2,repid)
        ii = ii + 3
      end do
    end if

    ! output settings
    !
    if (main_rank .and. rpath%first_replica) &
     write(MsgOut, '(/,"Iter. ",i5)') niter
    if (rpath%eneout_period > 0) then
      if (mod(niter, rpath%eneout_period) == 0) then
        if (replica_main_rank) write(output%logunit, '("Iter. ",i5,/)') niter
        minimize%eneout_period = min_ene_period

      else
        minimize%eneout_period = 0

      end if
    end if

    ! optimize surrounding atoms
    !
    if (rpath%mep_partial_opt) then
      save_step = dynvars%step
      if (rpath%opt_micro) then
        call minimize_micro(output, molecule, enefunc, dynvars, &
                           minimize, pairlist, boundary, rpath, niter, ireplica)

      else
        if (minimize%method == MinimizeMethodSD) then
          call steepest_descent(output, molecule, enefunc, dynvars, &
                                   minimize, pairlist, boundary)

        else if (minimize%method == MinimizeMethodLBFGS) then
          call minimize_lbfgs(output, molecule, enefunc, dynvars,   &
                               minimize, pairlist, boundary)
        end if

      end if
      dynvars%step = save_step

    end if

    ! set the energy and force
    !
    call energy_and_force_replica(dynvars, rpath)

    ! Store the old force by micro iteration
    !  
    do i = 1, minimize%num_optatoms_micro
        iatom = minimize%optatom_micro_id(i)
        rpath%micro_force(1:3, i, ireplica) = dynvars%force(1:3, iatom)
    end do

    ! check convergence
    !   
    if (ireplica == nrep_per_proc) call check_convergence(rpath, conv)

    ! output variables and images
    !
    call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                      boundary, rpath)    
    if (repid == nrep_per_proc) call print_images(rpath)

    ! path update
    !
    call path_update_string_replica(molecule, rpath, dynvars)

    return

  end subroutine run_rpath_string_replica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    path_update_string
  !> @brief        evolve path (MEP version)
  !! @authors      YK, YM, YA, KY
  !! @param[in]    molecule : molecular info
  !! @param[inout] rpath    : RPATH information
  !! @param[inout] dynvars  : dynamics variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine path_update_string(molecule, rpath, dynvars)

    ! formal arguments
    type(s_molecule),         intent(in)    :: molecule
    type(s_rpath),    target, intent(inout) :: rpath
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    integer           :: i, ii, nreplica
    integer           :: repid
    integer           :: d, dimension
    real(wp)          :: delta
    real(wp), pointer :: mep_force(:,:)
    real(wp), pointer :: mep_coord(:,:)
    real(wp), pointer :: recv_buff(:,:)
    real(wp), pointer :: coord(:,:)

    real(wp),    allocatable :: path_reparm(:)
    real(wp),    allocatable :: nonuniform_mesh(:)
    real(wp),    allocatable :: coeff(:,:)
    real(wp)                 :: tmp


    nreplica  =  rpath%nreplica
    dimension =  rpath%dimension
    delta     =  rpath%delta
    mep_force => rpath%mep_force
    mep_coord => rpath%mep_coord
    recv_buff => rpath%recv_buff
    coord     => dynvars%coord

    repid = my_country_no + 1

    if (rpath%massWeightCoord) call trans_mass_weight_coord(molecule, rpath)

    ! evolve images
    !
    do i = 1, dimension
      mep_coord(i, repid) = mep_coord(i, repid) + delta * mep_force(i, repid)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(mep_coord(1,repid), rpath%dimension, mpi_wp_real,&
                       recv_buff(1,1),     rpath%dimension, mpi_wp_real,&
                       mpi_comm_airplane, ierror)
    mep_coord = recv_buff
#endif

    ! calc arc lengths
    !
    rpath%pathlength_prev = rpath%pathlength
    call calc_pathlength(rpath, .false.)

    ! non-uniform normalized mesh along path
    !
    allocate(nonuniform_mesh(nreplica))
    do i = 1, nreplica
      nonuniform_mesh(i) = rpath%mep_length(i) / rpath%pathlength
    end do

    ! interpolate images to make equi-distance distribution
    !
    allocate(path_reparm(nreplica), coeff(4,nreplica-1))
    do d = 1, dimension
      call cubic_spline_coeff(nreplica, nonuniform_mesh, mep_coord(d,:), coeff)

      do i = 1, nreplica
        tmp = real(i-1,wp) / real(nreplica-1,wp)
        call cubic_spline(nreplica, nonuniform_mesh, mep_coord(d,:), coeff, &
                          tmp, path_reparm(i))
      end do
      mep_coord(d,:) = path_reparm
    end do
    deallocate(path_reparm, coeff)
    deallocate(nonuniform_mesh)

    ! update lengths
    !
    call calc_pathlength(rpath, .false.)

#ifdef HAVE_MPI_GENESIS
    ! broadcast MEP coordinates
    call mpi_bcast(mep_coord, rpath%dimension*rpath%nreplica,&
                   mpi_wp_real, 0, mpi_comm_world, ierror)
#endif

    if (rpath%massWeightCoord) call backtrans_mass_weight_coord(molecule, rpath)

    ! update the coordinates 
    ii = 1
    do i = 1, rpath%mep_natoms
      coord(1:3, rpath%mepatom_id(i)) = rpath%mep_coord(ii:ii+2,repid)
      ii = ii + 3
    end do

    return

  end subroutine path_update_string

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    path_update_string_replica
  !> @brief        evolve path (MEP version)
  !! @authors      SI
  !! @param[in]    molecule : molecular info
  !! @param[inout] rpath    : RPATH information
  !! @param[inout] dynvars  : dynamics variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine path_update_string_replica(molecule, rpath, dynvars)

    ! formal arguments
    type(s_molecule),         intent(in)    :: molecule
    type(s_rpath),    target, intent(inout) :: rpath
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    integer           :: i, ii, nreplica
    integer           :: repid
    integer           :: d, dimension
    real(wp)          :: delta
    real(wp), pointer :: mep_force(:,:)
    real(wp), pointer :: mep_coord(:,:)
    real(wp), pointer :: recv_buff(:,:)
    real(wp), pointer :: coord(:,:)

    real(wp),    allocatable :: path_reparm(:)
    real(wp),    allocatable :: nonuniform_mesh(:)
    real(wp),    allocatable :: coeff(:,:)
    real(wp)                 :: tmp


    nreplica  =  rpath%nreplica
    dimension =  rpath%dimension
    delta     =  rpath%delta
    mep_force => rpath%mep_force
    mep_coord => rpath%mep_coord
    recv_buff => rpath%recv_buff
    coord     => dynvars%coord

    repid = my_replica_no

    if (rpath%massWeightCoord) &
     call trans_mass_weight_coord(molecule, rpath)

    ! evolve images
    !
    do i = 1, dimension
      mep_coord(i, repid) = mep_coord(i, repid) + delta * mep_force(i, repid)
    end do

      !do i=1, nrep_per_proc
      !  ii = nrep_per_proc*(my_country_no) + i
      !  temp_ene(i) = rpath%mep_energy(ii)
      !  temp_ene_prev(i) = rpath%mep_energy_prev(ii)
      !end do
    i = nrep_per_proc*(my_country_no) + 1
    ii = i + nrep_per_proc - 1
    if (.not. rpath%first_replica) then
#ifdef HAVE_MPI_GENESIS    
      call mpi_allgather(mep_coord(1,i), rpath%dimension*nrep_per_proc, mpi_wp_real,&
                         recv_buff(1,1), rpath%dimension*nrep_per_proc, mpi_wp_real,&
                         mpi_comm_airplane, ierror)    
      mep_coord = recv_buff
#endif
      if (main_rank) then
        ! calc arc lengths
        !
        rpath%pathlength_prev = rpath%pathlength
        call calc_pathlength(rpath, .false.)

        ! non-uniform normalized mesh along path
        !
        allocate(nonuniform_mesh(nreplica))
        do i = 1, nreplica
          nonuniform_mesh(i) = rpath%mep_length(i) / rpath%pathlength
        end do

        ! interpolate images to make equi-distance distribution
        !
        allocate(path_reparm(nreplica), coeff(4,nreplica-1))
        do d = 1, dimension
          call cubic_spline_coeff(nreplica, nonuniform_mesh, mep_coord(d,:), coeff)

          do i = 1, nreplica
            tmp = real(i-1,wp) / real(nreplica-1,wp)
            call cubic_spline(nreplica, nonuniform_mesh, mep_coord(d,:), coeff, &
                              tmp, path_reparm(i))
          end do
          mep_coord(d,:) = path_reparm
        end do
        deallocate(path_reparm, coeff)
        deallocate(nonuniform_mesh)

        ! update lengths
        !
        call calc_pathlength(rpath, .false.) 

        if (rpath%massWeightCoord) &
         call backtrans_mass_weight_coord_replica(molecule, rpath)
      end if

#ifdef HAVE_MPI_GENESIS
      ! broadcast MEP coordinates
      call mpi_bcast(mep_coord, rpath%dimension*rpath%nreplica,&
                     mpi_wp_real, 0, mpi_comm_world, ierror)
#endif

    end if

    return

  end subroutine path_update_string_replica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    minimize_micro
  !> @brief        Optimization using ESP charges
  !! @authors      KY, YA, SI
  !! @param[in]    output       : output info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[inout] dynvars      : dynvars info
  !! @param[in]    minimize     : minimize info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[in]    rpath        : rpath info
  !! @param[in]    niter        : iteration of MEP search
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine minimize_micro(output, molecule, enefunc, dynvars, minimize, &
                            pairlist, boundary, rpath, niter, ireplica)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    integer,                  intent(in)    :: niter
    integer,        optional, intent(in)    :: ireplica

    ! local variables
    integer                :: i
    integer                :: natom_micro
    integer, pointer       :: optatom_micro_id(:)
    real(wp)               :: energy0corr, e0
    real(wp), allocatable  :: coord0(:,:), force0corr(:,:)

    real(wp)    , pointer  :: coord(:,:), force(:,:)
    type(s_qmmm), pointer  :: qmmm


    ! pointers
    !
    coord => dynvars%coord
    force => dynvars%force
    qmmm => enefunc%qmmm

    ! energy and gradient correction terms for micro_iteration
    !
    natom_micro       =  minimize%num_optatoms_micro
    optatom_micro_id  => minimize%optatom_micro_id

    energy0corr = dynvars%energy%total
    allocate(force0corr(3,natom_micro), coord0(3,natom_micro))
    do i = 1, natom_micro
      coord0(1:3,i)     = coord(1:3,optatom_micro_id(i))
      force0corr(1:3,i) = force(1:3,optatom_micro_id(i))
      if (nrep_per_proc > 1) &
        force0corr(1:3,i) = rpath%micro_force(1:3,i,ireplica)
    end do

    ! compute the ESP/MM energy and force
    !
    qmmm%qm_classical = .true.
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        enefunc%nonb_limiter, &
                        dynvars%coord,        &
                        dynvars%trans,        &
                        dynvars%coord_pbc,    &
                        dynvars%energy,       &
                        dynvars%temporary,    &
                        dynvars%force,        &
                        dynvars%force_omp,    &
                        dynvars%virial,       &
                        dynvars%virial_extern)

    ! get the correction
    !
    energy0corr = energy0corr - dynvars%energy%total
    do i = 1, natom_micro
      force0corr(1:3,i) = force0corr(1:3,i) - force(1:3,optatom_micro_id(i))
    end do

    ! QM internal energy
    !
    !rpath%qm_energy(my_country_no + 1) = energy0corr

    ! minimize the structure using micro-iteration
    !
    call micro_iter(output, molecule, enefunc, dynvars, minimize, &
             pairlist, boundary, coord0, energy0corr, force0corr, niter)

    deallocate(coord0, force0corr)

    ! compute the QM/MM energy and force
    !
    enefunc%qmmm%qm_classical = .false.
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        .false.,           &
                        dynvars%coord,     &
                        dynvars%trans,     &
                        dynvars%coord_pbc, &
                        dynvars%energy,    &
                        dynvars%temporary, &
                        dynvars%force,     &
                        dynvars%force_omp, &
                        dynvars%virial,    &
                        dynvars%virial_extern)

    return

  end subroutine minimize_micro

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gradient_correction
  !> @brief        energy and gradient correction for micro iteration
  !! @authors      YA, SI
  !! @param[in]    minimize     : minimize info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[inout] dynvars      : dynvars info
  !! @param[in]    rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gradient_correction(minimize, molecule, enefunc, pairlist, &
                boundary, dynvars, rpath, energy0corr, coord0, force0corr)

    ! formal arguments
    type(s_minimize), target,   intent(in)    :: minimize
    type(s_molecule),           intent(inout) :: molecule
    type(s_enefunc), target,    intent(inout) :: enefunc
    type(s_pairlist),           intent(in)    :: pairlist
    type(s_boundary),           intent(in)    :: boundary
    type(s_dynvars), target,    intent(inout) :: dynvars
    type(s_rpath),              intent(inout) :: rpath
    real(wp),                   intent(out)   :: energy0corr
    real(wp),                   intent(out)   :: coord0(3,*)
    real(wp),                   intent(out)   :: force0corr(3,*)

    ! local variables
    integer               :: i, repid, natom_micro
    integer, pointer      :: optatom_micro_id(:)
    real(wp), pointer     :: coord(:,:), force(:,:)
    type(s_qmmm), pointer :: qmmm


    natom_micro = minimize%num_optatoms_micro

    ! Pointers
    !
    coord => dynvars%coord
    force => dynvars%force
    qmmm => enefunc%qmmm
    optatom_micro_id => minimize%optatom_micro_id

    ! energy and gradient correction terms for micro_iteration
    !
    energy0corr = dynvars%energy%total
    do i = 1, natom_micro
      coord0(1:3,i)     = coord(1:3,optatom_micro_id(i))
      force0corr(1:3,i) = force(1:3,optatom_micro_id(i))
    end do

    qmmm%qm_classical = .TRUE.
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        enefunc%nonb_limiter, &
                        dynvars%coord,        &
                        dynvars%trans,        &
                        dynvars%coord_pbc,    &
                        dynvars%energy,       &
                        dynvars%temporary,    &
                        dynvars%force,        &
                        dynvars%force_omp,    &
                        dynvars%virial,       &
                        dynvars%virial_extern)

    energy0corr = energy0corr - dynvars%energy%total
    do i = 1, natom_micro
      force0corr(1:3,i) = force0corr(1:3,i) - force(1:3,optatom_micro_id(i))
    end do

    ! QM internal energy
    !
    repid = my_country_no + 1
    if (nrep_per_proc > 1) repid = my_replica_no

    rpath%qm_energy(repid) = energy0corr

    return

  end subroutine gradient_correction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    energy_and_force
  !> @brief        copy the energy and force from dynvars to rpath
  !! @authors      YA, KY
  !! @param[inout] dynvars      : dynvars info
  !! @param[in]    rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine energy_and_force(dynvars, rpath)

    ! formal arguments
    type(s_dynvars),  intent(inout) :: dynvars
    type(s_rpath),    intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: iatom, i


    repid = my_country_no + 1

    ! energy
    !
    rpath%energy = dynvars%energy%total
    rpath%mep_energy_prev = rpath%mep_energy
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(rpath%energy,     1, mpi_wp_real, &
                       rpath%mep_energy, 1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif

    ! force
    !
    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      rpath%mep_force(i*3-2:i*3, repid) = dynvars%force(1:3, iatom)
    end do

    ! clear force for terminals
    !
    if (rpath%fix_terminal) then
      if ((repid == 1) .or. (repid == rpath%nreplica)) then
        rpath%mep_force(:, repid) = 0.0_wp
      end if
    end if

    return

  end subroutine energy_and_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    ene_force_mepmd
  !> @brief        copy the energy and force from dynvars to rpath
  !! @authors      KY
  !! @param[inout] dynvars     : dynvars info
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] rpath       : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine ene_force_mepmd(dynvars, enefunc, rpath)

    ! formal arguments
    type(s_dynvars),  intent(inout) :: dynvars
    type(s_enefunc),  intent(inout) :: enefunc
    type(s_rpath),    intent(inout) :: rpath

    ! local variables
    integer               :: repid
    integer               :: i, j, ij
    real(wp)              :: gs, vnorm
    real(wp), allocatable :: gsi(:), enesi(:), vec(:)


    repid = my_country_no + 1

    ! force
    !
    do i = 1, rpath%mep_natoms
      rpath%mep_force(i*3-2:i*3, repid)          &
                    = rpath%mepmd_qmforce(:,i)   &
                    + enefunc%mepmd_avforce(:,i)
    end do

    !dbg if (main_rank) then
    !dbg   write(MsgOut,*) 'debug: MM average force'
    !dbg   do i = 1, rpath%mep_natoms
    !dbg     write(MsgOut,'(i8, 3f12.4)') rpath%mepatom_id(i), enefunc%mepmd_avforce(:,i)
    !dbg   end do
    !dbg   write(MsgOut,*) 'debug: MEP force'
    !dbg   do i = 1, rpath%mep_natoms
    !dbg     write(MsgOut,'(i8, 3f12.4)') rpath%mepatom_id(i), rpath%mep_force(i*3-2:i*3, repid)
    !dbg   end do
    !dbg end if

    ! energy
    !

    allocate(gsi(rpath%nreplica), vec(rpath%mep_natoms*3))

    if (repid == 1) then
      vec = rpath%mep_coord(:,repid+1) - rpath%mep_coord(:,repid)
    else if ((repid == rpath%nreplica)) then
      vec = rpath%mep_coord(:,repid)   - rpath%mep_coord(:,repid-1)
    else
      vec = rpath%mep_coord(:,repid+1) - rpath%mep_coord(:,repid-1)
    end if

    vnorm = 0.0_wp
    do ij = 1, rpath%mep_natoms*3
      vnorm = vnorm + vec(ij)*vec(ij)
    end do
    vnorm = sqrt(vnorm)
    vec   = vec / vnorm

    gs = 0.0_wp
    do ij = 1, rpath%mep_natoms*3
      gs = gs - vec(ij)*rpath%mep_force(ij, repid)
    end do

    !dbg if (main_rank) then
    !dbg   write(MsgOut,*) 'debug: vec and gs'
    !dbg   write(MsgOut,'(3f12.4)') vec
    !dbg   write(MsgOut,'(f12.6)') gs
    !dbg end if

#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(gs,  1, mpi_wp_real, gsi, 1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif

    allocate(enesi(rpath%nreplica))
    enesi(1) = 0.0_wp
    do i = 2, rpath%nreplica
      enesi(i) = enesi(i-1) + (rpath%mep_length(i) - rpath%mep_length(i-1)) &
                            * (gsi(i-1) + gsi(i)) * 0.5_wp
    end do

    !dbg if (main_rank) then
    !dbg   write(MsgOut,*) 'debug: rpath component of gradient and energy'
    !dbg   do i = 1, rpath%nreplica
    !dbg      write(MsgOut,'(i8, 2f12.4)') i, gsi(i), enesi(i)
    !dbg   end do
    !dbg end if

    rpath%mep_energy_prev = rpath%mep_energy
    rpath%mep_energy = enesi

    deallocate(gsi, vec, enesi)

    ! clear force for terminals
    !
    if(rpath%fix_terminal) then
      if((repid == 1) .or. (repid == rpath%nreplica)) then
        rpath%mep_force(:, repid) = 0.0_wp
      end if
    end if

    return

  end subroutine ene_force_mepmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    energy_and_force_replica
  !> @brief        copy the energy and force from dynvars to rpath
  !! @authors      SI
  !! @param[inout] dynvars      : dynvars info
  !! @param[in]    rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine energy_and_force_replica(dynvars, rpath)

    ! formal arguments
    type(s_dynvars),  intent(inout) :: dynvars
    type(s_rpath),    intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: iatom, i, ii
    real(wp),           allocatable :: temp_ene(:), temp_ene_prev(:)


    repid = my_replica_no

    ! energy
    !
    rpath%mep_energy_prev(repid) = rpath%mep_energy(repid)
    rpath%mep_energy(repid) = dynvars%energy%total

#ifdef HAVE_MPI_GENESIS
    if (mod(repid,nrep_per_proc) == 0) then
      allocate(temp_ene(nrep_per_proc))
      allocate(temp_ene_prev(nrep_per_proc))
      do i=1, nrep_per_proc
        ii = nrep_per_proc*(my_country_no) + i
        temp_ene(i) = rpath%mep_energy(ii)
        temp_ene_prev(i) = rpath%mep_energy_prev(ii)
      end do
      call mpi_allgather(temp_ene, nrep_per_proc, mpi_wp_real, &
                         rpath%mep_energy, nrep_per_proc, mpi_wp_real, &
                         mpi_comm_airplane, ierror)
      call mpi_allgather(temp_ene_prev, nrep_per_proc, mpi_wp_real, &
                         rpath%mep_energy_prev, nrep_per_proc, mpi_wp_real, &
                         mpi_comm_airplane, ierror)      
      deallocate(temp_ene)
      deallocate(temp_ene_prev)
    end if
#endif

    ! force
    !
    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      rpath%mep_force(i*3-2:i*3, repid) = dynvars%force(1:3, iatom)
    end do

    ! clear force for terminals
    !
    if (rpath%fix_terminal) then
      if ((repid == 1) .or. (repid == rpath%nreplica)) then
        rpath%mep_force(:, repid) = 0.0_wp
      end if
    end if

    return

  end subroutine energy_and_force_replica  

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_convergence
  !> @brief        Check the convergence of the energy and pathlength
  !! @authors      YA, KY
  !! @param[in]    rpath        : rpath info
  !! @param[out]   conv         : if true, convergence achieved
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_convergence(rpath, conv)

    ! formal arguments
    type(s_rpath),      intent(inout) :: rpath
    logical,            intent(out)   :: conv

    ! local variables
    integer                         :: i, iatom
    logical                         :: my_conv, conv_e, conv_path
    integer                         :: repid
    real(wp)                        :: disp, delta, dmax
    real(wp),           allocatable :: delta_ene(:)


    ! energy convergence
    !
    allocate(delta_ene(rpath%nreplica))
    delta_ene = rpath%mep_energy - rpath%mep_energy_prev
    dmax = delta_ene(1)
    do i = 2, rpath%nreplica
      if (abs(dmax) < abs(delta_ene(i))) dmax = delta_ene(i)
    end do
    deallocate(delta_ene)

    if (abs(dmax) < rpath%tol_energy) then
      conv_e = .true.
    else
      conv_e = .false.
    end if

    ! pathlength convergence
    !
    if (abs(rpath%pathlength - rpath%pathlength_prev) < rpath%tol_path) then
      conv_path = .true.
    else
      conv_path = .false.
    end if

    if (conv_e .and. conv_path) then
      conv = .true.
    else
      conv = .false.
    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(conv, 1, mpi_logical, 0, mpi_comm_world, ierror)
#endif
    
    return

  end subroutine check_convergence

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_images
  !> @brief        Print the information of images 
  !! @authors      KY
  !! @param[in]    rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_images(rpath)

    ! formal arguments
    type(s_rpath),    intent(inout) :: rpath

    ! local variables
    integer                         :: i
    real(wp)                        :: dmax
    real(wp),           allocatable :: delta_ene(:)


    if (main_rank) then

      allocate(delta_ene(rpath%nreplica))
      delta_ene = rpath%mep_energy - rpath%mep_energy_prev

      dmax = delta_ene(1)
      do i = 2, rpath%nreplica
        if (abs(dmax) < abs(delta_ene(i))) dmax = delta_ene(i)
      end do
    
      
      write(MsgOut, '(/,10x, &
                        "         Path Length", &
                        "   Energy (kcal/mol)", &
                        "         Relative E.", &
                        "        Energy Conv.")')
      write(MsgOut, '(90("-"))')
      do i = 1, rpath%nreplica
        write(MsgOut, '("Image ",i3,x,4f20.4)') &
          i,                   &
          rpath%mep_length(i), &
          rpath%mep_energy(i), &
          rpath%mep_energy(i) - rpath%mep_energy(1), &
          delta_ene(i)
      end do
      write(MsgOut, '(90("-"))')

      deallocate(delta_ene)
      
      write(MsgOut, '(3X, "Energy Conv. (Max) = ", F20.8)') dmax
      write(MsgOut, '(3X, "Path length: current value / variation = ", &
                      F15.8, " /", F15.8)') &
                      rpath%pathlength, &
                      rpath%pathlength - rpath%pathlength_prev

    end if

    return

  end subroutine print_images

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_convergence_neb
  !> @brief        
  !! @authors      YA, KY
  !! @param[in]    rpath        : rpath info
  !! @param[in]    molecule     : molecular info
  !! @param[in]    enefunc      : enefunc info
  !! @param[in]    pairlist     : pairlist info
  !! @param[in]    boundary     : boundary info
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[inout] dynvars      : dynvars info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_convergence_neb(rpath, niter, conv)

    ! formal arguments
    type(s_rpath),       intent(inout) :: rpath
    integer,             intent(in)    :: niter
    logical,             intent(out)   :: conv

    ! local variables
    character(3)                  :: stat
    logical                       :: conv_image(rpath%nreplica)
    integer                       :: repid, i, image
    real(wp)                      :: disp, delta, absg, dmax, rmsg_max, maxg_max
    real(wp),         allocatable :: work(:), rmsg(:), maxg(:)


    repid = my_country_no + 1

    allocate(rmsg(rpath%nreplica))
    allocate(maxg(rpath%nreplica))
    allocate(work(rpath%nreplica))

    delta = rpath%energy - rpath%energy_prev
    rpath%energy_prev = rpath%energy
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(delta, 1, mpi_wp_real, work, 1, mpi_wp_real, &
      mpi_comm_airplane, ierror)
#endif
    ! RMSG and MaxG
    !
    rmsg(:) = 0.0_wp
    maxg(:) = 0.0_wp
    do image = 1, rpath%nreplica
      do i = 1, rpath%mep_natoms * 3
        rmsg(image) = rmsg(image) + rpath%mep_force(i,image)**2
        absg = abs(rpath%mep_force(i,image))
        if (absg > maxg(image)) maxg(image) = absg
      end do
      rmsg(image) = sqrt(rmsg(image) / real(rpath%mep_natoms*3,wp))
      conv_image(image) = (rmsg(image) < rpath%tol_rmsg) .and. (maxg(image) < rpath%tol_maxg)
    end do

    rmsg_max = maxval(rmsg(:))
    maxg_max = maxval(maxg(:))
    conv = (rmsg_max < rpath%tol_rmsg) .and. (maxg_max < rpath%tol_maxg)

    dmax = 0.0_wp
    do i = 1, rpath%nreplica
      if (abs(work(i)) > dmax) dmax = abs(work(i))
    end do

    ! Output
    !
    if (main_rank) then
      
      write(MsgOut, '(/,"Image ", &
                        "     Energy (kcal/mol)", &
                        "           Relative E.", &
                        "                DeltaE", &
                        "                  RMSG", &
                        "                  MaxG", &
                        "           Convergence")')
      write(MsgOut, '("------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------")')
      do i = 1, rpath%nreplica
        if (conv_image(i)) then
          stat = "yes"
        else
          stat = " no"
        end if
        write(MsgOut, '(i5,X,5F22.10,19x,a3)') i, &
                             rpath%mep_energy(i), &
       rpath%mep_energy(i) - rpath%mep_energy(1), &
                                         work(i), &
                                         rmsg(i), &
                                         maxg(i), &
                                         stat
      end do
      write(MsgOut, '("------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------", &
                      "----------------------")')

      write(MsgOut, '(3X, "Energy Conv. (Max) = ", F25.15)') dmax
      write(MsgOut, '(3X, "RMSG   Conv. (Max) = ", F25.15)') rmsg_max
      write(MsgOut, '(3X, "MAXG   Conv. (Max) = ", F25.15)') maxg_max

    end if

    ! Switch-on CINEB
    if (rpath%climbing_image .and. rmsg_max < rpath%tol_rmsg_cineb ) then
      rpath%do_cineb = .true.
    end if
    
    deallocate(work)
    deallocate(rmsg)
    deallocate(maxg)
    
    return

  end subroutine check_convergence_neb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trans_mass_weight_coord
  !> @brief        Transform to mass-weighted Cartesian coordinates
  !! @authors      YA, KY
  !! @param[in]    molecule     : molecular info
  !! @param[inout] rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trans_mass_weight_coord(molecule, rpath)

    ! formal arguments
    type(s_molecule),    intent(in)    :: molecule
    type(s_rpath),       intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: i, ii, iatom, offset
    real(wp)                        :: massfac


    repid = my_country_no + 1
    if (nrep_per_proc > 1) repid = my_replica_no

    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      massfac = sqrt(molecule%mass(iatom))

      offset = (i - 1) * 3
      rpath%mep_coord(offset+1:offset+3,repid) = &
        rpath%mep_coord(offset+1:offset+3,repid) * massfac
      rpath%mep_force(offset+1:offset+3,repid) = &
        rpath%mep_force(offset+1:offset+3,repid) / massfac
    end do

    return

  end subroutine trans_mass_weight_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trans_mass_weight_coord
  !> @brief        Transform to mass-weighted Cartesian coordinates
  !! @authors      YA, KY
  !! @param[in]    molecule     : molecular info
  !! @param[inout] rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trans_mass_weight_coord_replica(molecule, rpath)

    ! formal arguments
    type(s_molecule),    intent(in)    :: molecule
    type(s_rpath),       intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: i, ii, iatom, irep, offset
    real(wp)                        :: massfac


    do irep = 1, nrep_per_proc
      repid = nrep_per_proc*(my_country_no) + irep
      do i = 1, rpath%mep_natoms
        iatom = rpath%mepatom_id(i)
        massfac = sqrt(molecule%mass(iatom))
  
        offset = (i - 1) * 3
        rpath%mep_coord(offset+1:offset+3,repid) = &
          rpath%mep_coord(offset+1:offset+3,repid) * massfac
        rpath%mep_force(offset+1:offset+3,repid) = &
          rpath%mep_force(offset+1:offset+3,repid) / massfac
      end do
    end do

    return

  end subroutine trans_mass_weight_coord_replica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    backtrans_mass_weight_coord
  !> @brief        
  !! @authors      YA, KY
  !! @param[in]    molecule     : molecular info
  !! @param[inout] rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine backtrans_mass_weight_coord(molecule, rpath)

    ! formal arguments
    type(s_molecule),    intent(in)    :: molecule
    type(s_rpath),       intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: i, ii, iatom, offset
    real(wp)                        :: massfac


    repid = my_country_no + 1

    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      massfac = sqrt(molecule%mass(iatom))

      offset = (i - 1) * 3
      rpath%mep_coord(offset+1:offset+3,repid) = &
        rpath%mep_coord(offset+1:offset+3,repid) / massfac
      rpath%mep_force(offset+1:offset+3,repid) = &
        rpath%mep_force(offset+1:offset+3,repid) * massfac
    end do

    return

  end subroutine backtrans_mass_weight_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    backtrans_mass_weight_coord
  !> @brief        
  !! @authors      YA, KY
  !! @param[in]    molecule     : molecular info
  !! @param[inout] rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine backtrans_mass_weight_coord_replica(molecule, rpath)

    ! formal arguments
    type(s_molecule),    intent(in)    :: molecule
    type(s_rpath),       intent(inout) :: rpath

    ! local variables
    integer                         :: repid
    integer                         :: i, ii, iatom, irep, offset
    real(wp)                        :: massfac


    do irep = 1, rpath%nreplica
      do i = 1, rpath%mep_natoms
        iatom = rpath%mepatom_id(i)
        massfac = sqrt(molecule%mass(iatom))
  
        offset = (i - 1) * 3
        rpath%mep_coord(offset+1:offset+3,irep) = &
          rpath%mep_coord(offset+1:offset+3,irep) / massfac
        rpath%mep_force(offset+1:offset+3,irep) = &
          rpath%mep_force(offset+1:offset+3,irep) * massfac
      end do
    end do

    return

  end subroutine backtrans_mass_weight_coord_replica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_pathlength
  !> @brief        Calculate the arc length of current MEP
  !! @authors      KY
  !! @param[inout] rpath : rpath info
  !! @param[in]    gather_coord : gather mep_coord, if true
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_pathlength(rpath, gather_coord)

    ! formal arguments
    type(s_rpath), target, intent(inout) :: rpath
    logical, optional,     intent(in)    :: gather_coord

    ! local variables
    integer           :: i, n
    integer           :: repid
    logical           :: flag
    real(wp)          :: dd
    real(wp), pointer :: mep_coord(:,:)
    real(wp), pointer :: recv_buff(:,:)
    real(wp), pointer :: mep_length(:)


    mep_coord  => rpath%mep_coord
    recv_buff  => rpath%recv_buff
    mep_length => rpath%mep_length

    if (.not. present(gather_coord)) then
        flag = .true.
    else
        flag = gather_coord
    end if

    if (flag) then
#ifdef HAVE_MPI_GENESIS
      repid = my_country_no + 1
      call mpi_allgather(mep_coord(1,repid), rpath%dimension, mpi_wp_real,&
                         recv_buff(1,1),     rpath%dimension, mpi_wp_real,&
                         mpi_comm_airplane, ierror)
      mep_coord = recv_buff
#endif
    end if

    mep_length(1) = 0.0_wp
    do i = 2, rpath%nreplica

      dd = 0.0_wp
      do n = 1, rpath%dimension
        dd = dd + (mep_coord(n,i) - mep_coord(n,i-1))**2
      end do
      mep_length(i) = mep_length(i-1) + sqrt(dd)

    end do
    rpath%pathlength = mep_length(rpath%nreplica)

    return

  end subroutine calc_pathlength

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_pathlength_replica
  !> @brief        Calculate the arc length of current MEP
  !! @authors      KY, SI
  !! @param[inout] rpath : rpath info
  !! @param[in]    gather_coord : gather mep_coord, if true
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_pathlength_replica(rpath, gather_coord)

    ! formal arguments
    type(s_rpath), target, intent(inout) :: rpath
    logical, optional,     intent(in)    :: gather_coord

    ! local variables
    integer               :: i, n
    integer               :: repid
    real(wp)              :: dd
    real(wp), pointer     :: mep_coord(:,:)
    real(wp), pointer     :: recv_buff(:,:)
    real(wp), pointer     :: mep_length(:)


    mep_coord  => rpath%mep_coord
    recv_buff  => rpath%recv_buff
    mep_length => rpath%mep_length

    if (.not. present(gather_coord) .or. &
        (present(gather_coord) .and. gather_coord)) then
#ifdef HAVE_MPI_GENESIS
      repid = my_replica_no

      call mpi_allgather(mep_coord(1,(repid-nrep_per_proc)+1), &
                         (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                         recv_buff(1,1), &    
                         (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                         mpi_comm_airplane, ierror)
      mep_coord = recv_buff
#endif
    end if

    mep_length(1) = 0.0_wp
    do i = 2, rpath%nreplica

      dd = 0.0_wp
      do n = 1, rpath%dimension
        dd = dd + (mep_coord(n,i) - mep_coord(n,i-1))**2
      end do
      mep_length(i) = mep_length(i-1) + sqrt(dd)

    end do
    rpath%pathlength = mep_length(rpath%nreplica)

    return

  end subroutine calc_pathlength_replica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_neb_glbfgs
  !> @brief        global LBFGS-NEB
  !! @authors      YA, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_neb_glbfgs(output, molecule, enefunc, dynvars, &
                    minimize, dynamics, pairlist, boundary, rpath, conv, niter)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    logical,                  intent(out)   :: conv
    integer,                  intent(out)   :: niter

    ! local variables
    integer                :: i, j
    integer                :: replicaid
    integer                :: ii, iatom, ioff, image

    real(wp)    , pointer  :: coord(:,:)
    type(s_qmmm), pointer  :: qmmm

    logical                :: bnd, bnd_qmonly
    real(wp)               :: maxmove
    real(wp)               :: neb_energy, ddot
    real(wp), allocatable  :: tmp(:)
    integer, parameter     :: start_minimize = 2

    ! local variabels for LBFGS
    character(60)          :: csave, task
    logical                :: lsave(4)
    integer                :: no3, ncorr, nwork, dim
    integer                :: isave(44)
    integer                :: iprint = -1
    integer, allocatable   :: list_bound(:), iwa(:)
    real(wp)               :: dsave(29), factr, pgtol
    real(wp), allocatable  :: vec(:), gradient(:), wa(:), lower(:), upper(:)

    integer                :: min_ene_period, save_step


    ! save eneout_period for minimize
    !
    min_ene_period = minimize%eneout_period

    ! initialize
    !
    rpath%energy_prev     = 0.0_wp
    rpath%mep_force       = 0.0_wp
    rpath%mep_length      = 0.0_wp
    rpath%pathlength_prev = 0.0_wp

    ! Pointers
    !
    coord => dynvars%coord
    qmmm  => enefunc%qmmm

    ! My replica ID
    !
    replicaid = my_country_no + 1

    ! allocate work arrays for G-LBFGS
    !
    ncorr = rpath%ncorrection
    no3 = rpath%mep_natoms * 3
    dim = no3 * rpath%nreplica
    if (main_rank) then
      nwork = (2*dim + 11*ncorr + 8) * ncorr + 5*dim
      allocate(vec(dim))
      allocate(gradient(dim))
      allocate(lower(dim))
      allocate(upper(dim))
      allocate(list_bound(dim))
      allocate(iwa(dim*3))
      allocate(wa(nwork))
      allocate(tmp(no3))

      bnd        = rpath%lbfgs_bnd
      bnd_qmonly = rpath%lbfgs_bnd_qmonly
      maxmove    = rpath%lbfgs_bnd_maxmove

      list_bound(:) = 0
      if (bnd) then
        if (bnd_qmonly) then
          ! set boundary to QM atoms
          do i = 1, rpath%mep_natoms
            do j = 1, enefunc%qmmm%qm_natoms
              if (rpath%mepatom_id(i) == enefunc%qmmm%qmatom_id(j)) then
                list_bound(3*i-2:3*i) = 2
                exit
              end if
            end do
          end do
        else
          ! set boundary to all atoms
          list_bound = 2
        end if

        ioff = no3
        do image = 2, rpath%nreplica
          list_bound(ioff+1:ioff+no3) = list_bound(1:no3)
          ioff = ioff + no3
        end do

      end if

    end if

    ! run rpath
    !
    dynvars%step = 0
    if (main_rank) then
      write(MsgOut, *) "--- Start NEB iteration ---"
    end if

    task = 'START'

    niter = -1
    do
      if (task(1:2) .eq. 'FG') then
        niter = niter + 1
        dynvars%step = niter
        if (niter /= 0) rpath%first_iter = .false.

        ! output settings
        !
        if (main_rank) write(MsgOut, '(/,"Iter. ",i5)') niter
        if (rpath%eneout_period > 0) then
          if (mod(niter, rpath%eneout_period) == 0) then
            minimize%eneout_period = min_ene_period

          else
            minimize%eneout_period = 0

          end if
        end if

        ! optimize surrounding atoms
        !
        if (rpath%mep_partial_opt) then
          if (niter > start_minimize) then

            save_step = dynvars%step
            if (rpath%opt_micro) then
              call minimize_micro (output, molecule, enefunc, dynvars, &
                                   minimize, pairlist, boundary, rpath, niter)
            else
              
              if (minimize%method == MinimizeMethodSD) then
                call steepest_descent (output, molecule, enefunc, dynvars, &
                                       minimize, pairlist, boundary)
                
              else if (minimize%method == MinimizeMethodLBFGS) then
                call minimize_lbfgs (output, molecule, enefunc, dynvars, &
                                     minimize, pairlist, boundary)
              end if
            end if
            dynvars%step = save_step

          else
            ! compute the QM/MM energy and force
            !
            enefunc%qmmm%qm_classical = .false.
            call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                                .false.,           &
                                dynvars%coord,     &
                                dynvars%trans,     &
                                dynvars%coord_pbc, &
                                dynvars%energy,    &
                                dynvars%temporary, &
                                dynvars%force,     &
                                dynvars%force_omp, &
                                dynvars%virial,    &
                                dynvars%virial_extern)

          end if
          
        end if
        
        ! energy and NEB force
        !
        call energy_and_force(dynvars, rpath)
        call calc_neb_force(rpath)

        ! check convergence and output variables
        !
        call check_convergence(rpath, conv)
        if (conv .or. (niter == rpath%ncycle))  exit

        ! output rpath info (rstmep)
        !
        call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                          boundary, rpath)
        call print_images(rpath)

      end if

      ! path update
      !
      if (rpath%massWeightCoord) call trans_mass_weight_coord(molecule, rpath)
#ifdef HAVE_MPI_GENESIS
      call mpi_allgather(rpath%mep_coord(1,replicaid), rpath%dimension, mpi_wp_real,&
                         rpath%recv_buff(1,1),         rpath%dimension, mpi_wp_real,&
                         mpi_comm_airplane, ierror)
      rpath%mep_coord = rpath%recv_buff
      call mpi_allgather(rpath%mep_force(1,replicaid), rpath%dimension, mpi_wp_real,&
                         rpath%recv_buff(1,1),         rpath%dimension, mpi_wp_real,&
                         mpi_comm_airplane, ierror)
      rpath%mep_force = rpath%recv_buff
#endif

      if (main_rank) then
        
        ! global vector and gradient
        ioff = 0
        do image = 1, rpath%nreplica
          vec(ioff+1:ioff+no3)      =  rpath%mep_coord(1:no3,image)
          gradient(ioff+1:ioff+no3) = -rpath%mep_force(1:no3,image)
          ioff = ioff + no3
        end do

        ! global energy
        neb_energy = 0.0_wp
        do image = 1, rpath%nreplica
          neb_energy = neb_energy + rpath%mep_energy(image)
          if (image > 1) then
            tmp(:) = rpath%mep_coord(:,image) - rpath%mep_coord(:,image-1)
            !YA_20180119 neb_energy = neb_energy + &
            !              0.5_wp * ddot(no3, tmp, 1, tmp, 1) * rpath%k_spring
            neb_energy = neb_energy                         &
                       + 0.5_wp * ddot(no3, tmp, 1, tmp, 1) &
                       * rpath%k_spring * dble(rpath%nreplica - 1)
          end if
        end do

        ! boundary
        if (bnd) then
          if (bnd_qmonly) then
            ! set boundary to QM atoms
            do i = 1, dim
              if (list_bound(i) == 2) then
                upper(i) = vec(i) + maxmove
                lower(i) = vec(i) - maxmove
              end if
            end do
          else
            ! set boundary to all atoms
            do i = 1, dim
              upper(i) = vec(i) + maxmove
              lower(i) = vec(i) - maxmove
            end do
          end if
        end if

        factr = 0.0_wp
        pgtol = 0.0_wp
        if (rpath%verbose .and. main_rank) iprint = 1
        call setulb(dim, ncorr, vec, lower, upper, list_bound, neb_energy, &
                    gradient, factr, pgtol, wa, iwa, task, &
                    iprint, csave, lsave, isave, dsave, my_country_no)

        ioff = 0
        do image = 1, rpath%nreplica
          rpath%mep_coord(1:no3,image) = vec(ioff+1:ioff+no3)
          ioff = ioff + no3
        end do
      
      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(task, 60, mpi_character, 0, mpi_comm_world, ierror)
      call mpi_bcast(rpath%mep_coord(1,1), dim, mpi_wp_real, &
                     0, mpi_comm_world, ierror)
#endif

      ! update pathlength
      !
      rpath%pathlength_prev = rpath%pathlength
      call calc_pathlength(rpath, .false.)

      if (rpath%massWeightCoord) call backtrans_mass_weight_coord(molecule, rpath)

      ! update the coordinates 
      !
      ii = 0
      do i = 1, rpath%mep_natoms
        iatom = rpath%mepatom_id(i)
        coord(1:3,iatom) = rpath%mep_coord(ii+1:ii+3,replicaid)
        ii = ii + 3
      end do

    end do

    ! output final status
    !
    rpath%eneout = .true.
    rpath%crdout = .true.
    rpath%rstout = .true.
    call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                      boundary, rpath)
    rpath%eneout = .false.
    rpath%crdout = .false.
    rpath%rstout = .false.

    call print_images(rpath)

    ! free work arrays
    !
    if (main_rank) then
      deallocate(vec)
      deallocate(gradient)
      deallocate(lower)
      deallocate(upper)
      deallocate(list_bound)
      deallocate(iwa)
      deallocate(wa)
      deallocate(tmp)
    end if

    return

  end subroutine run_rpath_neb_glbfgs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_neb_glbfgs_replica
  !> @brief        global LBFGS-NEB if npreplica > nproc
  !! @authors      SI
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] rpath       : RPATH information
  !! @param[out]   conv        : convergence flag
  !! @param[in]    niter       : current iteration number
  !! @param[in]    ireplica    : replica number in mpi process  
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_neb_glbfgs_replica(output, molecule, enefunc, dynvars, &
                                          minimize, dynamics, pairlist, &
                                          boundary, rpath, conv, niter, &
                                          ireplica)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_rpath),            intent(inout) :: rpath
    logical,                  intent(out)   :: conv
    integer,                  intent(in)    :: niter
    integer,                  intent(in)    :: ireplica

    ! local variables
    integer                :: i, j
    integer                :: replicaid
    integer                :: ii, iatom, ioff, image

    real(wp)    , pointer  :: coord(:,:), force(:,:), force_omp(:,:,:)
    type(s_qmmm), pointer  :: qmmm

    logical                :: bnd, bnd_qmonly
    real(wp)               :: maxmove
    real(wp)               :: neb_energy, ddot
    real(wp), allocatable  :: tmp(:)
    integer, parameter     :: start_minimize = 2

    ! local variabels for LBFGS
    integer                :: no3, ncorr, nwork, dim
    integer                :: iprint = -1
    integer, allocatable   :: list_bound(:)
    real(wp)               :: factr, pgtol
    real(wp), allocatable  :: vec(:), gradient(:), lower(:), upper(:)

    integer                :: min_ene_period, save_step


    ! Pointers
    !
    coord => dynvars%coord
    force => dynvars%force
    force_omp => dynvars%force_omp
    qmmm  => enefunc%qmmm

    ! My replica ID
    !
    replicaid = my_replica_no

    ! initialize
    !
    ncorr = rpath%ncorrection
    no3   = rpath%mep_natoms * 3
    dim   = no3 * rpath%nreplica
    nwork = (2*dim + 11*ncorr + 8) * ncorr + 5*dim
    rpath%rstout     = .true.
    rpath%neb_output = .false.

    if ((niter == 0) .and. (ireplica == 1)) then
      ! save eneout_period for minimize
      !
      min_ene_period = minimize%eneout_period

      ! initialize
      !
      rpath%neb_cycle       = 0
      rpath%energy_prev     = 0.0_wp
      rpath%mep_force       = 0.0_wp
      rpath%mep_length      = 0.0_wp
      rpath%pathlength_prev = 0.0_wp   
    end if   

    if (main_rank) then
      ! allocate work arrays for G-LBFGS
      !
      allocate(vec(dim))
      allocate(gradient(dim))
      allocate(lower(dim))
      allocate(upper(dim))
      allocate(list_bound(dim))
      allocate(tmp(no3))

      bnd        = rpath%lbfgs_bnd
      bnd_qmonly = rpath%lbfgs_bnd_qmonly
      maxmove    = rpath%lbfgs_bnd_maxmove

      list_bound(:) = 0
      if (bnd) then
        if (bnd_qmonly) then
          ! set boundary to QM atoms
          do i = 1, rpath%mep_natoms
            do j = 1, enefunc%qmmm%qm_natoms
              if (rpath%mepatom_id(i) == enefunc%qmmm%qmatom_id(j)) then
                list_bound(3*i-2:3*i) = 2
                exit
              end if
            end do
          end do
        else
          ! set boundary to all atoms
          list_bound = 2
        end if

        ioff = no3
        do image = 2, rpath%nreplica
          list_bound(ioff+1:ioff+no3) = list_bound(1:no3)
          ioff = ioff + no3
        end do

      end if

    end if 

    ! Run RPATH
    !    
    if (niter == 0) then
      rpath%task = 'START'

      dynvars%step = 0
      if (main_rank .and. (ireplica == 1)) then
        write(MsgOut, *) "--- Start NEB iteration ---"
      end if  
    end if

    ! update the coordinates
    !
    if (niter /= 0) then
      ii = 1
      do i = 1, rpath%mep_natoms
        coord(1:3, rpath%mepatom_id(i)) = rpath%mep_coord(ii:ii+2,replicaid)
        ii = ii + 3
      end do
    end if

    if (rpath%task(1:2) .eq. 'FG') then

      dynvars%step = rpath%neb_cycle
      rpath%neb_output = .true.
      !if (ireplica == rpath%nreplica) rpath%neb_cycle = rpath%neb_cycle + 1

      ! output settings
      !
      if (main_rank .and. (ireplica == 1)) &
       write(MsgOut, '(/,"Iter. ",i5)') rpath%neb_cycle
      if (rpath%eneout_period > 0) then
        if (mod(rpath%neb_cycle, rpath%eneout_period) == 0) then
          minimize%eneout_period = min_ene_period

        else
          minimize%eneout_period = 0

        end if
      end if

      ! optimize surrounding atoms
      !
      if (rpath%mep_partial_opt) then
        if (rpath%neb_cycle > start_minimize) then
          save_step = dynvars%step
          if (rpath%opt_micro) then
            call minimize_micro (output, molecule, enefunc, dynvars, &
                                 minimize, pairlist, boundary, rpath, &
                                 rpath%neb_cycle, ireplica)
          else
              
            if (minimize%method == MinimizeMethodSD) then
              call steepest_descent (output, molecule, enefunc, dynvars, &
                                     minimize, pairlist, boundary)
                
            else if (minimize%method == MinimizeMethodLBFGS) then
              call minimize_lbfgs (output, molecule, enefunc, dynvars, &
                                   minimize, pairlist, boundary)
            end if
          end if
          dynvars%step = save_step

        else
          ! compute the QM/MM energy and force
          !
          enefunc%qmmm%qm_classical = .false.
          call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                              .false.,           &
                              dynvars%coord,     &
                              dynvars%trans,     &
                              dynvars%coord_pbc, &
                              dynvars%energy,    &
                              dynvars%temporary, &
                              dynvars%force,     &
                              dynvars%force_omp, &
                              dynvars%virial,    &
                              dynvars%virial_extern)
        end if
      end if
        
      ! energy and NEB force
      !
      call energy_and_force_replica(dynvars, rpath)

      ! Store the old force by micro iteration
      !
      do i = 1, minimize%num_optatoms_micro
        iatom = minimize%optatom_micro_id(i)
        rpath%micro_force(1:3, i, ireplica) = dynvars%force(1:3, iatom)
      end do

      if (ireplica == nrep_per_proc) call calc_neb_force_replica(rpath)

      ! check convergence and output variables
      !
      if (ireplica == nrep_per_proc) &
       call check_convergence(rpath, conv)
       
      ! output rpath info
      !
      call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                        boundary, rpath)
      if (ireplica == nrep_per_proc) then 
        call print_images(rpath)
        rpath%neb_cycle = rpath%neb_cycle + 1
      end if  

      rpath%neb_output = .false.
    end if

    ! path update
    !
    if (ireplica == nrep_per_proc .and. .not. conv) then
      if (rpath%massWeightCoord) call trans_mass_weight_coord_replica(molecule, rpath)
#ifdef HAVE_MPI_GENESIS
      call mpi_allgather(rpath%mep_coord(1,(replicaid-nrep_per_proc)+1), &
                         (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                         rpath%recv_buff(1,1), &    
                         (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                         mpi_comm_airplane, ierror)      
      rpath%mep_coord = rpath%recv_buff

      call mpi_allgather(rpath%mep_force(1,(replicaid-nrep_per_proc)+1), &
                         (rpath%dimension)*nrep_per_proc, mpi_wp_real,&
                         rpath%recv_buff(1,1), &
                         (rpath%dimension)*nrep_per_proc, mpi_wp_real,&
                         mpi_comm_airplane, ierror)        
      rpath%mep_force = rpath%recv_buff
#endif

      if (main_rank) then
        ! global vector and gradient
        ioff = 0
        do image = 1, rpath%nreplica
          vec(ioff+1:ioff+no3)      =  rpath%mep_coord(1:no3,image)
          gradient(ioff+1:ioff+no3) = -rpath%mep_force(1:no3,image)
          ioff = ioff + no3
        end do

        ! global energy
        neb_energy = 0.0_wp
        do image = 1, rpath%nreplica
          neb_energy = neb_energy + rpath%mep_energy(image)
          if (image > 1) then
            tmp(:) = rpath%mep_coord(:,image) - rpath%mep_coord(:,image-1)
            !YA_20180119 neb_energy = neb_energy + &
            !              0.5_wp * ddot(no3, tmp, 1, tmp, 1) * rpath%k_spring
            neb_energy = neb_energy                         &
                       + 0.5_wp * ddot(no3, tmp, 1, tmp, 1) &
                       * rpath%k_spring * dble(rpath%nreplica - 1)
          end if
        end do

        ! boundary
        if (bnd) then
          if (bnd_qmonly) then
            ! set boundary to QM atoms
            do i = 1, dim
              if (list_bound(i) == 2) then
                upper(i) = vec(i) + maxmove
                lower(i) = vec(i) - maxmove
              end if
            end do
          else
            ! set boundary to all atoms
            do i = 1, dim
              upper(i) = vec(i) + maxmove
              lower(i) = vec(i) - maxmove
            end do
          end if
        end if

        factr = 0.0_wp
        pgtol = 0.0_wp
        if (rpath%verbose .and. main_rank) iprint = 1
        call setulb(dim, ncorr, vec, lower, upper, list_bound, neb_energy, &
                    gradient, factr, pgtol, rpath%wa(:,my_replica_no), &
                    rpath%iwa(:,my_replica_no), rpath%task, iprint, &
                    rpath%csave(my_replica_no), rpath%lsave(:,my_replica_no), &
                    rpath%isave(:,my_replica_no), &
                    rpath%dsave(:,my_replica_no), my_country_no)

        ioff = 0
        do image = 1, rpath%nreplica
          rpath%mep_coord(1:no3,image) = vec(ioff+1:ioff+no3)
          ioff = ioff + no3
        end do
      
      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(rpath%task, 60, mpi_character, 0, mpi_comm_world, ierror)
      call mpi_bcast(rpath%mep_coord(1,1), dim, mpi_wp_real, &
                       0, mpi_comm_world, ierror)
#endif

      ! update pathlength
      !
      rpath%pathlength_prev = rpath%pathlength
      call calc_pathlength(rpath, .false.)

      if (rpath%massWeightCoord) call backtrans_mass_weight_coord_replica(molecule, rpath)

    end if

    if (conv .or. rpath%neb_cycle == (rpath%ncycle+1)) then
      if (ireplica == nrep_per_proc) then
        ! free work arrays
        !
        deallocate(rpath%wa)
        deallocate(rpath%iwa)
        deallocate(rpath%isave)
        deallocate(rpath%dsave)
        deallocate(rpath%lsave)
        deallocate(rpath%csave)
      end if          
    end if

    ! free work arrays
    !
    if (main_rank) then
      deallocate(vec)
      deallocate(gradient)
      deallocate(lower)
      deallocate(upper)
      deallocate(list_bound)
      deallocate(tmp)
    end if

    ! output restart file
    !
    call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                      boundary, rpath)

    return

  end subroutine run_rpath_neb_glbfgs_replica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_neb_force
  !> @brief        compute NEB force
  !! @authors      YA, KY
  !! @param[in]    rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_neb_force(rpath)

    ! formal arguments
    type(s_rpath), target,    intent(inout) :: rpath

    ! local variables
    logical           :: isCI(rpath%nreplica)
    integer           :: natom, image, repid
    real(wp), parameter :: zero = 0.0_wp
    real(wp), pointer :: mep_energy(:)
    real(wp), pointer :: mep_coord(:,:), mep_force(:,:), recv_buff(:,:)
    real(wp), allocatable :: tmp1(:), tmp2(:), tau(:,:)
    real(wp)          :: norm, norm1, norm2, deltaE_f, deltaE_b, fac1, fac2, fac
    real(wp)          :: ddot
    real(wp), allocatable :: recvbuf(:,:)


    natom = rpath%mep_natoms
    repid = my_country_no + 1

    ! Pointers
    !
    mep_energy => rpath%mep_energy
    mep_coord  => rpath%mep_coord
    mep_force  => rpath%mep_force
    recv_buff  => rpath%recv_buff

    ! collect image information
    !
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(rpath%energy, 1, mpi_wp_real, mep_energy, 1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
    call mpi_allgather(mep_coord(1,repid), rpath%dimension, mpi_wp_real,&
                       recv_buff(1,1),     rpath%dimension, mpi_wp_real,&
                       mpi_comm_airplane, ierror)
    mep_coord = recv_buff

    call mpi_allgather(mep_force(1,repid), rpath%dimension, mpi_wp_real,&
                       recv_buff(1,1),     rpath%dimension, mpi_wp_real,&
                       mpi_comm_airplane, ierror)
    mep_force = recv_buff
#endif
    
    ! identify climbing image
    !
    isCI(:) = .FALSE.
    if (rpath%do_cineb) then
      do image = 2, rpath%nreplica - 1
        if (mep_energy(image) > mep_energy(image-1) .and. &
          mep_energy(image) > mep_energy(image+1)) then
          isCI(image) = .TRUE.
          if (main_rank) then
            write(MsgOut, '(" Climbing image: ",i0)') image
          end if
        end if
      end do
    end if

    ! tangential vector
    !
    allocate(tau(rpath%dimension, rpath%nreplica))
    allocate(tmp1(rpath%dimension))
    allocate(tmp2(rpath%dimension))

    tau(:,:) = zero
    do image = 2, rpath%nreplica - 1
      deltaE_f = mep_energy(image+1) - mep_energy(image)
      deltaE_b = mep_energy(image)   - mep_energy(image-1)
      tmp1(:) = mep_coord(:,image+1) - mep_coord(:,image)
      tmp2(:) = mep_coord(:,image)   - mep_coord(:,image-1)
      if (deltaE_b > zero .and. deltaE_f > zero) then
        tau(:,image) = tmp1(:)
      else if (deltaE_b < zero .and. deltaE_f < zero) then
        tau(:,image) = tmp2(:)
      else
        fac1 = max(abs(deltaE_b),abs(deltaE_f))
        fac2 = min(abs(deltaE_b),abs(deltaE_f))
        if (mep_energy(image+1) > mep_energy(image-1)) then
          tau(:,image) = fac1 * tmp1(:) + fac2 * tmp2(:)
        else
          tau(:,image) = fac2 * tmp1(:) + fac1 * tmp2(:)
        end if
      end if
      ! normalization
      norm = ddot(rpath%dimension, tau(1,image), 1, tau(1,image), 1)
      tau(:,image) = tau(:,image) / sqrt(norm)
    end do

    ! NEB force: perpeudicular component
    !
    do image = 2, rpath%nreplica - 1
      norm = ddot(rpath%dimension, tau(1,image), 1, rpath%mep_force(1,image),1)
      if (isCI(image)) then
        mep_force(:,image) = mep_force(:,image) - 2.0_wp * norm * tau(:,image)
      else
        mep_force(:,image) = mep_force(:,image) - norm * tau(:,image)
      end if
    end do

    ! NEB force: spring force
    !
    do image = 2, rpath%nreplica - 1
      if (isCI(image)) cycle
      tmp1(:) = mep_coord(:,image+1) - mep_coord(:,image)
      tmp2(:) = mep_coord(:,image)   - mep_coord(:,image-1)
      norm1 = ddot(rpath%dimension, tmp1, 1, tmp1, 1)
      norm2 = ddot(rpath%dimension, tmp2, 1, tmp2, 1)
      fac = rpath%k_spring * (sqrt(norm1) - sqrt(norm2))
      mep_force(:,image) = mep_force(:,image) + fac * tau(:,image)      
    end do

    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(tau)

    ! clear force for terminals
    !
    if (rpath%fix_terminal) then
      mep_force(:,1) = zero
      mep_force(:,rpath%nreplica) = zero
    end if

    return

  end subroutine calc_neb_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_neb_force_replica
  !> @brief        compute NEB force if npreplica > nproc
  !! @authors      SI
  !! @param[in]    rpath        : rpath info
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_neb_force_replica(rpath)

    ! formal arguments
    type(s_rpath), target,    intent(inout) :: rpath

    ! local variables
    logical           :: isCI(rpath%nreplica)
    integer           :: natom, image, repid
    real(wp), parameter :: zero = 0.0_wp
    real(wp), pointer :: mep_energy(:)
    real(wp), pointer :: mep_coord(:,:), mep_force(:,:), recv_buff(:,:)
    real(wp), allocatable :: tmp1(:), tmp2(:), tau(:,:)
    real(wp)          :: norm, norm1, norm2, deltaE_f, deltaE_b, fac1, fac2, fac
    real(wp)          :: ddot
    real(wp), allocatable :: recvbuf(:)


    natom = rpath%mep_natoms
    repid = my_replica_no

    ! Pointers
    !
    mep_energy => rpath%mep_energy
    mep_coord  => rpath%mep_coord
    mep_force  => rpath%mep_force
    recv_buff  => rpath%recv_buff

    ! collect image information
    !  
#ifdef HAVE_MPI_GENESIS
    allocate(recvbuf(rpath%nreplica))
    !! Energy
    call mpi_allgather(mep_energy((repid-nrep_per_proc)+1), &
                       nrep_per_proc, mpi_wp_real, recvbuf(1), &
                       nrep_per_proc, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
    mep_energy = recvbuf

    !! Coord
    call mpi_allgather(mep_coord(1,(repid-nrep_per_proc)+1), &
                       (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                       recv_buff(1,1), &    
                       (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                       mpi_comm_airplane, ierror)        
    mep_coord = recv_buff

    !! Force
    call mpi_allgather(mep_force(1,(repid-nrep_per_proc)+1), &
                       (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                       recv_buff(1,1), &    
                       (rpath%dimension)*nrep_per_proc, mpi_wp_real, &
                       mpi_comm_airplane, ierror)   
    mep_force = recv_buff
    deallocate(recvbuf)
#endif
    
    ! identify climbing image
    !
    isCI(:) = .FALSE.
    if (rpath%do_cineb) then
      do image = 2, rpath%nreplica - 1
        if (mep_energy(image) > mep_energy(image-1) .and. &
          mep_energy(image) > mep_energy(image+1)) then
          isCI(image) = .TRUE.
          if (main_rank) then
            write(MsgOut, '(" Climbing image: ",i0)') image
          end if
        end if
      end do
    end if

    ! tangential vector
    !
    allocate(tau(rpath%dimension, rpath%nreplica))
    allocate(tmp1(rpath%dimension))
    allocate(tmp2(rpath%dimension))

    tau(:,:) = zero
    do image = 2, rpath%nreplica - 1
      deltaE_f = mep_energy(image+1) - mep_energy(image)
      deltaE_b = mep_energy(image)   - mep_energy(image-1)
      tmp1(:) = mep_coord(:,image+1) - mep_coord(:,image)
      tmp2(:) = mep_coord(:,image)   - mep_coord(:,image-1)
      if (deltaE_b > zero .and. deltaE_f > zero) then
        tau(:,image) = tmp1(:)
      else if (deltaE_b < zero .and. deltaE_f < zero) then
        tau(:,image) = tmp2(:)
      else
        fac1 = max(abs(deltaE_b),abs(deltaE_f))
        fac2 = min(abs(deltaE_b),abs(deltaE_f))
        if (mep_energy(image+1) > mep_energy(image-1)) then
          tau(:,image) = fac1 * tmp1(:) + fac2 * tmp2(:)
        else
          tau(:,image) = fac2 * tmp1(:) + fac1 * tmp2(:)
        end if
      end if
      ! normalization
      norm = ddot(rpath%dimension, tau(1,image), 1, tau(1,image), 1)
      tau(:,image) = tau(:,image) / sqrt(norm)
    end do

    ! NEB force: perpeudicular component
    !
    do image = 2, rpath%nreplica - 1
      norm = ddot(rpath%dimension, tau(1,image), 1, rpath%mep_force(1,image),1)
      if (isCI(image)) then
        mep_force(:,image) = mep_force(:,image) - 2.0_wp * norm * tau(:,image)
      else
        mep_force(:,image) = mep_force(:,image) - norm * tau(:,image)
      end if
    end do

    ! NEB force: spring force
    !
    do image = 2, rpath%nreplica - 1
      if (isCI(image)) cycle
      tmp1(:) = mep_coord(:,image+1) - mep_coord(:,image)
      tmp2(:) = mep_coord(:,image)   - mep_coord(:,image-1)
      norm1 = ddot(rpath%dimension, tmp1, 1, tmp1, 1)
      norm2 = ddot(rpath%dimension, tmp2, 1, tmp2, 1)
      fac = rpath%k_spring * (sqrt(norm1) - sqrt(norm2))
      mep_force(:,image) = mep_force(:,image) + fac * tau(:,image)      
    end do

    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(tau)

    ! clear force for terminals
    !
    if (rpath%fix_terminal) then
      mep_force(:,1) = zero
      mep_force(:,rpath%nreplica) = zero
    end if

    return

  end subroutine calc_neb_force_replica

end module at_rpath_mep_mod

