!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_rpath_fep_mod
!> @brief   Free-energy perturbation along given MEP
!! @authors Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_rpath_fep_mod

  use at_energy_mod
  use at_md_leapfrog_mod
  use at_md_vverlet_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_dynvars_mod
  use at_ensemble_str_mod
  use at_rpath_str_mod
  use at_rpath_mep_mod
  use at_constraints_str_mod
  use at_boundary_str_mod
  use at_boundary_mod
  use at_pairlist_str_mod
  use at_pairlist_mod
  use at_enefunc_str_mod
  use at_output_str_mod
  use at_output_mod
  use fileio_rst_mod
  use fileio_rstmep_mod
  use molecules_str_mod
  use select_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_rpath_fep
  public  :: run_rpath_fep
  private :: calc_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath_fep
  !> @brief        setup RPATH for FEP
  !! @authors      YA, KY
  !! @param[in]    fep_period : the period of FEP calculation
  !! @param[in]    mepatm_select_index : index of MEP atoms
  !! @param[in]    sel_info   : selector input information
  !! @param[in]    rst        : restart information
  !! @param[in]    rstmep     : rstmep data
  !! @param[in]    molecule   : molecule information
  !! @param[in]    qmmm       : qmmm     information
  !! @param[in]    dynamics   : dynamics information
  !! @param[inout] dynvars    : dynamic variables
  !! @param[inout] rpath      : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath_fep(fep_period, mepatm_select_index, sel_info, rst, &
                             rstmep, molecule, qmmm, dynamics, dynvars, rpath)

    ! formal arguments
    integer,                 intent(in)    :: fep_period
    character(256),          intent(in)    :: mepatm_select_index
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_rst),     target, intent(in)    :: rst
    type(s_rstmep),          intent(in)    :: rstmep
    type(s_dynamics),target, intent(in)    :: dynamics
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm),            intent(inout) :: qmmm
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_rpath),   target, intent(inout) :: rpath

    ! local variables
    integer                  :: ns, i, j, k, m, ii, iatom
    integer                  :: shift, next, id, atomid
    integer                  :: replicaid, errorcode, my_rank, ierr


    ! setup FEP period
    !
    if (fep_period > 0) then
      if (mod(dynamics%nsteps, fep_period) /= 0) then
        call error_msg('Setup_Rpath> mod(nsteps,fep_period) must be zero')
      else
        rpath%equilibration_only = .false.
        rpath%fep_period = fep_period
        rpath%ncycle = dynamics%nsteps/rpath%fep_period
      end if
    else if (fep_period == 0) then
      rpath%equilibration_only = .true.
      rpath%fep_period = dynamics%nsteps
      rpath%ncycle = 1
    else
      call error_msg('Setup_Rpath> rpath_period should be non-negative integer')
    end if


    if (.not. qmmm%do_qmmm) then
      rpath%esp_energy = .false.
      rpath%esp_md     = .false.
    end if

    ! define atoms
    !
    call setup_mepatoms(mepatm_select_index, sel_info, molecule, rpath)

    if (qmmm%do_qmmm) then
      call setup_mepatoms_qmmm(molecule, qmmm, rpath)
    end if

    rpath%dimension = rpath%mep_natoms * 3
    if (qmmm%do_qmmm) then
      call alloc_rpath_fep(rpath, rpath%dimension, rpath%nreplica, &
                           qmmm%qm_natoms)
    else
      call alloc_rpath_fep(rpath, rpath%dimension, rpath%nreplica)
    end if

    ! set MEP variables
    !
    replicaid = my_country_no + 1

    if (.not. rstmep%restart) &
      call error_msg("Setup_Rpath_FEP> rstmep is needed for FEP.")
    if (rstmep%mep_natoms /= rpath%mep_natoms) &
      call error_msg("Setup_Rpath_FEP> The number of MEP atoms don't match &
                     &in the ctrlfile and rstmep.")
    rpath%mep_coord(:,replicaid) = rstmep%mep_coord

    if (qmmm%do_qmmm) then
      if (rpath%esp_energy .or. rpath%esp_md) then
        if (.not. allocated(rstmep%qm_charge)) &
          call error_msg("Setup_Rpath_FEP> QM charges needed for esp_energy &
                         &and/or esp_md are not found in rstmep!") 

        rpath%qm_charge(:,replicaid) = rstmep%qm_charge
        rpath%qm_energy(replicaid)   = rstmep%qm_energy
      end if

    end if 

#ifdef HAVE_MPI_GENESIS
    ! broadcast rpath information
    do i = 1, rpath%nreplica
      call mpi_bcast(rpath%mep_coord(1,i), rpath%dimension, mpi_wp_real, & 
                     nproc_country*(i-1), mpi_comm_world, ierr)
    end do
    if (qmmm%do_qmmm) then
      do i = 1, rpath%nreplica
        call mpi_bcast(rpath%qm_charge(1,i), qmmm%qm_natoms, mpi_wp_real, &
                     nproc_country*(i-1), mpi_comm_world, ierr)
        call mpi_bcast(rpath%qm_energy(i), 1, mpi_wp_real, &
                     nproc_country*(i-1), mpi_comm_world, ierr)
      end do
    end if
#endif

    ! initialize FEP
    !
    rpath%num_fep = 0
    rpath%pfunc_f = 0.0_wp
    rpath%pfunc_b = 0.0_wp
    if (rstmep%num_fep > 0) then
      rpath%num_fep = rstmep%num_fep
      rpath%pfunc_f = rstmep%pfunc(1)
      rpath%pfunc_b = rstmep%pfunc(2)
    end if

    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Rpath_FEP> Rpath information'
      write(MsgOut,'(A)') ''

      if (rstmep%num_fep > 0) then
        write(MsgOut, '(a)') ' '
        write(MsgOut,'(a)') 'Setup_Rpath_FEP> Restart retrieved from mep files.'
        write(MsgOut,'(a,i0)') 'Setup_Rpath_FEP> num_fep = ', rstmep%num_fep
        write(MsgOut, '(a)') ' '
      end if

      ! Punch out atom info.
      !
      write(MsgOut,'(a)') "Setup_Rpath_FEP> Atoms involved in MEP"
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
      write(MsgOut,'(a,i0)') "  number of atoms in MEP = ", rpath%mep_natoms
      write(MsgOut, '(a)') ' '

    end if

    return

  end subroutine setup_rpath_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_fep
  !> @brief        run rpath FEP
  !! @authors      YA, KY
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

  subroutine run_rpath_fep(output, molecule, enefunc, dynvars, dynamics,   &
                           pairlist, boundary, constraints, ensemble, rpath)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_rpath),            intent(inout) :: rpath

    ! local variables
    integer                   :: i, j, k
    integer                   :: n, niter
    integer                   :: iloop_start, iloop_end
    integer                   :: replicaid, errorcode, ierr
    integer                   :: ii, iatom
    real(wp)                  :: dene_qm
    real(wp),     allocatable :: work(:,:)


    ! Open output files
    !
    call open_output(output)
    DynvarsOut = output%logunit

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
      call update_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                           dynvars%coord_pbc, pairlist)
    end if

    if (rpath%esp_energy .or. rpath%esp_md) then
      replicaid = my_country_no + 1
      enefunc%qmmm%qm_charge = rpath%qm_charge(:,replicaid) 
    end if

    !---------
    ! MD step
    do n = 1, rpath%ncycle
      if (.not. rpath%equilibration_only) &
        rpath%num_fep = rpath%num_fep + 1

      dynamics%istart_step  = (n-1)*rpath%fep_period + 1
      dynamics%iend_step    =  n   *rpath%fep_period
      dynamics%initial_time = dynvars%time

      ! MD main loop
      if (rpath%esp_md) then
        enefunc%qmmm%qm_classical = .true.
      else
        enefunc%qmmm%qm_classical = .false.
        if(rpath%equilibration_only) then
          write(enefunc%qmmm%qmindex,'("equ",i0)') n
        else
          write(enefunc%qmmm%qmindex,'("fep",i0)') n
        end if
      end if
      if (dynamics%integrator == IntegratorLEAP) then
        call leapfrog_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
                               
      else if (dynamics%integrator == IntegratorVVER) then
        call vverlet_dynamics (output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      end if
      if (.not. rpath%esp_md) enefunc%qmmm%qmindex = ''

      if (rpath%equilibration_only) then
        ! output restart
        call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                          boundary, rpath)

      else
        ! calc FEP energies
        call calc_fep(molecule, enefunc, dynvars, rpath, pairlist, &
                      boundary, ensemble)

        ! output restart
        call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                          boundary, rpath)

        ! print summary
        allocate(work(rpath%nreplica,3))
     
        ! Forward
#ifdef HAVE_MPI_GENESIS
        call mpi_gather(rpath%dfene_f, 1, mpi_wp_real, work(1,1), 1, &
                        mpi_wp_real, 0, mpi_comm_airplane, ierror)
#endif
        work(1,2) = 0.0_wp
        do i = 2, rpath%nreplica
          work(i,2) = work(i-1,2) + work(i-1,1)
        end do

        !dbg  if(main_rank) then
        !dbg    do i = 1, rpath%nreplica
        !dbg      write(MsgOut,'(2F20.10)') work(i,1), work(i,2)
        !dbg    end do
        !dbg    write(MsgOut,*)
        !dbg end if

        ! Backward
#ifdef HAVE_MPI_GENESIS
        call mpi_gather(rpath%dfene_b, 1, mpi_wp_real, work(1,1), 1, &
                        mpi_wp_real, 0, mpi_comm_airplane, ierror)
#endif
        work(1,3) = 0.0_wp
        do i = 2, rpath%nreplica
          work(i,3) = work(i-1,3) - work(i,1)
        end do

        !dbg if(main_rank) then
        !dbg   do i = 1, rpath%nreplica
        !dbg     write(MsgOut,'(2F20.10)') work(i,1), work(i,3)
        !dbg   end do
        !dbg end if

        if (main_rank) then
          write(MsgOut, '("Rpath_FEP> FEP iteration: ",i8)') rpath%num_fep
          write(MsgOut, '(/," +++ QM/MM-FEP summary +++")')

          if (rpath%esp_energy) then
            write(MsgOut, '(/,"Image ", &
                              "   dEqm     (kcal/mol)", &
                              "   dFqmmm_f (kcal/mol)", &
                              "   dFqmmm_b (kcal/mol)"  &
                              "   dF_f     (kcal/mol)", &
                              "   dF_b     (kcal/mol)")')
            write(MsgOut, '(116("-"))')
            do i = 1, rpath%nreplica

              dene_qm = rpath%qm_energy(i) - rpath%qm_energy(1)
              write(MsgOut, '(i5,X,5F22.10)') i, dene_qm, work(i,2), work(i,3), &
                                      dene_qm + work(i,2), dene_qm + work(i,3)
            end do
            write(MsgOut, '(116("-"),/)')

          else
            write(MsgOut, '(/,"Image ", &
                              "   dF_f     (kcal/mol)", &
                              "   dF_b     (kcal/mol)")')
            write(MsgOut, '(50("-"))')
            do i = 1, rpath%nreplica

              write(MsgOut, '(i5,X,2F22.10)') i, work(i,2), work(i,3)
            end do
            write(MsgOut, '(50("-"),/)')

          end if

        end if
     
        deallocate(work)

        !dbg call mpi_barrier(mpi_comm_world, ierr)
        !dbg call mpi_abort(mpi_comm_world, errorcode, ierr)

      end if

    end do

    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_rpath_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_fep
  !> @brief        calculate the free-energy using FEP
  !! @authors      YA, KY
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] dynvars  : dynamical variables information
  !! @param[inout] rpath    : RPATH information
  !! @param[inout] pairlist : non-bond pair list information
  !! @param[inout] boundary : boundary conditions information
  !! @param[inout] ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_fep(molecule, enefunc, dynvars, rpath, pairlist, &
                      boundary, ensemble)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_rpath),            intent(inout) :: rpath
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_ensemble),         intent(inout) :: ensemble

    ! local variables
    type(s_qmmm), pointer :: qmmm
    real(wp),     pointer :: coord(:,:)
    real(wp), allocatable :: mep_coord_tmp(:,:)

    integer    :: replicaid
    integer    :: i, ii, iatom
    real(wp)   :: beta, dene, e1f, e1b, e0
    character(1) :: fcb

    
    coord => dynvars%coord
    qmmm  => enefunc%qmmm

    replicaid = my_country_no + 1

    allocate(mep_coord_tmp(3,rpath%mep_natoms))
    beta = 1.0_wp / (KBOLTZ * ensemble%temperature)

    ! Save current MEP atom coordinates
    !
    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      mep_coord_tmp(1:3,i) =  coord(1:3,iatom)
    end do

    ! Energy of neighboring image (forward)
    fcb = 'f'
    if (replicaid /= rpath%nreplica) then
      e1f = get_energy(replicaid+1)
    else
      e1f = 0.0_wp
    end if

    ! Energy of neighboring image (backward)
    fcb = 'b'
    if (replicaid /= 1) then
      e1b = get_energy(replicaid-1)
    else
      e1b = 0.0_wp
    end if

    ! Energy of this image
    fcb = 'c'
    e0 = get_energy(replicaid)

    ! Retrieve MEP atom coordinates
    !
    do i = 1, rpath%mep_natoms
      iatom = rpath%mepatom_id(i)
      coord(1:3,iatom) = mep_coord_tmp(1:3,i)
    end do
    deallocate(mep_coord_tmp)

    if (replicaid == 1) then
      rpath%pfunc_f = rpath%pfunc_f + exp(-beta * (e1f - e0))
      rpath%dfene_f = -log(rpath%pfunc_f / dble(rpath%num_fep)) / beta

    else if (replicaid == rpath%nreplica) then
      rpath%pfunc_b = rpath%pfunc_b + exp(-beta * (e1b - e0))
      rpath%dfene_b = -log(rpath%pfunc_b / dble(rpath%num_fep)) / beta

    else
      rpath%pfunc_f = rpath%pfunc_f + exp(-beta * (e1f - e0))
      rpath%dfene_f = -log(rpath%pfunc_f / dble(rpath%num_fep)) / beta

      rpath%pfunc_b = rpath%pfunc_b + exp(-beta * (e1b - e0))
      rpath%dfene_b = -log(rpath%pfunc_b / dble(rpath%num_fep)) / beta
    end if

    contains

    function get_energy(repID)
      real(wp) :: get_energy
      integer  :: repID
      logical  :: ene_only_org
      
      character(256) :: folder,basename
      character      :: num*5
      logical        :: savefile

      ii = 0
      do i = 1, rpath%mep_natoms
        iatom = rpath%mepatom_id(i)
        coord(1:3,iatom) = rpath%mep_coord(ii+1:ii+3,repID)
        ii = ii + 3
      end do
      call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
          dynvars%coord_pbc, pairlist)

      if (qmmm%do_qmmm) then
        ene_only_org  = qmmm%ene_only
        qmmm%ene_only = .true.

        if(rpath%esp_energy) then
          qmmm%qm_charge = rpath%qm_charge(:,repID) 
          qmmm%qm_classical = .true.

        else
          !folder   = 'qmmm_fep'
          !write(num,'(i5.5)') rpath%num_fep
          !basename = 'fep'//num//'_'//fcb

          !if (mod(rpath%num_fep,rpath%qmsave_period) == 0) then
          !  savefile = .true.
          !else
          !  savefile = .false.
          !end if
          !call set_runtime_qmmm(enefunc%qmmm, folder, basename, savefile)
          qmmm%qm_classical = .false.

        end if
      end if

      call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                          enefunc%nonb_limiter,  &
                          dynvars%coord,         &
                          dynvars%trans,         &
                          dynvars%coord_pbc,     &
                          dynvars%energy,        &
                          dynvars%temporary,     &
                          dynvars%force,         &
                          dynvars%force_omp,     &
                          dynvars%virial,        &
                          dynvars%virial_extern)

      get_energy = dynvars%energy%total - &
        sum(dynvars%energy%restraint(1:enefunc%num_restraintfuncs))

      if (qmmm%do_qmmm) qmmm%ene_only = ene_only_org

    end function get_energy

  end subroutine calc_fep

end module at_rpath_fep_mod

