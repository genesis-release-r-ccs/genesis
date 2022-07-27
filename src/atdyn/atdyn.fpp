!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  ATDYN
!! @brief   Molecular Dynamics Simulation of Biomolecules using 
!!          Atom Decomposition Scheme
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ), Chigusa Kobayashi (CK),
!!          Yasuhiro Matsunaga (YM), Takashi Imai (TI), Takao Yoda (TY),
!!          Yuji Sugita (YS), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

program atdyn

  use at_setup_mpi_mod
  use at_remd_mod
  use at_rpath_mod
  use at_morph_mod
  use at_minimize_mod
  use at_vibration_mod
  use at_dynamics_mod
  use at_setup_atdyn_mod
  use at_control_mod
  use at_energy_mod
  use at_energy_pme_mod
  use at_qmmm_mod
  use at_minimize_str_mod
  use at_vibration_str_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_ensemble_str_mod
  use at_output_str_mod
  use at_remd_str_mod
  use at_morph_str_mod
  use at_rpath_str_mod
  use at_constraints_str_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use hardwareinfo_mod
  use timers_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none

  integer             :: genesis_run_mode
  integer, parameter  :: GenesisMD    = 1
  integer, parameter  :: GenesisMIN   = 2
  integer, parameter  :: GenesisREMD  = 3
  integer, parameter  :: GenesisRPATH = 4
  integer, parameter  :: GenesisVIB   = 5
  integer, parameter  :: GenesisBD    = 6
  integer, parameter  :: GenesisMORPH = 7

  ! local variables
  character(MaxFilename)      :: ctrl_filename
  type(s_ctrl_data)           :: ctrl_data
  type(s_molecule)            :: molecule
  type(s_enefunc)             :: enefunc
  type(s_dynvars)             :: dynvars
  type(s_pairlist)            :: pairlist
  type(s_boundary)            :: boundary
  type(s_constraints)         :: constraints
  type(s_ensemble)            :: ensemble
  type(s_dynamics)            :: dynamics
  type(s_minimize)            :: minimize
  type(s_output)              :: output
  type(s_remd)                :: remd
  type(s_rpath)               :: rpath
  type(s_vibration)           :: vibration
  type(s_morph)               :: morph
  integer                     :: omp_get_max_threads


#ifdef HAVE_MPI_GENESIS
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world, my_world_rank, ierror)
  call mpi_comm_size(mpi_comm_world, nproc_world,   ierror)
  main_rank = (my_world_rank == 0)
#ifdef QSIMULATE
  nprow = 1
  npcol = nproc_world
  call sl_init(scalapack_context, nprow, npcol)
  call blacs_gridinfo(scalapack_context, nprow, npcol, my_prow, my_pcol)
#endif
#else
  my_world_rank = 0
  nproc_world   = 1
  main_rank     = .true.
#endif

#ifdef OMP
  nthread = omp_get_max_threads()
#else
  nthread = 1
#endif

#ifdef __PGI
   call error_msg('Atdyn> atdyn with PGI is not allowed at this moment')
#endif
 
  ! show usage
  !
  call usage(ctrl_filename)

  ! get run mode from control file
  !
  call get_genesis_mode(ctrl_filename, genesis_run_mode)

  ! run genesis
  !
  call atomic_decomposition_genesis(ctrl_filename, genesis_run_mode)

#ifdef HAVE_MPI_GENESIS
#ifdef QSIMULATE
  call blacs_exit(1)
#endif
  call mpi_finalize(ierror)
#endif

  stop

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_genesis_mode
  !> @brief        get genesis run mode
  !! @authors      TM
  !! @param[in]    ctrl_filename    : control file name
  !! @param[out]   genesis_run_mode : run MD, MIN, REMD, RPATH, VIB
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_genesis_mode(ctrl_filename, genesis_run_mode)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_filename
    integer,                 intent(inout) :: genesis_run_mode

    integer                                :: nthread


    if (find_ctrlfile_section(ctrl_filename, 'REMD')) then
      genesis_run_mode = GenesisREMD

    else if (find_ctrlfile_section(ctrl_filename, 'RPATH')) then
      genesis_run_mode = GenesisRPATH

    else if (find_ctrlfile_section(ctrl_filename, 'DYNAMICS')) then
      genesis_run_mode = GenesisMD

    else if (find_ctrlfile_section(ctrl_filename, 'MINIMIZE')) then
      genesis_run_mode = GenesisMIN

    else if (find_ctrlfile_section(ctrl_filename, 'VIBRATION')) then
      genesis_run_mode = GenesisVIB

    else if (find_ctrlfile_section(ctrl_filename, 'MORPH')) then
      genesis_run_mode = GenesisMORPH

    else
      call error_msg('Get_Genesis_Mode> ERROR: Unknown control file format.')

    end if


    return

  end subroutine get_genesis_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    atomic_decomposition_genesis
  !> @brief        run genesis using atomic decomposition scheme
  !! @authors      TM
  !! @param[in]    ctrl_filename   : control file name
  !! @param[in]    genesis_run_mod : run MD, MIN, REMD, RPATH
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine atomic_decomposition_genesis(ctrl_filename, genesis_run_mode)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_filename
    integer,                 intent(in)    :: genesis_run_mode

    ! local variables
    character(256)         :: folder, basename
    character(5)           :: num

    ! set timer
    !
    call timer(TimerTotal, TimerOn)

    ! [Step0] Architecture & Compiler information
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP0] Architecture and Compiler Information'
      write(MsgOut,'(A)') ' '

      call hw_information

    end if


    ! [Step1] Read control file
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP1] Read Control Parameters'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD)

      call control_md  (ctrl_filename, ctrl_data)

    case (GenesisMIN)

      call control_min (ctrl_filename, ctrl_data)

    case (GenesisREMD)

      call control_remd(ctrl_filename, ctrl_data)

    case (GenesisRPATH)

      call control_rpath(ctrl_filename, ctrl_data)

    case (GenesisVIB)

      call control_vib(ctrl_filename, ctrl_data)

    case (GenesisMORPH)

      call control_morph(ctrl_filename, ctrl_data)

    end select

    ! [Step2] Setup MPI
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2] Setup MPI'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD,GenesisMIN,GenesisMORPH)

      call setup_mpi_md  (ctrl_data%ene_info)

    case (GenesisREMD)

      call setup_mpi_remd(ctrl_data%ene_info, ctrl_data%rep_info)

    case (GenesisRPATH)

      call setup_mpi_rpath(ctrl_data%ene_info, ctrl_data%rpath_info)

    case (GenesisVIB)

      call setup_mpi_vib(ctrl_data%ene_info, ctrl_data%vib_info)

    end select


    ! [Step3] Set relevant variables and structures 
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP3] Set Relevant Variables and Structures'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD)

      call setup_atdyn_md  (ctrl_data, output, molecule, enefunc,     &
                            pairlist, dynvars, dynamics, constraints, &
                            ensemble, boundary)
    case (GenesisMIN)

      call setup_atdyn_min (ctrl_data, output, molecule, enefunc,     &
                            pairlist, dynvars, minimize, boundary)

    case (GenesisREMD)

      call setup_atdyn_remd(ctrl_data, output, molecule, enefunc,     &
                            pairlist, dynvars, dynamics, constraints, &
                            ensemble, boundary, remd)

    case (GenesisRPATH)

      call setup_atdyn_rpath(ctrl_data, output, molecule, enefunc,    &
                            pairlist, dynvars, dynamics, constraints, &
                            ensemble, boundary, rpath, minimize)


    case (GenesisVIB)

      call setup_atdyn_vib (ctrl_data, output, molecule, enefunc,     &
                            pairlist, dynvars, vibration, boundary)

    case (GenesisMORPH)

      call setup_atdyn_morph (ctrl_data, output, molecule, enefunc,   &
                            pairlist, dynvars, morph, boundary)
    end select


    ! [Step4] Compute single point energy for molecules
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP4] Compute Single Point Energy for Molecules'
      write(MsgOut,'(A)') ' '
    end if

    if (enefunc%qmmm%do_qmmm .and. &
        (genesis_run_mode == GenesisMIN .or. genesis_run_mode == GenesisVIB)) then
      if (enefunc%qmmm%qm_debug) then
        call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            enefunc%nonb_limiter,  &
                            dynvars%coord,     &
                            dynvars%trans,     &
                            dynvars%coord_pbc, &
                            dynvars%energy,    &
                            dynvars%temporary, &
                            dynvars%force,     &
                            dynvars%force_omp, &
                            dynvars%virial,    &
                            dynvars%virial_extern)

        call output_energy(dynvars%step, enefunc, dynvars%energy)

        ! reset count
        enefunc%qmmm%qm_count = 0
        if (enefunc%qmmm%is_qm_charge) then
          enefunc%qmmm%qm_classical = .true.
          call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                              enefunc%nonb_limiter,  &
                              dynvars%coord,     &
                              dynvars%trans,     &
                              dynvars%coord_pbc, &
                              dynvars%energy,    &
                              dynvars%temporary, &
                              dynvars%force,     &
                              dynvars%force_omp, &
                              dynvars%virial,    &
                              dynvars%virial_extern)

          call output_energy(dynvars%step, enefunc, dynvars%energy)
          enefunc%qmmm%qm_classical = .false.
        end if

      else
        if (main_rank) then
          write(MsgOut,'(A)') 'SKIPPED: Energy calculation is omitted when QMMM is performed'
          write(MsgOut,'(A)') ' '
        end if

      end if

    else
      call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                          enefunc%nonb_limiter,  &
                          dynvars%coord,     &
                          dynvars%trans,     &
                          dynvars%coord_pbc, &
                          dynvars%energy,    &
                          dynvars%temporary, &
                          dynvars%force,     &
                          dynvars%force_omp, &
                          dynvars%virial,    &
                          dynvars%virial_extern, &
                          constraints)

      call output_energy(dynvars%step, enefunc, dynvars%energy)

    end if

    ! [Step5] Perform MD/REMD/RPATH simulation or Energy minimization
    !
    call timer(TimerDynamics, TimerOn)

    select case (genesis_run_mode)

    case (GenesisMD)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Molecular Dynamics Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_md(output, molecule, enefunc, dynvars, dynamics, &
                  pairlist, boundary, constraints, ensemble)

    case (GenesisMIN)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Energy Minimization'
        write(MsgOut,'(A)') ' '
      end if

      call run_min(output, molecule, enefunc, dynvars, minimize, &
                   pairlist, boundary)

    case (GenesisREMD)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Replica-Exchange MD Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_remd(output, molecule, enefunc, dynvars, dynamics, &
                    pairlist, boundary, constraints, ensemble, remd)

    case (GenesisRPATH)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Replica Path MD Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_rpath(ctrl_data%inp_info, ctrl_data%out_info, &
                     ctrl_data%qmmm_info, ctrl_data%rpath_info, &
                     ctrl_data%sel_info, &
                     output, molecule, enefunc, &
                     dynvars, minimize, dynamics, pairlist, boundary, &
                     constraints, ensemble, rpath)  

    case (GenesisVIB)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Vibrational Analysis'
        write(MsgOut,'(A)') ' '
      end if

      call run_vib(molecule, enefunc, dynvars, vibration,  &
                   output, pairlist, boundary)

    case (GenesisMORPH)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Morphing Simulation'
        write(MsgOut,'(A)') ' '
      end if

       call run_morph(output, molecule, enefunc, dynvars, morph, &
                              boundary, pairlist)

    end select

    call timer(TimerDynamics, TimerOff)


    ! [Step6] Deallocate arrays
    !
    if (main_rank) then
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '[STEP6] Deallocate Arrays'
      write(MsgOut,'(A)') ' '
    end if

    call dealloc_pme
    call dealloc_remd_all       (remd)
    call dealloc_rpath_all       (rpath)
    call dealloc_constraints_all(constraints)
    call dealloc_boundary_all   (boundary)
    call dealloc_pairlist_all   (pairlist)
    call dealloc_energy_all     (dynvars%energy)
    call dealloc_dynvars_all    (dynvars)
    call dealloc_enefunc_all    (enefunc)
    call dealloc_molecules_all  (molecule)

    if(enefunc%qmmm%do_qmmm) call qmmm_finalize(enefunc%qmmm)

    call timer(TimerTotal, TimerOff)


    ! output process time
    !
    call output_time

    return

  end subroutine atomic_decomposition_genesis

end program atdyn
