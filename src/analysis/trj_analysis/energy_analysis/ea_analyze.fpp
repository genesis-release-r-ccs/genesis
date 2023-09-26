!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ea_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Motoshi Kamiya (MK)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ea_analyze_mod

  use ea_option_str_mod

  use at_enefunc_str_mod
  use at_pairlist_str_mod
  use at_pairlist_mod
  use at_dynvars_str_mod
  use at_boundary_str_mod
  use at_boundary_mod
  use at_energy_mod
  use at_remd_str_mod
  use at_ensemble_str_mod
  use at_remd_mod

  use fileio_trj_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: analyze

contains

  !======1=========2=========3=========4=========5=========6=========7=========8  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      MK
  !! @param[in]    output     : output information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] enefunc    : enefunc information
  !! @param[inout] pairlist   : pairlist information
  !! @param[inout] dynvars    : dynvars information
  !! @param[inout] boundary   : boundary information
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(output, molecule, enefunc, pairlist, dynvars, boundary, &
                     ensemble, remd, trj_list, trajectory, option)

    ! formal arguments
    type(s_output),     intent(in)    :: output
    type(s_molecule),   intent(inout) :: molecule
    type(s_enefunc),    intent(inout) :: enefunc
    type(s_pairlist),   intent(inout) :: pairlist
    type(s_dynvars),    intent(inout) :: dynvars
    type(s_boundary),   intent(inout) :: boundary
    type(s_ensemble),   intent(inout) :: ensemble
    type(s_remd),       intent(inout) :: remd
    type(s_trj_list),   intent(in)    :: trj_list
    type(s_trajectory), intent(inout) :: trajectory
    type(s_option),     intent(inout) :: option

    ! local variables
    type(s_trj_file) :: trj_in
    integer ::  num_trjfiles, natoms, ene_out, rem_inp
    integer :: ifile, istep, i, cur_step, rstep, ierr
    integer :: parmsetid, prev_parmsetid, nreplicas
    logical :: has_rem, has_rep_trj, mbar
    real(wp) :: tot
    real(wp), allocatable :: tot_mbar(:), tot_rest_comp(:)
    real(wp) :: c1, c2, c3, rc1, rc2, rc3, d12, d13
    real(wp) :: euu, euv, evv


    if ( option%check_only ) return

    ! open output file if specified
    !
    if ( output%enefile /= '' ) then
      if (main_rank) then
        call open_file(ene_out, output%enefile, IOFileOutputNew)
      end if
    else
      ene_out = MsgOut
    end if

    has_rem     = .false.
    has_rep_trj = .false.

    parmsetid      = 0
    prev_parmsetid = 0
    if ( option%remfile /= '' .and. .not. option%mbar ) then
      ! note: every processes read the file
      !
      call open_file(rem_inp, option%remfile, IOFileInput)
      if ( rem_inp /= 0 ) then
        has_rem     = .true.
        has_rep_trj = .true.
        read(rem_inp,*) rstep, prev_parmsetid
        !write(MsgOut,*) "Remfile found.", rstep, prev_parmsetid
        parmsetid = prev_parmsetid
      end if
      call assign_condition(option%tgt_parmid, remd, ensemble, &
                            molecule, enefunc)
    end if

    mbar      = option%mbar
    !nreplicas = remd%nreplicas(1)
    nreplicas = 1
    if ( option%mbar .or. option%rest_component ) has_rem = .true.

    allocate( tot_mbar( max( 1, nreplicas ) ) )
    if ( mbar .and. main_rank ) then
      write(MsgOut,'("nreplicas = ", i6 )') nreplicas
    end if

    natoms = molecule%num_atoms

    ! analysis loop
    !
    num_trjfiles = size(trj_list%md_steps)
    if (main_rank) then
      write(ene_out,'("# num trj files = ",i10)') num_trjfiles
      write(ene_out,'("# target rep id = ",i10)') option%tgt_parmid
      write(ene_out,'("#")')
    end if

    ! header
    if (main_rank) then
      write(ene_out,'("# Column 1: STEP")')
      if ( option%component .and. .not. has_rem ) then
        write(ene_out,'("# Column 2: TOTAL POTENTIAL ENERGY")')
        write(ene_out,'("# Column 3: BOND ENERGY")')
        write(ene_out,'("# Column 4: ANGLE ENERGY")')
        write(ene_out,'("# Column 5: UREY-BRADLEY ENERGY")')
        write(ene_out,'("# Column 6: DIHEDRAL ENERGY")')
        write(ene_out,'("# Column 7: IMPROPER ENERGY")')
        write(ene_out,'("# Column 8: CMAP ENERGY")')
        write(ene_out,'("# Column 9: ELECTROSTATIC ENERGY")')
        write(ene_out,'("# Column 10: LJ ENERGY")')
        write(ene_out,'("# Column 11: RESTRAINTS ENERGY")')
      else if ( mbar ) then
        do i = 1, nreplicas
          write(ene_out,'("# Column ",i3,": POTENTIAL ENERGY")') i + 1
        end do
      else if ( has_rem ) then
        write(ene_out,'("# Column 2: REPLICA PARAMETER INDEX")')
        write(ene_out,'("# Column 3: TOTAL POTENTIAL ENERGY")')
      else
        write(ene_out,'("# Column 2: TOTAL POTENTIAL ENERGY")')
      end if
      write(ene_out,'("# unit of output energy is kcal/mol")')
      write(ene_out,'("#")')
    end if

    do ifile = 1, num_trjfiles

      ! open trajectory file
      ! NOTE: every processes read the file!
      !
      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(ifile)

        cur_step = trj_list%mdout_periods(ifile) * istep

        ! read trajectory
        !   coordinates of one MD snapshot are saved in trajectory%coord
        !
        call read_trj(trj_in, trajectory)

        if ( mod(istep, trj_list%ana_periods(ifile)) /= 0 ) cycle

        ! copy coordinate and box size to dynvars
        !
        dynvars%coord(1:3,1:natoms) = trajectory%coord(1:3,1:natoms)

        if ( boundary%type /= BoundaryTypeNOBC ) then
          boundary%box_size_x = trajectory%pbc_box(1,1)
          boundary%box_size_y = trajectory%pbc_box(2,2)
          boundary%box_size_z = trajectory%pbc_box(3,3)
          call update_boundary(enefunc%table%table,  &
                               enefunc%pairlistdist, &
                               boundary)
        end if

        ! always needs to update pairlist
        !
        call update_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

        ! compute energy
        !
        if ( .not. mbar .and. .not. option%rest_component ) then
          call compute_energy(molecule, enefunc, pairlist, boundary, &
                              .true., .false.,                       &
                              dynvars%coord, dynvars%trans,          &
                              dynvars%coord_pbc, dynvars%energy,     &
                              dynvars%temporary, dynvars%force,      &
                              dynvars%force_omp,                     &
                              dynvars%virial, dynvars%virial_extern)

          if (main_rank .and. .not. has_rem) then
            call write_energy(0, ene_out, option%component, cur_step, &
                              enefunc, dynvars)
          end if
          tot = dynvars%energy%total
        end if

        if ( has_rem ) then

          if ( .not. option%mbar .and. has_rep_trj ) then
            do while ( cur_step > rstep )
              prev_parmsetid = parmsetid
              read(rem_inp,*,iostat=ierr) rstep, parmsetid
              if ( ierr /= 0 ) exit
            end do
          end if

          if (main_rank) then
            write(MsgOut,*) "step: ", cur_step, " parmsetid: ", prev_parmsetid
          end if

          if ( .not. mbar ) then

            call assign_condition(prev_parmsetid, remd, ensemble, &
                                  molecule, enefunc)

            ! compute energy
            !
            call compute_energy(molecule, enefunc, pairlist, boundary, &
                                .true., .false.,                       &
                                dynvars%coord, dynvars%trans,          &
                                dynvars%coord_pbc, dynvars%energy,     &
                                dynvars%temporary, dynvars%force,      &
                                dynvars%force_omp,                     &
                                dynvars%virial, dynvars%virial_extern)

            if (main_rank) then
              !call write_energy(prev_parmsetid, ene_out, option%component, &
              !                  cur_step, enefunc, dynvars)
              write(ene_out,'(2I16,F16.4)') &
                  cur_step, prev_parmsetid, dynvars%energy%total
            end if

            ! put back parameters
            !
            call assign_condition(option%tgt_parmid, remd, ensemble, &
                                  molecule, enefunc)

            !! check
            !call compute_energy(molecule, enefunc, pairlist, boundary, &
            !                    .true., .false.,                       &
            !                    dynvars%coord, dynvars%trans,          &
            !                    dynvars%coord_pbc, dynvars%energy,     &
            !                    dynvars%temporary, dynvars%force,      &
            !                    dynvars%force_omp,                     &
            !                    dynvars%virial, dynvars%virial_extern)

            !if (main_rank) then
            !  call write_energy(ene_out, option%component, cur_step, enefunc, &
            !                    dynvars)
            !end if
            
          else

            do i = 1, nreplicas

              call assign_condition(i, remd, ensemble, &
                                    molecule, enefunc)

              ! compute energy
              !
              call compute_energy(molecule, enefunc, pairlist, boundary, &
                                .true., .false.,                       &
                                dynvars%coord, dynvars%trans,          &
                                dynvars%coord_pbc, dynvars%energy,     &
                                dynvars%temporary, dynvars%force,      &
                                dynvars%force_omp,                     &
                                dynvars%virial, dynvars%virial_extern)

              !if (main_rank) then
              !  call write_energy(i, ene_out, option%component, &
              !                    cur_step, enefunc, dynvars)
              !end if

              tot_mbar(i) = dynvars%energy%total

            end do

            if (main_rank) then
              write(ene_out,'(I16,2x)',advance="NO") cur_step
              do i = 1, nreplicas
                write(ene_out,'(F16.4)',advance="NO") tot_mbar(i)
              end do
              write(ene_out,*)
            end if

          end if

        end if

      end do

    end do

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8  !
  !  Subroutine    write_energy
  !> @brief        write energy
  !! @authors      MK
  !! @param[in]    ene_out    : 
  !! @param[in]    component  :
  !! @param[in]    istep      :
  !! @param[inout] enefunc    : enefunc information
  !! @param[inout] dynvars    : dynvars information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_energy(repid, ene_out, component, istep, enefunc, dynvars)

    integer,         intent(in)    :: repid
    integer,         intent(in)    :: ene_out
    logical,         intent(in)    :: component
    integer,         intent(in)    :: istep
    type(s_enefunc), intent(inout) :: enefunc
    type(s_dynvars), intent(inout) :: dynvars

    real(wp) :: tot_rest
    integer :: i

    if ( component ) then
      tot_rest = sum(dynvars%energy%restraint(1:enefunc%num_restraintfuncs))
      write(ene_out,'(I3,I13,10F16.4)',advance="no") &
                             repid, &
                             istep, &
                             dynvars%energy%total, &
                             dynvars%energy%bond, &
                             dynvars%energy%angle, &
                             dynvars%energy%urey_bradley, &
                             dynvars%energy%dihedral, &
                             dynvars%energy%improper, &
                             dynvars%energy%cmap, &
                             dynvars%energy%electrostatic, &
                             dynvars%energy%van_der_waals, &
                             tot_rest
      do i = 1, enefunc%num_restraintfuncs
        write(ene_out,'(10F16.4)',advance="no") &
          dynvars%energy%restraint(i)
      end do
      write(ene_out,*)
    else
      write(ene_out,'(I16,F16.4)') istep, dynvars%energy%total
    end if

    return

  end subroutine write_energy

end module ea_analyze_mod
