!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   contact_analyze_mod
!> @brief   run contact analysis
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module contact_analyze_mod

  use contact_option_str_mod
  use sa_domain_mod
  use sa_boundary_mod
  use sa_tool_mod
  use sa_domain_str_mod
  use sa_boundary_str_mod
  use sa_ensemble_str_mod
  use sa_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use fileio_trj_mod
  use fileio_mod
  use molecules_str_mod
  use measure_mod
  use select_atoms_str_mod
  use sa_timers_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
  use string_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
  use omp_lib
#endif

  implicit none
  private

  ! subroutines
  !
  public  :: run_contact
  public  :: check_options

  ! structures for atom index
  !
  type, public :: s_atomidx
    integer, allocatable :: ig(:)  
    integer, allocatable :: is(:) 
  end type s_atomidx

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_contact 
  !> @brief        run analyzing trajectories
  !! @authors      IY
  !! @param[inout] molecule       : molecule information
  !! @param[in]    trj_list       : trajectory file list information
  !! @param[in]    output         : output information
  !! @param[inout] option         : option information
  !! @param[in]    contact_option : contact_option information
  !! @param[inout] trajectory     : trajectory information
  !! @param[inout] boundary       : boundary information
  !! @param[inout] domain         : domain information
  !! @param[in]    ensemble       : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_contact(molecule, trj_list, output, option, &
                         contact_option, trajectory, boundary, domain, ensemble)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option), target,  intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_contact_option),  intent(in)    :: contact_option

    ! local variables
    type(s_trj_file)         :: trj_in
    real(wp)                 :: cofm(3)
    real(wp)                 :: move(3)
    real(wp)                 :: csize_x, csize_y, csize_z, cvolume, acvolume
    real(wp)                 :: i_x, i_y, i_z, j_x, j_y, j_z
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: dist_ij, global_mindist, local_mindist
    integer                  :: ifile, istep, totstep, nfile
    integer                  :: iatom, jatom 
    integer                  :: txt_out 
    integer                  :: alloc_step
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: natom_selec
    integer                  :: ithread
    integer                  :: id, my_id
    integer                  :: icount, iatom_count
    integer                  :: ig, is, jg, js
    integer                  :: natom_group_i, natom_group_j, natom_pair
    integer                  :: ngroup, igroup, jgroup
    integer                  :: icelg, jcelg
    integer                  :: global_minrank, local_minthread
    integer                  :: minatom_i, minatom_j
    integer                  :: inside, inside_i, inside_j
    integer                  :: ntalcontact, ntalcontact_allg
    integer                  :: irank

    type(s_parray),          allocatable   :: selec_atom(:)
    type(s_atomidx),         allocatable   :: selatoms_inbnd(:)
    real(wp),                allocatable   :: mindist_proc(:)
    real(wp),                allocatable   :: mindist_all(:)
    real(wp),                allocatable   :: mindist_ij_thread(:)
    integer,                 allocatable   :: minatom_i_proc(:)
    integer,                 allocatable   :: minatom_j_proc(:)
    integer,                 allocatable   :: minatom_i_all(:)
    integer,                 allocatable   :: minatom_j_all(:)
    integer,                 allocatable   :: ncontact_thread(:)
    integer,                 allocatable   :: ncontact_proc(:)
    integer,                 allocatable   :: ncontact_all(:)
    integer,                 allocatable   :: minatom_i_thread(:)
    integer,                 allocatable   :: minatom_j_thread(:)


    ! check cntrl parameter
    !
    totstep = 0
    if (main_rank) then
      call check_options(molecule, boundary, domain, trj_list, &
                         option, contact_option, totstep)
    end if

    nfile  = size(trj_list%md_steps)
    ngroup = size(option%analysis_atoms)

    alloc_stat   = 0
    dealloc_stat = 0

    !allocate output arrays
    !

    allocate(selec_atom(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(selatoms_inbnd(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    allocate(mindist_ij_thread(nthread), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(minatom_i_thread(nthread), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(minatom_j_thread(nthread), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    mindist_ij_thread(1:nthread)  = 0.0_wp
    minatom_i_thread(1:nthread)   = 0
    minatom_j_thread(1:nthread)   = 0

    if (contact_option%mode == ModeMindist) then

      allocate(mindist_all(nproc_country), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      allocate(mindist_proc(nproc_country), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      mindist_all(1:nproc_country)  = 0
      mindist_proc(1:nproc_country) = 0

      allocate(minatom_i_all(nproc_country), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      allocate(minatom_j_all(nproc_country), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      allocate(minatom_i_proc(nproc_country), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      allocate(minatom_j_proc(nproc_country), stat =alloc_stat)

      minatom_i_all(1:nproc_country) = 0
      minatom_j_all(1:nproc_country) = 0
      minatom_i_proc(1:nproc_country) = 0
      minatom_j_proc(1:nproc_country) = 0

    end if

    if (contact_option%mode == ModeNumber) then

      allocate(ncontact_all(nproc_country), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      allocate(ncontact_proc(nproc_country), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc
      allocate(ncontact_thread(nthread), stat =alloc_stat)
      if (alloc_stat /= 0) call error_msg_alloc

      ncontact_all(1:nproc_country)  = 0
      ncontact_proc(1:nproc_country) = 0
      ncontact_thread(1:nthread)     = 0

    end if

    acvolume = 0

    ! open output file
    !
    if (main_rank) then
      if (output%txtfile .ne. '') then
        call open_file(txt_out, output%txtfile, IOFileOutputReplace)
      else
        txt_out = MsgOut
      end if
    end if

    natom_selec = 0
    do igroup = 1, ngroup
      selec_atom(igroup)%idx => option%analysis_atoms(igroup)%idx
      natom_selec = natom_selec + size(option%analysis_atoms(igroup)%idx)
    end do

    domain%num_atom_selec   = natom_selec
    domain%num_group        = ngroup

    if (main_rank) then
      if (totstep == 0) then

        if (contact_option%mode == ModeMindist) then
          write(txt_out, '(A)') "# This file was created with contact_analysis"
          write(txt_out, '(A)') "#"
          write(txt_out, '(A)') "# The closest atom pair (k,l) between group i and j"
          write(txt_out, '(A)') "# Column 1: number of analysis frame"
          write(txt_out, '(A)') "# Column 2: target group i"
          write(txt_out, '(A)') "# Column 3: target group j"
          write(txt_out, '(A)') "# Column 4: index of atom k in the group i"
          write(txt_out, '(A)') "# Column 5: index of atom l in the group j"
          write(txt_out, '(A)') "# Column 6: distance between atoms k and l (A)"
          write(txt_out, '(A)') "# ###: minimum distance is out of the range"
          write(txt_out, '(A)') ""
        end if

        if (contact_option%mode == ModeNumber) then
          write(txt_out, '(A)') "# This file was created with contact_analysis"
          write(txt_out, '(A)') "#"
          write(txt_out, '(A)') "# Number of Contact between group i and j"
          write(txt_out, '(A)') "# Column 1: number of analysis frame"
          write(txt_out, '(A)') "# Column 2: target group i"
          write(txt_out, '(A)') "# Column 3: target group j"
          write(txt_out, '(A)') "# Column 4: Number of contact between i and j"
          write(txt_out, '(A)') ""
        end if

      end if
    end if

    ! set timer
    !
    call timer(TimerAnalysis, TimerOn)


    ! analysis loop
    !
    do ifile = 1, nfile

      ! open trajectory file
      !
      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(ifile)

        ! read trajectory
        !   coordinates of one MD snapshot are saved in trajectory%coord)
        !
        call read_trj(trj_in, trajectory)

        if (mod(istep, trj_list%ana_periods(ifile)) == 0) then

          totstep = totstep + 1

          ! get volume of a cell
          bsize_x = boundary%box_size_x
          bsize_y = boundary%box_size_y
          bsize_z = boundary%box_size_z

          csize_x = bsize_x/boundary%num_cells_x
          csize_y = bsize_y/boundary%num_cells_y
          csize_z = bsize_z/boundary%num_cells_z

          cvolume = csize_x*csize_y*csize_z
          acvolume = acvolume + cvolume

          ! refine atomcoodinate
          do iatom = 1, molecule%num_atoms
            molecule%atom_coord(1,iatom) = trajectory%coord(1,iatom)
            molecule%atom_coord(2,iatom) = trajectory%coord(2,iatom)
            molecule%atom_coord(3,iatom) = trajectory%coord(3,iatom)
          end do

          if (contact_option%recenter /= 0) then
            call get_cofm(molecule, selec_atom, contact_option%recenter, cofm)
            move(1) = -cofm(1)
            move(2) = -cofm(2)
            move(3) = -cofm(3)
            call shift_molecule(move, molecule)
          end if

          ! wrap molecule
          !
          call wrap_molecule(option, trajectory, ensemble, &
                             boundary, molecule)

          ! setup domains
          !
          call setup_atom(option, ensemble, &
                          trajectory, selec_atom, molecule, boundary, domain)

          ! check cntrl parameter not to include image atoms
          ! for NPT trajectory
          !
          if (main_rank) then
            if (boundary%type == BoundaryTypePBC) then
              if (option%determine_box == DetermineBoxTrajectory) then
                if (ensemble%ensemble == EnsembleNPT  .or. &
                    ensemble%ensemble == EnsembleNPAT .or. &
                    ensemble%ensemble == EnsembleNPgT ) then
                  if (main_rank) then
                    call check_options(molecule, boundary, domain, trj_list, &
                         option, contact_option, totstep)
                  end if
                end if
              end if
            end if
          end if

          ! main analysis starts
          !
          ! determin the solute atom inside the boundary
          !
          do alloc_step = 1, 2
            iatom_count = 0
            do igroup = 1, ngroup
              icount = 0
              do iatom = 1, size(selec_atom(igroup)%idx)
                iatom_count = iatom_count +1
                ig = selec_atom(igroup)%idx(iatom)
                icelg = domain%selec_atom2cell(iatom_count)

                if(icelg > 0)then
                  inside = domain%cell_g2l(icelg)
                  inside = inside + domain%cell_g2b(icelg)
                else
                  inside = 0
                endif

                if (inside > 0) then
                  icount = icount +1
                  if (alloc_step == 2) then
                    selatoms_inbnd(igroup)%ig(icount)=ig
                    selatoms_inbnd(igroup)%is(icount)=iatom_count
                  end if
                end if

              end do ! iatom

              if (alloc_step == 1) then

                if (allocated(selatoms_inbnd(igroup)%ig))then
                  deallocate(selatoms_inbnd(igroup)%ig, stat=dealloc_stat)
                  deallocate(selatoms_inbnd(igroup)%is, stat=dealloc_stat)
                end if
                allocate(selatoms_inbnd(igroup)%ig(icount),  stat =alloc_stat)
                allocate(selatoms_inbnd(igroup)%is(icount), stat =alloc_stat)

                selatoms_inbnd(igroup)%ig(1:icount) = 0
                selatoms_inbnd(igroup)%is(1:icount) = 0
              end if

            end do ! igroup

          end do ! alloc_step


          ! determin the minimum atom pair or contact number between group ij
          !
          ntalcontact_allg = 0
          do igroup = 1 , ngroup

            natom_group_i   = size(selatoms_inbnd(igroup)%ig)

            do jgroup = 1 , ngroup

              if (jgroup > igroup) then
                natom_group_j   = size(selatoms_inbnd(jgroup)%ig)

                ! initialize
                !
                natom_pair = natom_group_i * natom_group_j
                mindist_ij_thread(1:nthread) = 8999998.0
                minatom_i_thread(1:nthread)  = 0
                minatom_j_thread(1:nthread)  = 0

                if (contact_option%mode == ModeNumber) then
                  ncontact_all(1:nproc_country)  = 0
                  ncontact_proc(1:nproc_country) = 0
                  ncontact_thread(1:nthread)     = 0
                end if

                if (igroup /= jgroup .and. natom_pair > 0) then

                  !$omp parallel default(shared)                           &
                  !$omp private(id, iatom, my_id, i_x, i_y, i_z, ig, is)   &
                  !$omp private(icelg, inside_i, jatom, jg, j_x, j_y, j_z) &
                  !$omp private(js, jcelg, inside_j, dist_ij)
                  !
#ifdef OMP
                  id  = omp_get_thread_num()
#else
                  id  = 0
#endif
                  my_id = id
                  do iatom = my_id+1, natom_group_i, nthread
                    ig = selatoms_inbnd(igroup)%ig(iatom)

                    i_x = molecule%atom_coord(1,ig)
                    i_y = molecule%atom_coord(2,ig)
                    i_z = molecule%atom_coord(3,ig)

                    is = selatoms_inbnd(igroup)%is(iatom)
                    icelg = domain%selec_atom2cell(is)
                    inside_i = domain%cell_g2l(icelg)

                    !iatom should be in domain (not boundary)
                    !
                    if (inside_i > 0) then 
                      do jatom = 1 , natom_group_j
                        jg = selatoms_inbnd(jgroup)%ig(jatom)

                        j_x = molecule%atom_coord(1,jg)
                        j_y = molecule%atom_coord(2,jg)
                        j_z = molecule%atom_coord(3,jg)

                        js = selatoms_inbnd(jgroup)%is(jatom)
                        jcelg = domain%selec_atom2cell(js)
                        inside_j = domain%cell_g2b(jcelg)

                        ! translate jatom to the closest position from iatom
                        !
                        j_x       = j_x - bsize_x*anint((j_x-i_x)/bsize_x)
                        j_y       = j_y - bsize_y*anint((j_y-i_y)/bsize_y)
                        j_z       = j_z - bsize_z*anint((j_z-i_z)/bsize_z)

                        dist_ij = (i_x - j_x)*(i_x - j_x)+ &
                                  (i_y - j_y)*(i_y - j_y)+ &
                                  (i_z - j_z)*(i_z - j_z)

                        if (dist_ij < mindist_ij_thread(id+1) .and. &
                            dist_ij > 0.0) then
                          mindist_ij_thread(id+1) = dist_ij
                          minatom_i_thread(id+1)  = ig
                          minatom_j_thread(id+1)  = jg
                        end if

                        if(contact_option%mode == ModeNumber) then
                          if (sqrt(dist_ij) < contact_option%ana_range) then
                            ncontact_thread(id+1) = ncontact_thread(id+1) + 1
                          end if
                        end if

                      end do !jatom
                    end if ! inside g2l>0
                  end do !iatom
                  !$omp end parallel

                end if !natom_pair > 0

                local_mindist   = minval(mindist_ij_thread,1)
                local_minthread = minloc(mindist_ij_thread,1)

                local_mindist = sqrt(local_mindist)

                ! combine and output the results
                ! 
                if (contact_option%mode == ModeMindist) then

                  mindist_proc(my_city_rank+1)   = local_mindist 
                  minatom_i_proc(my_city_rank+1) = &
                       minatom_i_thread(local_minthread)
                  minatom_j_proc(my_city_rank+1) = &
                       minatom_j_thread(local_minthread)

#ifdef HAVE_MPI_GENESIS
                  call mpi_allgather( &
                       mindist_proc(my_city_rank+1), 1, mpi_wp_real, &
                       mindist_all(1), 1, mpi_wp_real, mpi_comm_country, ierror)

                  call mpi_allgather( &
                       minatom_i_proc(my_city_rank+1), 1, mpi_integer, &
                       minatom_i_all(1), 1, mpi_integer, mpi_comm_country,ierror)

                  call mpi_allgather( &
                       minatom_j_proc(my_city_rank+1), 1, mpi_integer, &
                       minatom_j_all(1), 1, mpi_integer, mpi_comm_country,ierror)

                  ! find the process which gives the global minimum
                  global_mindist      = minval(mindist_all,1)
                  global_minrank      = minloc(mindist_all,1)

                  minatom_i           = minatom_i_all(global_minrank)
                  minatom_j           = minatom_j_all(global_minrank)
#else
                  global_mindist      = mindist_proc(my_city_rank+1)
                  minatom_i           = minatom_i_proc(my_city_rank+1)
                  minatom_j           = minatom_j_proc(my_city_rank+1)
#endif

                  if (main_rank) then
                    if (global_mindist <= contact_option%ana_range) then
                    write(txt_out,'(I8,A,I8,A,I8,A,I8,A,I8,A,ES25.16E3)') &
                            totstep," ",igroup, " ", &
                            jgroup, " ", minatom_i, " ", minatom_j," ", &  
                            global_mindist 
                    else
                      write(txt_out,'(I8,A,I8,A,I8,A,A)') &
                            totstep," ",igroup, " ", jgroup, " ", &
                            "     ###     ###     ###"
                    end if
                  end if

                end if ! mode = mindist

                if (contact_option%mode == ModeNumber) then

                  ! accumulate ncontact over threads
                  !
                  do ithread = 1, nthread
                    ncontact_proc(my_city_rank+1) =      &
                         ncontact_proc(my_city_rank+1) + & 
                         ncontact_thread(ithread)
                  end do

#ifdef HAVE_MPI_GENESIS
                  call mpi_allgather( &
                       ncontact_proc(my_city_rank+1), 1, mpi_integer, &
                       ncontact_all(1), 1, mpi_integer,mpi_comm_country,ierror)
#else
                  ncontact_all(1) = ncontact_proc(my_city_rank+1)
#endif

                  ! calculate total number of contact
                  !
                  ntalcontact = 0
                  if (main_rank) then
                    do irank = 1, nproc_country
                      ntalcontact = ntalcontact + ncontact_all(irank)
                    end do
                    ntalcontact_allg = ntalcontact_allg + ntalcontact
                    write(txt_out,'(I8,A,I8,A,I8,A,I8)') &
                          totstep," ",igroup, " ", jgroup, " ",ntalcontact
                  end if ! main_rank
                end if ! mode = number

              end if ! jgroup > igroup
            end do ! jgroup

          end do ! igroup

          if (contact_option%mode == ModeNumber) then
            if (main_rank) then
              write(txt_out,'(I8,A,I8)')totstep,"      total= ",ntalcontact_allg 
              write(txt_out,'(A)') " "
            end if

          end if ! mode= number

          ! if (main_rank) then
          !   write(MsgOut,*) "Analysis for ", totstep, "-th step finished"
          ! end if

        end if ! istep = ana_period
      end do ! istep

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do

    if (main_rank) then
      if (output%txtfile .ne. '') then
        call close_file(txt_out)
      end if
    end if

    deallocate(selatoms_inbnd, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(selec_atom, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    deallocate(mindist_ij_thread, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(minatom_i_thread, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(minatom_j_thread, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    if (contact_option%mode == ModeMindist) then
      deallocate(mindist_proc, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(mindist_all, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(minatom_i_proc, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(minatom_j_proc, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(minatom_i_all, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(minatom_j_all, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    end if

    if (contact_option%mode == ModeNumber) then
      deallocate(ncontact_proc, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(ncontact_all, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
      deallocate(ncontact_thread, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    end if

    ! close output file
    !

    call timer(TimerAnalysis, TimerOff)

    return

  end subroutine run_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_options 
  !> @brief        check options in control file 
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_options(molecule, boundary, domain, trj_list, option,&
                           contact_option, totstep)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_option),          intent(in)    :: option
    type(s_contact_option),  intent(in)    :: contact_option
    integer,                 intent(in)    :: totstep

    ! local variables
    real(wp)                 :: max_csize
    real(wp)                 :: csize_xyz(3)
    real(wp)                 :: min_bsize
    real(wp)                 :: bsize_xyz(3)
    logical                  :: find_error


    csize_xyz(1)   =  boundary%cell_size_x
    csize_xyz(2)   =  boundary%cell_size_y
    csize_xyz(3)   =  boundary%cell_size_z
    max_csize      =  maxval(csize_xyz)
    bsize_xyz(1)   =  boundary%box_size_x
    bsize_xyz(2)   =  boundary%box_size_y
    bsize_xyz(3)   =  boundary%box_size_z
    min_bsize      =  minval(bsize_xyz)

    find_error  = .false.

    ! check options before the analysis loop
    !
    if (totstep == 0) then

      if (option%buffer < contact_option%ana_range) then
        write(MsgOut,'(a)') &
          'Check_Options> [SPANA_OPTION]:buffer should be larger than [CONTACT_OPTION]:range'
        find_error  = .true.
      end if

      if (contact_option%ana_range == 0.0_wp) then
        write(MsgOut,'(a)') &
          'Check_Options> [CONTACT_OPTION]:range is not specified'
        find_error  = .true.
      end if

      if (boundary%type == BoundaryTypeNOBC) then
        if (option%determine_box == DetermineBoxTrajectory) then
          write(MsgOut,'(a)') &
            'Check_Options> [SPANA_OPTION]:TRAJECTORY is available only in PBC'
          find_error  = .true.
        end if
      else if (boundary%type == BoundaryTypePBC) then
        if (option%determine_box == DetermineboxMax) then
          write(MsgOut,'(a)') &
            'Check_Options> [SPANA_OPTION]:MAX is available only in NOBC'
          find_error  = .true.
        end if
      end if

      if (option%wrap) then
        if (trj_list%trj_type /= 2) then
          write(MsgOut,'(a)') &
            'Check_Options> cannot wrap atoms when [TRAJECTORY]:trj_type is not COOR+BOX'
          find_error  = .true.
        end if
      end if

      if (boundary%num_domain(1) /= 0 .and. &
          boundary%num_domain(2) /= 0 .and. &
          boundary%num_domain(3) /= 0) then

        if (product(boundary%num_domain) /= nproc_country) then
          write(MsgOut,'(a)') &
            'Check_Options> # of process is not domain_x * domain_y * domain_z '
          find_error  = .true.
        end if
      end if
    end if ! totstep ==0

    ! check buffer length no to include self image
    ! this is better to repeat for every analysis frame in the case of NPT
    !
    if (option%buffer+max_csize >= min_bsize/2.0_wp) then
      write(MsgOut,'(a)') &
        'Check_Options> [SPANA_OPTION]:buffer is too large. There is a risk to include the self image atom'
      find_error  = .true.
    end if

    if (find_error) then
      write(MsgOut,'(a, I8)')'Check_Options> error detected at step ',totstep 
      call error_msg(&
        'Check_Options> analysis terminated because of option error')
    end if

    return

  end subroutine check_options

end module contact_analyze_mod
