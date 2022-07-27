!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fa_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module fa_analyze_mod

  use fa_option_str_mod
  use fitting_mod
  use fileio_trj_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: assign_mass
  private :: sum_vcv_matrix
  private :: out_vcv_matrix

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, trajectory, fitting, option, &
                     input, output)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_pdb)              :: pdb_in
    type(s_trj_file)         :: trj_in
    real(wp)                 :: nstru_r
    integer                  :: natom, ntraj, nsela, nself, nvcv, nstru
    integer                  :: iatom, itraj, isel, istep
    integer                  :: i, j, alloc_stat

    real(wp), allocatable  :: sqrt_mass (:)
    real(wp), allocatable  :: ave_coord (:,:)
    real(wp), allocatable  :: aft_coord (:,:)
    real(wp), allocatable  :: trj_coord (:,:)
    real(wp), allocatable  :: vcv_vector(:)
    real(wp), allocatable  :: vcv_matrix(:,:)


    if (option%check_only) &
      return

    natom = molecule%num_atoms
    ntraj = size(trj_list%md_steps)
    nsela = size(option%analysis_atom%idx)
    nself = size(fitting%fitting_atom%idx)
    if (option%vcv_matrix == VcvMatrixGlobal) then
      nvcv = nsela * 3
    else
      nvcv = nself * 3
    end if

    ! allocate memory
    !
    allocate(sqrt_mass (natom),   &
             ave_coord (3,natom), &
             aft_coord (3,natom), &
             trj_coord (3,natom), &
             vcv_vector(nvcv),    &
             vcv_matrix(nvcv, nvcv), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! prepare data
    !

    ! setup mass
    if (fitting%mass_weight) then
      call assign_mass(molecule)
    end if

    if (fitting%mass_weight) then
      do iatom = 1, natom
        sqrt_mass(iatom) = sqrt(molecule%mass(iatom))
      end do
    else
      do iatom = 1, natom
        sqrt_mass(iatom) = 1.0_wp
      end do
    end if

    ! setup average coordintes
    !
    call input_pdb(input%pdb_avefile, pdb_in)

    if (nsela /= size(pdb_in%atom_name)) &
      call error_msg('Analyze> different count between pdb_avefile atom '// &
                     'and analysis_atom selection.')

    do isel = 1, nsela
      iatom = option%analysis_atom%idx(isel)
      ave_coord(1:3,iatom) = pdb_in%atom_coord(1:3,isel) * sqrt_mass(iatom)
    end do

    call dealloc_pdb_all(pdb_in)

    ! setup fitted average coordinates
    !
    call input_pdb(input%pdb_aftfile, pdb_in)

    if (nself /= size(pdb_in%atom_name)) &
      call error_msg('Analyze> different count between pdb_aftfile atom '// &
                     'and fitting_atom selection.')

    do isel = 1, nself
      iatom = fitting%fitting_atom%idx(isel)
      aft_coord(1:3,iatom) = pdb_in%atom_coord(1:3,isel) * sqrt_mass(iatom)
    end do

    call dealloc_pdb_all(pdb_in)

    ! setup vcv matrix
    !
    do i = 1, nvcv
      vcv_vector(i) = 0.0_wp
      do j = 1, nvcv
        vcv_matrix(i, j) = 0.0_wp
      end do
    end do

    ! analysis loop
    !
    nstru = 0
    do itraj = 1, ntraj

      call open_trj(trj_in, &
                    trj_list%filenames(itraj), &
                    trj_list%trj_format,       &
                    trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(itraj)

        call read_trj(trj_in, trajectory)

        if (mod(istep, trj_list%ana_periods(itraj)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru

          ! geometrical fitting
          !
          call run_fitting(fitting,             &
                           molecule%atom_coord, &
                           trajectory%coord,    &
                           trajectory%coord)

          ! mass-weighted fitting
          !
          do iatom = 1, natom
            trj_coord(1:3,iatom) = trajectory%coord(1:3,iatom)*sqrt_mass(iatom)
          end do

          call run_fitting(fitting,   &
                           aft_coord, &
                           trj_coord, &
                           trj_coord)

          if (option%vcv_matrix == VcvMatrixGlobal) then

            call sum_vcv_matrix(option%analysis_atom%idx, &
                                trj_coord, &
                                ave_coord, &
                                vcv_vector, vcv_matrix)
          else

            call sum_vcv_matrix(fitting%fitting_atom%idx, &
                                trj_coord, &
                                aft_coord, &
                                vcv_vector, vcv_matrix)
          end if

        end if

      end do

      call close_trj(trj_in)

    end do


    ! make the varience-covarience matrix
    !
    nstru_r = real(nstru, wp)
    
    do i = 1, nvcv
      do j = 1, nvcv
        vcv_matrix(i,j) = vcv_matrix(i,j) / nstru_r
      end do
    end do


    ! output RMS deviation and Varience-covarience matrix
    !
    write(MsgOut,*) 'INFO> output RMSF and VCV Matrix'
    write(MsgOut,*) ' '

    if (option%vcv_matrix == VcvMatrixGlobal) then

      call out_vcv_matrix(output,   &
                          molecule, &
                          option%analysis_atom%idx, &
                          sqrt_mass, &
                          vcv_matrix)

    else

      call out_vcv_matrix(output,   &
                          molecule, &
                          fitting%fitting_atom%idx, &
                          sqrt_mass, &
                          vcv_matrix)

    end if


    ! Output summary
    !
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [rmsfile] ' // trim(output%rmsfile)
    write(MsgOut,'(A)') '    Column 1: Index of the analysis atom'
    write(MsgOut,'(A)') '    Column 2: Residue number of the analysis atom'
    write(MsgOut,'(A)') '    Column 3: Atom name of the analysis atom'
    write(MsgOut,'(A)') '    Column 4: Root-mean-square fluctuation (RMSF)(angstrom)'
    write(MsgOut,'(A)') '    Column 5: B-factor ([(8*PI^2)/3]*(RMSF)^2)'
    write(MsgOut,'(A)') ''


    ! deallocate memory
    !
    deallocate(vcv_matrix, vcv_vector, &
               trj_coord, aft_coord, ave_coord, sqrt_mass)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_mass
  !> @brief        assign mass
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_mass(molecule)

    ! parameters
    real(wp),                parameter     :: MassH   =  1.008000_wp
    real(wp),                parameter     :: MassC   = 12.011000_wp
    real(wp),                parameter     :: MassN   = 14.007000_wp
    real(wp),                parameter     :: MassO   = 15.999000_wp
    real(wp),                parameter     :: MassS   = 32.060000_wp
    real(wp),                parameter     :: MassP   = 30.974000_wp
    real(wp),                parameter     :: MassMG  = 24.305000_wp
    real(wp),                parameter     :: MassZN  = 65.370000_wp

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i


    write(MsgOut,'(a)'),      'WARNING: atom mass is not assigned.'
    write(MsgOut,'(a)'),      '   uses default mass.'
    write(MsgOut,'(a,f9.6)')  '      1) H   : ', MassH
    write(MsgOut,'(a,f9.6)')  '      2) C   : ', MassC
    write(MsgOut,'(a,f9.6)')  '      3) N   : ', MassN
    write(MsgOut,'(a,f9.6)')  '      4) O   : ', MassO
    write(MsgOut,'(a,f9.6)')  '      5) S   : ', MassS
    write(MsgOut,'(a,f9.6)')  '      6) P   : ', MassP
    write(MsgOut,'(a,f9.6)')  '      7) MG  : ', MassMG
    write(MsgOut,'(a,f9.6)')  '      8) ZN  : ', MassZN


    do i = 1, molecule%num_atoms

      if (molecule%atom_name(i)(1:1) == 'H') then
        molecule%mass(i) = MassH
      else if (molecule%atom_name(i)(1:1) == 'C') then
        molecule%mass(i) = MassC
      else if (molecule%atom_name(i)(1:1) == 'N') then
        molecule%mass(i) = MassN
      else if (molecule%atom_name(i)(1:1) == 'O') then
        molecule%mass(i) = MassO
      else if (molecule%atom_name(i)(1:1) == 'S') then
        molecule%mass(i) = MassS
      else if (molecule%atom_name(i)(1:1) == 'P') then
        molecule%mass(i) = MassP
      else if (molecule%atom_name(i)(1:2) == 'MG') then
        molecule%mass(i) = MassMG
      else if (molecule%atom_name(i)(1:2) == 'ZN') then
        molecule%mass(i) = MassZN
      else
        write(MsgOut,'(a,a)') 'Assign_Mass> Unknown atom :', &
             molecule%atom_name(i)
      end if

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine assign_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    sum_vcv_matrix
  !> @brief        sum vcv matrix
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sum_vcv_matrix(idx, trj_coord, ave_coord, vcv_vector, vcv_matrix)

    ! formal arguments
    integer,                 intent(in)    :: idx(:)
    real(wp),                intent(in)    :: trj_coord(:,:)
    real(wp),                intent(in)    :: ave_coord(:,:)
    real(wp),                intent(inout) :: vcv_vector(:)
    real(wp),                intent(inout) :: vcv_matrix(:,:)

    ! local variables
    integer                  :: i, j, iatom, nsel


    nsel = size(idx)

    do i = 1, nsel
      iatom = idx(i)
      vcv_vector(3*(i-1)+1) = trj_coord(1,iatom) - ave_coord(1,iatom)
      vcv_vector(3*(i-1)+2) = trj_coord(2,iatom) - ave_coord(2,iatom)
      vcv_vector(3*(i-1)+3) = trj_coord(3,iatom) - ave_coord(3,iatom)
    end do

    do i = 1, nsel * 3
      do j = 1, nsel * 3
        vcv_matrix(i,j) = vcv_matrix(i,j) + vcv_vector(i)*vcv_vector(j)
      end do
    end do

    return

  end subroutine sum_vcv_matrix

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_vcv_matrix
  !> @brief        output vcv matrix
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_vcv_matrix(output, molecule, idx, sqrt_mass, vcv_matrix)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_molecule),        intent(in)    :: molecule
    integer,                 intent(in)    :: idx(:)
    real(wp),                intent(in)    :: sqrt_mass(:)
    real(wp),                intent(inout) :: vcv_matrix(:,:)

    ! local variables
    real(wp)                 :: ttl_mass, chk_msf, bfct
    integer                  :: pca_out, rms_out, vcv_out, crs_out
    integer                  :: iatom, jatom, ires, jres, isel, jsel
    integer                  :: nsel, nvcv, nres, maxres, minres
    integer                  :: i, j, k, alloc_stat

    integer,     allocatable :: resno(:)
    real(wp),    allocatable :: rmsf(:), vcvmat2(:,:), crsmat2(:,:)


    call open_file(pca_out, output%pcafile, IOFileOutputNew)
    call open_file(rms_out, output%rmsfile, IOFileOutputNew)
    call open_file(vcv_out, output%vcvfile, IOFileOutputNew)
    call open_file(crs_out, output%crsfile, IOFileOutputNew)

    nsel = size(idx)

    allocate(rmsf(nsel), resno(nsel), stat=alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_alloc


    ! output Varience-covarience matrix

    write(MsgOut,*) 'Out_Vcvmat> output Varience-Covarience Matrix'
    write(MsgOut,*) ' '

    ttl_mass = 0.0_wp
    do isel = 1, nsel
      iatom = idx(isel)
      ttl_mass = ttl_mass + sqrt_mass(iatom)**2
    end do

    nvcv = nsel * 3
    chk_msf = 0.0_wp

    write(pca_out,*) nsel
    do i = 1, nvcv
      do j = 1, nvcv
        write(pca_out,*) vcv_matrix(i,j)
      end do
      chk_msf = chk_msf + vcv_matrix(i,i)
    end do

    write(MsgOut,*) 'Out_Vcvmat> rmsf = ', chk_msf, sqrt(chk_msf/ttl_mass)
    write(MsgOut,*) ' '
  
    ! output Rms.fluctuation

    write(MsgOut,*) 'Out_Vcvmat> output RMS fluctuation'
    write(MsgOut,*) ' '

    do i = 1, nsel
      rmsf(i) = 0.0_wp
    end do

    do i = 1, nvcv
      isel = (i - 1) / 3 + 1
      rmsf(isel) = rmsf(isel) + vcv_matrix(i,i)

      if (mod(i,3) == 0) then
        iatom = idx(isel)
        rmsf(isel) = rmsf(isel) / sqrt_mass(iatom)**2

        bfct = rmsf(isel) * (8.0_wp / 3.0_wp) * PI**2
        write(rms_out, '(2(i5,1x),a4,1x,2(f8.3,1x))') &
              isel,                       &
              molecule%residue_no(iatom), &
              molecule%atom_name(iatom),  &
              sqrt(rmsf(isel)), bfct
      end if
    end do

    ! output Varience-covarience matrix and Cross-correlation matrix for plot
    
    write(MsgOut,*) 'Out_Vcvmat> output Varience-Covarience matrix'
    write(MsgOut,*) ' '

    nres   = 0   
    maxres = 0
    minres = 999999

    do isel = 1, nsel
      iatom = idx(isel)
      if (molecule%atom_name(iatom) == 'CA') then
        nres = nres + 1
        resno(nres) = molecule%residue_no(iatom)
        maxres = max(resno(nres), maxres)
        minres = min(resno(nres), minres)
      end if
    end do

    if (minres > 1) minres = 1

    allocate(vcvmat2(minres:maxres,minres:maxres), &
             crsmat2(minres:maxres,minres:maxres), &
             stat=alloc_stat)

    if (alloc_stat /= 0) &
      call error_msg_alloc

    vcvmat2(:,:) = 0.0_wp
    crsmat2(:,:) = 0.0_wp

    do isel = 1, nsel
      do jsel = 1, nsel

        iatom = idx(isel)  
        jatom = idx(jsel)  

        if((molecule%atom_name(iatom) == 'CA').and. &
           (molecule%atom_name(jatom) == 'CA')) then 

          ires = molecule%residue_no(iatom)
          jres = molecule%residue_no(jatom)

          do k = 1, 3

            i = 3 * (isel - 1) + k
            j = 3 * (jsel - 1) + k

            vcvmat2(ires,jres) = vcvmat2(ires,jres) + &
                 vcv_matrix(i,j) / (sqrt_mass(iatom)*sqrt_mass(jatom))
          end do

          crsmat2(ires,jres) = &
               vcvmat2(ires,jres) / sqrt(rmsf(isel) * rmsf(jsel))

        end if
      end do
    end do

    do ires = 1, nres
      write(vcv_out,'(1000(f8.3,1x))') &
           (vcvmat2(resno(ires),resno(jres)), jres=1, nres)
      write(crs_out,'(1000(f8.3,1x))') &
           (crsmat2(resno(ires),resno(jres)), jres=1, nres)
    end do


    deallocate(rmsf, vcvmat2, crsmat2, resno)

    call close_file(pca_out)
    call close_file(rms_out)
    call close_file(vcv_out)
    call close_file(crs_out)

    return

  end subroutine out_vcv_matrix

end module fa_analyze_mod
