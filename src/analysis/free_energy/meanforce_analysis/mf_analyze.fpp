!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   mf_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module mf_analyze_mod

  use mf_option_str_mod
  use fileio_trj_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use fitting_mod
  use fitting_str_mod
  use molecules_str_mod
  use fileio_trj_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: check_dcdfile
  private :: get_replicate_name1

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    input    : input information
  !! @param[in]    output   : output information
  !! @param[in]    fitting  : fitting information
  !! @param[in]    option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, input, output, fitting, option)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(in)    :: option

    ! local variables
    type(s_trj_file)         :: file
    type(s_trajectory)       :: trajectory
    real(wp)                 :: tmp
    integer                  :: nimage, nstep, nsel, natom
    integer                  :: iimage, istep, isel
    integer                  :: pfile, ofile, i
    character(MaxFilename)   :: filename

    real(wp), allocatable    :: path(:,:), f(:,:), path_prv(:,:), f_prv(:,:)
    real(wp), allocatable    :: F_k(:,:), F_k_prv(:,:), Fdat(:)


    ! check check only
    !

    if (option%check_only) &
      return


    ! check nstep and allocate trajectory
    !

    call check_dcdfile(input%dcdfile, nstep, natom)


    ! allocate memory
    !

    nimage = option%nimage
    nsel   = size(option%cv_atom%idx)

    allocate(path    (3,nsel), &
             path_prv(3,nsel), &
             f       (3,nsel), &
             f_prv   (3,nsel), &
             F_k     (3,nsel), &
             F_k_prv (3,nsel), &
             Fdat    (nimage))

    call alloc_trajectory(trajectory, natom)

    do isel = 1, nsel
      F_k(1:3,isel) = 0.0_wp
    end do


    ! analyze
    !

    do iimage = 1, nimage

      do isel = 1, nsel
        f(:,isel) = 0.0_wp
      end do

      ! read path file

      filename = get_replicate_name1(input%pathfile, iimage)

      call open_file(pfile, filename, IOFileInput)
      read(pfile,*) tmp, ((path(i,isel),i=1,3),isel=1,nsel)
      call close_file(pfile)

      ! read dcd file

      filename = get_replicate_name1(input%dcdfile, iimage)
      write(MsgOut,'(a,a)') '  read and analyze trajectory: ',trim(filename)

      call open_trj(file,           &
                    filename,       &
                    TrjFormatDCD,   &
                    TrjTypeCoorBox, &
                    IOFileInput)

      do istep = 1, nstep

        call read_trj(file, trajectory)

        call run_fitting(fitting,             &
                         molecule%atom_coord, &
                         trajectory%coord,    &
                         trajectory%coord)

        do isel = 1, nsel
          f(1:3,isel) = f(1:3,isel) + &
                (trajectory%coord(1:3,option%cv_atom%idx(isel))-path(1:3,isel))
        end do

      end do

      call close_trj(file)

      ! calc mean-force
      do isel = 1, nsel
        f(1:3,isel) = &
             2.0_wp * option%force_constant * (f(1:3,isel) / real(nstep,wp))
      end do

      ! free energy profile
      if (iimage > 1) then
        do isel = 1, nsel
          F_k(1:3,isel) = F_k_prv(1:3,isel) + &
               (path(1:3,isel) - path_prv(1:3,isel))*(-0.5_wp) * &
                  (f(1:3,isel) + f_prv(1:3,isel))
        end do
      end if

      do isel = 1, nsel
        F_k_prv(1:3,isel)  = F_k(1:3,isel)
        path_prv(1:3,isel) = path(1:3,isel)
        f_prv(1:3,isel)    = f(1:3,isel)
      end do

      Fdat(iimage) = 0.0_wp
      do isel = 1, nsel
        do i = 1,3
          Fdat(iimage) = Fdat(iimage) + F_k(i,isel)
        end do
      end do

    end do


    ! output result
    !
    call open_file(ofile, output%fenefile, IOFileOutputReplace)

    do iimage = 1, nimage
      write(ofile,'(es25.16e3)') Fdat(iimage)
    end do

    call close_file(ofile)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [fenefile] ' // trim(output%fenefile)
    write(MsgOut,'(A)') '    Column 1: Free energy profile (kcal/mol)'
    write(MsgOut,'(A)') ''


    ! deallocate memory
    !
    call dealloc_trajectory(trajectory)

    deallocate(path, path_prv, f, f_prv, F_k, Fdat)

    return

  end subroutine analyze

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

end module mf_analyze_mod
