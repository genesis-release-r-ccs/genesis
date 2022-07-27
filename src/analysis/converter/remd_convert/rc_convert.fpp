!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rc_convert_mod
!> @brief   convert remd trajectory files
!! @authors Takaharu Mori (TM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rc_convert_mod

  use constants_mod
  use rc_option_str_mod
  use pbc_correct_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_trj_mod
  use fileio_mod
  use string_mod
  use measure_mod
  use messages_mod
 
  implicit none
  private

  ! subroutines
  public  :: convert
  private :: setup_convert
  private :: centering
  private :: get_replicate_name

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert remd trajectory files
  !! @authors      TM
  !! @param[in]    input      : input information
  !! @param[in]    output     : output information
  !! @param[in]    option     : option information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert(input, output, option, molecule, trajectory, fitting)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting

    ! local variables
    integer                       :: i, j, k, m, n
    integer                       :: num_conv_id, nlines
    integer                       :: num_open_files
    integer                       :: istep, parmid, unit_no, imdstep
    logical                       :: conv_crd, conv_log, conv_ene, do_skip
    logical                       :: done
    integer                       :: ierr
    character(MaxFilename)        :: filename
    character(2000)               :: line
    type(s_selatoms)              :: trjout_atom_trj
    integer,          allocatable :: conv_id(:), open_id(:), rem_unit_in(:)
    integer,          allocatable :: log_unit_in(:), log_unit_out(:)
    integer,          allocatable :: ene_unit_in(:), ene_unit_out(:)
    integer,          allocatable :: param_index(:)
    logical,          allocatable :: do_open_file(:)
    type(s_trj_file), allocatable :: trj_in(:), trj_out(:)


    ! check check only
    !
    if (option%check_only) &
      return


    ! check remfile, dcdfile, logfile, enefile existing
    !
    call setup_convert(input, output, option, conv_crd, conv_log, conv_ene)


    ! make convert id lists
    !
    if (.not. allocated(option%convert_ids)) then
      allocate(conv_id(option%num_replicas))
      do i = 1, option%num_replicas
        conv_id(i) = i
      end do
    else
      allocate(conv_id(size(option%convert_ids)))
      conv_id(:) = option%convert_ids(:)
      do i = 1, size(conv_id)
        if (conv_id(i) > option%num_replicas) &
          call error_msg('Convert> Convert IDs is out-of-range.')
      end do
    end if
    num_conv_id = size(conv_id)


    ! other setup
    !
    allocate(open_id     (option%num_replicas), &
             do_open_file(option%num_replicas), &
             trj_in      (option%num_replicas), &
             trj_out     (option%num_replicas), &
             rem_unit_in (option%num_replicas), &
             log_unit_in (option%num_replicas), &
             log_unit_out(option%num_replicas), &
             ene_unit_in (option%num_replicas), &
             ene_unit_out(option%num_replicas), &
             param_index (option%num_replicas))


    ! make lists to be opened
    !
    if (option%convert_type == ConvertTypeParameter) then

      do_open_file(:) = .false.
      do i = 1, option%num_replicas

        nlines = option%nsteps/option%exchange_period
        filename = get_replicate_name(input%remfile,i)
        call open_file(unit_no, filename, IOFileInput)

        do j = 1, nlines
          read(unit_no,*) istep, parmid
          do k = 1, num_conv_id
            m = conv_id(k)
            if (parmid == m) then
              do_open_file(i) = .true.
            end if
          end do
          if (do_open_file(i)) exit
        end do

        call close_file(unit_no)
      end do

    else if (option%convert_type == ConvertTypeReplica) then

      do_open_file(:) = .false.
      do i = 1, option%num_replicas

        do k = 1, num_conv_id
          m = conv_id(k)
          if (m == i) then
            do_open_file(i) = .true.
          end if
        end do

      end do

    end if


    ! print open file index lists
    !
    write(MsgOut,'(A)')   'Convert> Replicas that contain "convert_ids"'
    write(MsgOut,'(A,$)') '  Replicas: '
    do i = 1, option%num_replicas
      if (do_open_file(i)) then
        write(MsgOut,'(I5,$)') i
      end if
    end do
    write(MsgOut,'(A)') ''


    ! open input files
    !
    trj_in(:)%unit_no = 0
    rem_unit_in(:) = 0
    log_unit_in(:) = 0
    ene_unit_in(:) = 0

    do i = 1, option%num_replicas
       if (do_open_file(i)) then

         ! assign file unit no for remfiles
         if (option%convert_type == ConvertTypeParameter) then
           filename = get_replicate_name(input%remfile,i)
           call open_file(unit_no, filename, IOFileInput)
           rem_unit_in(i) = unit_no
         end if

         ! assign file unit no for dcdfiles
         if (conv_crd) then
           filename = get_replicate_name(input%dcdfile,i)
           call open_trj(trj_in(i), filename, TrjFormatDCD, &
                         TrjTypeCoorBox, IOFileInput)
         end if

         ! assign file unit no for logfiles
         if (conv_log) then
           filename = get_replicate_name(input%logfile,i)
           call open_file(unit_no, filename, IOFileInput)
           log_unit_in(i) = unit_no
         end if

         ! assign file unit no for enefiles of gREST
         if (conv_ene) then
           filename = get_replicate_name(input%enefile,i)
           call open_file(unit_no, filename, IOFileInput)
           ene_unit_in(i) = unit_no
         end if

       end if
    end do


    ! open output files
    !
    trj_out(:)%unit_no = 0
    log_unit_out(:)    = 0
    ene_unit_out(:)    = 0

    do k = 1, num_conv_id
      m = conv_id(k)

      if (conv_crd) then
        filename = get_replicate_name(output%trjfile,m)
        call open_trj(trj_out(k), filename, option%trjout_format, &
                      option%trjout_type, IOFileOutputNew)
      end if

      if (conv_log) then
        filename = get_replicate_name(output%logfile,m)
        call open_file(unit_no, filename, IOFileOutputNew)
        log_unit_out(k) = unit_no
      end if

      if (conv_ene) then
        filename = get_replicate_name(output%enefile,m)
        call open_file(unit_no, filename, IOFileOutputNew)
        ene_unit_out(k) = unit_no
      end if

    end do


    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Convert> Assign Input/Output file unit no'
    write(MsgOut,*)     ' Input  remfiles :', rem_unit_in(:)
    if (conv_crd) &
      write(MsgOut,*)   ' Input  dcdfiles :', trj_in(:)%unit_no
    if (conv_log) &
      write(MsgOut,*)   ' Input  logfiles :', log_unit_in(:)
    if (conv_ene) &
      write(MsgOut,*)   ' Input  logfiles :', ene_unit_in(:)
    if (conv_crd) &
      write(MsgOut,*)   ' Output trjfiles :', trj_out(:)%unit_no
    if (conv_log) &
      write(MsgOut,*)   ' Output logfiles :', log_unit_out(:)
    if (conv_ene) &
      write(MsgOut,*)   ' Output logfiles :', ene_unit_out(:)
    write(MsgOut,'(A)') ''


    ! read 0 step parameter index
    !
    if (option%convert_type == ConvertTypeParameter) then

      param_index(:) = 0
      do i = 1, option%num_replicas
        if (do_open_file(i)) then
          read(rem_unit_in(i),*) istep, param_index(i)
        end if
      end do

    else if (option%convert_type == ConvertTypeReplica) then

      param_index(:) = 0
      do i = 1, option%num_replicas
        if (do_open_file(i)) then
          param_index(i) = i
        end if
      end do

    end if


    ! read(i)/write(j) logfile header
    !
    if (conv_log) then

      do i = 1, option%num_replicas
        if (do_open_file(i)) then

          ! check length of the header
          read (log_unit_in(i),'(A)') line
          if (line(1:17) == 'Generate_Velocity') then
            n = 7
          else
            n = 1
          end if

          do_skip = .true.
          do j = 1, num_conv_id
            m = conv_id(j)
            if (m == param_index(i)) then
              write(log_unit_out(j),'(A)') trim(line)
              do k = 1, n
                read (log_unit_in (i),'(A)') line
                write(log_unit_out(j),'(A)') trim(line)
              end do
              do_skip = .false.
            end if
          end do

          if (do_skip) then
            do k = 1, n
              read (log_unit_in(i),'(A)') line
            end do
          end if

        end if
      end do

    end if


    ! main loop
    !
    DO imdstep = 1, option%nsteps

      ! read dcdfile(i) / write trjfile(j)
      !
      if (conv_crd) then

        if (mod(imdstep,option%crdout_period) == 0) then

          write(MsgOut,'(A,I10,A)') 'Convert> ', imdstep, ' step  read dcdfile > write trjfile'
          write(MsgOut,'(A,$)')     '         '

          do i = 1, option%num_replicas
            if (do_open_file(i)) then

              do_skip = .true.
              do j = 1, num_conv_id
                m = conv_id(j)
                if (m == param_index(i)) then
                  ! read trj
                  call read_trj(trj_in(i), trajectory)

                  ! selection
                  call reselect_atom(molecule,               &
                                     option%trjout_atom_exp, &
                                     trajectory%coord,       &
                                     option%trjout_atom,     &
                                     trjout_atom_trj)

                  ! centering
                  call centering(molecule, trajectory%coord, &
                                 option)

                  ! pbc correct
                  call run_pbc_correct(option%pbcc_mode,     &
                                       molecule,             &
                                       trajectory)

                  ! fitting
                  if (fitting%mass_weight) then
                    call run_fitting(fitting,                &
                                     molecule%atom_coord,    &
                                     trajectory%coord,       &
                                     trajectory%coord,       &
                                     molecule%mass)

                  else
                    call run_fitting(fitting,                &
                                     molecule%atom_coord,    &
                                     trajectory%coord,       &
                                     trajectory%coord)
                  end if

                  ! write trj(j)
                  if (mod(imdstep,option%trjout_period) == 0) then
                    call write_trj(trj_out(j), trajectory,     &
                                 trjout_atom_trj, molecule)
                    write(MsgOut,'(I5,A,I5,$)') trj_in(i)%unit_no,'>',trj_out(j)%unit_no
                  end if
                  do_skip = .false.

                end if
              end do

              if (do_skip) then
                ! read trj only
                call read_trj(trj_in(i), trajectory)
              end if

            end if
          end do

          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') ''

        end if

      end if


      ! read logfile(i) / write logfile(j)
      !
      if (conv_log) then

        if (mod(imdstep,option%eneout_period) == 0) then

          write(MsgOut,'(A,I10,A)') 'Convert> ', imdstep, ' step  read logfile > write logfile'
          write(MsgOut,'(A,$)')     '         '

          do i = 1, option%num_replicas
            if (do_open_file(i)) then

              do_skip = .true.
              do j = 1, num_conv_id
                m = conv_id(j)
                if (m == param_index(i)) then

                  done = .false.
                  do while (.not. done)

                    read(log_unit_in(i),'(A)',iostat=ierr) line

                    if (ierr < 0) then
                      write(0,*) trim(line)
                      write(MsgOut,'(A,I10,A)') 'Convert> Something is wrong around ', imdstep, ' step in logfile'
                      write(MsgOut,'(A,I10,A,I10)') 'Convert> error: ', ierr, &
                                                    ' unit: ', log_unit_in(i)
                      call error_msg('Convert> Convert is failed due to unrecognized lines in logfile')
                    end if

                    if (line(1:5) == 'INFO:') then
                      if (mod(imdstep,option%logout_period) == 0) then
                        write(log_unit_out(j),'(A)') trim(line)
                        write(log_unit_out(j),'(A)') ''
                      endif
                      done = .true.
                    end if

                  end do

                  do_skip = .false.
                  write(MsgOut,'(I5,A,I5,$)') log_unit_in(i),'>',log_unit_out(j)
                end if
              end do

              if (do_skip) then

                done = .false.
                do while (.not. done)
                  read (log_unit_in(i),'(A)') line
                  if (line(1:5) == 'INFO:') then
                    read (log_unit_in(i),'(A)') line
                    done = .true.
                  end if
                end do

              end if

            end if
          end do

          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') ''

        end if

      end if


      ! read enefile(i) / write enefile(j) for gREST
      !
      if (conv_ene) then

        if (mod(imdstep,option%crdout_period) == 0) then

          write(MsgOut,'(A,I10,A)') 'Convert> ', imdstep, ' step  read enefile > write enefile'
          write(MsgOut,'(A,$)')     '         '

          do i = 1, option%num_replicas
            if (do_open_file(i)) then

              do_skip = .true.
              do j = 1, num_conv_id
                m = conv_id(j)
                if (m == param_index(i)) then

                  done = .false.
                  do while (.not. done)

                    read(ene_unit_in(i),'(A)',iostat=ierr) line

                    if (ierr < 0) then
                      write(0,*) trim(line)
                      write(MsgOut,'(A,I10,A)') 'Convert> Something is wrong around ', imdstep, ' step in enefile'
                      write(MsgOut,'(A,I10,A,I10)') 'Convert> error: ', ierr, &
                                                    ' unit: ', ene_unit_in(i)
                      call error_msg('Convert> Convert is failed due to unrecognized lines in enefile')
                    end if


                    if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
                      write(ene_unit_out(j),'(A)') trim(line)
                      done = .true.
                    end if

                  end do

                  do_skip = .false.
                  write(MsgOut,'(I5,A,I5,$)') ene_unit_in(i),'>',ene_unit_out(j)
                end if
              end do

              if (do_skip) then

                done = .false.
                do while (.not. done)
                  read (ene_unit_in(i),'(A)') line
                  if (line(1:5) == 'INFO:') then
                    read (ene_unit_in(i),'(A)') line
                    done = .true.
                  end if
                end do

              end if

            end if
          end do

          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') ''

        end if

      end if


      ! read remfile(i)
      !
      if (mod(imdstep,option%exchange_period) == 0 .and. &
          imdstep /= option%nsteps) then

        if (option%convert_type == ConvertTypeParameter) then

          do i = 1, option%num_replicas
            if (do_open_file(i)) then
              read(rem_unit_in(i),*) istep, param_index(i)
            end if
          end do

        else if (option%convert_type == ConvertTypeReplica) then

          do i = 1, option%num_replicas
            if (do_open_file(i)) then
              param_index(i) = i
            end if
          end do

        end if

      end if

    END DO


    ! close files
    !
    do i = 1, option%num_replicas
      if (do_open_file(i)) then
        if (option%convert_type == ConvertTypeParameter) &
          call close_file(rem_unit_in(i))
        if (conv_crd) call close_trj(trj_in(i))
        if (conv_log) call close_file(log_unit_in(i))
        if (conv_ene) call close_file(ene_unit_in(i))
      end if
    end do

    do j = 1, num_conv_id
      if (conv_crd) call close_trj(trj_out(j))
      if (conv_log) call close_file(log_unit_out(j))
      if (conv_ene) call close_file(ene_unit_out(j))
    end do


    return

  end subroutine convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_convert
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_convert(input, output, option, conv_crd, conv_log, conv_ene)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option
    logical,                 intent(out)   :: conv_crd
    logical,                 intent(out)   :: conv_log
    logical,                 intent(out)   :: conv_ene


    conv_crd = .false.
    conv_log = .false.
    conv_ene = .false.

    ! check remfile, dcdfile, logfile existing
    !
    if (option%convert_type == ConvertTypeParameter) then
      if (input%remfile == '') then
        call error_msg('Setup_Convert> Error: remfile is not specified in [INPUT]')
      end if
    end if

    if ((input%dcdfile == '' .and. output%trjfile /= '') .or. &
        (input%dcdfile /= '' .and. output%trjfile == '')) then
      call error_msg('Setup_Convert> dcdfile in [INPUT] or trjfile in [OUTPUT] is not specified')
    else if (input%dcdfile == '' .and. output%trjfile == '') then
      conv_crd = .false.
    else
      conv_crd = .true.
    end if

    if ((input%logfile == '' .and. output%logfile /= '') .or. &
        (input%logfile /= '' .and. output%logfile == '')) then
      call error_msg('Setup_Convert> logfile in [INPUT] or logfile in [OUTPUT] is not specified')
    else if (input%logfile == '' .and. output%logfile == '') then
      conv_log = .false.
    else
      conv_log = .true.
    end if

    if ((input%enefile == '' .and. output%enefile /= '') .or. &
        (input%enefile /= '' .and. output%enefile == '')) then
      call error_msg('Setup_Convert> enefile in [INPUT] or enefile in [OUTPUT] is not specified')
    else if (input%enefile == '' .and. output%enefile == '') then
      conv_ene = .false.
    else
      conv_ene = .true.
    end if

    if (conv_crd .and. option%crdout_period == 0) &
      call error_msg('Setup_Convert> error: crdout_period in [OPTION] = 0')

    if (conv_log .and. option%eneout_period == 0) &
      call error_msg('Setup_Convert> error: eneout_period in [OPTION] = 0')

    if (conv_ene .and. option%crdout_period == 0) &
      call error_msg('Setup_Convert> error: crdout_period in [OPTION] = 0')

    if (option%nsteps == 0) &
      call error_msg('Setup_Convert> error: nsteps in [OPTION] = 0')

    if (option%num_replicas == 0) &
      call error_msg('Setup_Convert> error: num_replicas in [OPTION] = 0')

    if (.not. conv_crd .and. .not. conv_log .and. .not. conv_ene) then
      write(MsgOut,'(A)') 'Setup_Convert> Nothing was done, because dcdfile, trjfile, logfile, and enefile are empty'
      write(MsgOut,'(A)') ''
      return
    end if

    return

  end subroutine setup_convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    centering
  !> @brief        move the COM of the fitting target to the origin
  !! @authors      DM
  !! @param[in]    molecule   : molecule information
  !! @param[inout] coord      : atom coordinates
  !! @param[in]    fitting    : fitting information
  !! @param[in]    option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine centering(molecule, coord, option)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    real(wp),         intent(inout) :: coord(:,:)
    type(s_option),   intent(in)    :: option

    ! local variables
    integer  :: iatm, natm
    real(wp) :: com(3)

    if (.not. option%centering) return

    natm = size(molecule%atom_no)
    com  = compute_com(coord, molecule%mass, option%centering_atom%idx)

    do iatm = 1, natm
      coord(:,iatm) = coord(:,iatm) - com(:) + option%center_coord(:)
    end do

    return

  end subroutine centering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_replicate_name
  !> @brief        get replicate name
  !! @authors      NT
  !! @param[in]    filename   : filename
  !! @param[in]    replica_no : replica number
  !! @return       replicate name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name(filename, replica_no)

    ! return
    character(Maxfilename)   :: get_replicate_name

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: replica_no

    ! local variables
    integer                  :: bl, br
    character(MaxLine)       :: bl_str, br_str
    character(20)            :: fmt
    integer                  :: i, ndigit, comp


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0) then
      get_replicate_name = filename
      call error_msg('Get_Replicate_Name> Error: "{}" is not found in the file name')
    end if

    if (bl > br) &
      call error_msg('Get_Replicate_Name> Syntax error.')

    if (bl > 1) then
      bl_str = filename(1:bl-1)
    else
      bl_str = ''
    end if

    if (br < len_trim(filename)) then
      br_str = filename(br+1:)
    else
      br_str = ''
    end if
    
    do i = 1, 100
      comp = 10**i
      if (replica_no < comp) then
        ndigit = i
        exit
      end if
    end do
    
    write(fmt,'(a,i0,a,i0,a)') '(a,i', ndigit, '.', ndigit, ',a)'
    write(get_replicate_name,fmt) trim(bl_str), replica_no, trim(br_str)

    return

  end function get_replicate_name

end module rc_convert_mod
