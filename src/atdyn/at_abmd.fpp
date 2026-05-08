!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_abmd_mod
!> @brief   ABMD method
!           Refs: 
!! @authors Junichiro Koga (JK), Shingo Ito (SI)
! 
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

module at_abmd_mod
  use timers_mod
  use at_dynamics_str_mod
  use at_enefunc_str_mod
  use at_enefunc_reac_coord_mod
  use at_enefunc_abmd_mod
  !use at_domain_str_mod 
  use at_remd_mod 
  use at_remd_str_mod
  use at_bias_mod 
  use at_bias_str_mod  
  use fileio_rst_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use bspline_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  type, public :: s_abmd_info
    logical  :: abmd
    logical  :: wt    
    real(wp) :: tauf
    real(wp) :: wt_temperature    
    real(wp) :: u0
    real(wp) :: sigma
    logical  :: update_biasing_potential
    logical  :: biasing_potential_force
    logical  :: stop_when_out_of_range
    logical  :: reduce_bias_potential
    integer  :: update_period
    integer  :: monitor_period
    integer  :: monitor_period_spline
  end type s_abmd_info

  ! public subroutines
  public :: show_ctrl_abmd
  public :: read_ctrl_abmd
  public :: setup_abmd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_abmd
  !> @brief        show ABMD section usage
  !! @authors      JK
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "abmd"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_abmd(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode

    if (show_all) then
      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[ABMD]'
        write(MsgOut,'(A)') 'abmd          = no    # Adaptively-biased MD (ABMD) [yes, no]'
        write(MsgOut,'(A)') 'wt            = no    # Well-Tempered (WT) [yes, no]'

      end select
    else

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[ABMD]'
        write(MsgOut,'(A)') 'abmd          = no    # Adaptively-biased MD (ABMD) [yes, no]'
        write(MsgOut,'(A)') 'wt            = no    # Well-Tempered (WT) [yes, no]'

      end select

    end if

    return
  end subroutine show_ctrl_abmd
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_abmd
  !> @brief        read ABMD section in the control file
  !! @authors      JK
  !! @param[in]    handle    : unit number of control file
  !! @param[out]   abmd_info : ABMD section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_abmd(handle, abmd_info, remd_info)

    ! parameters
    character(*),            parameter     :: Section = 'ABMD'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_abmd_info),       intent(inout) :: abmd_info
    type(s_rep_info), optional, intent(in) :: remd_info
    abmd_info%abmd                     = .false.
    abmd_info%wt                       = .false.
    abmd_info%tauf                     = 0.0_wp
    abmd_info%wt_temperature           = 0.0_wp
    abmd_info%u0                       = 0.0_dp
    abmd_info%sigma                    = 0.5_dp
    abmd_info%update_biasing_potential = .true.
    abmd_info%stop_when_out_of_range   = .true.
    abmd_info%reduce_bias_potential    = .true.
    if(present(remd_info)) then
      abmd_info%reduce_bias_potential  = remd_info%exchange_period == 0
    endif
    abmd_info%update_period            = 1
    abmd_info%monitor_period           = 1
    abmd_info%monitor_period_spline    = 1
    call begin_ctrlfile_section(handle, Section)
    !call read_ctrlfile_logical (handle, Section,                  &
    !                   'abmd', abmd_info%abmd)
    !call read_ctrlfile_logical (handle, Section,                  &
    !                   'wt', abmd_info%wt)
    call read_ctrlfile_real    (handle, Section, 'tau',           &
                                abmd_info%tauf)
    call read_ctrlfile_real    (handle, Section, 'wt_temperature', &
                                abmd_info%wt_temperature)
    call read_ctrlfile_integer (handle, Section, 'monitor_period', &
                                abmd_info%monitor_period)                                                                  
    !! For debug
    call read_ctrlfile_real    (handle, Section, 'u0',            &
                                abmd_info%u0)
    call read_ctrlfile_real    (handle, Section, 'sigma',         &
                                abmd_info%sigma)
    call read_ctrlfile_integer (handle, Section, 'update_period', &
                                abmd_info%update_period)    
    call read_ctrlfile_integer (handle, Section, 'monitor_period_spline', &
                                abmd_info%monitor_period_spline)
    call read_ctrlfile_logical (handle, Section, 'reduce_bias_potential', &
                                abmd_info%reduce_bias_potential)
    if(abmd_info%u0>0) abmd_info%update_biasing_potential = .false.
    call read_ctrlfile_logical (handle, Section, 'update_biasing_potential', &
                                abmd_info%update_biasing_potential)
    call read_ctrlfile_logical (handle, Section, 'stop_when_out_of_range', &
                                abmd_info%stop_when_out_of_range)
    call end_ctrlfile_section(handle)
    
    if (abmd_info%tauf > EPS)           abmd_info%abmd = .true.
    if (abmd_info%wt_temperature > EPS) abmd_info%wt   = .true.
    
    if (main_rank) then
      if (abmd_info%abmd) then
        write(MsgOut,'(a)')       'Read_Ctrl_ABMD> Parameters of ABMD'
        write(MsgOut,'(a)')       '  ABMD enabled'
        write(MsgOut,'(a,f10.2)') '  tau             = ', abmd_info%tauf
        if(abmd_info%wt) then
          write(MsgOut,'(a)')       '  Well-Tempered method enabled'
          write(MsgOut,'(a,f10.2)') '  wt_temperature  = ', abmd_info%wt_temperature
        end if
        write(MsgOut,'(a)')
        !write(MsgOut,'(a,i8)')    '  update_period   = ', abmd_info%update_period
        write(MsgOut,'(a,i8)')    '  monitor_period  = ', abmd_info%monitor_period
        !write(MsgOut,'(a,i8)')    '  monitor_period(spline coeffs)  = ', abmd_info%monitor_period_spline
        if(abmd_info%u0>0) then
          write(MsgOut,'(a,f10.2)') '  u0              = ', abmd_info%u0
          write(MsgOut,'(a,f10.2)') '  sigma           = ', abmd_info%sigma
        endif
        if(present(remd_info)) then
        write(MsgOut,'(a,l2)') '  reduce_bias_potential = ', &
                           abmd_info%reduce_bias_potential
        endif
      end if
    end if

    !if (main_rank) then 
    !  if (.not. abmd_info%abmd) then
    !    write(MsgOut,'(a)')       'Read_Ctrl_ABMD> ABMD is not enabled'
    !  else
    !    write(MsgOut,'(a)')       'Read_Ctrl_ABMD> Parameters of ABMD'
    !    write(MsgOut,'(a)')       '  ABMD enabled'
    !    if(abmd_info%wt) then
    !      write(MsgOut,'(a)')       '  Well-Tempered method enabled'
    !      write(MsgOut,'(a,f10.2)') '  wt_temperature  = ', abmd_info%wt_temperature
    !    else
    !      write(MsgOut,'(a)')       '  Well-Tempered method disabled'
    !    end if
    !    write(MsgOut,'(a,f10.2)') '  tau             = ', abmd_info%tauf
    !    write(MsgOut,'(a)')
    !    write(MsgOut,'(a,i8)')    '  update_period   = ', abmd_info%update_period
    !    write(MsgOut,'(a,i8)')    '  monitor_period  = ', abmd_info%monitor_period
    !    write(MsgOut,'(a,i8)')    '  monitor_period(spline coeffs)  = ', abmd_info%monitor_period_spline
    !    if(abmd_info%u0>0) then
    !      write(MsgOut,'(a,f10.2)') '  u0              = ', abmd_info%u0
    !      write(MsgOut,'(a,f10.2)') '  sigma           = ', abmd_info%sigma
    !    endif
    !    if(present(remd_info)) then
    !    write(MsgOut,'(a,l2)') '  reduce_bias_potential = ', &
    !                       abmd_info%reduce_bias_potential
    !    endif
    !  endif
    !endif

  end subroutine read_ctrl_abmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_abmd
  !> @brief        setup ABMD
  !! @authors      JK, SI
  !! @param[in]    abmd_info : ABMD section control parameters information
  !! @param[in]    spline    : information of spline coeffs read from file
  !! @param[inout] dynamics  : dynamics information
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[inout] remd      : REMD information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_abmd(abmd_info, rst, dynamics, enefunc, remd, bias)

    ! formal arguments
    type(s_abmd_info),       intent(in) :: abmd_info
    type(s_rst),             intent(in)    :: rst
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_enefunc),target,  intent(inout) :: enefunc
    type(s_remd),  optional, intent(inout) :: remd
    type(s_bias),  optional, intent(inout) :: bias 

    ! local variables
    type(s_enefunc_abmd),       pointer :: abmd
    type(s_enefunc_reac_coord), pointer :: reac_coord
    integer                             :: i, j, ixi, inext, max_nwalkers
    integer                             :: inum, jnum, dockid, num_groups, max_atoms
    real(kind=wp)                       :: xi, xi0, sigma, width
    real(kind=dp)                       :: utmp
    integer, pointer                    :: indices(:)
    logical                             :: match

    enefunc%abmd_use = abmd_info%abmd

    if (.not. enefunc%abmd_use) then
      return
    endif

    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_ABMD> ABMD information'
      write(MsgOut,'(A)') ''
    end if

    if (enefunc%reac_coord%nreac_coord .eq. 0) then
      call error_msg('Setup_Abmd> Reaction coordinate has not been specified.')
    end if

    abmd                          => enefunc%abmd
    reac_coord                    => enefunc%reac_coord
    abmd%wt                       =  abmd_info%wt
    abmd%wt_temperature           =  abmd_info%wt_temperature    
    abmd%tauf                     =  abmd_info%tauf
    abmd%update_period            =  abmd_info%update_period
    abmd%monitor_period           =  abmd_info%monitor_period
    !abmd%monitor_period_spline    =  abmd_info%monitor_period_spline
    abmd%monitor_period_spline    =  dynamics%rstout_period
    abmd%stop_when_out_of_range   =  abmd_info%stop_when_out_of_range
    abmd%update_biasing_potential =  abmd_info%update_biasing_potential
    abmd%reduce_bias_potential    =  abmd_info%reduce_bias_potential

    if(present(remd)) then
      remd%reduce_bias_potential  =  abmd%reduce_bias_potential
    endif
    
    !! Copy parameters from enefunc%reac_coord to enefunc%abmd
    if ( present (bias) ) then
      dockid                = my_country_no + 1
      max_nwalkers          = bias%nwalkers
      jnum                  = bias%dockid2repsetid(dockid)
      inum                  = bias%repsetid2paramid(jnum,1)
      abmd%nreac_coord      = reac_coord%nreac_func
      abmd%nreac_coord_comb = 1
      do i = 1, abmd%nreac_coord
        abmd%nreac_coord_comb = &
        abmd%nreac_coord_comb * (reac_coord%nxi(bias%ifunc(inum,i))+2)
      enddo

      call alloc_enefunc(enefunc, EneFuncABMD, &
                         enefunc%reac_coord%nreac_coord, &
                         maxval(enefunc%reac_coord%nxi), &
                         abmd%nreac_coord_comb, &
                         bias%ndocks)
    else
      abmd%nreac_coord      = reac_coord%nreac_coord
      abmd%nreac_coord_comb = reac_coord%nreac_coord_comb

      call alloc_enefunc(enefunc, EneFuncABMD, &
                         enefunc%reac_coord%nreac_coord, &
                         maxval(enefunc%reac_coord%nxi), &
                         abmd%nreac_coord_comb)
    end if
    abmd%bondslist             = 0
    abmd%bondslist_to_atomlist = 0

    if (present (bias)) then
      abmd%dockid2repsetid  = bias%dockid2repsetid
      abmd%repsetid2paramid = bias%repsetid2paramid
      abmd%ifunc            = bias%ifunc
      do i = 1, abmd%nreac_coord
        abmd%nxi(i)         = reac_coord%nxi(abmd%ifunc(inum,i))
        abmd%dxi(i)         = reac_coord%dxi(abmd%ifunc(inum,i))
        abmd%invdxi(i)      = reac_coord%invdxi(abmd%ifunc(inum,i))
        abmd%ximin(i)       = reac_coord%ximin(abmd%ifunc(inum,i))
        abmd%ximax(i)       = reac_coord%ximax(abmd%ifunc(inum,i))
        abmd%kind(i)        = reac_coord%kind(abmd%ifunc(inum,i))
        abmd%periodicity(i) = reac_coord%periodicity(abmd%ifunc(inum,i))
        abmd%bondslist(:,:,i) = reac_coord%bondslist(:,:,bias%ifunc(inum,i))
        abmd%bondslist_to_atomlist(:,i) = &
         reac_coord%bondslist_to_atomlist(:,bias%ifunc(inum,i))
      end do
    else
      do i = 1, abmd%nreac_coord
        abmd%ifunc(1,i)   = i
      end do
      abmd%nxi            = reac_coord%nxi
      abmd%dxi            = reac_coord%dxi
      abmd%invdxi         = reac_coord%invdxi
      abmd%ximin          = reac_coord%ximin
      abmd%ximax          = reac_coord%ximax
      abmd%kind           = reac_coord%kind
      abmd%periodicity    = reac_coord%periodicity    
      abmd%bondslist      = reac_coord%bondslist
      abmd%bondslist_to_atomlist = reac_coord%bondslist_to_atomlist 
    end if
    !! Read restart files 
    if (abmd_info%u0>0) then ! debug mode
      indices                       => reac_coord%indices
      do i = 1, reac_coord%nreac_coord_comb
        call bspline_resolve_indices(i, enefunc%reac_coord%nreac_coord, &
                                        enefunc%reac_coord%nxi, indices)
        utmp = 1.0_dp
        do j = 1, reac_coord%nreac_coord
          if (indices(j)==-1 .or. indices(j) == reac_coord%nxi(j)) cycle
          ixi                = indices(j)
          width              = reac_coord%ximax(j) - reac_coord%ximin(j)
          xi0                = reac_coord%ximin(j) + width * 0.5_wp
          xi                 = reac_coord%ximin(j) + ixi * reac_coord%dxi(j)
          sigma              = width * abmd_info%sigma
          utmp               = utmp  * gauss(xi, xi0, sigma)
        enddo
        enefunc%abmd%uk(i) = abmd_info%u0 * utmp * KBOLTZ
      enddo
      do i = 1, reac_coord%nreac_coord_comb
        call bspline_resolve_indices(i, enefunc%reac_coord%nreac_coord, &
                                        enefunc%reac_coord%nxi, indices)
        do j = 1, reac_coord%nreac_coord
          if (indices(j) == -1 ) then
            inext = i + 2 * bspline_get_next_index_offset(&
                            j, reac_coord%nreac_coord, reac_coord%nxi)
            enefunc%abmd%uk(i) = enefunc%abmd%uk(inext)
          else if (indices(j) == reac_coord%nxi(j)) then
            inext = i - 2 * bspline_get_next_index_offset(&
                            j, reac_coord%nreac_coord, reac_coord%nxi)
            enefunc%abmd%uk(i) = enefunc%abmd%uk(inext)
          endif
        enddo
      enddo
    endif
    ! obtain b-spline coefficients from bias potential read from file
    if (allocated(rst%abmd_nxi) .and. allocated(rst%abmd_dxi) .and. &
        allocated(rst%abmd_ubuf)) then
      ! resolve bspline coeffs on the previous mesh
      if(rst%abmd_nreac_coord /= abmd%nreac_coord) then
	 call error_msg('Setup_Abmd> Number of reaction coordinates '//&
	                'does not match the bias file.')
      endif
      match = rst%abmd_nreac_coord_comb == abmd%nreac_coord_comb
      if(match) then
        do i = 1, abmd%nreac_coord
          if(rst%abmd_nxi(i) /= abmd%nxi(i)) then
            match = .false.
            exit
          endif
        enddo
      endif
     
      if(match) then
        ! if previous mesh == present mesh, then just copy the bspline coeffs
        abmd%uk(:) = rst%abmd_ubuf(:)
      else
        ! bspline-interpolate on to the new mesh
        call bspline_interpolation(rst%abmd_nreac_coord_comb,     &
                                   abmd%nreac_coord_comb, &
                                   abmd%nreac_coord,      &
                                   rst%abmd_nxi, abmd%nxi,  &
                                   abmd%ximin,            &
                                   rst%abmd_dxi, abmd%dxi,  &
                                   rst%abmd_ubuf, abmd%uk)

        ! resolve bspline coeffs on the present mesh
        call bspline_solve_for_coeffs(abmd%nreac_coord_comb,            &
                                      abmd%nreac_coord, abmd%nxi, &
                                      abmd%uk, positive = .true.)
      endif
    endif
    
    return

    contains

    function gauss(x, x0, sig) result(y)
      implicit none
      real(wp), intent(in)  :: x, x0, sig
      real(dp) :: y
      y = exp(-((x-x0)/sig)**2)
    end function gauss

  end subroutine setup_abmd

end module at_abmd_mod
