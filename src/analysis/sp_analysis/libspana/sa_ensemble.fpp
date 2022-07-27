!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_ensemble_mod
!> @brief   module for analysis ensembles
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_ensemble_mod

  use sa_ensemble_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_ens_info
    integer    :: ensemble = EnsembleNVE
  end type s_ens_info

  ! subroutines
  public  :: show_ctrl_ensemble
  public  :: read_ctrl_ensemble
  public  :: setup_ensemble

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_ensemble_
  !> @brief        show control parameters in OPTION section
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_ensemble

    write(MsgOut,'(A)') '[ENSEMBLE]'
    write(MsgOut,'(A)') 'ensemble     = NVE             # (NVT/NPT/NVE)'
    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_ensemble

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_ensemble
  !> @brief        read control parameters in OPTION section
  !! @authors      IY
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   ens_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_ensemble(handle, ens_info)

    ! parameters
    character(*),            parameter     :: Section = 'Ensemble'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_ens_info),        intent(inout) :: ens_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)


    call read_ctrlfile_type(handle, Section, 'ensemble', &
                            ens_info%ensemble, EnsembleTypes)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if(main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_ENSEMBLE> Parameters of ENSEMBLE'
      write(MsgOut,'(A20,A)') '  ensemble        = ', &
                                 trim(EnsembleTypes(ens_info%ensemble))
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_ensemble

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_ensemble
  !> @brief        setup ensemble information
  !! @authors      IY
  !! @param[in]    ens_info : ENSEMBLE section control parameters information
  !! @param[inout] ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_ensemble(ens_info, ensemble)

    ! formal arguments
    type(s_ens_info),        intent(in)    :: ens_info
    type(s_ensemble),        intent(inout) :: ensemble


    ensemble%ensemble = ens_info%ensemble

    return

  end subroutine setup_ensemble

end module sa_ensemble_mod
