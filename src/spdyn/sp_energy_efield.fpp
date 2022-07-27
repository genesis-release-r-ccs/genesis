!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_efield_mod
!> @brief   calculate energy with external electric field
!! @authors Jaewoon Jung (JJ), Yuji Sugia (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_efield_mod

  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_efield

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_efield
  !> @brief        calculate energy and force with external electric field
  !! @authors      JJ, YS
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] efield  : bond energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_efield(domain, enefunc, coord, force, virial, &
                                   efield)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: efield(nthread)

    ! local variables
    real(wp)                 :: efield_temp, work(3), center(3), viri(3)
    integer                  :: i, ix
    integer                  :: id, omp_get_thread_num

    integer,         pointer :: natom(:), ncell_local
    real(wp),        pointer :: charge(:,:)
    real(wp)                 :: field(3)


    ncell_local => domain%num_cell_local
    natom       => domain%num_atom
    charge      => domain%charge

    field(1:3)  = enefunc%efield(1:3)*23.06_wp 
    if (enefunc%efield_normal) &
    field(1:3)  = field(1:3)*domain%system_size_ini(1:3)/domain%system_size(1:3)

    center(1:3) = domain%system_size(1:3) / 2.0_wp

    ! calculate electric field energy
    !
    !$omp parallel default(shared)                      &
    !$omp private(id, i, ix, efield_temp, work)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      efield_temp = 0.0_wp
      viri(1:3) = 0.0_wp

      do ix = 1, natom(i)

        work(1) = field(1)*charge(ix,i)
        work(2) = field(2)*charge(ix,i)
        work(3) = field(3)*charge(ix,i)
        efield_temp = efield_temp - work(1)*(coord(1,ix,i)+center(1))
        efield_temp = efield_temp - work(2)*(coord(2,ix,i)+center(2))
        efield_temp = efield_temp - work(3)*(coord(3,ix,i)+center(3))
        force(1,ix,i,id+1) = force(1,ix,i,id+1) + work(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) + work(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) + work(3)
        viri(1) = viri(1) + coord(1,ix,i)*work(1)
        viri(2) = viri(2) + coord(2,ix,i)*work(2)
        viri(3) = viri(3) + coord(3,ix,i)*work(3)

      end do

      if (enefunc%efield_virial) then
        virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) + viri(3)
      end if
      efield(id+1) = efield(id+1) + efield_temp

    end do

    !$omp end parallel 
    
    return

  end subroutine compute_energy_efield

end module sp_energy_efield_mod
