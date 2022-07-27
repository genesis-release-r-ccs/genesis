!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_minfo_mod
!> @brief   read/write data for vibrational analysis
!! @authors Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_minfo_mod

  use atom_libs_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use molecules_str_mod

  implicit none
  private

  ! structures
  type, public :: s_minfo

    logical  :: atom_data   = .false.
    logical  :: elec_data   = .false.
    logical  :: vib_data    = .false.

    integer  :: nat = 0
    integer  :: nat_sub = 0
    integer, allocatable :: vibatom_id(:)
    integer, allocatable :: subatom_id(:)

    real(wp) :: energy
    real(wp), allocatable :: gradient(:), hessian(:)
    real(wp), allocatable :: dipole(:), dipole_derv(:,:)

    real(wp), allocatable :: omega(:), vec(:,:)

  end type s_minfo

  ! subroutines
  public  :: init_minfo_atom
  public  :: init_minfo_elec
  public  :: init_minfo_vib
  public  :: output_minfo
  public  :: read_minfo

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_minfo_atom
  !> @brief        initialize minfo for Atomic data
  !! @authors      KY
  !! @param[in]    nat     : number of atoms
  !! @param[out]   minfo   : minfo data
  !! @param[in]    nat_sub : number of sub_atoms (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_minfo_atom(nat, minfo, nat_sub)

    ! formal arguments
    integer      ,           intent(in)    :: nat
    type(s_minfo),           intent(inout) :: minfo
    integer      , optional, intent(in)    :: nat_sub


    minfo%nat = nat
    allocate(minfo%vibatom_id(nat))

    if(present(nat_sub)) then
      minfo%nat_sub = nat_sub
      allocate(minfo%subatom_id(nat_sub))
    end if

    minfo%atom_data   = .true.
    return

  end subroutine init_minfo_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_minfo_elec
  !> @brief        initialize minfo for Electronic data
  !! @authors      KY
  !! @param[in]    nat     : number of atoms
  !! @param[out]   minfo   : minfo data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_minfo_elec(nat, minfo)

    ! formal arguments
    integer      ,           intent(in)    :: nat
    type(s_minfo),           intent(inout) :: minfo

    ! local variables
    integer :: nat3


    minfo%nat    = nat
    minfo%energy = 0.0_wp
    nat3 = nat*3
    allocate(minfo%gradient(nat3))
    allocate(minfo%hessian(nat3*(nat3+1)/2))
    allocate(minfo%dipole(3))
    allocate(minfo%dipole_derv(3,nat3))

    minfo%gradient    = 0.0_wp
    minfo%hessian     = 0.0_wp
    minfo%dipole      = 0.0_wp
    minfo%dipole_derv = 0.0_wp
    minfo%elec_data   = .true.

    return

  end subroutine init_minfo_elec

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_minfo_vib
  !> @brief        initialize minfo for Vibrational data
  !! @authors      KY
  !! @param[in]    nfree   : number of vibrational degrees of freedom
  !! @param[out]   minfo   : minfo data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_minfo_vib(nfree, minfo)

    ! formal arguments
    integer      ,           intent(in)    :: nfree
    type(s_minfo),           intent(inout) :: minfo


    allocate(minfo%omega(nfree))
    allocate(minfo%vec(nfree, nfree))

    minfo%vib_data    = .true.

    return

  end subroutine init_minfo_vib

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_minfo
  !> @brief        output minfo data to a file
  !! @authors      KY
  !! @param[in]    fname    : name of the file
  !! @param[in]    minfo    : minfo data
  !! @param[inout] molecule : molecule data
  !! @param[in]    coord    : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_minfo(fname, minfo, molecule, coord)

    ! formal arguments
    character(*),            intent(in)     :: fname
    type(s_minfo),           intent(in)     :: minfo
    type(s_molecule),        intent(inout), optional :: molecule
    real(wp),                intent(in),    optional :: coord(:,:)

    ! local variables
    integer                  :: ifile
    integer                  :: i, j, k, iatom, ncolumn, isize
    integer                  :: atomic_no
    integer                  :: nat
    real(wp)                 :: mm_mass

    integer                  :: nd
    integer, allocatable     :: domain_nat(:), domain_nf(:)
    integer, allocatable     :: domain_idx(:,:)


    ifile   = 10
    ncolumn = 5
    nat     = minfo%nat

    call open_file(ifile, fname, IOFileOutputReplace)
    write(ifile,'(''# minfo File version 2:'')')
    write(ifile,'(''#'')')

    if (minfo%atom_data .and. present(molecule) .and. present(coord)) then
      write(ifile,'(''[ Atomic Data ]'')')
      write(ifile,'(i5)') nat
      do i = 1, nat
        iatom = minfo%vibatom_id(i)
        mm_mass = molecule%mass(iatom)
        !call qm_atomic_number(mm_mass, atomic_no)
        !write(MsgOut,'(2i8,a6)') i, iatom, molecule%atom_cls_name(iatom)
        atomic_no = atomic_number(mm_mass, molecule%atom_cls_name(iatom))
        write(ifile,'(a6,'', '',i4,'', '',f12.4,'', '',2(f17.10,'', ''),f17.10)') &
          trim(molecule%atom_name(iatom)),  &
          atomic_no,                  &
          mm_mass,                    &
          coord(:,iatom) / CONV_UNIT_LEN
      end do
      if (minfo%nat_sub > 0) then 
        write(ifile,'(i5)') minfo%nat_sub
        do i = 1, minfo%nat_sub
          iatom = minfo%subatom_id(i)
          mm_mass = molecule%mass(iatom)
          !call qm_atomic_number(mm_mass, atomic_no)
          !write(MsgOut,'(2i8,a6)') i, iatom, molecule%atom_cls_name(iatom)
          atomic_no = atomic_number(mm_mass, molecule%atom_cls_name(iatom))
          write(ifile,'(a6,'', '',i4,'', '',f12.4,'', '',2(f17.10,'', ''),f17.10)') &
            trim(molecule%atom_name(iatom)),  &
            atomic_no,                  &
            mm_mass,                    &
            coord(:,iatom) / CONV_UNIT_LEN
        end do
      end if
      write(ifile,*)
    end if

    if (minfo%elec_data) then
      write(ifile,'(''[ Electronic Data ]'')')
      write(ifile,'(''Energy'')')
      write(ifile,'(f25.14)') minfo%energy

      if (any(minfo%gradient /= 0.0_wp)) then
        write(ifile,'(''Gradient'')')
        isize = nat*3
        write(ifile,'(i5)') isize
        k = 1
        do i = 1, nat
          do j = 1, 3
            write(ifile,'(es15.8,$)') minfo%gradient(k)
            if (mod(k,ncolumn) == 0 .or. k == isize) then
              write(ifile,*)
            else
              write(ifile,'('', '',$)')
            end if
            k = k + 1
          end do
        end do
      end if

      if (any(minfo%hessian /= 0.0_wp)) then
        write(ifile,'(''Hessian'')')
        isize = nat*3*(nat*3+1)/2
        write(ifile,'(i0)') isize
        do k = 1, isize
          write(ifile,'(es15.8,$)') minfo%hessian(k)
          if (mod(k,ncolumn) == 0 .or. k == isize) then
            write(ifile,*)
          else
            write(ifile,'('', '',$)')
          end if
        end do
      end if

      if (any(minfo%dipole /= 0.0_wp)) then
        write(ifile,'(''Dipole Moment'')')
        isize = 3
        write(ifile,'(i5)') isize
        do k = 1, isize
          write(ifile,'(es15.8,$)') minfo%dipole(k)
          if(mod(k,ncolumn) == 0 .or. k == isize) then
            write(ifile,*)
          else
            write(ifile,'('', '',$)')
          end if
        end do
      end if

      if (any(minfo%dipole_derv /= 0.0_wp)) then
        write(ifile,'(''Dipole Derivative'')')
        isize = 3*nat*3
        write(ifile,'(i5)') isize
        k = 1
        do i = 1, nat*3
          do j = 1, 3
            write(ifile,'(es15.8,$)') minfo%dipole_derv(j, i)
            if(mod(k,ncolumn) == 0 .or. k == isize) then
              write(ifile,*)
            else
              write(ifile,'('', '',$)')
            end if
            k = k + 1
          end do
        end do
      end if
      write(ifile,*)

    end if

    if (minfo%vib_data) then
      nd = 1
      allocate(domain_nat(nd), domain_nf(nd), domain_idx(nat,nd))
      domain_nat(1) = nat
      domain_nf(1)  = nat*3
      do i = 1, nat
        domain_idx(i,1) = i
      end do

      write(ifile,'(''[ Vibrational Data ]'')')
      write(ifile,'(''Number of Domain'')')
      write(ifile,'(i4)') nd
      do i = 1, nd
        write(ifile,'(''Domain '',i4)') nd
        write(ifile,'(''Atom Index'')')
        isize = domain_nat(i)
        write(ifile,'(i4)') isize
        do j = 1, isize
          write(ifile,'(i15,$)') domain_idx(j,i)
          if (mod(j,ncolumn) == 0 .or. j == isize) then
            write(ifile,*)
          else
            write(ifile,'('', '',$)')
          end if
        end do

        write(ifile,'(''Local Normal modes'')')
        write(ifile,'(''Vibrational Frequency'')')
        isize = domain_nf(i)
        write(ifile,'(i4)') isize
        do j = 1, isize
          write(ifile,'(es15.8,$)') minfo%omega(j)
          if (mod(j,ncolumn) == 0 .or. j == isize) then
            write(ifile,*)
          else
            write(ifile,'('', '',$)')
          end if
        end do

        write(ifile,'(''Vibrational vector'')')
        do j = 1, domain_nf(i)
          write(ifile,'(''Mode '',i6)') j
          isize = domain_nat(i)*3
          write(ifile,'(i4)') isize
          do k = 1, isize
            write(ifile,'(es15.8,$)') minfo%vec(k,j)
            if (mod(k,ncolumn) == 0 .or. k == isize) then
              write(ifile,*)
            else
              write(ifile,'('', '',$)')
            end if
          end do
        end do

      end do
      write(ifile,*)

    end if

    call close_file(ifile)

    return

  end subroutine output_minfo

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_minfo
  !> @brief        read minfo data to a file
  !! @authors      KY
  !! @param[in]    fname    : name of the file
  !! @param[out]   minfo    : minfo data
  !! @param[in]    molecule : molecule data
  !! @param[in]    coord    : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_minfo(fname, minfo, molecule, coord)

    ! formal arguments
    character(*),            intent(in)    :: fname
    type(s_minfo),           intent(out)   :: minfo
    type(s_molecule),        intent(inout),optional :: molecule
    real(wp),                intent(in),   optional :: coord(:,:)

    ! local variables
    integer        :: ifile, isize
    integer        :: i, j, k
    character(120) :: line
    real(wp)       :: vibatm_mass, vibatm_coord(3), dd
    integer        :: atom_num
    character(6)   :: atom_name
    logical        :: match

    real(wp), parameter :: thresh = 1.0e-03


    call open_file(ifile, fname, IOFileInput)
    read(ifile,*)
    read(ifile,*)

    do while(.true.)
      read(ifile,'(a)',end=10) line

      if (index(line,'Atomic Data') > 0 .and. present(molecule) &
         .and. present(coord)) then
        minfo%atom_data = .true.
        read(ifile,*) minfo%nat
        allocate(minfo%vibatom_id(minfo%nat))
        do i = 1, minfo%nat
          read(ifile,*) atom_name, atom_num, vibatm_mass, vibatm_coord
          vibatm_coord = vibatm_coord * CONV_UNIT_LEN

          do j = 1, molecule%num_atoms
            if (trim(atom_name) == trim(molecule%atom_name(j))) then 
              match = .true.
              do k = 1, 3
                dd = abs(coord(k,j) - vibatm_coord(k))/vibatm_coord(k)
                if (dd > thresh) then
                  match = .false.
                  exit
                end if
              end do

              if (match) then
                minfo%vibatom_id(i) = j
                exit
              end if
            end if
          end do

        end do

        read(ifile,*,err=5) minfo%nat_sub
        allocate(minfo%subatom_id(minfo%nat_sub))
        do i = 1, minfo%nat_sub
          read(ifile,*) atom_name, atom_num, vibatm_mass, vibatm_coord
          vibatm_coord = vibatm_coord * CONV_UNIT_LEN

          do j = 1, molecule%num_atoms
            if (trim(atom_name) == trim(molecule%atom_name(j))) then 
              match = .true.
              do k = 1, 3
                dd = abs((coord(k,j) - vibatm_coord(k))/vibatm_coord(k))
                if (dd > thresh) then
                  match = .false.
                  exit
                end if
              end do

              if(match) then
                minfo%subatom_id(i) = j
                exit
              end if
            end if
          end do
        end do

        5 continue

      end if

      if(index(line,'Electronic Data') > 0) then
        minfo%elec_data = .true.

        read(ifile,'(a)') line
        read(ifile,'(f20.14)') minfo%energy

        do while(.true.)
          read(ifile,'(a)',end=10) line

          if (index(line,'Gradient') > 0) then
            read(ifile,*) isize
            allocate(minfo%gradient(isize))
            read(ifile,*) minfo%gradient(1:isize)
  
          else if (index(line,'Hessian') > 0) then
            read(ifile,*) isize
            allocate(minfo%hessian(isize))
            read(ifile,*) minfo%hessian(1:isize)
  
          else if (index(line,'Dipole Moment') > 0) then
            read(ifile,*) isize
            allocate(minfo%dipole(isize))
            read(ifile,*) minfo%dipole(1:isize)

          else if (index(line,'Dipole Derivative') > 0) then
            read(ifile,*) isize
            allocate(minfo%dipole_derv(3,isize/3))
            read(ifile,*) minfo%dipole_derv(1:3,1:isize/3)

          else
            exit
  
          end if
        end do

      end if

      ! TODO
      !if(index(line,'Vibrational Data') > 0) then
      !  minfo%vib_data = .false.
      !end if

    end do

    10 continue

    call close_file(ifile)

    return

  end subroutine read_minfo

end module fileio_minfo_mod
