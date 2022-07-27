!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_gpu_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_gpu_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: update_pairlist_pbc_gpu
  public  :: update_pairlist_pbc_check_gpu

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_gpu
  !> @brief        update pairlist in each domain with periodic boundary
  !!               condition for GPU
  !! @authors      JJ, Naruse, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_gpu(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

#ifdef USE_GPU
    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: cutoffdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: atmcls(:,:), atmcls_pbc(:)
    real(wp)                  :: dij1, dij2, dij3
    real(wp)                  :: rtmp1, rtmp2, rtmp3

    ! for GPU
    integer                   :: univ_ij
    integer                   :: num_all, num_target, num_calc
    integer                   :: max_atoms, max_load
    integer,      allocatable :: index(:), count(:)
    integer                   :: id_load, iix, iiy
    integer                   :: ix_natom, iy_natom
    integer                   :: idx, base_idx
    integer                   :: load_ij, load_ji
    integer                   :: all_iter_32x1, all_iter_16x2, all_iter_8x4
    integer                   :: act_iter
    real(8)                   :: rate_iter
    integer(1)                :: mark
    integer                   :: max_natom, univ_mask2_size
    integer(1),       pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer(1),       pointer :: univ_ix_list(:,:), univ_iy_list(:,:)
    integer,          pointer :: univ_ix_natom(:), univ_iy_natom(:)
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)

    integer                   :: ncell_local, ncell_bound
    integer                   :: ret_shape(1), ncell_max
    integer                   :: num
    integer                   :: univ_ix_num(MaxAtom), univ_iy_num(MaxAtom)

    integer                   :: pack_univ_mask2_size
    integer                   :: return_size


    ncell               => domain%num_cell_local
    nboundary           => domain%num_cell_boundary
    natom               => domain%num_atom
    start_atom          => domain%start_atom
    coord               => domain%coord
    trans1              => domain%trans_vec
    trans2              => domain%translated
    cell_move           => domain%cell_move
    system_size         => domain%system_size
    cell_pairlist1      => domain%cell_pairlist1
    atmcls_pbc          => domain%atmcls_pbc
    charge              => domain%charge
    atmcls              => domain%atom_cls_no
    exclusion_mask      => enefunc%exclusion_mask
    exclusion_mask1     => enefunc%exclusion_mask1

    pairdist2           =  pairlist%pairlistdist * pairlist%pairlistdist
    cutoffdist2         =  enefunc%cutoffdist * enefunc%cutoffdist

    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_list        => pairlist%univ_iy_list
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    ncell_local         =  domain%num_cell_local
    ncell_bound         =  domain%num_cell_boundary

    ! allocation of mask
    !
    max_natom = 0
    do i = 1, ncell+nboundary
      if (max_natom < natom(i)) then
        max_natom = natom(i)
      end if
    end do
    start_atom(1:ncell+nboundary) = 0
    ij = natom(1)
    do i = 2, ncell+nboundary
      start_atom(i) = start_atom(i-1) + natom(i-1)
      ij = ij + natom(i)
    end do
    domain%num_atom_domain = ij

    !$omp parallel do default(shared) &
    !$omp private(i, ix, k)
    do i = 1, ncell+nboundary
      k = start_atom(i)
      do ix = 1, natom(i)
        trans2(     k+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(  ij+k+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(2*ij+k+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
        trans2(3*ij+k+ix,1,1) = charge(ix,i) 
        atmcls_pbc(k+ix)      = atmcls(ix,i)
      end do
    end do
    !$omp end parallel do

#ifdef DEBUG
    if (max_natom > MaxAtom) then
      call error_msg('Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
    end if
#endif
    max_load = 16 * ((max_natom+7)/8 * (max_natom+3)/4) + max_natom
    allocate(index(0:max_load), count(0:max_load) )

    univ_mask2_size = max_natom * max_natom
    if (pairlist%univ_mask2_size < univ_mask2_size) then
      if (allocated(pairlist%univ_mask2)) then
        deallocate(pairlist%univ_mask2)
        call unset_pinned_memory(pairlist%pack_univ_mask2)
        deallocate(pairlist%pack_univ_mask2)
      end if
      max_natom = max_natom + max_natom/10  ! to reduce number of memory allocation
      univ_mask2_size = max_natom * max_natom
      allocate(pairlist%univ_mask2(univ_mask2_size,univ_ncell_near))
      pairlist%univ_mask2_size = univ_mask2_size

      pack_univ_mask2_size = (univ_mask2_size * univ_ncell_near + 7) / 8
      allocate(pairlist%pack_univ_mask2(pack_univ_mask2_size))
      call set_pinned_memory(pairlist%pack_univ_mask2, pack_univ_mask2_size)
      pairlist%pack_univ_mask2_size = pack_univ_mask2_size
    end if
    univ_natom_max = max_natom

    !
    ! ******** Initialization of GPU ********
    !
    call gpu_init()

    !
    ! ******** make pairlist on GPU ********
    !
    ret_shape = shape(natom)
    ncell_max = ret_shape(1)

    call gpu_launch_build_pairlist( &
         trans2, atmcls_pbc, cell_move, natom, start_atom,                   &
         univ_cell_pairlist1, univ_ix_list, univ_iy_list, univ_ix_natom,     &
         univ_iy_natom, domain%num_atom_domain, MaxAtom, ncell_local,        &
         ncell_bound, ncell_max, univ_maxcell, univ_maxcell1, pairdist2,     &
         cutoffdist2, system_size(1), system_size(2), system_size(3))

    ! Initialization of mask on CPU
    !
    !$omp parallel default(shared) private(id, univ_ij, ij, i, j, ix, iy, idx)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do univ_ij = id+1, ncell, nthread
      i = univ_cell_pairlist1(1,univ_ij)
      do ix = 1, natom(i)
        do iy = 1, ix
          idx = 1 + (iy-1) + (ix-1)*max_natom
          pairlist%univ_mask2(idx, univ_ij) = exclusion_mask1(ix,iy,i)
          idx = 1 + (ix-1) + (iy-1)*max_natom
          pairlist%univ_mask2(idx, univ_ij) = exclusion_mask1(ix,iy,i)
        end do
      end do
    end do
    do univ_ij =ncell+id+1, univ_ncell_near, nthread
      ij = univ_ij - ncell
      i = univ_cell_pairlist1(1,univ_ij)
      j = univ_cell_pairlist1(2,univ_ij)
      do ix = 1, natom(i)
        do iy = 1, natom(j)
          idx = 1 + (iy-1) + (ix-1)*max_natom
          pairlist%univ_mask2(idx, univ_ij) = exclusion_mask(iy,ix,ij)
        end do
      end do
    end do
    !$omp end parallel

    call gpu_allocate_packdata(pairlist%pack_univ_mask2, &
                               pairlist%univ_mask2_size, univ_ncell_near)
    call send_mask_data(pairlist%univ_mask2, pairlist%univ_mask2_size, &
                        pairlist%pack_univ_mask2_size,          &
                        univ_ncell_near,                        &
                        pairlist%pack_univ_mask2)

    call gpu_wait_build_pairlist(univ_ix_list, univ_iy_list, univ_ix_natom, &
                                 univ_iy_natom, pairlist%univ_mask2,        &
                                 pairlist%univ_mask2_size, univ_ncell_near )

    ! Sort pairlist
    !
    !$omp parallel default(shared) &
    !$omp private(univ_ij, ix_natom, iy_natom, base_idx, idx, id_load)
    !$omp do schedule(static,1)
    do univ_ij =1, univ_maxcell
      ix_natom = univ_ix_natom(univ_ij)
      iy_natom = univ_iy_natom(univ_ij)

      id_load = 16*(((ix_natom+7)/8) * ((iy_natom+3)/4)) + (ix_natom+7)/8;

      if (univ_ij <= ncell) then
        ! self cell-pair
        id_load = max_load - max_natom + (ix_natom + iy_natom)/2
      else if (univ_ij <= univ_ncell_near) then
        ! near
        if (id_load > max_load - max_natom) then
          id_load = max_load - max_natom
        end if
      else
        ! far
        if (id_load > max_load - max_natom - 1) then
          id_load = max_load - max_natom - 1
        end if
      end if
      pairlist%univ_ij_load(univ_ij) = id_load
    end do
    !$omp end parallel

    count(:) = 0
    do univ_ij = 1, univ_maxcell
      id_load = pairlist%univ_ij_load(univ_ij)
      count(id_load) = count(id_load) + 1
    end do

    pairlist%univ_ncell_nonzero = univ_maxcell - count(0)

    index(max_load) = 0
    do id_load = max_load-1, 0, -1
      index(id_load) = index(id_load+1) + count(id_load+1)
    end do

    do univ_ij = 1, univ_maxcell
      id_load = pairlist%univ_ij_load(univ_ij)
      index(id_load) = index(id_load) + 1
      pairlist%univ_ij_sort_list(index(id_load)) = univ_ij
    end do

    deallocate(index, count)

#endif

    return

  end subroutine update_pairlist_pbc_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_check_gpu
  !> @brief        update pairlist in each domain with periodic boundary
  !!               condition for GPU
  !! @authors      JJ, Naruse, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_check_gpu(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain
    type(s_pairlist), target, intent(inout) :: pairlist


    call update_pairlist_pbc_gpu(enefunc, domain, pairlist)

    return

  end subroutine update_pairlist_pbc_check_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    send_mask_data
  !> @brief        Compress and Send mask data
  !! @authors      JJ, Sakamoto
  !! @param[in]    univ_mask2           : mask data
  !! @param[in]    univ_mask2_size      : size of mask data
  !! @param[in]    pack_univ_mask2_size : size of compressed mask2 data
  !! @param[in]    univ_ncell_near      : univ_ncell_near
  !! @param[inout] pack_univ_mask2      : compressed mask2 data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine send_mask_data(univ_mask2, univ_mask2_size, &
                            pack_univ_mask2_size, &
                            univ_ncell_near, pack_univ_mask2)

    integer(1), intent(in)    :: univ_mask2(*)
    integer,    intent(in)    :: univ_mask2_size
    integer,    intent(in)    :: pack_univ_mask2_size
    integer,    intent(in)    :: univ_ncell_near
    integer(1), intent(inout) :: pack_univ_mask2(*)

    integer                   :: all_size
    integer                   :: div_index


    all_size = univ_mask2_size*univ_ncell_near
    div_index = int(all_size/2/8)*8

#ifdef USE_GPU
    call pack_mask_data(univ_mask2,           1, div_index, pack_univ_mask2)
    call gpu_copy_mask2(pack_univ_mask2, univ_mask2_size, univ_ncell_near, &
                        0, div_index/8-1, 1)

    call pack_mask_data(univ_mask2, div_index+1,  all_size, pack_univ_mask2)
    call gpu_copy_mask2(pack_univ_mask2, univ_mask2_size, univ_ncell_near, &
                        div_index/8, pack_univ_mask2_size-1, 2)
#endif

    return

  end subroutine send_mask_data

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_mask_data
  !> @brief        Compress mask2 data
  !! @authors      JJ, Sakamoto
  !! @param[in]    mask2       : mask data
  !! @param[in]    start_index : Start index
  !! @param[in]    end_index   : End index
  !! @param[inout] pack_mask   : Compressed mask2 data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pack_mask_data(mask2, start_index, end_index, pack_mask)

    ! formal arguments
    integer(1), intent(in)    :: mask2(*)
    integer,    intent(in)    :: start_index
    integer,    intent(in)    :: end_index
    integer(1), intent(inout) :: pack_mask(*)

    integer                        :: i,j,k
    integer                        :: masksize
    integer                        :: packblock_size
    integer                        :: mod_size
    integer                        :: packblock_offset

    integer(kind=1)                :: packed
    integer                        :: offset_index
    integer                        :: base_index
    integer(kind=1)                :: setval


    ! Size of data to pack
    masksize = end_index - start_index + 1
    ! Size of data after packing
    packblock_size = masksize/8
    packblock_offset = start_index/8

    ! Loop of packblock
    !$omp parallel do default(shared) private(packed,setval,offset_index)
    do k = 0, packblock_size - 1
      setval = 1
      packed = 0
      do j = 0, 7
        offset_index = k * 8 + j + start_index
        if (mask2(offset_index) == 1) packed = IOR(packed, setval)
        setval = ISHFT(setval, 1)
      end do
      pack_mask(packblock_offset + k + 1) = packed
    end do
    !$omp end parallel do

    ! Check if there is a remainder
    mod_size = mod(masksize, 8)
    if (mod_size > 0) then
      setval = 1
      packed = 0
      do j = 0, mod_size - 1
        offset_index = packblock_size * 8 + j + start_index
        if (mask2(offset_index) == 1) packed = IOR(packed, setval)
        setval = ISHFT(setval, 1)
      end do
      pack_mask(packblock_offset + packblock_size + 1) = packed
    end if

    return

  end subroutine pack_mask_data

end module sp_pairlist_gpu_mod
