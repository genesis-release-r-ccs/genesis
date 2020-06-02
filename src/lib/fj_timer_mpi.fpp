#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#define N_TIM 500
#define IUNIT 97
!cc#define USE_MPI_WTIME
#ifdef __INTEL_COMPILER
#define gettod clockx
#endif
#ifdef INTEL
#define gettod clockx
#endif

module fj_timer_mod

  use mpi

  implicit none

  ! timer data
  real(8), save :: timers(N_TIM), tsta(N_TIM), ttot(3)
  integer, save :: ncount(N_TIM), ncnt_e(N_TIM)
  integer, save :: iregion_id, in_region
  ! 出力済みタイマー情報退避領域
  real(8), save :: timers_old(N_TIM)
  integer, save :: ncount_old(N_TIM), ncnt_e_old(N_TIM)

end module fj_timer_mod


    subroutine timer_init

      use fj_timer_mod
      integer i, id_env
      character*256 env
      real*4 rdum, rarray(2)

      external lnblnk
      integer  lnblnk
      external etime
      real*4   etime

      ! real, user, sys of total
      rdum = etime(rarray)


!#ifdef PKTIMER
      call gettod(ttot(1))
!#endif

      ttot(2) = rarray(1)    ! user
      ttot(3) = rarray(2)    ! sys

      ! clear counters
      do i = 1, N_TIM
        timers(i) = 0.0d0
        ncount(i) = 0
        ncnt_e(i) = 0
      enddo

      ! check timer region setting
      iregion_id = 0  ! default: region setting is off
      in_region  = 1
      call getenv('FJ_TIMER_REGION', env)
      if (lnblnk(env) .gt. 0) then
        read(env,*) id_env
        if (id_env .gt. 0) then
          iregion_id = id_env
          in_region  = 0
        endif
      endif


      return
    end subroutine timer_init


!-----------------------------------------------------------------------
    !  タイマー情報に出力（全体を出力する）
    subroutine timer_fin

      use fj_timer_mod
      integer myrank, nproc, ierr
      character(20) filename

      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

      if (myrank .eq. 0) then
!       filename = 'timer_'//char(nproc)
        write(filename,'(A,i0)') 'timer_',nproc
        if (IUNIT .ne. 6) open(IUNIT, file=filename)
      endif

      call timer_out_file( N_TIM, timers, ncount, ncnt_e, 1 )

      if (myrank .eq. 0) then
        if (IUNIT .ne. 6) close(IUNIT)
      endif

      return
    end subroutine timer_fin

!-----------------------------------------------------------------------
    !  タイマー情報に出力（前回との差分を出力するる）
    subroutine timer_out_range ( istep_no )

      use fj_timer_mod
      integer istep_no
      integer myrank, ierr

      integer iflg_done
      data iflg_done  / 0 /

      double precision timers_range(N_TIM)
      integer ncount_range(N_TIM), ncnt_e_range (N_TIM)

      character*128  file_name
      integer i

      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

      if (myrank .eq. 0) then
        write( file_name,"('timer.out_stepNo',i6.6)" ) istep_no
        if (IUNIT .ne. 6) open(IUNIT, file=file_name)
      endif

      if( iflg_done .eq. 0 )  then
        iflg_done = 1
        timers_old(:) = 0.0
        ncount_old(:) = 0
        ncnt_e_old(:) = 0
      endif

      do i = 1, N_TIM
        timers_range(i) = timers(i) - timers_old(i)
        ncount_range(i) = ncount(i) - ncount_old(i)
        ncnt_e_range(i) = ncnt_e(i) - ncnt_e_old(i)
        timers_old(i) = timers(i)
        ncount_old(i) = ncount(i)
        ncnt_e_old(i) = ncnt_e(i)
      enddo

      call timer_out_file( N_TIM, timers_range, &
                          ncount_range, ncnt_e_range, 0 )

      if (myrank .eq. 0) then
        if (IUNIT .ne. 6) close(IUNIT)
      endif

      ! ファイル出力のオーバヘッドがあるため、同期を取る
      ! timer_sta()が呼ばれてtimer_end()が呼ばれていないものは
      ! オーバヘッド分処理時間がかかるように見えるので注意 
      call mpi_barrier(mpi_comm_world, ierr )

      return
    end subroutine timer_out_range

!-----------------------------------------------------------------------
      !  タイマー情報をファイルに出力
      !  ファイルのopen/closeは外側で行うこと
    subroutine timer_out_file( &
               n_tim_, timers_, ncount_, ncnt_e_, &
               count_chk )
      use fj_timer_mod
      integer n_tim_
      double precision timers_(n_tim_)
      integer ncount_(n_tim_), ncnt_e_(n_tim_)
           ! 個数カウントチェックの有無 =1 check  =0 no check
           ! 実行途中でtimer_out_fileを呼びだすと
           ! ncount, ncnt_eが不一致となるため
      integer count_chk

      integer ign_zero
      parameter (ign_zero = 0)
      double precision g
      character*64 name(n_tim_)

      integer myrank, nprocs, ierr, status(MPI_STATUS_SIZE)
      integer i, m, n_err, n_err_g
      integer pid
      real*4 rdum, rarray(2)

      external etime, getpid
      real*4   etime
      integer getpid

      ! set timer name
      do i = 1, n_tim_
         name(i) = ''
      enddo

      name(  1)     = 'genesis'     ! Total_Calc_Time tion_genesis : Timc(8)
      name(  3)     = 'run_md'      ! Total_Calc_Time run_md : Timc(7)
      name(  4)     = 'mpi_barrier in Main Loop'
      name(  5)     = 'Main Loop'

  !**************************************
  ! 評価区間
  !**************************************
  
  !--- Calculation -------------------------------
      ! name( 11)     = 'Nonb15F compute_force_nonbond_table_charmm()'      ! Total_Calc_Time Nonb15F : Timc(1)
      ! name( 12)     = 'Nonb15F compute_force_nonbond_table_charmm_npt()'  ! Total_Calc_Time Nonb15F : Timc(1)
      name( 11)     = 'Nonb15F'      ! Total_Calc_Time Nonb15F : Timc(1)
      name( 12)     = 'Nonb15F_novirial'  ! Total_Calc_Time Nonb15F : Timc(1)
      name( 13)     = 'PairList'    ! Total_Calc_Time PairList : Timc(5)

      name( 21)     = 'Recip PRE'   ! Total_Calc_Time Recip PRE  : Timc(3)
      name( 22)     = 'Recip POST'  ! Total_Calc_Time Recip POST : Timc(4)
      name( 23)     = 'FFT3D'       ! Total_Calc_Time Recip FFT  : Timc(2)
      name( 24)     = 'Constraint'  ! Total_Calc_Time constraint : Timc(6)

      name( 31)     = 'bussi_v1'
      name( 32)     = 'bussi_v2'

  !--- Communication -------------------------------
      name( 41)     = 'communicate_coor'       ! Total_Trans_Time_coor         : Timt(3)
      name( 42)     = 'communicate_force'      ! Total_Trans_Time_force        : Timt(4)
      name( 43)     = 'FFT3D(Allgather)'       ! Total_Trans_Time_FFT_allgather: Timt(1)
      name( 44)     = 'FFT3D(Alltoall)'        ! Total_Trans_Time_FFT_alltoall : Timt(2)
      name( 45)     = 'tb_bcast'               ! Total_Trans_Time_tb_bcast     : Timt(5)
      name( 46)     = 'tb_allreduce'           ! Total_Trans_Time_tb_allreduce : Timt(6)
      name( 47)     = 'trans_pre'              ! Total_Trans_Time_pre          : Timt(7)
      name( 48)     = 'trans_post'             ! Total_Trans_Time_post         : Timt(8)

  !--- Barrier -------------------------------
      name( 51)     = 'barrier in FFT'         ! Total_barrier_Time in FFT     : Timb(1)
      name( 52)     = 'barrier in coor'        ! Total_barrier_Time in coor    : Timb(2)
      name( 53)     = 'barrier in force'       ! Total_barrier_Time in force   : Timb(3)
      name( 54)     = 'barrier in barosta'     ! Total_barrier_Time in barosta : Timb(4)
      name( 55)     = 'barrier in prepost'     ! Total_barrier_Time in prepost : Timb(5)

!#ifdef FJ_TIMER_DETAIL

  !**************************************
  ! 評価区間内の詳細タイマー
  !**************************************

  !--- Recip PRE -------------------------------
      name(111)     = 'Recip PRE 1'
      name(112)     = 'Recip PRE 2'
      name(113)     = 'Recip PRE 3'
      name(114)     = 'Recip PRE 4'
      name(115)     = 'Recip PRE 5'
      name(116)     = 'Recip PRE 6'

  !--- Recip PRE -------------------------------
      name(121)     = 'Recip POST 1'
      name(122)     = 'Recip POST 2'

  !--- FFT3D -------------------------------
      name(130)     = 'FFT3D pme_recip'
      ! file name: fft3d_1dalltoall.fpp 
      name(131)     = 'FFT3D fft3d_1d_alltoall 1'
      name(132)     = 'FFT3D fft3d_1d_alltoall 2'
      name(133)     = 'FFT3D fft3d_1d_alltoall 3'
      name(134)     = 'FFT3D fft3d_1d_alltoall 4'
      name(135)     = 'FFT3D bfft3d_1d_alltoall 1'
      name(136)     = 'FFT3D bfft3d_1d_alltoall 2'
      name(137)     = 'FFT3D bfft3d_1d_alltoall 3'
      name(138)     = 'FFT3D bfft3d_1d_alltoall 4'
      ! file name: fft3d_2dalltoall.fpp 
      name(139)     = 'FFT3D fft3d_2d_alltoall 1'
      name(140)     = 'FFT3D fft3d_2d_alltoall 2'
      name(141)     = 'FFT3D fft3d_2d_alltoall 3'
      name(142)     = 'FFT3D bfft3d_2d_alltoall 1'
      name(143)     = 'FFT3D bfft3d_2d_alltoall 2'
      name(144)     = 'FFT3D bfft3d_2d_alltoall 3'
      ! file name: fft3d_slab.fpp 
      name(145)     = 'FFT3D fft3d_slab 1'
      name(146)     = 'FFT3D fft3d_slab 2'
      name(147)     = 'FFT3D bfft3d_slab 1'
      name(148)     = 'FFT3D bfft3d_slab 2'
      name(149)     = 'FFT3D bfft3d_slab 3'

  !--- Constraints -------------------------------
      name(171)     = 'Constraint ConstraintModeLEAP  compute_settle()'
      name(172)     = 'Constraint ConstraintModeLEAP  compute_shake()'
      name(173)     = 'Constraint ConstraintModeVVER1 compute_rattle_fast_vv1()'
      name(174)     = 'Constraint ConstraintModeVVER1 compute_rattle_vv1()'
      name(175)     = 'Constraint ConstraintModeVVER2 compute_rattle_fast_vv2()'
      name(176)     = 'Constraint ConstraintModeVVER2 compute_rattle_vv2()'
!#endif

  !**************************************
  ! 評価区間外　高コスト部
  !**************************************

  !--- sp_energy.fpp ------------------
      name(201)     = 'compute_energy_short: clear force etc.'
      name(202)     = 'compute_energy_short: gather values'
      name(204)     = 'compute_energy_allatom clear force,force_long'
      name(205)     = 'compute_energy_allatom clear force_omp,force_pbc'
      name(206)     = 'compute_energy_allatom omp loop3'

  !--- sp_energy_dihedrals.fpp ------------------
      name(210)     = 'compute_energy_dihed'

  !--- sp_energy_angles.fpp ------------------
      name(212)     = 'compute_energy_angle'

  !--- sp_enefunc_charmm.fpp ------------------
      name(214)     = 'count_nonb_excl_omp_loop1'

  !--- sp_energy_table_linear_bondcorr.fpp ------------------
      name(220)     = 'pme_bond_corr_linear_general'

  !--- sp_energy_table_linear_nowater.fpp ------------------
      name(222)     = 'compute_energy_nonbond_table_linear'

  !--- sp_energy_table_linear.fpp ------------------
      name(224)     = 'compute_energy_nonbond14_table_linear_charmm'

!#ifdef FJ_TIMER_2

  !**************************************
  ! MainLoop配下に階層で入れているタイマー
  !**************************************

  !--- run_md sp_dynamics.fppのタイマー -----------
      name(281)     = 'run_md: open_output()'
      name(282)     = 'run_md: leapfrog_dynamics()'
      name(283)     = 'run_md: vverlet_dynamics()'
      name(284)     = 'run_md: vverlet_respa_dynamics()'
      name(285)     = 'run_md: close_output()'

  !--- sp_md_respa.fpp ------------------
      name(291)     = 'vverlet_respa_dynamics: initial_velocity()'
      name(292)     = 'vverlet_respa_dynamics: stop_trans_rotation()'
      name(293)     = 'vverlet_respa_dynamics: initial_vverlet()'

  !   sp_md_respa.fpp MainLoop 配下

      ! name(299)     = 'clear mpi_tran,mpi_bari'    ! very small
      name(300)     = 'set coord_ref,vel_ref'
      name(301)     = 'integrate_vv1'
      name(302)     = 'domain_interaction_update_md'
      name(303)     = 'communicate_coor'
      name(304)     = 'compute_energy_short'
      name(305)     = 'communicate_force'
      name(306)     = 'compute_energy'
      name(307)     = 'communicate_force'
      name(308)     = 'integrate_vv2'
      name(309)     = 'add force short and long'
      name(310)     = 'output_md'
      name(311)     = 'output_prst_md'
      ! name(312)     = 'set timer value'     ! very small

  !--- sp_energy.fpp ------------------
      name(321)     = 'compute_energy_short:init_energy'
      name(322)     = 'compute_energy_short:clear force etc.'
      name(323)     = 'compute_energy_short:compute_energy_nonbond_pme_short'
      name(324)     = 'compute_energy_short:compute_energy_bond'
      name(325)     = 'compute_energy_short:compute_energy_angle'
      name(326)     = 'compute_energy_short:compute_energy_dihed_localres'
      name(327)     = 'compute_energy_short:compute_energy_dihed'
      name(328)     = 'compute_energy_short:compute_energy_improp'
      name(329)     = 'compute_energy_short:compute_energy_cmap'
      name(330)     = 'compute_energy_short:compute_energy_improp_cos'
      name(331)     = 'compute_energy_short:compute_energy_rb_dihed'
      name(332)     = 'compute_energy_short:compute_energy_nonbond14_table_linear'
      name(333)     = 'compute_energy_short:pme_bond_corr_linear'
      name(334)     = 'compute_energy_short:compute_energy_restraints_short'
      name(335)     = 'compute_energy_short:compute_stats'
      name(336)     = 'compute_energy_short:gather values'
      name(337)     = 'compute_energy_short:mpi_allreduce'

      name(351)     = 'compute_energy:compute_energy_allatom'
      name(352)     = 'compute_energy:compute_energy_go'

      name(361)     = 'compute_energy_allatom:init_energy'
      name(362)     = 'compute_energy_allatom:clear force,force_long'
      name(363)     = 'compute_energy_allatom:clear force_omp,force_pbc'
      name(364)     = 'compute_energy_allatom:clear virial_cell etc.'
      name(365)     = 'compute_energy_allatom:compute_energy_nonbond_pme'
      name(366)     = 'compute_energy_allatom:compute_energy_bond'
      name(367)     = 'compute_energy_allatom:compute_energy_angle'
      name(368)     = 'compute_energy_allatom:compute_energy_dihed_localres'
      name(369)     = 'compute_energy_allatom:compute_energy_dihed'
      name(370)     = 'compute_energy_allatom:compute_energy_improp'
      name(371)     = 'compute_energy_allatom:compute_energy_cmap'
      name(372)     = 'compute_energy_allatom:compute_energy_improp_cos'
      name(373)     = 'compute_energy_allatom:compute_energy_rb_dihed'
      name(374)     = 'compute_energy_allatom:compute_energy_nonbond14_table_linear'
      name(375)     = 'compute_energy_allatom:pme_bond_corr_linear'
      name(376)     = 'compute_energy_allatom:set virial_omp'
      name(377)     = 'compute_energy_allatom:compute_energy_restraints'
      name(378)     = 'compute_energy_allatom:compute_stats'
      name(379)     = 'compute_energy_allatom:comm_reciprocal_to_real_force'
      name(380)     = 'compute_energy_allatom:set force'
      name(381)     = 'compute_energy_allatom:mpi_allreduce'
      name(382)     = 'compute_energy_allatom:mpi_bcast'
      name(383)     = 'compute_energy_allatom:comm_real_to_reciprocal_force'
      name(384)     = 'compute_energy_allatom:reduce_ene'

  !--- sp_energy_nonbonds.fpp ------------------
      ! compute_energy_nonbond_pme_short
      name(401)     = 'nonbond_pme_short:compute_force_nonbond_table_linear'

      ! compute_energy_nonbond_pme
      name(403)     = 'nonbond_pme:compute_energy_nonbond_table_linear'
      name(404)     = 'nonbond_pme:compute_force_nonbond_table_linear'
      name(405)     = 'nonbond_pme:pme_pre'
      name(406)     = 'nonbond_pme:mpi_barrier'
      name(407)     = 'nonbond_pme:pme_recip'

  !--- sp_energy_table_linear_nowater.fpp ------------------
      ! compute_force_nonbond_table_linear
      name(411)     = 'nonbond_table_linear:compute_force_nonbond_table_charmm_npt'
      name(412)     = 'nonbond_table_linear:compute_force_nonbond_table_charmm'
!#endif


      ! user and system cpu time
      rdum = etime(rarray)


!#ifdef PKTIMER
      call gettod(g)
!#endif
      ttot(1) = g - ttot(1)
      ttot(1) = ttot(1) * 1.0d-6
      do i = 1, n_tim_
        timers_(i) = timers_(i) * 1.0d-6
      enddo

      ttot(2) = rarray(1)  ! user
      ttot(3) = rarray(2)  ! sys

      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      pid = getpid()

      ! check count of timer call
      if( count_chk .ne. 0 )  then
        n_err = 0
        do i = 1, n_tim_
          if (ncount_(i) .ne. ncnt_e_(i)) then
            !write(*,*) 'Error: timer_fin: not pair timer_(sta|end), ', &
            write(IUNIT,*) 'Error: rank=',myrank,' timer_fin: not pair timer_(sta|end), ', &
                'ID = ', i
            n_err = n_err + 1
          endif
        enddo
        !call MPI_ALLREDUCE(n_err, n_err_g, 1, MPI_INTEGER, MPI_MAX, &
        !                MPI_COMM_WORLD, ierr)
        !if (n_err_g .gt. 0) then
        !  return
        !endif
      endif

      !if (myrank .eq. 0) then
      !  if (IUNIT .ne. 6) open(IUNIT, file='timer.out')
      !endif

      do m = 0, nprocs - 1

        if (myrank .eq. 0) then
          if (m .ne. 0) then
            call MPI_RECV(ttot, 3, MPI_DOUBLE_PRECISION, m, 0, &
                         MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(ncount_, n_tim_, MPI_INTEGER, m, 0, &
                         MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(timers_, n_tim_, MPI_DOUBLE_PRECISION, m, 0, &
                         MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(pid, 1, MPI_INTEGER, m, 0, &
                         MPI_COMM_WORLD, status, ierr)
          endif

          write(IUNIT,'("# rank = ",i6,2x,"pid = ",i8,2x,"real,user,sys = ",3f13.3)') m, pid, ttot
          do i = 1, n_tim_
            ! if (ign_zero .ne. 1 .or. ncount_(i) .ne. 0) then
            if (ign_zero .ne. 1 ) then
              if (name(i) .ne. '') then
                    write(IUNIT,'(" ",i4,1x,a,1x,f13.3,1x,i8)') &
                  i, name(i), timers_(i), ncount_(i)
              endif
            endif
          enddo
          write(IUNIT,'()')

        else if (myrank .eq. m) then   ! myrank .ne. 0
          call MPI_SEND(ttot, 3, MPI_DOUBLE_PRECISION, 0, 0, &
                     MPI_COMM_WORLD, ierr)
          call MPI_SEND(ncount_, n_tim_, MPI_INTEGER, 0, 0, &
                     MPI_COMM_WORLD, ierr)
          call MPI_SEND(timers_, n_tim_, MPI_DOUBLE_PRECISION, 0, 0, &
                     MPI_COMM_WORLD, ierr)
          call MPI_SEND(pid, 1, MPI_INTEGER, 0, 0, &
                     MPI_COMM_WORLD, ierr)
        endif

!-- mpi_barrier for MRQ error on high parallel ---------
        if( m/10000.ge.1 .and. mod(m,10000).eq.0 ) then
!          write(*,*) '***** mpi_barrier start : rank=',myrank,'    m=',m
          call mpi_barrier(mpi_comm_world, ierr )
!          write(*,*) '***** mpi_barrier  end  : rank=',myrank,'    m=',m
        endif

      enddo  ! do m = 0, nprocs - 1

      !if (myrank .eq. 0) then
      !  write(97,'()')
      !  if (97 .ne. 6) close(97)
      !endif

      return
    end subroutine timer_out_file

!-----------------------------------------------------------------------
    subroutine timer_sta(id)
      use fj_timer_mod
      integer id
      double precision g




!!$omp master
      if (id .eq. iregion_id) then
        in_region = 1
      endif

      if (in_region .eq. 1) then


!#ifdef PKTIMER
        call gettod(g)
!#endif
        tsta(id) = g

        ncount(id) = ncount(id) + 1
!!!#ifdef _FPCOLL_PA
!!!        write(cid,'(I8.8)') id
!!!        call start_collection(cid)
!!!#endif
      endif
!!$omp end master

      return
    end subroutine timer_sta

!-----------------------------------------------------------------------
    subroutine timer_end(id)
      use fj_timer_mod
      integer id
      double precision g

!!$omp master
      if (in_region .eq. 1) then

!#ifdef PKTIMER
        call gettod(g)
!#endif
        timers(id) = timers(id) + (g - tsta(id))

        ncnt_e(id) = ncnt_e(id) + 1
      endif

      if (id .eq. iregion_id) then
        in_region = 0
      endif
!!$omp end master

      return
    end subroutine timer_end

!-----------------------------------------------------------------------
    subroutine timer_bar(id)
      use fj_timer_mod
      integer id
      double precision g
      integer ierr

!$omp master
      if (in_region .eq. 1) then

!#ifdef PKTIMER
        call gettod(g)
!#endif
        tsta(id) = g

        ncount(id) = ncount(id) + 1
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!#ifdef PKTIMER
        call gettod(g)
!#endif
        timers(id) = timers(id) + (g - tsta(id))

        ncnt_e(id) = ncnt_e(id) + 1
      endif
!$omp end master

      return
    end subroutine timer_bar

!-----------------------------------------------------------------------
    subroutine get_time( id, sec )
      use fj_timer_mod
      integer id
      double precision sec

      sec = timers(id) * 1.0d-6

      return
    end subroutine get_time

#if 0
!-----------------------------------------------------------------------
    subroutine timer_add(id, add_time)
      integer id
      double precision add_time

!$omp master
      if (in_region .eq. 1) then
        ncount(id) = ncount(id) + 1
        timers(id) = timers(id) + add_time
        ncnt_e(id) = ncnt_e(id) + 1
      endif
!$omp end master

      return
    end subroutine timer_add
#endif

!-----------------------------------------------------------------------
    subroutine timer_set(id, set_time)
      use fj_timer_mod
      integer id
      double precision set_time

!$omp master
      if (in_region .eq. 1) then
        ncount(id) = ncount(id) + 1
        timers(id) = set_time
        ncnt_e(id) = ncnt_e(id) + 1
      endif
!$omp end master

      return
    end subroutine timer_set

