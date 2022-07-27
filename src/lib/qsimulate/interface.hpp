#ifdef QSIMULATE
  interface
    function qsimulate_interface(mpicomm, myrank, input,             &
                                 natoms, atoms, coord, charges,      &
                                 force, qmcharges,                   &
                                 born_radii, error) bind(C)
      use constants_mod
      use iso_c_binding

      type(C_ptr) :: qsimulate_interface
      integer :: mpicomm, myrank
      character(kind=c_char) :: input(*)
      integer(c_int), value :: natoms
      type(c_ptr), value :: atoms
      type(c_ptr), value :: coord
      type(c_ptr), value :: charges
      type(c_ptr), value :: force
      type(c_ptr), value :: qmcharges
      type(c_ptr), value :: born_radii
      logical  :: error
    end function qsimulate_interface
  end interface
#endif
