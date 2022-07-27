!C
!C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!C
!C     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
!C                BY
!C         DAISUKE TAKAHASHI
!C         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
!C         UNIVERSITY OF TSUKUBA
!C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!C         E-MAIL: daisuke@cs.tsukuba.ac.jp
!C
!C
!C     HEADER FILE FOR PARAMETERS
!C
!C     FORTRAN77 SOURCE PROGRAM
!C
!C     WRITTEN BY DAISUKE TAKAHASHI
!C
!C The maximum supported number of processors is 65536.
      PARAMETER (MAXNPU=65536)
!C The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
!C The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
      PARAMETER (NDA4=256)
!C The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
!C      PARAMETER (NBLK=8)  (for PentiumIII and Athlon)
!C      PARAMETER (NBLK=16) (for Pentium4, Athlon XP, Opteron, Itanium
!C                           and Itanium2)
!C The parameter NP is a padding parameter to avoid cache conflicts in
!C the FFT routines.
      PARAMETER (NP=8)
!C      PARAMETER (NP=2) (for PentiumIII)
!C      PARAMETER (NP=4) (for Athlon, Athlon XP, Opteron and Itanium)
!C      PARAMETER (NP=8) (for Pentium4 and Itanium2)
!C Size of L2 cache
      PARAMETER (L2SIZE=2097152)
