!C
!C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!C
!C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!C                BY
!C         DAISUKE TAKAHASHI
!C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!C         UNIVERSITY OF TSUKUBA
!C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!C         E-MAIL: daisuke@cs.tsukuba.ac.jp
!C
!C
!C     3-D REAL-TO-COMPLEX FFT ROUTINE
!C
!C     FORTRAN77 SOURCE PROGRAM
!C
!C     CALL DZFFT3D(A,NX,NY,NZ,IOPT,B)
!C
!C     A(NX,NY,NZ) IS REAL INPUT VECTOR (REAL(WP))
!C     A(NX/2+1,NY,NZ) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C     B(NX/2+1,NY,NZ) IS WORK VECTOR (COMPLEX(WP))
!C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
!C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
!C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
!C       ------------------------------------
!C         NX = (2**IP) * (3**IQ) * (5**IR)
!C         NY = (2**JP) * (3**JQ) * (5**JR)
!C         NZ = (2**KP) * (3**KQ) * (5**KR)
!C       ------------------------------------
!C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
!C          = -1 FOR FORWARD TRANSFORM
!C
!C     WRITTEN BY DAISUKE TAKAHASHI
!C

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#if !defined(_SINGLE) && !defined(_MIXED)

      SUBROUTINE DZFFT3D(A,NX,NY,NZ,IOPT,B)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
      INCLUDE 'ffte_param.h'
      COMPLEX(WP) A(*),B(*)
      COMPLEX(WP) C((NDA3+NP)*NBLK),D(NDA3)
      COMPLEX(WP) WX(NDA3),WY(NDA3),WZ(NDA3)
      DIMENSION LNX(3),LNY(3),LNZ(3)
      SAVE WX,WY,WZ
!C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)
!C
      IF (IOPT .EQ. 0) THEN
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL(WZ,NZ)
        RETURN
      END IF
!C
!$OMP PARALLEL PRIVATE(C,D)
      CALL DZFFT3D0(A,A,B,C,C,C,D,WX,WY,WZ,NX,NY,NZ,LNX,LNY,LNZ)
!$OMP END PARALLEL
      RETURN
      END
      SUBROUTINE DZFFT3D0(DA,A,B,CX,CY,CZ,D,WX,WY,WZ,NX,NY,NZ, &
                         LNX,LNY,LNZ)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
      INCLUDE 'ffte_param.h'
      COMPLEX(WP) A(NX/2+1,NY,*),B(NX/2+1,NY,*)
      COMPLEX(WP) CX(*),CY(NY+NP,*),CZ(NZ+NP,*),D(*)
      COMPLEX(WP) WX(*),WY(*),WZ(*)
      DIMENSION DA(NX,NY,*)
      DIMENSION LNX(*),LNY(*),LNZ(*)
!C
!$OMP DO
      DO 150 K=1,NZ
        IF (MOD(NY,2) .EQ. 0) THEN
          DO 30 J=1,NY,2
            DO 10 I=1,NX
              CX(I)=DCMPLX(DA(I,J,K),DA(I,J+1,K))
   10       CONTINUE
            CALL FFT235(CX,D,WX,NX,LNX)
            B(1,J,K)=DBLE(CX(1))
            B(1,J+1,K)=DIMAG(CX(1))
!!DIR$ VECTOR ALIGNED
            DO 20 I=2,NX/2+1
              B(I,J,K)=0.5D0*(CX(I)+DCONJG(CX(NX-I+2)))
              B(I,J+1,K)=(0.0D0,-0.5D0)*(CX(I)-DCONJG(CX(NX-I+2)))
   20       CONTINUE
   30     CONTINUE
        ELSE
          DO 60 J=1,NY-1,2
            DO 40 I=1,NX
              CX(I)=DCMPLX(DA(I,J,K),DA(I,J+1,K))
   40       CONTINUE
            CALL FFT235(CX,D,WX,NX,LNX)
            B(1,J,K)=DBLE(CX(1))
            B(1,J+1,K)=DIMAG(CX(1))
!!DIR$ VECTOR ALIGNED
            DO 50 I=2,NX/2+1
              B(I,J,K)=0.5D0*(CX(I)+DCONJG(CX(NX-I+2)))
              B(I,J+1,K)=(0.0D0,-0.5D0)*(CX(I)-DCONJG(CX(NX-I+2)))
   50       CONTINUE
   60     CONTINUE
          DO 70 I=1,NX
            CX(I)=DCMPLX(DA(I,NY,K),0.0D0)
   70     CONTINUE
          CALL FFT235(CX,D,WX,NX,LNX)
!!DIR$ VECTOR ALIGNED
          DO 80 I=1,NX/2+1
            B(I,NY,K)=CX(I)
   80     CONTINUE
        END IF
        DO 140 II=1,NX/2+1,NBLK
          DO 100 I=II,MIN0(II+NBLK-1,NX/2+1)
!!DIR$ VECTOR ALIGNED
            DO 90 J=1,NY
              CY(J,I-II+1)=B(I,J,K)
   90       CONTINUE
  100     CONTINUE
          DO 110 I=II,MIN0(II+NBLK-1,NX/2+1)
            CALL FFT235(CY(1,I-II+1),D,WY,NY,LNY)
  110     CONTINUE
          DO 130 J=1,NY
!!DIR$ VECTOR ALIGNED
            DO 120 I=II,MIN0(II+NBLK-1,NX/2+1)
              B(I,J,K)=CY(J,I-II+1)
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
!$OMP DO
      DO 220 J=1,NY
        DO 210 II=1,NX/2+1,NBLK
          DO 170 I=II,MIN0(II+NBLK-1,NX/2+1)
!!DIR$ VECTOR ALIGNED
            DO 160 K=1,NZ
              CZ(K,I-II+1)=B(I,J,K)
  160       CONTINUE
  170     CONTINUE
          DO 180 I=II,MIN0(II+NBLK-1,NX/2+1)
            CALL FFT235(CZ(1,I-II+1),D,WZ,NZ,LNZ)
  180     CONTINUE
          DO 200 K=1,NZ
!!DIR$ VECTOR ALIGNED
            DO 190 I=II,MIN0(II+NBLK-1,NX/2+1)
              A(I,J,K)=CZ(K,I-II+1)
  190       CONTINUE
  200     CONTINUE
  210   CONTINUE
  220 CONTINUE
      RETURN
      END
#endif
