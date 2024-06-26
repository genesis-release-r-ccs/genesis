!C
!C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!C
!C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!C                BY
!!C         DAISUKE TAKAHASHI
!C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!C         UNIVERSITY OF TSUKUBA
!C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!C         E-MAIL: daisuke@cs.tsukuba.ac.jp
!C
!C
!C     PARALLEL 3-D REAL-TO-COMPLEX FFT ROUTINE
!C
!C     FORTRAN77 + MPI SOURCE PROGRAM
!C
!C     CALL PDZFFT3D(A,B,NX,NY,NZ,ICOMM,ME,NPU,IOPT)
!C
!C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
!C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
!C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
!C       ------------------------------------
!C         NX = (2**IP) * (3**IQ) * (5**IR)
!C         NY = (2**JP) * (3**JQ) * (5**JR)
!C         NZ = (2**KP) * (3**KQ) * (5**KR)
!C       ------------------------------------
!C     ICOMM IS THE COMMUNICATOR (INTEGER*4)
!C     ME IS THE RANK (INTEGER*4)
!C     NPU IS THE NUMBER OF PROCESSORS (INTEGER*4)
!C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
!C     IOPT = -1 FOR FORWARD TRANSFORM WHERE
!C              A(NX,NY,NZ/NPU) IS REAL INPUT VECTOR (REAL(WP))
!C!HPF$ DISTRIBUTE A(*,*,BLOCK)
!C              A(NX/2+1,NY,NZ/NPU) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(*,*,BLOCK)
!C              B(NX/2+1,NY,NZ/NPU) IS WORK VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(*,*,BLOCK)
!C     IOPT = -2 FOR FORWARD TRANSFORM WHERE
!C              A(NX,NY,NZ/NPU) IS REAL INPUT VECTOR (REAL(WP))
!C!HPF$ DISTRIBUTE A(*,*,BLOCK)
!C     ME = 0   A((NX/2)/NPU+1,NY,NZ) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C     ME > 0   A((NX/2)/NPU,NY,NZ) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(BLOCK,*,*)
!C              B(NX/2+1,NY,NZ/NPU) IS WORK VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(*,*,BLOCK)
!C
!C     WRITTEN BY DAISUKE TAKAHASHI
!C

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

      SUBROUTINE PDZFFT3D(A,B,NX,NY,NZ,ICOMM,ME,NPU,IOPT)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
      INCLUDE 'ffte_param.h'
      DIMENSION A(*)
      COMPLEX(WP) B(*)
      COMPLEX(WP) C((NDA3+NP)*NBLK),D(NDA3)
      COMPLEX(WP) WX(NDA3),WY(NDA3),WZ(NDA3)
      SAVE WX,WY,WZ
!C
      IF (IOPT .EQ. 0) THEN
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL(WZ,NZ)
        RETURN
      END IF
!C
!$omp PARALLEL PRIVATE(C,D)
      CALL PDZFFT3D0(A,A,A,B,B,C,C,C,D,WX,WY,WZ,NX,NY,NZ,ICOMM,ME,NPU, &
                     IOPT)
!$omp END PARALLEL
      RETURN
      END
      SUBROUTINE PDZFFT3D0(DA,A,AXYZ,B,BXYZ,CX,CY,CZ,D,WX,WY,WZ,  &
                           NX,NY,NZ,ICOMM,ME,NPU,IOPT)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
#ifdef HAVE_MPI_GENESIS
      INCLUDE 'mpif.h'
#endif
      INCLUDE 'ffte_param.h'
      COMPLEX(WP) A(*),AXYZ(NX/2+1,NY,*)
      COMPLEX(WP) B(*),BXYZ(NX/2+1,NY,*)
      COMPLEX(WP) CX(*),CY(NY+NP,*),CZ(NZ+NP,*),D(*)
      COMPLEX(WP) WX(*),WY(*),WZ(*)
      DIMENSION DA(NX,NY,*)
      DIMENSION ISCNT(MAXNPU),ISDSP(MAXNPU),IRCNT(MAXNPU),IRDSP(MAXNPU)
      DIMENSION LNX(3),LNY(3),LNZ(3)
!C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)
!C
      NNX=NX/NPU
      NNZ=NZ/NPU
!C
      ISCNT(1)=(NNX/2+1)*NY*NNZ
      ISDSP(1)=0
      DO 10 I=2,NPU
        ISCNT(I)=(NNX/2)*NY*NNZ
        ISDSP(I)=ISDSP(I-1)+ISCNT(I-1)
   10 CONTINUE
      IF (ME .EQ. 0) THEN
        IRCNT(1)=(NNX/2+1)*NY*NNZ
        IRDSP(1)=0
        DO 20 I=2,NPU
          IRCNT(I)=(NNX/2+1)*NY*NNZ
          IRDSP(I)=IRDSP(I-1)+IRCNT(I-1)
   20   CONTINUE
      ELSE
        IRCNT(1)=(NNX/2)*NY*NNZ
        IRDSP(1)=0
        DO 30 I=2,NPU
          IRCNT(I)=(NNX/2)*NY*NNZ
          IRDSP(I)=IRDSP(I-1)+IRCNT(I-1)
   30   CONTINUE
      END IF
!C
!$omp DO
      DO 120 K=1,NNZ
        IF (MOD(NY,2) .EQ. 0) THEN
          DO 60 J=1,NY,2
            DO 40 I=1,NX
              CX(I)=DCMPLX(DA(I,J,K),DA(I,J+1,K))
   40       CONTINUE
            CALL FFT235(CX,D,WX,NX,LNX)
            BXYZ(1,J,K)=REAL(CX(1),WP)
            BXYZ(1,J+1,K)=IMAG(CX(1))
!!DIR$ VECTOR ALIGNED
            DO 50 I=2,NX/2+1
              BXYZ(I,J,K)=0.5_WP*(CX(I)+CONJG(CX(NX-I+2)))
              BXYZ(I,J+1,K)=(0.0_WP,-0.5_WP)*(CX(I)-CONJG(CX(NX-I+2)))
   50       CONTINUE
   60     CONTINUE
        ELSE
          DO 90 J=1,NY-1,2
            DO 70 I=1,NX
              CX(I)=DCMPLX(DA(I,J,K),DA(I,J+1,K))
   70       CONTINUE
            CALL FFT235(CX,D,WX,NX,LNX)
            BXYZ(1,J,K)=REAL(CX(1),WP)
            BXYZ(1,J+1,K)=IMAG(CX(1))
!!DIR$ VECTOR ALIGNED
            DO 80 I=2,NX/2+1
              BXYZ(I,J,K)=0.5_WP*(CX(I)+CONJG(CX(NX-I+2)))
              BXYZ(I,J+1,K)=(0.0_WP,-0.5_WP)*(CX(I)-CONJG(CX(NX-I+2)))
   80       CONTINUE
   90     CONTINUE
          DO 100 I=1,NX
            CX(I)=DCMPLX(DA(I,NY,K),0.0_WP)
  100     CONTINUE
          CALL FFT235(CX,D,WX,NX,LNX)
!!DIR$ VECTOR ALIGNED
          DO 110 I=1,NX/2+1
            BXYZ(I,NY,K)=CX(I)
  110     CONTINUE
        END IF
  120 CONTINUE
!$omp DO
      DO 280 K=1,NNZ
        DO 190 II=1,NNX/2+1,NBLK
          DO 150 JJ=1,NY,NBLK
            DO 140 I=II,MIN0(II+NBLK-1,NNX/2+1)
!!DIR$ VECTOR ALIGNED
              DO 130 J=JJ,MIN0(JJ+NBLK-1,NY)
                CY(J,I-II+1)=BXYZ(I,J,K)
  130         CONTINUE
  140       CONTINUE
  150     CONTINUE
          DO 160 I=II,MIN0(II+NBLK-1,NNX/2+1)
            CALL FFT235(CY(1,I-II+1),D,WY,NY,LNY)
  160     CONTINUE
          DO 180 J=1,NY
!!DIR$ VECTOR ALIGNED
            DO 170 I=II,MIN0(II+NBLK-1,NNX/2+1)
              A(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY)=CY(J,I-II+1)
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
        DO 270 L=2,NPU
          DO 260 II=1,NNX/2,NBLK
            DO 220 JJ=1,NY,NBLK
              DO 210 I=II,MIN0(II+NBLK-1,NNX/2)
!!DIR$ VECTOR ALIGNED
                DO 200 J=JJ,MIN0(JJ+NBLK-1,NY)
                  CY(J,I-II+1)=BXYZ(I+(L-2)*(NNX/2)+(NNX/2+1),J,K)
  200           CONTINUE
  210         CONTINUE
  220       CONTINUE
            DO 230 I=II,MIN0(II+NBLK-1,NNX/2)
              CALL FFT235(CY(1,I-II+1),D,WY,NY,LNY)
  230       CONTINUE
            DO 250 J=1,NY
!!DIR$ VECTOR ALIGNED
              DO 240 I=II,MIN0(II+NBLK-1,NNX/2)
                A(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY                  &
                  +((L-2)*(NNX/2)+(NNX/2+1))*NY*NNZ)=CY(J,I-II+1)
  240         CONTINUE
  250       CONTINUE
  260     CONTINUE
  270   CONTINUE
  280 CONTINUE
!$omp BARRIER
!$omp MASTER
#ifdef HAVE_MPI_GENESIS
      CALL MPI_ALLTOALLV(A,ISCNT,ISDSP,MPI_WP_COMPLEX, &
                         B,IRCNT,IRDSP,MPI_WP_COMPLEX, &
                         ICOMM,IERR)
#else
      DO I=1,ISCNT(1)
         B(I) = A(I)
      ENDDO
#endif 
!$omp END MASTER
!$omp BARRIER
      IF (ME .EQ. 0) THEN
!$omp DO
        DO 360 J=1,NY
          DO 350 II=1,NNX/2+1,NBLK
            DO 310 L=1,NPU
              DO 300 I=II,MIN0(II+NBLK-1,NNX/2+1)
!!DIR$ VECTOR ALIGNED
                DO 290 K=1,NNZ
                  CZ(K+(L-1)*NNZ,I-II+1)                   &
                 =B(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY   &
                    +(L-1)*(NNX/2+1)*NY*NNZ)
  290           CONTINUE
  300         CONTINUE
  310       CONTINUE
            DO 320 I=II,MIN0(II+NBLK-1,NNX/2+1)
              CALL FFT235(CZ(1,I-II+1),D,WZ,NZ,LNZ)
  320       CONTINUE
            DO 340 K=1,NZ
!!DIR$ VECTOR ALIGNED
              DO 330 I=II,MIN0(II+NBLK-1,NNX/2+1)
                A(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY)=CZ(K,I-II+1)
  330         CONTINUE
  340       CONTINUE
  350     CONTINUE
  360   CONTINUE
      ELSE
!$omp DO
        DO 440 J=1,NY
          DO 430 II=1,NNX/2,NBLK
            DO 390 L=1,NPU
              DO 380 I=II,MIN0(II+NBLK-1,NNX/2)
!!DIR$ VECTOR ALIGNED
                DO 370 K=1,NNZ
                  CZ(K+(L-1)*NNZ,I-II+1)              &
                 =B(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY  &
                    +(L-1)*(NNX/2)*NY*NNZ)
  370           CONTINUE
  380         CONTINUE
  390       CONTINUE
            DO 400 I=II,MIN0(II+NBLK-1,NNX/2)
              CALL FFT235(CZ(1,I-II+1),D,WZ,NZ,LNZ)
  400       CONTINUE
            DO 420 K=1,NZ
!!DIR$ VECTOR ALIGNED
              DO 410 I=II,MIN0(II+NBLK-1,NNX/2)
                A(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY)=CZ(K,I-II+1)
  410         CONTINUE
  420       CONTINUE
  430     CONTINUE
  440   CONTINUE
      END IF
      IF (IOPT .EQ. -2) RETURN
!$omp BARRIER
!$omp MASTER
#ifdef HAVE_MPI_GENESIS
      CALL MPI_ALLTOALLV(A,IRCNT,IRDSP,MPI_WP_COMPLEX, &
                         B,ISCNT,ISDSP,MPI_WP_COMPLEX, &
                         ICOMM,IERR)
#else
      DO I=1,ISCNT(1)
         B(I) = A(I)
      ENDDO
#endif
!$omp END MASTER
!$omp BARRIER
!$omp DO
      DO 490 K=1,NNZ
        DO 480 J=1,NY
!!DIR$ VECTOR ALIGNED
          DO 450 I=1,NNX/2+1
            AXYZ(I,J,K)=B(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY)
  450     CONTINUE
          DO 470 L=2,NPU
!!DIR$ VECTOR ALIGNED
            DO 460 I=1,NNX/2
              AXYZ(I+((L-2)*(NNX/2)+(NNX/2+1)),J,K)      &
             =B(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY         &
                +((L-2)*(NNX/2)+(NNX/2+1))*NY*NNZ)
  460       CONTINUE
  470     CONTINUE
  480   CONTINUE
  490 CONTINUE
      RETURN
      END

