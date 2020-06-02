!C
!C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!C
!C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2009, ALL RIGHTS RESERVED
!C                BY
!C         DAISUKE TAKAHASHI
!C         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
!C         UNIVERSITY OF TSUKUBA
!C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!C         E-MAIL: daisuke@cs.tsukuba.ac.jp
!C
!C
!C     PARALLEL VOLUMETRIC 3-D COMPLEX FFT ROUTINE
!C
!C     FORTRAN77 + MPI SOURCE PROGRAM
!C
!C     CALL PZFFT3DV(A,B,NX,NY,NZ,ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
!C
!C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
!C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
!C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
!C       ------------------------------------
!C         NX = (2**IP) * (3**IQ) * (5**IR)
!C         NY = (2**JP) * (3**JQ) * (5**JR)
!C         NZ = (2**KP) * (3**KQ) * (5**KR)
!C       ------------------------------------
!C     ICOMMY IS THE COMMUNICATOR IN THE Y-DIRECTION (INTEGER*4)
!C     ICOMMZ IS THE COMMUNICATOR IN THE Z-DIRECTION (INTEGER*4)
!C     NPUY IS THE NUMBER OF PROCESSORS IN THE Y-DIRECTION (INTEGER*4)
!C     NPUZ IS THE NUMBER OF PROCESSORS IN THE Z-DIRECTION (INTEGER*4)
!C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
!C     IOPT = -1 FOR FORWARD TRANSFORM WHERE
!C              A(NX,NY/NPUY,NZ/NPUZ) IS COMPLEX INPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(*,BLOCK,BLOCK)
!C              B(NX,NY/NPUY,NZ/NPUZ) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(*,BLOCK,BLOCK)
!C     IOPT = +1 FOR INVERSE TRANSFORM WHERE
!C              A(NX,NY/NPUY,NZ/NPUZ) IS COMPLEX INPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(*,BLOCK,BLOCK)
!C              B(NX,NY/NPUY,NZ/NPUZ) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(*,BLOCK,BLOCK)
!C     IOPT = -2 FOR FORWARD TRANSFORM WHERE
!C              A(NX,NY/NPUY,NZ/NPUZ) IS COMPLEX INPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(*,BLOCK,BLOCK)
!C              B(NX/NPUY,NY/NPUZ,NZ) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(BLOCK,BLOCK,*)
!C     IOPT = +2 FOR INVERSE TRANSFORM WHERE
!C              A(NX/NPUY,NY/NPUZ,NZ) IS COMPLEX INPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(BLOCK,BLOCK,*)
!C              B(NX,NY/NPUY,NZ/NPUZ) IS COMPLEX OUTPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(*,BLOCK,BLOCK)
!C
!C     WRITTEN BY DAISUKE TAKAHASHI
!C

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

      SUBROUTINE PZFFT3DV(A,B,NX,NY,NZ,ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
      INCLUDE 'ffte_param.h'
      COMPLEX(WP) A(*),B(*)
      COMPLEX(WP) C(NDA3)
      COMPLEX(WP) WX(NDA3),WY(NDA3),WZ(NDA3)
      SAVE WX,WY,WZ

      NN=NX*(NY/NPUY)*(NZ/NPUZ)

      IF (IOPT .EQ. 0) THEN
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL(WZ,NZ)
        RETURN
      END IF

      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
!$omp PARALLEL DO
        DO 10 I=1,NN
          A(I)=CONJG(A(I))
   10   CONTINUE
      END IF

      IF (IOPT .EQ. -1 .OR. IOPT .EQ. -2) THEN
!$omp PARALLEL PRIVATE(C)
        CALL PZFFT3DVF(A,A,A,A,A,A,B,B,B,B,B,B,C,WX,WY,WZ,NX,NY,NZ, &
                       ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
!$omp END PARALLEL
      ELSE
!$omp PARALLEL PRIVATE(C)
        CALL PZFFT3DVB(A,A,A,A,A,A,B,B,B,B,B,B,C,WX,WY,WZ,NX,NY,NZ, &
                       ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
!$omp END PARALLEL
      END IF

      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
!JJ     DN=1.0_WP/(REAL(NX,WP)*REAL(NY,WP)*REAL(NZ,WP))
!$omp PARALLEL DO
        DO 20 I=1,NN
!JJ       B(I)=CONJG(B(I))*DN
          B(I)=CONJG(B(I))
   20   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE PZFFT3DVF(A,AXYZP,AXZYP,AYZXP,AZXY,AZXY2,B,BXPYZ,BXZYP, &
                           BYZX,BYZX2,BZXYP,C,WX,WY,WZ,NX,NY,NZ,         &
                           ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
      INCLUDE 'ffte_param.h'
      COMPLEX(WP) A(NX,NY/NPUY,*),AXYZP(NX/NPUY,NY/NPUZ,NZ/NPUZ,*),&
                  AXZYP(NX/NPUY,NZ/NPUZ,NY/NPUY,*), &
                  AYZXP(NY/NPUY,NZ/NPUZ,NX/NPUY,*), &
                  AZXY(NZ,NX/NPUY,*),AZXY2(NZ/NPUZ,NX/NPUY,*)
      COMPLEX(WP) B(NX/NPUY,NY/NPUZ,*),BXPYZ(NX/NPUY,NPUY,NY/NPUY,*),&
                  BXZYP(NX/NPUY,NZ/NPUZ,NY/NPUZ,*),&
                  BYZX(NY,NZ/NPUZ,*),BYZX2(NY/NPUY,NZ/NPUZ,*),&
                  BZXYP(NZ/NPUZ,NX/NPUY,NY/NPUZ,*)
      COMPLEX(WP) C(*)
      COMPLEX(WP) WX(*),WY(*),WZ(*)
      DIMENSION LNX(3),LNY(3),LNZ(3)

      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)

      NN=NX*(NY/NPUY)*(NZ/NPUZ)

!$omp DO
      DO 60 KK=1,NZ/NPUZ,NBLK
        DO 20 K=KK,MIN0(KK+NBLK-1,NZ/NPUZ)
          DO 10 J=1,NY/NPUY
            CALL FFT235(A(1,J,K),C,WX,NX,LNX)
   10     CONTINUE
   20   CONTINUE
        DO 50 I=1,NX
          DO 40 K=KK,MIN0(KK+NBLK-1,NZ/NPUZ)
            DO 30 J=1,NY/NPUY
              BYZX2(J,K,I)=A(I,J,K)
   30       CONTINUE
   40     CONTINUE
   50   CONTINUE
   60 CONTINUE
!$omp END DO
!$omp SINGLE
#ifdef MPI
      CALL MPI_ALLTOALL(B,NN/NPUY,MPI_WP_COMPLEX, &
                        A,NN/NPUY,MPI_WP_COMPLEX, &
                        ICOMMY,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            A(K,J,I) = B(K,J,I)
          ENDDO
        ENDDO
      ENDDO
#endif 
!$omp END SINGLE
!$omp DO
      DO 100 I=1,NX/NPUY
        DO 90 K=1,NZ/NPUZ
          DO 80 L=1,NPUY
            DO 70 J=1,NY/NPUY
              BYZX(J+(L-1)*(NY/NPUY),K,I)=AYZXP(J,K,I,L)
   70       CONTINUE
   80     CONTINUE
          CALL FFT235(BYZX(1,K,I),C,WY,NY,LNY)
   90   CONTINUE
  100 CONTINUE
!$omp DO
      DO 160 JJ=1,NY,NBLK
        DO 150 II=1,NX/NPUY,NBLK
          DO 140 KK=1,NZ/NPUZ,NBLK
            DO 130 J=JJ,MIN0(JJ+NBLK-1,NY)
              DO 120 I=II,MIN0(II+NBLK-1,NX/NPUY)
                DO 110 K=KK,MIN0(KK+NBLK-1,NZ/NPUZ)
                  AZXY2(K,I,J)=BYZX(J,K,I)
  110           CONTINUE
  120         CONTINUE
  130       CONTINUE
  140     CONTINUE
  150   CONTINUE
  160 CONTINUE
!$omp SINGLE
#ifdef MPI
      CALL MPI_ALLTOALL(A,NN/NPUZ,MPI_WP_COMPLEX, &
                        B,NN/NPUZ,MPI_WP_COMPLEX, &
                        ICOMMZ,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            B(K,J,I) = A(K,J,I) 
          ENDDO
        ENDDO
      ENDDO
#endif 
!$omp END SINGLE
!$omp DO
      DO 200 J=1,NY/NPUZ
        DO 190 I=1,NX/NPUY
          DO 180 L=1,NPUZ
            DO 170 K=1,NZ/NPUZ
              AZXY(K+(L-1)*(NZ/NPUZ),I,J)=BZXYP(K,I,J,L)
  170       CONTINUE
  180     CONTINUE
          CALL FFT235(AZXY(1,I,J),C,WZ,NZ,LNZ)
  190   CONTINUE
  200 CONTINUE
!$omp DO
      DO 260 KK=1,NZ,NBLK
        DO 250 JJ=1,NY/NPUZ,NBLK
          DO 240 II=1,NX/NPUY,NBLK
            DO 230 K=KK,MIN0(KK+NBLK-1,NZ)
              DO 220 J=JJ,MIN0(JJ+NBLK-1,NY/NPUZ)
                DO 210 I=II,MIN0(II+NBLK-1,NX/NPUY)
                  B(I,J,K)=AZXY(K,I,J)
  210           CONTINUE
  220         CONTINUE
  230       CONTINUE
  240     CONTINUE
  250   CONTINUE
  260 CONTINUE
      IF (IOPT .EQ. -2) RETURN
!$omp SINGLE
#ifdef MPI
      CALL MPI_ALLTOALL(B,NN/NPUZ,MPI_WP_COMPLEX, &
                        A,NN/NPUZ,MPI_WP_COMPLEX, &
                        ICOMMZ,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            A(K,J,I) = B(K,J,I)
          ENDDO
        ENDDO
      ENDDO
#endif 
!$omp END SINGLE
!$omp DO
      DO 300 L=1,NPUZ
        DO 290 J=1,NY/NPUZ
          DO 280 K=1,NZ/NPUZ
            DO 270 I=1,NX/NPUY
              BXZYP(I,K,J,L)=AXYZP(I,J,K,L)
  270       CONTINUE
  280     CONTINUE
  290   CONTINUE
  300 CONTINUE
!$omp SINGLE
#ifdef MPI
      CALL MPI_ALLTOALL(B,NN/NPUY,MPI_WP_COMPLEX, &
                        A,NN/NPUY,MPI_WP_COMPLEX, &
                        ICOMMY,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            A(K,J,I) = B(K,J,I)
          ENDDO
        ENDDO
      ENDDO
#endif 
!$omp END SINGLE
!$omp DO
      DO 340 K=1,NZ/NPUZ
        DO 330 J=1,NY/NPUY
          DO 320 L=1,NPUY
            DO 310 I=1,NX/NPUY
              BXPYZ(I,L,J,K)=AXZYP(I,K,J,L)
  310       CONTINUE
  320     CONTINUE
  330   CONTINUE
  340 CONTINUE
      RETURN
      END
      SUBROUTINE PZFFT3DVB(A,AXPYZ,AXZYP,AYZX,AYZX2,AZXYP,B,BXYZP,BXZYP, &
                           BYZXP,BZXY,BZXY2,C,WX,WY,WZ,NX,NY,NZ, &
                           ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
      INCLUDE 'ffte_param.h'
      COMPLEX(WP) A(NX/NPUY,NY/NPUZ,*),AXPYZ(NX/NPUY,NPUY,NY/NPUY,*), &
                  AXZYP(NX/NPUY,NZ/NPUZ,NY/NPUZ,*),AYZX(NY,NZ/NPUZ,*),&
                  AYZX2(NY/NPUY,NZ/NPUZ,*), &
                  AZXYP(NZ/NPUZ,NX/NPUY,NY/NPUZ,*)
      COMPLEX(WP) B(NX,NY/NPUY,*),BXYZP(NX/NPUY,NY/NPUZ,NZ/NPUZ,*), &
                  BXZYP(NX/NPUY,NZ/NPUZ,NY/NPUY,*), &
                  BYZXP(NY/NPUY,NZ/NPUZ,NX/NPUY,*), &
                  BZXY(NZ,NX/NPUY,*),BZXY2(NZ/NPUZ,NX/NPUY,*)
      COMPLEX(WP) C(*)
      COMPLEX(WP) WX(*),WY(*),WZ(*)
      DIMENSION LNX(3),LNY(3),LNZ(3)

      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)

      NN=(NX/NPUY)*(NY/NPUZ)*NZ

      IF (IOPT .EQ. 1) THEN
!$omp DO
        DO 40 L=1,NPUY
          DO 30 J=1,NY/NPUY
            DO 20 K=1,NZ/NPUZ
              DO 10 I=1,NX/NPUY
                BXZYP(I,K,J,L)=AXPYZ(I,L,J,K)
   10         CONTINUE
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
!$omp SINGLE
#ifdef MPI
        CALL MPI_ALLTOALL(B,NN/NPUY,MPI_WP_COMPLEX, &
                          A,NN/NPUY,MPI_WP_COMPLEX, &
                          ICOMMY,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            A(K,J,I) = B(K,J,I)
          ENDDO
        ENDDO
      ENDDO
#endif 
!$omp END SINGLE
!$omp DO
        DO 80 L=1,NPUZ
          DO 70 K=1,NZ/NPUZ
            DO 60 J=1,NY/NPUZ
              DO 50 I=1,NX/NPUY
                BXYZP(I,J,K,L)=AXZYP(I,K,J,L)
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE
   80   CONTINUE
!$omp SINGLE
#ifdef MPI
        CALL MPI_ALLTOALL(B,NN/NPUZ,MPI_WP_COMPLEX, &
                          A,NN/NPUZ,MPI_WP_COMPLEX, &
                          ICOMMZ,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            A(K,J,I) = B(K,J,I)
          ENDDO
        ENDDO
      ENDDO
#endif 

!$omp END SINGLE
      END IF
!$omp DO
      DO 140 JJ=1,NY/NPUZ,NBLK
        DO 130 II=1,NX/NPUY,NBLK
          DO 120 KK=1,NZ,NBLK
            DO 110 J=JJ,MIN0(JJ+NBLK-1,NY/NPUZ)
              DO 100 I=II,MIN0(II+NBLK-1,NX/NPUY)
                DO 90 K=KK,MIN0(KK+NBLK-1,NZ)
                  BZXY(K,I,J)=A(I,J,K)
   90           CONTINUE
  100         CONTINUE
  110       CONTINUE
  120     CONTINUE
  130   CONTINUE
  140 CONTINUE
!$omp DO
      DO 180 J=1,NY/NPUZ
        DO 170 I=1,NX/NPUY
          CALL FFT235(BZXY(1,I,J),C,WZ,NZ,LNZ)
          DO 160 L=1,NPUZ
            DO 150 K=1,NZ/NPUZ
              AZXYP(K,I,J,L)=BZXY(K+(L-1)*(NZ/NPUZ),I,J)
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
  180 CONTINUE
!$omp SINGLE
#ifdef MPI
      CALL MPI_ALLTOALL(A,NN/NPUZ,MPI_WP_COMPLEX, &
                        B,NN/NPUZ,MPI_WP_COMPLEX, &
                        ICOMMZ,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            B(K,J,I) = A(K,J,I)
          ENDDO
        ENDDO
      ENDDO
#endif 
!$omp END SINGLE
!$omp DO
      DO 240 II=1,NX/NPUY,NBLK
        DO 230 KK=1,NZ/NPUZ,NBLK
          DO 220 JJ=1,NY,NBLK
            DO 210 I=II,MIN0(II+NBLK-1,NX/NPUY)
              DO 200 K=KK,MIN0(KK+NBLK-1,NZ/NPUZ)
                DO 190 J=JJ,MIN0(JJ+NBLK-1,NY)
                  AYZX(J,K,I)=BZXY2(K,I,J)
  190           CONTINUE
  200         CONTINUE
  210       CONTINUE
  220     CONTINUE
  230   CONTINUE
  240 CONTINUE
!$omp DO
      DO 280 I=1,NX/NPUY
        DO 270 K=1,NZ/NPUZ
          CALL FFT235(AYZX(1,K,I),C,WY,NY,LNY)
          DO 260 L=1,NPUY
            DO 250 J=1,NY/NPUY
              BYZXP(J,K,I,L)=AYZX(J+(L-1)*(NY/NPUY),K,I)
  250       CONTINUE
  260     CONTINUE
  270   CONTINUE
  280 CONTINUE
!$omp SINGLE
#ifdef MPI
      CALL MPI_ALLTOALL(B,NN/NPUY,MPI_WP_COMPLEX, &
                        A,NN/NPUY,MPI_WP_COMPLEX, &
                        ICOMMY,IERR)
#else
      DO I = 1,NZ
        DO J = 1,NY
          DO K = 1,NX
            A(K,J,I) = B(K,J,I)
          ENDDO
        ENDDO
      ENDDO
#endif 
!$omp END SINGLE
!$omp DO
      DO 360 KK=1,NZ/NPUZ,NBLK
        DO 350 JJ=1,NY/NPUY,NBLK
          DO 320 II=1,NX,NBLK
            DO 310 K=KK,MIN0(KK+NBLK-1,NZ/NPUZ)
              DO 300 J=JJ,MIN0(JJ+NBLK-1,NY/NPUY)
                DO 290 I=II,MIN0(II+NBLK-1,NX)
                  B(I,J,K)=AYZX2(J,K,I)
  290           CONTINUE
  300         CONTINUE
  310       CONTINUE
  320     CONTINUE
          DO 340 K=KK,MIN0(KK+NBLK-1,NZ/NPUZ)
            DO 330 J=JJ,MIN0(JJ+NBLK-1,NY/NPUY)
              CALL FFT235(B(1,J,K),C,WX,NX,LNX)
  330       CONTINUE
  340     CONTINUE
  350   CONTINUE
  360 CONTINUE
      RETURN
      END
