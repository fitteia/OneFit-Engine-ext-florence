      PROGRAM TEST_CHEEVHQR
C     Tests of the production eigensolver CHEEVHQR (Householder
C     tridiagonalization CHTRIDR + implicit QL TRIDQLR in
C     ../Florence_f.f, which this program is linked with): a diagonal
C     and a non-diagonal known-answer case, then a sweep of random
C     Hermitian matrices against LAPACK's ZHEEV (eigenvalues, and the
C     reconstruction H*v = w*v of every eigenvector).
      IMPLICIT NONE
      LOGICAL ALLPASS
      ALLPASS = .TRUE.
      CALL CASE_DIAGONAL(ALLPASS)
      CALL CASE_NONDIAGONAL(ALLPASS)
      CALL SWEEP_VS_ZHEEV(ALLPASS)
      IF (ALLPASS) THEN
         WRITE(*,*) 'PASS: CHEEVHQR known-answer + sweep tests'
      ELSE
         WRITE(*,*) 'FAIL: CHEEVHQR known-answer + sweep tests'
         STOP 1
      ENDIF
      END

      SUBROUTINE CASE_DIAGONAL(ALLPASS)
      IMPLICIT NONE
      LOGICAL ALLPASS
      INTEGER N, I, J, INFO
      PARAMETER (N=3)
      DOUBLE PRECISION AR(N,N), AI(N,N), W(N), EXPECT(N), ERR
      DO 20 J = 1, N
         DO 10 I = 1, N
            AR(I,J) = 0.0D0
            AI(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
      AR(1,1) = 2.0D0
      AR(2,2) = 3.0D0
      AR(3,3) = 5.0D0
      EXPECT(1) = 2.0D0
      EXPECT(2) = 3.0D0
      EXPECT(3) = 5.0D0
      CALL CHEEVHQR(N, AR, AI, N, W, INFO)
      ERR = 0.0D0
      DO 30 I = 1, N
         ERR = MAX(ERR, DABS(W(I)-EXPECT(I)))
   30 CONTINUE
      WRITE(*,*) 'CASE_DIAGONAL info=', INFO, ' err=', ERR
      IF (INFO.NE.0 .OR. ERR.GT.1.0D-12) THEN
         WRITE(*,*) '  CASE_DIAGONAL FAIL'
         ALLPASS = .FALSE.
      ELSE
         WRITE(*,*) '  CASE_DIAGONAL PASS'
      ENDIF
      RETURN
      END

      SUBROUTINE CASE_NONDIAGONAL(ALLPASS)
      IMPLICIT NONE
      LOGICAL ALLPASS
      INTEGER N, INFO
      PARAMETER (N=2)
      DOUBLE PRECISION AR(N,N), AI(N,N), W(N)
      DOUBLE PRECISION SQ2, ERR
      SQ2 = DSQRT(2.0D0)
      AR(1,1) = 0.0D0
      AI(1,1) = 0.0D0
      AR(2,2) = 0.0D0
      AI(2,2) = 0.0D0
      AR(1,2) = 1.0D0
      AI(1,2) = -1.0D0
      AR(2,1) = 1.0D0
      AI(2,1) = 1.0D0
      CALL CHEEVHQR(N, AR, AI, N, W, INFO)
      ERR = MAX(DABS(W(1)-(-SQ2)), DABS(W(2)-SQ2))
      WRITE(*,*) 'CASE_NONDIAGONAL info=', INFO, ' err=', ERR
      IF (INFO.NE.0 .OR. ERR.GT.1.0D-12) THEN
         WRITE(*,*) '  CASE_NONDIAGONAL FAIL'
         ALLPASS = .FALSE.
      ELSE
         WRITE(*,*) '  CASE_NONDIAGONAL PASS'
      ENDIF
      RETURN
      END

      SUBROUTINE SWEEP_VS_ZHEEV(ALLPASS)
      IMPLICIT NONE
      LOGICAL ALLPASS
      INTEGER NMAX, NTRIAL, N, ITRIAL, I, J, K, INFO, INFOZ
      PARAMETER (NMAX=30, NTRIAL=30)
      DOUBLE PRECISION AR(NMAX,NMAX), AI(NMAX,NMAX)
      DOUBLE PRECISION ARSAVE(NMAX,NMAX), AISAVE(NMAX,NMAX)
      DOUBLE PRECISION ARREC(NMAX,NMAX), AIREC(NMAX,NMAX)
      COMPLEX*16 ZA(NMAX,NMAX), WORK(4*NMAX)
      DOUBLE PRECISION RWORK(3*NMAX-2)
      DOUBLE PRECISION W(NMAX), WZ(NMAX)
      DOUBLE PRECISION EIGERR, RECERR, WORSTEIG, WORSTREC
      DOUBLE PRECISION X, Y

      WORSTEIG = 0.0D0
      WORSTREC = 0.0D0

      DO 200 N = 2, NMAX
         DO 190 ITRIAL = 1, NTRIAL
            DO 110 J = 1, N
               CALL RANDOM_NUMBER(X)
               ARSAVE(J,J) = 2.0D0*X-1.0D0
               AISAVE(J,J) = 0.0D0
               DO 100 I = J+1, N
                  CALL RANDOM_NUMBER(X)
                  CALL RANDOM_NUMBER(Y)
                  ARSAVE(I,J) = 2.0D0*X-1.0D0
                  AISAVE(I,J) = 2.0D0*Y-1.0D0
                  ARSAVE(J,I) = ARSAVE(I,J)
                  AISAVE(J,I) = -AISAVE(I,J)
  100          CONTINUE
  110       CONTINUE

            DO 130 J = 1, N
               DO 120 I = 1, N
                  AR(I,J) = ARSAVE(I,J)
                  AI(I,J) = AISAVE(I,J)
                  ZA(I,J) = DCMPLX(ARSAVE(I,J),AISAVE(I,J))
  120          CONTINUE
  130       CONTINUE

            CALL CHEEVHQR(N, AR, AI, NMAX, W, INFO)
            CALL ZHEEV('V','L',N,ZA,NMAX,WZ,WORK,4*NMAX,RWORK,INFOZ)

            IF (INFO.NE.0 .OR. INFOZ.NE.0) THEN
               WRITE(*,*) 'SWEEP: solver failure N=',N,
     *                    ' trial=',ITRIAL,' info=',INFO,
     *                    ' infoz=',INFOZ
               ALLPASS = .FALSE.
               GOTO 190
            ENDIF

            EIGERR = 0.0D0
            DO 140 I = 1, N
               EIGERR = MAX(EIGERR, DABS(W(I)-WZ(I)))
  140       CONTINUE
            IF (EIGERR.GT.WORSTEIG) WORSTEIG = EIGERR

C           Reconstruction: AR/AI now hold eigenvectors (real/imag).
C           AREC = V * diag(W) * V^H should equal the original matrix.
            DO 165 J = 1, N
               DO 160 I = 1, N
                  ARREC(I,J) = 0.0D0
                  AIREC(I,J) = 0.0D0
                  DO 150 K = 1, N
C                    term = V(I,K)*W(K)*conj(V(J,K))
C                    V(I,K)=(AR(I,K),AI(I,K)), conj(V(J,K)) =
C                    (AR(J,K),-AI(J,K))
                     ARREC(I,J) = ARREC(I,J) + W(K)*
     *                    (AR(I,K)*AR(J,K) + AI(I,K)*AI(J,K))
                     AIREC(I,J) = AIREC(I,J) + W(K)*
     *                    (AI(I,K)*AR(J,K) - AR(I,K)*AI(J,K))
  150             CONTINUE
  160          CONTINUE
  165       CONTINUE

            RECERR = 0.0D0
            DO 180 J = 1, N
               DO 170 I = 1, N
                  RECERR = MAX(RECERR, DABS(ARREC(I,J)-ARSAVE(I,J)))
                  RECERR = MAX(RECERR, DABS(AIREC(I,J)-AISAVE(I,J)))
  170          CONTINUE
  180       CONTINUE
            IF (RECERR.GT.WORSTREC) WORSTREC = RECERR

            IF (EIGERR.GT.1.0D-9 .OR. RECERR.GT.1.0D-9) THEN
               WRITE(*,*) 'SWEEP FAIL N=',N,' trial=',ITRIAL,
     *                    ' eigerr=',EIGERR,' recerr=',RECERR
               ALLPASS = .FALSE.
            ENDIF
  190    CONTINUE
  200 CONTINUE

      WRITE(*,*) 'SWEEP worst eigerr=', WORSTEIG,
     *           ' worst recerr=', WORSTREC
      IF (WORSTEIG.LE.1.0D-9 .AND. WORSTREC.LE.1.0D-9) THEN
         WRITE(*,*) '  SWEEP PASS'
      ELSE
         WRITE(*,*) '  SWEEP FAIL'
      ENDIF
      RETURN
      END
