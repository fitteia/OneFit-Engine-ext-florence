C     Known-answer tests for OFHEIG, Florence's complex Hermitian
C     eigensolver (../Florence_f.f): eigenvalues, the phase convention
C     of the eigenvectors, and IFAIL on return (0 after success with
C     any value on entry; the failure code, without stopping, when
C     IFAIL is not 0 on entry), the same results for a matrix scaled
C     from 1D-300 to 1D300, and code 3 for a NaN or infinite element.
      PROGRAM TESTEIGENSOLVER
      IMPLICIT NONE
      INTEGER N, MBRANC
      PARAMETER (MBRANC=3)
      DOUBLE PRECISION AR(MBRANC,MBRANC), AI(MBRANC,MBRANC)
      DOUBLE PRECISION WR(MBRANC), VR(MBRANC,MBRANC), VI(MBRANC,MBRANC)
      DOUBLE PRECISION WK1(MBRANC), WK2(MBRANC), WK3(MBRANC)
      INTEGER IFAIL, I, J, K
      DOUBLE PRECISION BR(MBRANC,MBRANC), BI(MBRANC,MBRANC)
      DOUBLE PRECISION W0(MBRANC)
      DOUBLE PRECISION V0R(MBRANC,MBRANC), V0I(MBRANC,MBRANC)
      DOUBLE PRECISION S(6), ERR, X
      DATA S /1.0D-300, 1.0D-200, 1.0D-100, 1.0D100, 1.0D200, 1.0D300/
      DOUBLE PRECISION TOL
      PARAMETER (TOL=1.0D-12)
      LOGICAL OK
      OK = .TRUE.

C     Case 1: diagonal matrix, trivially known eigenvalues 2,3,5.
      N = 3
      AR(1,1)=2.0D0
      AR(1,2)=0.0D0
      AR(1,3)=0.0D0
      AR(2,1)=0.0D0
      AR(2,2)=3.0D0
      AR(2,3)=0.0D0
      AR(3,1)=0.0D0
      AR(3,2)=0.0D0
      AR(3,3)=5.0D0
      AI(1,1)=0.0D0
      AI(1,2)=0.0D0
      AI(1,3)=0.0D0
      AI(2,1)=0.0D0
      AI(2,2)=0.0D0
      AI(2,3)=0.0D0
      AI(3,1)=0.0D0
      AI(3,2)=0.0D0
      AI(3,3)=0.0D0
      IFAIL = 0
      CALL OFHEIG(AR,MBRANC,AI,MBRANC,N,WR,VR,MBRANC,VI,MBRANC,
     $            WK1,WK2,WK3,IFAIL)
      IF (IFAIL.NE.0) THEN
         PRINT *, 'CASE 1: unexpected IFAIL=', IFAIL
         OK = .FALSE.
      END IF
      IF (ABS(WR(1)-2.0D0).GT.TOL .OR. ABS(WR(2)-3.0D0).GT.TOL .OR.
     $    ABS(WR(3)-5.0D0).GT.TOL) THEN
         PRINT *, 'CASE 1 FAIL: eigenvalues', WR(1), WR(2), WR(3)
         OK = .FALSE.
      ELSE
         PRINT *, 'CASE 1 PASS'
      END IF

C     Case 2: [[0,1-i],[1+i,0]], known eigenvalues +-sqrt(2), and
C     known (up to the routine's own phase normalization) eigenvectors
C     - see the commit that added this test for the by-hand derivation.
      N = 2
      AR(1,1)=0.0D0
      AI(1,1)=0.0D0
      AR(2,1)=1.0D0
      AI(2,1)=1.0D0
      AR(2,2)=0.0D0
      AI(2,2)=0.0D0
      IFAIL = 0
      CALL OFHEIG(AR,MBRANC,AI,MBRANC,N,WR,VR,MBRANC,VI,MBRANC,
     $            WK1,WK2,WK3,IFAIL)
      IF (IFAIL.NE.0) THEN
         PRINT *, 'CASE 2: unexpected IFAIL=', IFAIL
         OK = .FALSE.
      END IF
      IF (ABS(WR(1)-(-SQRT(2.0D0))).GT.TOL .OR.
     $    ABS(WR(2)-SQRT(2.0D0)).GT.TOL) THEN
         PRINT *, 'CASE 2 FAIL: eigenvalues', WR(1), WR(2)
         OK = .FALSE.
      ELSE IF (ABS(VR(1,2)-0.70710678118654752D0).GT.TOL .OR.
     $         ABS(VR(2,2)-0.5D0).GT.TOL .OR.
     $         ABS(VI(2,2)-0.5D0).GT.TOL) THEN
         PRINT *, 'CASE 2 FAIL: eigenvector 2', VR(1,2), VR(2,2),
     $            VI(2,2)
         OK = .FALSE.
      ELSE
         PRINT *, 'CASE 2 PASS'
      END IF

C     Case 3: case 2 again with IFAIL = 1 on entry: success must
C     still return IFAIL = 0.
      IFAIL = 1
      CALL OFHEIG(AR,MBRANC,AI,MBRANC,N,WR,VR,MBRANC,VI,MBRANC,
     $            WK1,WK2,WK3,IFAIL)
      IF (IFAIL.NE.0 .OR. ABS(WR(2)-SQRT(2.0D0)).GT.TOL) THEN
         PRINT *, 'CASE 3 FAIL: IFAIL=', IFAIL, ' WR(2)=', WR(2)
         OK = .FALSE.
      ELSE
         PRINT *, 'CASE 3 PASS'
      END IF

C     Case 4: a non-zero diagonal element of AI is an error (code 2);
C     with IFAIL = 1 on entry it is returned, not a stop.
      AI(2,2)=1.0D0
      IFAIL = 1
      CALL OFHEIG(AR,MBRANC,AI,MBRANC,N,WR,VR,MBRANC,VI,MBRANC,
     $            WK1,WK2,WK3,IFAIL)
      IF (IFAIL.NE.2) THEN
         PRINT *, 'CASE 4 FAIL: IFAIL=', IFAIL, ' (expected 2)'
         OK = .FALSE.
      ELSE
         PRINT *, 'CASE 4 PASS'
      END IF

C     Case 5: a full complex Hermitian matrix (lower triangle) scaled
C     by S(K): the eigenvalues scale with it, the eigenvectors do not
C     change (relative to the unscaled solution, 1D-12). Without the
C     scaling in OFHEIG, 1D-200 lost the matrix (sums of squares
C     underflow; Codex audit 2026-10-30).
      N = 3
      DO 120 J = 1, N
         DO 110 I = 1, N
            BR(I,J) = 0.0D0
            BI(I,J) = 0.0D0
  110    CONTINUE
  120 CONTINUE
      BR(1,1) = 2.0D0
      BR(2,2) = 3.0D0
      BR(3,3) = 5.0D0
      BR(2,1) = 1.0D0
      BI(2,1) = 0.5D0
      BR(3,1) = 0.25D0
      BI(3,1) = -0.75D0
      BR(3,2) = 0.5D0
      BI(3,2) = 1.0D0
      IFAIL = 1
      CALL OFHEIG(BR,MBRANC,BI,MBRANC,N,W0,V0R,MBRANC,V0I,MBRANC,
     $            WK1,WK2,WK3,IFAIL)
      IF (IFAIL.NE.0) THEN
         PRINT *, 'CASE 5 FAIL: unscaled IFAIL=', IFAIL
         OK = .FALSE.
      END IF
      DO 160 K = 1, 6
         DO 140 J = 1, N
            DO 130 I = 1, N
               AR(I,J) = BR(I,J)*S(K)
               AI(I,J) = BI(I,J)*S(K)
  130       CONTINUE
  140    CONTINUE
         IFAIL = 1
         CALL OFHEIG(AR,MBRANC,AI,MBRANC,N,WR,VR,MBRANC,VI,MBRANC,
     $               WK1,WK2,WK3,IFAIL)
         ERR = 0.0D0
         DO 150 J = 1, N
            ERR = MAX(ERR, ABS(WR(J)/S(K)-W0(J))/ABS(W0(N)))
            DO 145 I = 1, N
               ERR = MAX(ERR, ABS(VR(I,J)-V0R(I,J)),
     $                   ABS(VI(I,J)-V0I(I,J)))
  145       CONTINUE
  150    CONTINUE
         IF (IFAIL.NE.0 .OR. .NOT.(ERR.LE.TOL)) THEN
            PRINT *, 'CASE 5 FAIL: scale', S(K), ' IFAIL=', IFAIL,
     $               ' deviation', ERR
            OK = .FALSE.
         END IF
  160 CONTINUE
      PRINT *, 'CASE 5 DONE'

C     Case 6: a NaN, then an infinite element is an error (code 3).
      DO 180 K = 1, 2
         DO 175 J = 1, N
            DO 170 I = 1, N
               AR(I,J) = BR(I,J)
               AI(I,J) = BI(I,J)
  170       CONTINUE
  175    CONTINUE
         X = 0.0D0
         IF (K.EQ.1) THEN
            AR(3,2) = X/X
         ELSE
            X = HUGE(X)
            AI(2,1) = X*2.0D0
         END IF
         IFAIL = 1
         CALL OFHEIG(AR,MBRANC,AI,MBRANC,N,WR,VR,MBRANC,VI,MBRANC,
     $               WK1,WK2,WK3,IFAIL)
         IF (IFAIL.NE.3) THEN
            PRINT *, 'CASE 6 FAIL: case', K, ' IFAIL=', IFAIL,
     $               ' (expected 3)'
            OK = .FALSE.
         END IF
  180 CONTINUE
      PRINT *, 'CASE 6 DONE'

      IF (OK) THEN
         PRINT *, 'PASS: OFHEIG eigensolver known-answer tests'
      ELSE
         PRINT *, 'FAIL: OFHEIG eigensolver known-answer tests'
         STOP 1
      END IF
      END
