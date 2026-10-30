C     Known-answer tests for OFHEIG, Florence's complex Hermitian
C     eigensolver (../Florence_f.f): eigenvalues, the phase convention
C     of the eigenvectors, and IFAIL on return (0 after success with
C     any value on entry; the failure code, without stopping, when
C     IFAIL is not 0 on entry).
      PROGRAM TESTEIGENSOLVER
      IMPLICIT NONE
      INTEGER N, MBRANC
      PARAMETER (MBRANC=3)
      DOUBLE PRECISION AR(MBRANC,MBRANC), AI(MBRANC,MBRANC)
      DOUBLE PRECISION WR(MBRANC), VR(MBRANC,MBRANC), VI(MBRANC,MBRANC)
      DOUBLE PRECISION WK1(MBRANC), WK2(MBRANC), WK3(MBRANC)
      INTEGER IFAIL
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

      IF (OK) THEN
         PRINT *, 'PASS: OFHEIG eigensolver known-answer tests'
      ELSE
         PRINT *, 'FAIL: OFHEIG eigensolver known-answer tests'
         STOP 1
      END IF
      END
