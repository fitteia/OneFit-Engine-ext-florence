C     OFHEIG with IFAIL = 0 on entry stops the program on a failure:
C     with a message on unit 0 and exit status 1, so a fit does not end
C     as if it had succeeded (run_tests.sh checks both; reaching the
C     PRINT below is a failure). Here: a non-zero diagonal element of
C     AI (code 2).
      PROGRAM TESTSOLVERFATAL
      IMPLICIT NONE
      DOUBLE PRECISION AR(2,2), AI(2,2), WR(2), VR(2,2), VI(2,2)
      DOUBLE PRECISION WK1(2), WK2(2), WK3(2)
      INTEGER IFAIL
      AR(1,1) = 1.0D0
      AR(2,1) = 0.0D0
      AR(2,2) = 2.0D0
      AI(1,1) = 1.0D0
      AI(2,1) = 0.0D0
      AI(2,2) = 0.0D0
      IFAIL = 0
      CALL OFHEIG(AR,2,AI,2,2,WR,VR,2,VI,2,WK1,WK2,WK3,IFAIL)
      PRINT *, 'returned with IFAIL', IFAIL, ' instead of stopping'
      END
