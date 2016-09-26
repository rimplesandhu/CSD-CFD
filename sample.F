C     -----------------------------------------------------
C     Author: Rimple Sandhu, Carleton University
C     Email: rimple_sandhu@outlook.com
C     -----------------------------------------------------

      PROGRAM SAMPLE

      USE CSD_module
      USE COMMON
      IMPLICIT NONE

C     Declare variables
      REAL(KIND=8) :: DT, U, T
      INTEGER      :: STATUS, NA, NS, NST, i, j, IFLAG_RBF
      INTEGER      :: IFLAG_IMPLICIT, NSIM_EXPLICIT, NSIM_IMPLICIT
      INTEGER      :: IFLAG_SPLINEB, IFLAG_SPLINE, IFLAG_EIGEN
      REAL(KIND=8), DIMENSION(8250)   :: testP
      REAL(KIND=8), DIMENSION(8250,3) :: testX0_A, testXT_A
      REAL(KIND=8), ALLOCATABLE:: testX0_S(:,:), testXT_S(:,:), F_S(:)

C     Define the simulation parameters  
      NA = 8250
      DT = 0.0002d0
      U = 0.901*340.3
      IFLAG_SPLINE  = 1
      IFLAG_SPLINEB = 1
      IFLAG_EIGEN   = 1
      IFLAG_RBF     = 0


C     Read aerodynamic nodal force from data_for_testing folder
      CALL MYIO_RVEC(NA, testP, 'READ',
     +      'data_for_testing/pressure1D_mode1to4.dat')

C     Read undeformed CFD boundary grid
      CALL MYIO_RMAT(NA, 3, testX0_A, 'READ',
     +           'data_for_testing/undefAe.dat')

C     Read initial deformed CFD boundary grid
      CALL MYIO_RMAT(NA, 3, testXT_A, 'READ',
     +           'data_for_testing/defAe_mode1to4.dat')

C     Preprocessing part A, read the NS, not X0_S
      CALL PREPROC_A(NA, testX0_A, NS, NST, U, STATUS)
      CALL MYCHECK(STATUS, 'CSD::PREPROC_A: failed')

C     Allocate to the matrices
      ALLOCATE(testX0_S(NS,3))
      ALLOCATE(testXT_S(NS,3))
      ALLOCATE(F_S(NS))
      
C     Read X0_S, calculate spline mat, perform eigen analysis 
      CALL PREPROC_B(NA, NS, testX0_S,IFLAG_SPLINE,IFLAG_EIGEN,STATUS)
      CALL MYCHECK(STATUS, 'CSD::PREPROC_B: failed')

C     Generate initial deflection for CSD grid 
      CALL INITIAl_X(NA, testXT_A, NS, NST, IFLAG_SPLINEB, STATUS)
      CALL MYCHECK(STATUS, 'CSD:INITIAL_X: failed')

C     do this if static defomation from elmer needed
C     Copy intial_static.str to restart.str
C      STATUS = SYSTEM('cd CSD_module/output/restart;
C     +                 cp initial_static.str restart.str')

C     Subsequent Implicit/explicit simulations
      NSIM_IMPLICIT = 1
      NSIM_EXPLICIT = 2000
C      testP = 0
      F_S = 0
      T = 0.0
      DO i=1,NSIM_EXPLICIT
        DO IFLAG_IMPLICIT=1,NSIM_IMPLICIT
         CALL CSD2CFD(NS, testX0_S, testXT_S, F_S, 
     +                NA, testX0_A, testXT_A, testP, 
     +                NST, DT, T, IFLAG_IMPLICIT-1, IFLAG_RBF, STATUS)
        END DO
        T = T + DT
      END DO

      END PROGRAM SAMPLE
