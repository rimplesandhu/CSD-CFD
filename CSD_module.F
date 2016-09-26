C     -----------------------------------------------------
C     Author: Rimple Sandhu, Carleton University
C     Email: rimple_sandhu@outlook.com
C     -----------------------------------------------------

      MODULE CSD_module
      USE COMMON
      IMPLICIT none
      CONTAINS

C     -----------------------------------------------------
C     DESC: Preprocessing for CFD2CSD: Part A
      SUBROUTINE PREPROC_A(NA, tX0_A, NS, NST, U, STATUS)
C     -----------------------------------------------------
C     IN--NA: No. of nodes on the CFD boundary
C     IN--tX0_A(NA,3): Undeformed CFD boundary node co-ordinates
C     OUT--NS: No. of nodes on the CSD boundary
C     OUT--NST: No. of total nodes in CSD mesh
C     IN--U: Free-stream velocity
C     OUT--STATUS: 1 = ran successfuly, otherwise 0
C     -----------------------------------------------------

C     Declare variables
      IMPLICIT None
      INTEGER                       :: STATUS, NA, NS, NST
      REAL(KIND=8)                  :: U
      REAL(KIND=8), DIMENSION(NA,3) :: X0_A, tX0_A
      REAL(KIND=8), DIMENSION(3,3)  :: tMat
      
C     Initialization
 101  FORMAT ('CSD::PREPROC_A: ', A)  
 102  FORMAT (I15 / I15 / I15 / ES15.8)   
      STATUS = 0
      tMat = 0.0D0
      tMat(1,1) = 1.0D0
      tMat(2,3) = 1.0D0
      tMat(3,2) = 1.0D0

C     Erase previous simulation data, excluding spline matrix
      WRITE(*,101,advance='no') 'Writing CFD/CSD grid...'
      STATUS = SYSTEM('cd CSD_module; make -s clean')
      CALL MYCHECK(STATUS, 'CSD::PREPROC_A: Make clean failed')

C     Write undeformed CFD grid data
      X0_A = MATMUL(tX0_A, tMat)
      CALL MYIO_RMAT(NA, 3, X0_A, 'WRITE', 'CSD_module/data/X0_A.dat')
      WRITE(*,*) 'Done.'

C     Extract CSD boundary grid data
      CALL EXTRACT_NS(NS, NST)

C     Write CFD simulation parameters, to be used by CSD
      OPEN (unit = 2, file = "CSD_module/data/PARAM_SIM.dat",
     +      form='formatted', action='write')
      WRITE (2,102) NST, NS, NA, U
      CLOSE(2)
      WRITE(*,101) 'Writing parameters in PARAM_SIM.dat...Done.'

C     Compile the code that assign the aeroloads to elmer 
      STATUS = SYSTEM('cd CSD_module; make -s compile')
      CALL MYCHECK(STATUS,'Make compile failed.')
      WRITE(*,101) 'Compiling AERO_PRESSURE.f90...Done.'
      END SUBROUTINE PREPROC_A

C     -----------------------------------------------------
C     DESC: Extract boundary nodes for CSD
      SUBROUTINE EXTRACT_NS(NS, NST)
C     -----------------------------------------------------
C     OUT--NS: No. of nodes on the CSD boundary
C     OUT--NST: No. of total nodes in CSD mesh
C     -----------------------------------------------------

C     Declare variables
      IMPLICIT None
      INTEGER       :: STATUS, NS, NST, i, iNode, j
      REAL(KIND=8)  :: xNode, yNode, zNode, PAR_S

C     Initialization   
 101  FORMAT ('CSD::EXTRACT_NS: ', A, I6)    
      NS = 0
      PAR_S = 0.0075

C     Get the total number of nodes
      OPEN (unit = 2, file = "CSD_module/output/mesh/mesh.header",
     +      form='formatted', action='read')
      READ (2,*) NST
      CLOSE(2)
      WRITE(*,101) 'No. of total nodes: ', NST

C     Extract the boundary nodes
      OPEN (unit = 3, file = "CSD_module/data/X0_S.dat",
     +      form='formatted', action='write')
      OPEN (unit = 4, file = "CSD_module/data/IND_S.dat",
     +      form='formatted', action='write')
      OPEN (unit = 5, file = "CSD_module/data/X0_ST.dat",
     +      form='formatted', action='write')
      OPEN (unit = 2, file = "CSD_module/output/mesh/mesh.nodes",
     +      form='formatted', action='read')
      DO i=1,NST
        READ (2,*) iNode, j, xNode, yNode, zNode
        WRITE(5,fmt = '(3ES24.12)') xNode, yNode, zNode
        IF  (ABS(zNode - PAR_S) < 1e-06 
     +  .OR. ABS(zNode + PAR_S) < 1e-06) THEN
          NS = NS + 1
          WRITE(3,fmt = '(3ES24.12)') xNode, yNode, zNode
          WRITE(4,*) iNode
        END IF
      END DO
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(5)
      WRITE(*,101) 'No. of boundary nodes: ', NS
      END SUBROUTINE EXTRACT_NS


C     -----------------------------------------------------
C     DESC: Preprocessing for CFD2CSD: Part B
      SUBROUTINE PREPROC_B(NA, NS, tX0_S, 
     +                    IFLAG_SPLINE, IFLAG_EIGEN, STATUS)
C     -----------------------------------------------------
C     IN--NA: N0. of nodes on the CFD boundary
C     IN--NS: N0. of nodes on the CSD boundary
C     OUT--tX0_S(NS,3): undeformed CSD boundary node co-ordinates
C     IN--IFLAG_SPLINE: 1 = calculate new spline matrix
C     IN--IFLAG_EIGEN: 1 = perform eigen value analysis
C     OUT--STATUS: 1 = ran successfuly, 0 otherwise
C     -----------------------------------------------------

C     Declare variables
      IMPLICIT None
      INTEGER                       :: STATUS, NS, NA, i, j
      INTEGER                       :: IFLAG_SPLINE, IFLAG_EIGEN
      REAL(KIND=8)                  :: rbfR
      REAL(KIND=8), DIMENSION(3,3)  :: tMat
      REAL(KIND=8), DIMENSION(NS,3) :: X0_S, tX0_S
      REAL(KIND=8), DIMENSION(NA,3) :: X0_A
      REAL(KIND=8), DIMENSION(NA,NS):: MAT_H
      
C     Initialization
 101  FORMAT ('CSD::PREPROC_B: ', A)
      rbfR = 20
      STATUS = 0
      tMat = 0.0D0
      tMat(1,1) = 1.0D0
      tMat(2,3) = 1.0D0
      tMat(3,2) = 1.0D0
      
C     Read the CSD boundary grid data
      CALL MYIO_RMAT(NA, 3, X0_A, 'READ','CSD_module/data/X0_A.dat')
      CALL MYIO_RMAT(NS, 3, X0_S, 'READ','CSD_module/data/X0_S.dat')
      tX0_S = MATMUL(X0_S,tMat)

C     Calculate and write the new spline matrix
      IF (IFLAG_SPLINE==1) THEN

C       Erase previous spline matrix
        WRITE(*,101,advance='no') 'Calculating new spline matrix...'
        STATUS = SYSTEM('cd CSD_module; make -s cleanSPLINE')
        CALL MYCHECK(STATUS,'CSD::PREPROC_B: Make cleanSPLINE failed')

C       Calculate and write the new spline matrix
        CALL RBF_HMAT(NA, NS, X0_A, X0_S, rbfR, MAT_H, STATUS)
        CALL MYIO_RMAT(NA,NS,MAT_H,'WRITE','CSD_module/data/MAT_H.dat')
        WRITE(*,*) 'Done.'
      ELSE
        WRITE(*,101) 'Loading available spline matrix...Done.'
      END IF

C     Perform eigen value analysis
      IF (IFLAG_EIGEN == 1) THEN
        WRITE(*,101,advance='no') 'Performing Eigen analysis...'
        STATUS = SYSTEM('cd CSD_module; make -s eigen')
        CALL MYCHECK(STATUS, 'CSD::PREPROC_B: Make eigen failed')
        WRITE(*,*) 'Done.'
      END IF  
      END SUBROUTINE PREPROC_B


C     -----------------------------------------------------
C     DESC: Generate initial deflection in CSD format
      SUBROUTINE INITIAl_X(NA, tXT_A, NS, NST, 
     +                     IFLAG_SPLINEB, STATUS)
C     -----------------------------------------------------
C     IN--NA: No. of nodes on the CFD boundary
C     IN--tXT_A(NA,3): deformed CFD boundary grid  
C     IN--NS: No. of nodes on the CSD boundary
C     IN--NST: No. of total nodes on the CSD mesh
C     IN--IFLAG_SPLINEB: 1=when CSD grid changed, 0 otherwise
C     OUT--STATUS: 1 = ran successfuly, 0 otherwise
C     -----------------------------------------------------

C     Declare variables
      IMPLICIT None
      INTEGER       :: STATUS, NA, NS, i, j, NST, IFLAG_SPLINEB
      REAL(KIND=8)  :: rbfR
      CHARACTER(len=100) :: row
      REAL(KIND=8), DIMENSION(NA,3)  :: X0_A, XT_A, tXT_A, DEF_A
      REAL(KIND=8), DIMENSION(NS,3)  :: DEF_S, X0_S
      REAL(KIND=8), DIMENSION(NST,3) :: X0_ST, DEF_ST
      REAL(KIND=8), DIMENSION(NA,NS) :: MAT_H
      REAL(KIND=8), DIMENSION(NS,NS) :: tempA, tempB
      REAL(KIND=8), DIMENSION(NS,NA) :: tempC
      REAL(KIND=8), DIMENSION(3,3)   :: tMat
      REAL(KIND=8), DIMENSION(NST,NS):: MAT_HB

C     Initialization
 101  FORMAT ('CSD::INITIAL_X: ', A)
      STATUS = 0
      rbfR = 20
      tMat = 0.0D0
      tMat(1,1) = 1.0D0
      tMat(2,3) = 1.0D0
      tMat(3,2) = 1.0D0

C     Load the data
      WRITE(*,101,advance='no') 'Generating initial.str...'
      CALL MYIO_RMAT(NA, NS, MAT_H, 'READ','CSD_module/data/MAT_H.dat')
      CALL MYIO_RMAT(NA,  3,  X0_A, 'READ','CSD_module/data/X0_A.dat')
      CALL MYIO_RMAT(NS,  3,  X0_S, 'READ','CSD_module/data/X0_S.dat')
      CALL MYIO_RMAT(NST, 3, X0_ST, 'READ','CSD_module/data/X0_ST.dat')

C     Calculate deformations (CFD Boundary)
      XT_A = MATMUL(tXT_A, tMat)
      DEF_A = XT_A - X0_A

C     Calculate deformations (CSD boundary)
      tempA = MATMUL(TRANSPOSE(MAT_H),MAT_H)
      CALL MAT_INVERSEB(tempA, tempB, NS)
      tempC = MATMUL(tempB, TRANSPOSE(MAT_H))
      DEF_S = MATMUL(tempC,DEF_A)

C     Only when CSD grid changed
      IF (IFLAG_SPLINEB == 1) THEN

C       Remove the previous version
        STATUS = SYSTEM('cd CSD_module; make -s cleanSPLINEB')
        CALL MYCHECK(STATUS,'CSD:INITIAL_X: Make cleanSPLINEB failed')

C       Calculate new spline matrix for CSD only
        CALL RBF_HMAT(NST, NS, X0_ST, X0_S, rbfR, MAT_HB, STATUS)
        CALL MYIO_RMAT(NST, NS, MAT_HB, 'WRITE',
     +                'CSD_module/data/MAT_HB.dat')
      ELSE
C       Read the version available
        CALL MYIO_RMAT(NST, NS, MAT_HB, 'READ',
     +                 'CSD_module/data/MAT_HB.dat')
      END IF

C     Calculate the deflected shape for new one
      DEF_ST = MATMUL(MAT_HB, DEF_S)

C     Write to initial.str
      OPEN (unit = 1, file = "CSD_module/output/restart/zero.str",
     +      form='formatted', action='read', status='old')
      OPEN (unit = 2, file = "CSD_module/output/restart/initial.str",
     +      form='formatted', action='write', status='unknown')
      DO j=1,NST+16
        READ(1,'(a)') row
        WRITE(2,'(a)')  row  
      END DO
      DO j=1,NST
        READ(1,'(a)') row
        WRITE(2,*)  DEF_ST(j,1) 
      END DO
      DO j=1,2
        READ(1,'(a)') row
        WRITE(2,'(a)')  row  
      END DO
      DO j=1,NST
        READ(1,'(a)') row
        WRITE(2,*)  DEF_ST(j,2) 
      END DO
      DO j=1,2
        READ(1,'(a)') row
        WRITE(2,'(a)')  row  
      END DO
      DO j=1,NST
        READ(1,'(a)') row
        WRITE(2,*)  DEF_ST(j,3) 
      END DO
      DO j=1,3*NST+6
        READ(1,'(a)') row
        WRITE(2,'(a)')  row
      END DO
      CLOSE(1)
      CLOSE(2)
      WRITE(*,*) 'Done.'

C     Copy intial.str to out.str and start.str
      STATUS = SYSTEM('cd CSD_module/output/restart;
     +                 cp initial.str restart.str')
      CALL MYCHECK(STATUS, 'CSD:INITIAL_X: Copying intial.str failed.')

      END SUBROUTINE INITIAL_X



C     -----------------------------------------------------
C     DESC: Transient simulation
      SUBROUTINE CSD2CFD(NS, tX0_S, tXT_S, F_S, NA, tX0_A, tXT_A,
     +                   F_A, NST, DT, T, IFLAG_IMPLICIT, 
     +                   IFLAG_RBF, STATUS)
C     -----------------------------------------------------
C     IN--NS: No. of nodes on the CSD boundary
C     IN--tX0_S(NS,3): undeformed CSD boundary nodes 
C     OUT--tXT_S(NS,3): deformed CSD boundary nodes
C     IN--F_S: Nodal forces at CSD boundary
C     IN--NA: No. of nodes on the aerodynamic boundary
C     IN--tX0_A(NS,3): undeformed CFD boundary nodes 
C     OUT--tXT_A(NS,3): deformed CFD boundary nodes
C     IN--F_A: Nodal Force on CFD boundary grid
C     IN--NST: No. of total nodes on the CSD mesh
C     IN--DT: Time integration step
C     IN--T: current time IN SIMULATION
C     IN--IFLAG_IMPLICIT = 1: Don't rewrite the start.str
C     IN--IFLAG_IMPLICIT = 0; Re-write the start.str
C     IN--IFLAG_RBF = 1; DO RBF in CSD module, in:F_A, out:tXT_A
C     IN--IFLAG_RBF = 0; DO RBF in CFD, in:F_S, out:tXT_S
C     OUT--STATUS: 1 = ran successfully, 0 otherwise.
C     -----------------------------------------------------

C     Declare the varibles
      IMPLICIT None
      INTEGER                       :: STATUS, NA, NS, NST, i, j
      INTEGER                       :: IFLAG_IMPLICIT, iNode, ITER
      INTEGER                       :: IFLAG_RBF, ja, jb, jc
      REAL(KIND=8)                  :: DT, T, xNode, yNode, zNode
      REAL(KIND=8)                  :: TIP_LE, TIP_TR, TIP_MID
      CHARACTER                     :: FILEN*100, tTEXT*100, tROW*100
      INTEGER,      DIMENSION(NS,3) :: IND_S
      REAL(KIND=8), DIMENSION(NS,3) :: tX0_S, tXT_S, XT_S
      REAL(KIND=8), DIMENSION(NA,3) :: tX0_A, tXT_A, XT_A
      REAL(KIND=8), DIMENSION(NA,NS):: MAT_H
      REAL(KIND=8), DIMENSION(NA)   :: F_A
      REAL(KIND=8), DIMENSION(NS)   :: F_S
      REAL(KIND=8), DIMENSION(3,3)  :: tMat

C     Format specifiers
 101  FORMAT ('CSD::CSD2CFD: ', A)
 102  FORMAT ('CSD::CSD2CFD: ', A, ES20.10)
 103  FORMAT (ES15.8 / ES15.8) 

C     Initialization
      tMat = 0.0D0
      tMat(1,1) = 1.0D0
      tMat(2,3) = 1.0D0
      tMat(3,2) = 1.0D0
      XT_S = 0.0
      ITER = 10000 + INT(T/DT)
      WRITE(*,*) '' 

C     Update the time par and initial CSD def for explicit sim
      IF (IFLAG_IMPLICIT == 0) THEN
C       Write the time parameters used by Elmer
        OPEN (unit = 2, file = "CSD_module/data/PARAM_TIME.dat",
     +        form='formatted', action='write')
        WRITE (2,103) DT, T
        CLOSE(2)

C       Write the time parameters to transient.sif file
        OPEN (unit = 2, file = "CSD_module/src/TRANSIENT_TEMP.SIF",
     +        form='formatted', action='read')
        OPEN (unit = 3, file = "CSD_module/src/TRANSIENT.SIF",
     +        form='formatted', action='write')
        READ(2,'(a)') tROW
        WRITE(3,'(A,ES20.12)')  '$dt = ', DT  
        READ(2,'(a)') tROW
        WRITE(3,'(A,ES20.12)')  '$T = ', T 
        DO j=1,140
          READ(2,'(a)') tROW
          WRITE(3,'(a)')  tROW  
        END DO
        CLOSE(2)
        CLOSE(3)

C       Change the starting deformation for explicit simulation
        STATUS = SYSTEM('cd CSD_module/output/restart;
     +                   cp restart.str start.str')
        CALL MYCHECK(STATUS, 'CSD:CSD2CFD: Copying out.str failed.')
        WRITE(*,101) 'EXPLICIT: ****** START *******'
      ELSE
        WRITE(*,101) 'IMPLICIT: ****** START *******'
      END IF

C     Calculate nodal forces F_S on CSD boundary
      IF (IFLAG_RBF == 1) THEN
        WRITE(*,101,advance='no') 'Calculating F_S from F_A...'
        CALL MYIO_RMAT(NA, NS,MAT_H,'READ','CSD_module/data/MAT_H.dat')
        CALL MYIO_RVEC(NA, F_A, 'WRITE', 'CSD_module/data/F_A.dat')
        F_S = MATMUL(TRANSPOSE(MAT_H),F_A)
        WRITE(*,*) 'Done.'
        WRITE(*,102) 'SUM(F_A)= ', SUM(F_A)
      ELSE
        WRITE(*,101) 'Using F_S available from CFD...Done.'
      END IF
      CALL MYIO_RVEC(NS, F_S, 'WRITE', 'CSD_module/data/F_S.dat')

C     CHECK: summation of structural forces
      WRITE(*,102) 'SUM(F_S)= ', SUM(F_S)
      WRITE(*,102) '       T = ', T
      WRITE(*,102) '      DT = ', DT

C     Run the CSD transient simulation
      STATUS = SYSTEM('cd CSD_module; make -s transient')

C     Extract CSD nodal deformations from raw elmer output
      CALL MYIO_IVEC(NS,IND_S,'READ','CSD_module/data/IND_S.dat')
      OPEN (unit=2,file = "CSD_module/output/transient/defSTraw.dat",
     +      form='formatted', action='read')
      i = 0
      DO WHILE (i .LT. NS)
        READ (2,*) ja, jb, jc, iNode, xNode, yNode, zNode
        CALL LOCATE_IN_VECTOR(NS, IND_S, iNode, j)
        IF (j .NE. 0) THEN
          XT_S(j,1) = xNode
          XT_S(j,2) = yNode
          XT_S(j,3) = zNode 
          i = i + 1 
        END IF
      END DO
      CLOSE(2)

C     Extract the tip responce at LE, TR and MID 
      TIP_LE  = XT_S(4,3)
      TIP_TR  = XT_S(3,3)
      TIP_MID = XT_S(42,3)
      OPEN (unit = 2, file = "CSD_module/data/XT_S_TIP.dat", 
     +      form='formatted', action='write', position = 'append')
      WRITE(2,*) T, TIP_LE, TIP_TR, TIP_MID
      CLOSE(2)

C     Calculate new CFD boundary grid deformation
      IF (IFLAG_RBF == 1) THEN
        WRITE(*,101,advance='no') 'Calculating tXT_A from XT_S...'
        XT_A = MATMUL(MAT_H,XT_S)
        tXT_A = MATMUL(XT_A,tMat)
        CALL MYIO_RMAT(NA, 3,XT_A,'WRITE', 'CSD_module/data/XT_A.dat')
        WRITE(*,*) 'Done.'
      ELSE
        tXT_S = MATMUL(XT_S,tMat)
        WRITE(*,101) 'Returning tXT_S to CFD...Done.'
      END IF
      CALL MYIO_RMAT(NS, 3, XT_S, 'WRITE', 'CSD_module/data/XT_S.dat')

C     Write the vtu files for every explicit iteration(for paraview)
      IF (IFLAG_IMPLICIT == 0) THEN       
        WRITE(FILEN,'(A,I5,A)') 'cd CSD_module/output/transient;
     +  cp transient.St.0001.vtu transient.CP.', ITER,'.vtu'
        STATUS = SYSTEM(FILEN)
        CALL MYCHECK(STATUS, 'CSD: Copying vtu file failed.')
        WRITE(*,101) 'EXPLICIT: ******* END ********'
      ELSE
        WRITE(*,101) 'IMPLICIT: ******* END ********'
      END IF
      WRITE(*,*) ''
      END SUBROUTINE CSD2CFD

C     -----------------------------------------------------
C     DESC: Calculation of spline matrix
      SUBROUTINE RBF_HMAT(N1, N2, X1, X2, R, H, STATUS)
C     -----------------------------------------------------
C     IN--N1: length of variable X1
C     IN--N2: length of variable X2
C     IN--X1(N1,3): Co-ordinates for variable 1
C     IN--X2(N2,3): Co-ordinates for variable 2
C     IN--R: Radius for rbf
C     OUT--H: New spline matrix
C     OUT--STATUS: 1= ran succesfully, else 0
C     -----------------------------------------------------

C     Declare variables
      IMPLICIT None
      INTEGER       :: STATUS, N1, N2, i, j
      REAL(KIND=8)  :: R
      REAL(KIND=8), DIMENSION(3) :: XRES
      REAL(KIND=8), DIMENSION(N1,3) :: X1
      REAL(KIND=8), DIMENSION(N2,3) :: X2
      REAL(KIND=8), DIMENSION(N1,N2):: H, M12  
      REAL(KIND=8), DIMENSION(N2,N2):: M22  
      REAL(KIND=8), DIMENSION(N2,4) :: P2 
      REAL(KIND=8), DIMENSION(N1,4) :: P1
      REAL(KIND=8), DIMENSION(N1,N2+4)   :: T1
      REAL(KIND=8), DIMENSION(N2+4,N2+4) :: T2, T2INV
      REAL(KIND=8), DIMENSION(N1,N2+4) :: T3
  
C     Initialize the arrays  
      STATUS = 0
      M12 = 0.0
      M22 = 0.0
      H = 0.0
      P1 = 1.0
      P2 = 1.0
      T1 = 0.0
      T2 = 0.0

C     Actual calculations
      DO i=1,N1
        P1(i,2:4) = X1(i,:)
        DO j=1,N2
          XRES = X1(i,:) - X2(j,:)
          CALL RBF_BECKERT_WENDLAND(3, XRES, R, M12(i,j)) 
        END DO
      END DO
      DO i=1,N2
        P2(i,2:4) = X2(i,:)
        DO j=1,N2
          XRES = X2(i,:) - X2(j,:)
          CALL RBF_BECKERT_WENDLAND(3, XRES, R, M22(i,j)) 
        END DO
      END DO

C     Array creation and adjustment
      T1(:,1:N2) = M12
      T1(:,N2+1:N2+4) = P1
      T2(1:N2,1:N2) = M22
      T2(N2+1:N2+4,1:N2) = TRANSPOSE(P2)
      T2(1:N2,N2+1:N2+4) = P2
      CALL MAT_INVERSEB(T2, T2INV, N2+4)
      T3 = MATMUL(T1,T2INV)
      H = T3(:,1:N2)

      END SUBROUTINE RBF_HMAT


C     -----------------------------------------------------
C     DESC: Radial Basis Function
      SUBROUTINE RBF_BECKERT_WENDLAND(N, X, R, RES)
C     -----------------------------------------------------
C     IN--N: length of vector X
C     IN--X: vector whose rbf value needs to be calculated
C     IN--R: radius of radial basis function
C     OUT--RES: the rbf value at given X
C     -----------------------------------------------------

C     Declare variables
      IMPLICIT None
      INTEGER       :: N, i, j
      REAL(KIND=8)  :: R, Xnorm, RES
      REAL(KIND=8), DIMENSION(N) :: X

C     Calculate norm of vector X
      CALL MYNORM2(N, X, Xnorm)

C     The beckert and wendland basis function
      IF (Xnorm < R) THEN
        RES = ((1.0 - Xnorm/R)**4)*(4.0*Xnorm/R + 1.0)
      ELSE
        RES = 4.0*Xnorm/R + 1.0
      END IF
      END SUBROUTINE RBF_BECKERT_WENDLAND


      END MODULE CSD_module

