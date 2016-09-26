C     -----------------------------------------------------
C     Author: Rimple Sandhu, Carleton University
C     Email: rimple_sandhu@outlook.com
C     -----------------------------------------------------


      FUNCTION getFORCE(MODEL, iNode, curT ) RESULT(iF_S)
      
      USE COMMON
      USE DefUtils
      IMPLICIT None

      TYPE(Model_t)              :: MODEL
      TYPE(Solver_t),    POINTER :: SOLVER 
      TYPE(Element_t),   POINTER :: curE

      INTEGER       :: NST, NS, NA,  NN_BELEM
      INTEGER       :: iNode, i, j, k,  NE_ACTIVE
      INTEGER       :: NN_ELEM, NN_BULK, NE_BOUND
      LOGICAL       :: firstVisit = .TRUE.
      REAL(KIND=dp) :: iF_S, curT, U
      REAL(KIND=dp), ALLOCATABLE :: F_S(:)
      INTEGER,       ALLOCATABLE :: IND_S(:)

      SAVE NST, NS, NA, U, F_S, IND_S
      SAVE NN_ELEM, NN_BULK, NE_BOUND, NE_ACTIVE, NN_BELEM
      SAVE firstVisit, k
      
 102  FORMAT (I15 / I15 / I15 / ES15.8) 

      IF (firstVisit) THEN

C       Pointer to the solver accesing this function
        SOLVER       => CurrentModel % Solver
        WRITE(*,*) ''
        WRITE(*,'(A)') 'SANDHU:  ************ START *************'

C       Number of bulk nodes
        NN_BULK  = SOLVER % mesh % NumberOfNodes
        WRITE(*,'(A,i6)') 'SANDHU:        # bulk nodes:', NN_BULK

C       Number of nodes in each element
        NN_ELEM   = SOLVER % mesh % MaxElementNodes
        WRITE(*,'(A,i6)') 'SANDHU:  # nodes in element:', NN_ELEM 

C       Number of nodes in boundary element        
        NN_BELEM = GetElementNOFNodes(CurrentModel % CurrentElement)
        WRITE(*,'(A,i6)') 'SANDHU:# nodes in B element:', NN_BELEM

C       Number of boundary elements
        NE_BOUND  = SOLVER % mesh % NumberOfBoundaryElements
        WRITE(*,'(A,i6)') 'SANDHU: # boundary elements:', NE_BOUND

C       Number of active elements for this boundary condition
        NE_ACTIVE = SOLVER % NumberOfActiveElements
        WRITE(*,'(A,i6)') 'SANDHU:   # active elements:', NE_ACTIVE
        
C       Read number of structural nodes from external file
        OPEN (unit = 2, file = "data/PARAM_SIM.dat",
     +        form='formatted', action='read')
        READ (2,102) NST, NS, NA, U
        CLOSE(2)

C       Allocation of arrays and vectors
        ALLOCATE(F_S(NS))
        ALLOCATE(IND_S(NS))
        CALL MYIO_RVEC(NS, F_S, 'READ', 'data/F_S.dat')
        CALL MYIO_IVEC(NS, IND_S, 'READ', 'data/IND_S.dat')

        firstVisit = .FALSE.
        k = 0 
      END IF

C     Calculate the pressure at the node
      CALL LOCATE_IN_VECTOR(NS, IND_S, iNode, j)
      IF (j .NE. 0) THEN 
        iF_S = F_S(j)
        k = k + 1
        !WRITE(*,*) k, iNode, curT, iF_S
      ELSE 
        iF_S = 0.0D0
      END IF

C     Get the co-ordinates of triangle      
      !curE => CurrentModel % CurrentElement
      !indE = curE % NodeIndexes
      !DO i=1,NN_BELEM
      !  eX(i,1) = MODEL % Nodes % x(indE(i)) 
      !  eX(i,2) = MODEL % Nodes % y(indE(i)) 
      !  eX(i,3) = MODEL % Nodes % z(indE(i))
      !END DO

      END FUNCTION getFORCE


C     ========================================================
C     Lessons learned
C     - permindexes are same for all deflcections 1, 2, 3
C     - deflection ofr all nodes are same at every node iteration
