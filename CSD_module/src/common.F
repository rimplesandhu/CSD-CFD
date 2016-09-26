C     -----------------------------------------------------
C     Author: Rimple Sandhu, Carleton University
C     Email: rimple_sandhu@outlook.com
C     -----------------------------------------------------

      MODULE COMMON
      IMPLICIT NONE
      CONTAINS

C     -----------------------------------------------------
C     DESC: Find out if the entry is present in the vector
      SUBROUTINE LOCATE_IN_VECTOR(inN, inVEC, inVAL, outIND)
C     -----------------------------------------------------
C     IN--inN: size of vector inVEC
C     IN--inVEC: vector that needs to be checked for value
C     IN--inVAL: the value that needs to be checked
C     OUT--outIND: index of value if found, else 0
C     -----------------------------------------------------

      IMPLICIT NONE
      INTEGER               :: ij, inN, inVAL, outIND
      INTEGER, DIMENSION(inN) :: inVEC

      outIND = 0
      DO ij = 1, inN
        IF (inVEC(ij) .eq. inVAL) THEN
          outIND = ij
        END IF 
      END DO
      END SUBROUTINE LOCATE_IN_VECTOR 
       

C     -----------------------------------------------------     
C     DESC: Calculate norm2
      SUBROUTINE MYNORM2(N, X, NRM)
C     -----------------------------------------------------
C     IN--N: size of vector X
C     IN--X: vector whose norm is to be found
C     OUT--NRM: Calculated norm
C     -----------------------------------------------------

      IMPLICIT NONE
      INTEGER      :: N, i
      REAL(KIND=8) :: NRM
      REAL(KIND=8), DIMENSION(N) :: X

      NRM = 0.0
      DO i=1,N
        NRM = NRM + X(i)*X(i)
      END DO
      NRM = SQRT(NRM)
      END SUBROUTINE MYNORM2

C     -----------------------------------------------------
C     DESC: Inverse of matrix
      SUBROUTINE MAT_INVERSE(a,c,n)
C     -----------------------------------------------------
      IMPLICIT NONE
      INTEGER     :: n, i, j, k
      REAL(KIND=8), DIMENSION(n,n):: a, c, L, U
      REAL(KIND=8), DIMENSION(n):: b, d, x
      REAL(KIND=8):: coeff

C     step 0: initialization for matrices L and U and b
C     Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

C     step 1: forward elimination
      DO k=1, n-1
        DO i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          DO j=k+1,n
            a(i,j) = a(i,j)-coeff*a(k,j)
          END DO
        END DO
      END DO

C     Step 2: prepare L and U matrices 
C     L matrix is a matrix of the elimination coefficient
C     + the diagonal elements are 1.0
      DO i=1,n
        L(i,i) = 1.0
      END DO
C     U matrix is the upper triangular part of A
      DO j=1,n
        DO i=1,j
          U(i,j) = a(i,j)
        END DO
      END DO

C     Step 3: compute columns of the inverse matrix C
      DO k=1,n
        b(k)=1.0
        d(1) = b(1)
C       Step 3a: Solve Ld=b using the forward substitution
        DO i=2,n
          d(i)=b(i)
          DO j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          END DO
        END DO
C       Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        DO i = n-1,1,-1
          x(i) = d(i)
          DO j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          END DO
          x(i) = x(i)/u(i,i)
        END DO
C     Step 3c: fill the solutions x(n) into column k of C
        DO i=1,n
          c(i,k) = x(i)
        END DO
        b(k)=0.0
      END DO
      END SUBROUTINE MAT_INVERSE


C     -----------------------------------------------------
C     DESC: Check the error and return a text
      SUBROUTINE MYCHECK(STATUS, ETEXT)
C     -----------------------------------------------------

      IMPLICIT NONE
      INTEGER   :: STATUS
      CHARACTER :: ETEXT*(*)

      IF (STATUS .NE. 0) THEN
        WRITE(*,*) 
        WRITE(*,'(2A)') 'CSD::ERROR:', ETEXT
        CALL EXIT(STATUS)
      END IF
      END SUBROUTINE MYCHECK


C     -----------------------------------------------------
C     DESC: Read and write the real matrix
      SUBROUTINE MYIO_RMAT(NR, NC, X, IO, IOF)
C     -----------------------------------------------------
      IMPLICIT NONE
      INTEGER    :: NR, NC, i, j
      CHARACTER  :: IO*(*), IOF*(*)
      REAL(KIND=8), DIMENSION(NR, NC) :: X
      i = 0
      IF (IO .EQ. 'WRITE') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='write')
        DO i=1,NR
          DO j=1,NC
            WRITE(2,'(ES24.12)',advance='no') X(i,j)
          END DO
        WRITE(2,*)
        END DO
        CLOSE(2)
      ELSE IF (IO .EQ. 'READ') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='read')
        DO i=1,NR
          DO j=1,NC
            READ(2,'(ES24.12)',advance='no') X(i,j)
          END DO
        READ(2,*)
        END DO
        CLOSE(2)
      ELSE
        WRITE(*,*) 
        WRITE(*,'(A)') 'CSD::MYIO_RMAT: Only READ or WRITE!'
        CALL EXIT(i)
      END IF
      END SUBROUTINE MYIO_RMAT

C     -----------------------------------------------------
C     DESC: Read and write the integer matrix
      SUBROUTINE MYIO_IMAT(NR, NC, X, IO, IOF)
C     -----------------------------------------------------
      IMPLICIT NONE
      INTEGER    :: NR, NC, i, j
      CHARACTER  :: IO*(*), IOF*(*)
      INTEGER, DIMENSION(NR, NC) :: X
      i = 0
      IF (IO .EQ. 'WRITE') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='write')
        DO i=1,NR
          DO j=1,NC
            WRITE(2,'(I16)',advance='no') X(i,j)
          END DO
        WRITE(2,*)
        END DO
        CLOSE(2)
      ELSE IF (IO .EQ. 'READ') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='read')
        DO i=1,NR
          DO j=1,NC
            READ(2,'(I16)',advance='no') X(i,j)
          END DO
        READ(2,*)
        END DO
        CLOSE(2)
      ELSE
        WRITE(*,*) 
        WRITE(*,'(A)') 'CSD::MYIO_IMAT: Only READ or WRITE!'
        CALL EXIT(i)
      END IF
      END SUBROUTINE MYIO_IMAT


C     -----------------------------------------------------
C     DESC: Read and write the integer vector
      SUBROUTINE MYIO_IVEC(NR, X, IO, IOF)
C     -----------------------------------------------------
      IMPLICIT NONE
      INTEGER    :: NR, i, j
      CHARACTER  :: IO*(*), IOF*(*)
      INTEGER, DIMENSION(NR) :: X
      i = 0
      IF (IO .EQ. 'WRITE') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='write')
        DO i=1,NR
            WRITE(2,'(I16)') X(i)
        END DO
        CLOSE(2)
      ELSE IF (IO .EQ. 'READ') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='read')
        DO i=1,NR
          READ(2,'(I16)') X(i)
        END DO
        CLOSE(2)
      ELSE
        WRITE(*,*) 
        WRITE(*,'(A)') 'CSD::MYIO_IVEC: Only READ or WRITE!'
        CALL EXIT(i)
      END IF
      END SUBROUTINE MYIO_IVEC

C     -----------------------------------------------------
C     DESC: Read and write the real vector
      SUBROUTINE MYIO_RVEC(NR, X, IO, IOF)
C     -----------------------------------------------------
      IMPLICIT NONE
      INTEGER    :: NR, i, j
      CHARACTER  :: IO*(*), IOF*(*)
      REAL(KIND=8), DIMENSION(NR) :: X
      i = 0
      IF (IO .EQ. 'WRITE') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='write')
        DO i=1,NR
          WRITE(2,'(ES24.12)') X(i)
        END DO
        CLOSE(2)
      ELSE IF (IO .EQ. 'READ') THEN
        OPEN (unit = 2, file = IOF, form='formatted', action='read')
        DO i=1,NR
          READ(2,'(ES24.12)') X(i)
        END DO
        CLOSE(2)
      ELSE
        WRITE(*,*) 
        WRITE(*,'(A)') 'CSD::MYIO_RVEC: Only READ or WRITE!'
        CALL EXIT(i)
      END IF
      END SUBROUTINE MYIO_RVEC


C     -----------------------------------------------------
      SUBROUTINE MAT_INVERSEB(a,b,n)
C     DESC: Calculate the inverse using Gauss-Jordan elimination
C     inverse of matrix a(n,n) is stored in the matrix b(n,n)
C     -----------------------------------------------------

      implicit none
      integer :: i,j,k,l,m,n,irow
      real*8:: big,a(n,n),b(n,n),dum
      real*8,allocatable:: a_tmp(:,:)

C     save the original matrix
      allocate(a_tmp(n,n))
      do i = 1,n
      do j = 1,n
      a_tmp(i,j) = a(i,j)
      end do
      end do

C     build the identity matrix
      do i = 1,n
      do j = 1,n
      b(i,j) = 0.0
      end do
      b(i,i) = 1.0
      end do

C     this is the big loop over all the columns of a(n,n)
      do i = 1,n 
C     in case the entry a(i,i) is zero, 
C     we need to find a good pivot; this pivot is chosen as 
C     the largest value on the column i from a(j,i) with j = 1,n
      big = abs(a(i,i))
      do j = i,n
      if (abs(a(j,i)).gt.big) then
      big = abs(a(j,i))
      irow = j
      end if
      end do
C     interchange lines i with irow for both a() and b() matrices
      if (big.gt.abs(a(i,i))) then
      do k = 1,n
      dum = a(i,k) ! matrix a()
      a(i,k) = a(irow,k)
      a(irow,k) = dum
      dum = b(i,k) ! matrix b()
      b(i,k) = b(irow,k)
      b(irow,k) = dum
      end do
      end if
C     divide all entries in line i from a(i,j) by the value a(i,i);
C     same operation for the identity matrix
      dum = a(i,i)
      do j = 1,n
      a(i,j) = a(i,j)/dum
      b(i,j) = b(i,j)/dum
      end do
C     make zero all entries in the column a(j,i); 
C     same operation for indent()
      do j = i+1,n
      dum = a(j,i)
      do k = 1,n
      a(j,k) = a(j,k) - dum*a(i,k)
      b(j,k) = b(j,k) - dum*b(i,k)
      end do
      end do
      end do
C     substract appropiate multiple of row j from row j-1
      do i = 1,n-1
      do j = i+1,n
      dum = a(i,j)
      do l = 1,n
      a(i,l) = a(i,l)-dum*a(j,l)
      b(i,l) = b(i,l)-dum*b(j,l)
      end do
      end do
      end do

C     reset the original matrix
      do i = 1,n
      do j = 1,n
      a(i,j) = a_tmp(i,j)
      end do
      end do
      deallocate(a_tmp)

      END SUBROUTINE MAT_INVERSEB

      END MODULE COMMON	
