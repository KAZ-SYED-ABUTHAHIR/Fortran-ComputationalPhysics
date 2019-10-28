! PROGRAM TO EVALUATE INTEGRAL OF THE FUNCTION 1.0/(1+X^2) USING THE TRAPEZOIDAL METHOD
! THE INTEGRAL ANALYTICALLY EVALUATES TO PI DIVIDED BY FOUR

PROGRAM TRAPEZOIDAL
    IMPLICIT NONE
    REAL*8      :: H,A,B,X
    REAL*8      :: FUNC
    REAL*8      :: INTEGRAL_SUM
    INTEGER     :: N,I

    INTEGRAL_SUM = 0.0D0

    WRITE(*,*) "ENTER THE NUMBER OF GRID POINTS: "
    READ(*,*)  N

    WRITE(*,*) "ENTER THE LOWER LIMIT OF INTEGRATION: "
    READ(*,*)  A

    WRITE(*,*) "ENTER THE UPPER LIMIT OF INTEGRATION: "
    READ(*,*)  B

    H = (B-A)/REAL(N)

    INTEGRAL_SUM = (FUNC(A)+FUNC(B))/2.0D0
    DO I=1,N-1,1
        X = A + I*H
        INTEGRAL_SUM = INTEGRAL_SUM + FUNC(X)
    END DO
    INTEGRAL_SUM = INTEGRAL_SUM*H

    PRINT*, "THE INTEGRAL BY THE TRAPEZOIDAL METHOD IS ",INTEGRAL_SUM
    PRINT*, "FOUR TIMES THE INTEGRAL BY THE TRAPEZOIDAL METHOD IS ",4.0D0*INTEGRAL_SUM
END PROGRAM TRAPEZOIDAL


FUNCTION FUNC(X) RESULT(RES)
    IMPLICIT NONE
    REAL*8              :: RES
    REAL*8, INTENT(IN)  :: X

    RES = 1.0D0/(1+X**2)
END FUNCTION FUNC
