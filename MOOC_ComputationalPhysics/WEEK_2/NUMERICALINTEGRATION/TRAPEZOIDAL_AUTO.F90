! PROGRAM TO EVALUATE INTEGRAL OF THE FUNCTION 4.0/(1+X^2) USING THE TRAPEZOIDAL METHOD
! THE INTEGRAL ANALYTICALLY EVALUATES TO PI DIVIDED BY FOUR

PROGRAM TRAPEZOIDAL_AUTO
    IMPLICIT NONE
    REAL*8      :: H,A,B,X
    REAL*8      :: FUNC
    REAL*8      :: INTEGRAL_SUM,PI,ERROR
    INTEGER     :: N,I
    
    OPEN(UNIT = 1, FILE = "TRAPEZOIDAL_PI.DAT")
    
    PI = 3.1415926535897932
    INTEGRAL_SUM = 0.0D0
    
    N = 1
    
    
    WRITE(*,*) "ENTER THE LOWER LIMIT OF INTEGRATION: "
    READ(*,*)  A
    
    WRITE(*,*) "ENTER THE UPPER LIMIT OF INTEGRATION: "
    READ(*,*)  B
    
    DO 
    
        N = N*5
        
        H = (B-A)/REAL(N)
        
        INTEGRAL_SUM = (FUNC(A)+FUNC(B))/2.0D0
        DO I=1,N-1,1
            X = A + I*H
            INTEGRAL_SUM = INTEGRAL_SUM + FUNC(X)
        END DO
        INTEGRAL_SUM = INTEGRAL_SUM*H
        ERROR = PI - INTEGRAL_SUM
        
        WRITE(1,*) "INTEGRAL = ",INTEGRAL_SUM," N = ", N, " ERROR = ", ERROR
        IF (N .GE. 1000000000) EXIT
        
    END DO
END PROGRAM TRAPEZOIDAL_AUTO


FUNCTION FUNC(X) RESULT(RES)
    IMPLICIT NONE
    REAL*8              :: RES
    REAL*8, INTENT(IN)  :: X

    RES = 4.0D0/(1+X**2)
END FUNCTION FUNC