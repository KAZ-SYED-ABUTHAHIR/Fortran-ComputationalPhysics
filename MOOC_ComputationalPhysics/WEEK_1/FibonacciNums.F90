PROGRAM FibonacciNums
    IMPLICIT NONE
    INTEGER :: X_N,X_NMINUS1,X_NMINUS2,N = 3
    REAL :: GoldenRatioEst
    X_NMINUS2 = 1
    X_NMINUS1 = 1
    PRINT*, "Fibonacci Number ",1," is ",X_NMINUS2
    PRINT*, "Fibonacci Number ",2," is ",X_NMINUS1
    GoldenRatioEst = (X_NMINUS2/X_NMINUS1)
    PRINT '(A30 f9.4)', "The Golden Ratio Estimation is ",GoldenRatioEst
    DO  WHILE (N < 47)
        X_N = X_NMINUS1 + X_NMINUS2
        PRINT*, "Fibonacci Number ",N," is ",X_N
        GoldenRatioEst = (X_NMINUS1/REAL(X_N))
        PRINT '(A30 f9.4)', "The Golden Ratio Estimation is ",GoldenRatioEst
        X_NMINUS2 = X_NMINUS1
        X_NMINUS1 = X_N
        N = N+1        
    END DO
END PROGRAM FibonacciNums