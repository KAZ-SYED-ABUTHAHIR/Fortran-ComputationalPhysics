PROGRAM WEEK1_Q6_Q7
    IMPLICIT NONE
    INTEGER :: X_N,X_NMINUS1,X_NMINUS2,N = 3,SUM_FIBO = 0
    
    X_NMINUS2 = 1
    X_NMINUS1 = 1
    
    PRINT*, "Fibonacci Number ",1," is ",X_NMINUS2
    PRINT*, "Fibonacci Number ",2," is ",X_NMINUS1
    
    SUM_FIBO = SUM_FIBO + X_NMINUS1 + X_NMINUS2
    
    DO  WHILE (N <= 40)
        X_N = X_NMINUS1 + X_NMINUS2
        PRINT*, "Fibonacci Number ",N," is ",X_N,"The running sum is ",SUM_FIBO
        
        X_NMINUS2 = X_NMINUS1
        X_NMINUS1 = X_N
        SUM_FIBO = SUM_FIBO + X_N
        N = N+1        
    END DO
    PRINT*,"______________________________________________"
    PRINT*,"The Sum of 40 Fibonacci Numbers is ",SUM_FIBO
END PROGRAM WEEK1_Q6_Q7