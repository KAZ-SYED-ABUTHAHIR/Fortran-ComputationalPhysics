!Program to Study Importance Sampling Method
Program ImportanceSampling
    IMPLICIT NONE
    INTEGER             ::      i,n,j,num
    REAL*8              ::      r,x
    REAL*8              ::      Integrand,WeightFunc,RandomTransform
    REAL*8              ::      Integral,IS_Integral
    REAL*8, PARAMETER   ::      e = 2.71828183d0, Exact = 0.746824d0
    REAL*8              ::      Error,IS_Error

    num = 8
    n = 1
    Write(*,"(A20,A30,A39,A24,A34)") "Numeber of Samples","Integral Value","Integral Value(Importance Sampled)","Error" &
                                    ,"Error(Importance Sampling)"
    DO j = 1,num
        Integral    = 0.0d0
        IS_Integral = 0.0d0
        n = n*10
        DO i = 1,n 
           CALL RANDOM_NUMBER(r)
           Integral  = Integral + Integrand(r)
           r = r * ((e-1.0d0)/e)
           x = RandomTransform(r)
           IS_Integral = IS_Integral + Integrand(x)/(WeightFunc(x)*e/(e-1))
        END DO
        Integral    = Integral / Real(n)
        IS_Integral = IS_Integral / Real(n)
        Error = Exact - Integral
        IS_Error = Exact - IS_Integral
        Write(*,"(I20,F30.8,F39.8,F24.6,F34.6)") n,Integral,IS_Integral,Error,IS_Error
    END DO
End Program ImportanceSampling 

Function Integrand(x) Result(Result)
    REAL*8              :: x
    REAL*8              :: Result

    Result = EXP(-x**2)
End Function Integrand

Function WeightFunc(x) Result(Result)
    REAL*8              :: x
    REAL*8              :: Result

    Result = EXP(-x)
End Function WeightFunc

Function RandomTransform(u) Result(Result)
    REAL*8              :: u
    REAL*8              :: Result

    Result = -LOG(1.0d0 - u)
End Function RandomTransform
