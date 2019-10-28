PROGRAM  FormatSpec
    IMPLICIT NONE
    INTEGER:: i = 12345
    REAL:: r = 3.14159265358979
    REAL, DIMENSION(5)::  v = (/ 1.1,1.2,1.9,1.3,1.5 /)   
    PRINT '(i5)',i
    PRINT '(f10.8)',r
    PRINT '(5f8.3)',v
    PRINT '(e10.3)',r
    PRINT '(i5 f15.8)',i,r
END PROGRAM FormatSpec