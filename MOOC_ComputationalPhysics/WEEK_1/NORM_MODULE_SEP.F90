MODULE NORM_MODULE_SEP 
IMPLICIT NONE
CONTAINS
    FUNCTION Norm(u,v,w) RESULT(res)
        IMPLICIT NONE
        REAL, INTENT(IN) :: u,v,w
        REAL :: res
        res = SQRT(u**2+v**2+w**2)
    END FUNCTION Norm
END MODULE NORM_MODULE_SEP