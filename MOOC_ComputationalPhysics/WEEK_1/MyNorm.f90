PROGRAM MyNorm
    IMPLICIT NONE
    REAL :: u,v,w ! Components of a Vector
    REAL :: a     ! Variable to Store Norm of the Vector
    REAL :: Norm  ! Function Norm Returns a Real Value
    
    PRINT*, "Enter the Components of Vector: "
    READ*, u,v,w
    a = Norm(u,v,w)
    PRINT*, "Norm of the Vector is :",a
END PROGRAM MyNorm

FUNCTION Norm(u,v,w)
    IMPLICIT NONE
    REAL :: Norm
    REAL :: u,v,w

    Norm = SQRT(u**2+v**2+w**2)
END FUNCTION Norm