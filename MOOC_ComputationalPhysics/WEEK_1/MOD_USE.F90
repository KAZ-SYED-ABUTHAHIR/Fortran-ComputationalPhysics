USE NORM_MODULE_SEP
IMPLICIT NONE
    REAL :: a,x,y,z

    PRINT*, "Enter the vector components: "
    READ*, x,y,z

    a = NORM(x,y,z)
    PRINT*, "The Norm of the Vector is : ",a
END