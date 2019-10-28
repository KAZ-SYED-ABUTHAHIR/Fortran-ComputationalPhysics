! Program to sort an Array using bubble sort

PROGRAM Sort

    IMPLICIT NONE
    INTEGER, DIMENSION(:), ALLOCATABLE :: NUMBERS
    INTEGER :: i,j,n,temp
    
    PRINT*, "Enter the Size: "
    READ*, n
    ALLOCATE(NUMBERS(n))
    
    PRINT '(a4 i3 a25)', "Type",n, " numbers, Each on a line: "
    DO i=1,n
        READ*,NUMBERS(i)
    END DO
    
    DO i=1,n-1
        DO j=i+1,n
            IF (NUMBERS(i) > NUMBERS(j)) THEN
                temp = NUMBERS(i)
                NUMBERS(i) = NUMBERS(j)
                NUMBERS(j) = temp
            END IF
        END DO
    END DO
    
    PRINT*, "The Sorted List is: "
    DO i=1,n
        PRINT*,NUMBERS(i)
    END DO
    
    DEALLOCATE(NUMBERS)
END PROGRAM Sort