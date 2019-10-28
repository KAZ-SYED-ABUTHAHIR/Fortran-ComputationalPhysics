PROGRAM SortSub
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
    
    CALL Sort(NUMBERS,n);
    
    PRINT*, "The Sorted List is: "
    DO i=1,n
        PRINT*,NUMBERS(i)
    END DO
    
    DEALLOCATE(NUMBERS) !Clear the memory allocated for NUMBERS array
END PROGRAM SortSub



SUBROUTINE Sort(array,length)
    IMPLICIT NONE
    INTEGER :: length, array(length),i,j,temp
    
    DO i=1,length-1
        DO j=i+1,length
            IF (array(i) > array(j)) THEN
                temp = array(i)
                array(i) = array(j)
                array(j) = temp
            END IF
        END DO
    END DO

END SUBROUTINE Sort