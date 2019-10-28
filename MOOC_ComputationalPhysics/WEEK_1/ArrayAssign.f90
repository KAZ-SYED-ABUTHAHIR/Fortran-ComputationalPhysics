PROGRAM ArrayAssign
    implicit none
    integer :: i
    real, dimension(5) :: v = (/ 1.1,1.2,1.3,1.4,1.5 /)
    print*, v
    do i = 1,5
        v(i) = 1+i*0.1
    end do
    print*, v
END PROGRAM ArrayAssign