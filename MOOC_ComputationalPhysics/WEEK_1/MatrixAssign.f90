PROGRAM MatrixAssign
    integer :: i,j
    real, dimension(2,2) :: A
    do i = 1,2
        do j=1,2
            A(i,j) = i*real(j) + 2.0
        end do
    end do
    print*, A
END PROGRAM MatrixAssign