PROGRAM AllocateArray
    implicit none
    real,dimension(:),allocatable :: vector
    integer :: n,i
    print*, "Enter the Vector Dimension: "
    read*, n
    allocate(vector(n))
    do i = 1,n
        vector(i) = acos(-1.0)**i
    end do
    print*,vector
    deallocate(vector)
END PROGRAM AllocateArray