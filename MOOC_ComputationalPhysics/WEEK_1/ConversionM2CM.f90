program conversion
    integer :: in_m,out_cm
    print*, "Enter Length in Meter:"
    read '(i5)', in_m
    out_cm = in_m*100
    print *, "Length in CM:"
    print '(i5)',out_cm
end program conversion