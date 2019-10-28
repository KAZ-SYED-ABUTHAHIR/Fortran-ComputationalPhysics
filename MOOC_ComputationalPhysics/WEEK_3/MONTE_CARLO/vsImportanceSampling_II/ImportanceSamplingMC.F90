program ImportanceSampling
    use Random_Seed_Mod
    implicit none
    integer     :: i,j,n
    real*8      :: r,x
    real*8      :: g_of_x,g_of_x_crude 
    real*8      :: integral,weight,func,integral_crude
    real*8      :: sumsqr,sumsqr_crude,variance,variance_crude

    sumsqr = 0.0d0 ; sumsqr_crude = 0.0d0 ; variance_crude = 0.0d0 ; variance = 0.0d0
    integral = 0.0d0 ; weight = 0.0d0 ; func = 0.0d0 ; integral_crude = 0.0d0

    n = 1
    do j = 1,8
        n = n * 10
        call setup_random()
        do i = 1,n
            call RANDOM_NUMBER(r)
            g_of_x_crude =  r**(-1.0d0/3.0d0) + (r/10.0d0)
            integral_crude = integral_crude + g_of_x_crude
            sumsqr_crude = sumsqr_crude + g_of_x_crude * g_of_x_crude
            x = r**(1.50d0)
            weight = (2.0d0/3.0d0)*x**(-1.0d0/3.0d0)
            func   = x**(-1.0d0/3.0d0) + (x/10.0d0)
            g_of_x  = func / weight
            integral = integral + g_of_x
            sumsqr = sumsqr + g_of_x * g_of_x
        end do
        integral = integral / real(n)
        integral_crude = integral_crude / real(n)
        sumsqr = sumsqr / real(n) 
        variance = SQRT(ABS((sumsqr - (integral*integral))/real(n-1)))
        sumsqr_crude = sumsqr_crude / real(n)
        variance_crude = SQRT(ABS((sumsqr_crude - (integral_crude*integral_crude))/real(n-1)))
        !write(*,*) "n=",n,"Integral Crude=",integral_crude,"Integral=",integral, &
        !           "Variance Crude=",variance_crude,"Variance=",variance
        write(*,*) n,integral_crude,integral,variance_crude,variance,(variance_crude/variance)
    end do
end program ImportanceSampling