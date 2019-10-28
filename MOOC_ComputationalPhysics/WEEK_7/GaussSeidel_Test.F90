Program Jacobi_GaussSeidel
    ! Program to solve [y''-5y'+10y = 10x] Using Jacobi and Gauss-Seidel Methods
    Use DE_SUBS
    Implicit None
    Interface
        Function yi(xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi) Result(R)
            Implicit None
            Real*8          :: xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi,h
            Real*8          :: R
        End Function yi
    End Interface

    Integer                 :: i,j,NoP
    Integer,Parameter       :: MAX_ITERS = 50000
    Real*8,Parameter        :: XL = 0.d0,XU = 0.9d0,YL = 0.d0,YU = 50.0d0, h = 0.01d0, CONVERGE = 0.000001d0
    Real*8,Allocatable      :: x(:),y(:)
    Real*8                  :: coeff_yiminus1,coeff_yiplus1,coeff_xi
    Real*8                  :: Slope    ! Slope to Init Y Grid points to a Line Function
    Logical                 :: isFinished  = .False.
    
    NoP = (XU - XL)/h
    Allocate(x(NoP)) ; Allocate(y(NoP)) ; 

    coeff_yiminus1 = (1.0d0 + 2.50d0*h)/(2.0d0 - 10.0d0*h*h)
    coeff_yiplus1  = (1.0d0 - 2.50d0*h)/(2.0d0 - 10.0d0*h*h)
    coeff_xi       = (-10.0d0*h*h)/(2.0d0 - 10.0d0*h*h)

    Write(*,*) "Beginning Gauss-Seidel Method..."
    Write(*,*) "XL    = ",XL,"XU    = ",XU
    Write(*,*) "Y(XL) = ",YL,"Y(XU) = ",YU
    Write(*,*) "Number of Grid Points : ",NoP

    Slope = (YU - YL)/(XU - XL)
    ! Init Grid Points to a Line Function 
    Do i = 1,NoP,1
        x(i) = h*(i-1) ; y(i) = Slope*x(i) 
    End Do 

    Do j = 1,MAX_ITERS,1
        ! yi is the Finite Difference Function Defined after the Program
        Call GaussSeidelIterator(NoP,x,y,CONVERGE,yi,coeff_yiminus1,coeff_yiplus1,coeff_xi,isFinished)
        If (isFinished) Exit 
    End Do
    Write(*,*) "Gauss-Seidel Finte Difference Method finished in ",j," Iterations"

    Open(UNIT = 11,file = "GaussSeidelData.csv",Status = 'Unknown')
    Do j = 1,NoP,1
        Write(11,*) x(j),",",y(j)
    End Do
    Close(11)
End Program Jacobi_GaussSeidel

Function yi(xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi) Result(R)
    Implicit None
    Real*8          :: xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi
    Real*8          :: R
    R = coeff_yiminus1 * yiminus1 + coeff_yiplus1 * yiplus1 + coeff_xi * xi
End Function yi