! Program to compute Ground State of Symmetric Finite Potential Well of Depth 'V0'
! Width 'b' using expansion on a plane wave basis set(Period 'a') and diagonalization 
! Using atomic units :: hbar^2/2m -----> 1
! Requires LAPAK dsyev routine

Program FiniteWell
    Implicit None
    Integer,Parameter       :: dp = SELECTED_REAL_KIND(14,200)
    Real(dp),Parameter      :: PI = 3.14159265358979_dp
    Integer                 :: N,NPW
    Real(dp)                :: a,b,V0 ! Potential Well Parameters
    Real(dp),Allocatable    :: K_n(:),E(:),H(:,:),Work(:)
    Real(dp)                :: x,dx,Norm,Probability
    Complex(dp)             :: f
    Integer                 :: i,j,NR,lWork,Info

    ! Input Data

    ! V = -V0    for |x| < b/2   and  V = 0    for |x| > b/2
    Write(*,"('Enter Parameters for Potential Well V0(Depth) and b(Breadth):', $)")
    Read(*,*) V0,b
    If(V0 <= 0.0_dp .OR. b <= 0.0_dp) STOP 'Wrong Input Parameters !'
    Write(*,"('  V0,b = ',2f10.4)") V0,b 

    ! Plane waves between -a/2 and +a/2. K_i = -i*(2*PI/a) to -i*(2*PI/a) i = 0,1,2,3,...,n
    100 Write(*,"('Enter Parameters for Plane Waves a(Period) and n(Number of Basis Waves):', $)")
    Read(*,*) a,n
    If(a < b .OR. n <= 0 ) STOP "Wrong Input Parameters !"
    Write(*,"('a,n = ',f8.4,i6)") a,n
    NPW = 2*n + 1 ! Number of Plane Waves in the Basis

    Allocate(K_n(NPW),E(NPW),H(NPW,NPW),Work(3*NPW))

    ! Assign values of K_n = 0,+1,-1,+2,-2,...
    K_n(1) = 0.0_dp
    Do i = 2,NPW-1,2
        K_n(i)   = +(2.0_dp*PI/a)*(i/2)
        K_n(i+1) = -(2.0_dp*PI/a)*(i/2)
    End Do

    ! Cleaning up Hamitonian Matrix
    H(:,:) = 0.0_dp

    ! Assigning Matrix Elements of the Hamiltonian Based on Expanded Plane Wave Basis
    Do i = 1,NPW,1
        Do j = 1,NPW,1
            If(i == j) Then 
                H(i,j) = K_n(i)*K_n(i) - V0*(b/a)
            Else
                H(i,j) = -(V0/a) * sin(b*(K_n(j)-K_n(i))/2.0_dp)/((K_n(j)-K_n(i))/2.0_dp)
            End If
            !Print '(2i4,f12.6)',i,j,H(i,j)
        End Do 
    End Do
    
    ! Solution Coefficients are stored into H(j,i)
    ! j = Basis Function Index, i = Eigenvaue Index
    lWork = 3*NPW
    !   The 'Leading Dimension of Array' H is  its' First Dimension
    !   For Dynamically Allocated Matrices see the 'Allocate' command

    ! IMPORTANT Call to LAPAK+BLAS Subroutine DSyEV - Douple Precsion Symmentric Eigen Value
    Call DSyEV('V','U',NPW,H,NPW,E,Work,lWork,Info) 

    If( Info .NE. 0) STOP 'H - Matrix Diagonalization Failed !'

    ! Write First 5 Eigen Values 
    Write(*,"('Lowest 5 Eigen Values: ',3f12.6)") (E(i),i=1,5)

    ! Write to Output File GroundStatePW.csv the Ground Eigen State
    Open(Unit = 77,File = 'GroundStatePW.csv',Status = 'Unknown',Form = 'Formatted')
    dx = 0.01_dp
    NR = NINT((a/2.0_dp)/dx)
    Norm = 0.0_dp
    Do i = -NR,NR,1
        x = i*dx
        f = 0.d0
        Do j = 1,NPW,1
            f  = f + H(j,1)*Exp((0.0,1.0)*K_n(j)*x)/SQRT(a)
        End Do
        Probability = f * CONJG(f)
        Norm = Norm + Probability * dx
        Write(77,"(f12.6,3f10.6)") x,Probability,f 
    End Do

    ! Verify Normalization if desired.
    Write(*,"('Norm = ',f12.6)") Norm
    Close(77)
    Deallocate(H,Work,E,K_n)
    GoTo 100

End Program FiniteWell
