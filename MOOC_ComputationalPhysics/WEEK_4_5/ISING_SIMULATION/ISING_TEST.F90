! Program to simulate a 3D Ising Lattice
Program Ising_Test
    Use ISING 
    Implicit None
    Integer             :: i,j,k,time                           ! Lattice and Iteration Loop indices 
    Integer             :: L,Iter,N                             ! Lattice Length,Number of Iterations and Size of the Lattice
    Real*8              :: r                                    ! Random Number 
    Real*8 ,Allocatable :: Lattice(:,:,:)                       ! 3D Ising Lattice to be simulated
    Real*8               :: T , J_Ising = 1.0d0               ! Temperature and Units
    Real*8              :: E,M,Ei,Ef,dE,U,M_Per_Spin,E_Per_Spin ! Energy and Magnetic Moment
    Real*8              :: Energy_Avg,Moment_Avg
    Real*8              :: B = 0.0d0

    E = 0.0d0
    M = 0.0d0
    Energy_Avg = 0.0d0 ; Moment_Avg = 0.0d0

    Open(UNIT = 10,FILE = "Energy_Per_Spin.dat")
    Open(Unit = 11,FILE = "MagneticMoment_Per_Spin.dat")

    Write(*,*) "Enter the Temperature: "
    Read(*,*) T 

    Write(*,*) "Enter the Lattice Dimension(1D): "
    Read(*,*) L

    Write(*,*) "Enter the Number of Iterations: "
    Read(*,*) Iter 

    Allocate(Lattice(L,L,L))
    N = L * L * L

    !Lattice = 1.0d0
    Call RandomizeLattice(Lattice,L)

    ! Calculation of Initial Magnetization and Energy
    Call CalcMoment_Energy(Lattice,L,J_Ising,B,M,E)
    M_Per_Spin = M / Real(N) ; E_Per_Spin = E / Real(N)
    Write(*,*) "Total Initial Magnetization          : ",M
    Write(*,*) "Total Initial Energy                 : ",E
    Write(*,*) "Total Initial Magnetization per Spin : ",M_Per_Spin
    Write(*,*) "Total Initial Energy per Spin        : ",E_Per_Spin

    ! Ising Iteration Begins
    Do time = 1,Iter
        ! Choose a Random Lattice Site
        Call RANDOM_NUMBER(r); i = int(float(L)*r)+1
        Call RANDOM_NUMBER(r); j = int(float(L)*r)+1
        Call RANDOM_NUMBER(r); k = int(float(L)*r)+1

        Call EnergyOfSite(Lattice,L,i,j,k,J_Ising,B,Ei)
        Call FlipSite(Lattice,L,i,j,k)
        Call EnergyOfSite(Lattice,L,i,j,k,J_Ising,B,Ef)
        dE = Ef - Ei

        If (dE < 0) Then
            ! Spin Flip Acceppted and Instantaneous Energy and Magnetic Moment Calculated
            E = E + dE
            M = M + 2.0d0*Lattice(i,j,k)
        Else
            ! Compute Boltzman Factor
            U = Exp(-dE/T)
            Call RANDOM_NUMBER(r)
            If (r < U) Then
                ! Spin Flip Acceppted(with a probability of Boltzman Factor) and Instantaneous Energy and Magnetic Moment Calculated
                E = E + dE
                M = M + 2.0d0*Lattice(i,j,k)
            Else
                ! Spin Flip Not Accepted. Revert it back
                Call FlipSite(Lattice,L,i,j,k)
            End If 
        End If
        Energy_Avg = Energy_Avg + E ; Moment_Avg = Moment_Avg + M
        M_Per_Spin = M / Real(N) ; E_Per_Spin = E / Real(N)
        Write(10,*) time,E_Per_Spin
        Write(11,*) time,M_Per_Spin
    End Do       
    Write(*,*) "Magnetization Average(Per Spin)   = ",M/Real(N) 
    Write(*,*) "Energy Average(Per Spin)          = ",E/Real(N)
    Deallocate(Lattice)
    Close(10) ;  Close(11)  ! Close Data Files
End Program Ising_Test

