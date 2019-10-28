! Program to simulate a 3D Ising Lattice
Program Ising
    Implicit None
    Integer             :: i,j,k,time                           ! Lattice and Iteration Loop indices 
    Integer             :: L,Iter,N                             ! Lattice Length,Number of Iterations and Size of the Lattice
    Real                :: r                                    ! Random Number 
    Real,Allocatable    :: Lattice(:,:,:)                       ! 3D Ising Lattice to be simulated
    Real                :: T = 2.0, J_Ising = 1.0               ! Units
    Real                :: E,M,Ei,Ef,dE,U,M_Per_Spin,E_Per_Spin ! Energy and Magnetic Moment
    Integer             :: Trials 
    Real                :: Energy_Avg,Moment_Avg
    Integer,Parameter   :: Trials_Num = 20000


    E = 0.0
    M = 0.0
    Energy_Avg = 0.0 ; Moment_Avg = 0.0

    Open(UNIT = 10,FILE = "Energy_Per_Spin.dat")
    Open(Unit = 11,FILE = "MagneticMoment_Per_Spin.dat")

    ! Write(*,*) "Enter the Lattice Dimension(1D): "
    ! Read(*,*) L

    ! Write(*,*) "Enter the Number of Iterations: "
    ! Read(*,*) Iter 

    L = 5 ; Iter = 50000

    N = L * L * L
    Allocate(Lattice(L,L,L))

    !Lattice = 1.0
    Call RandomizeLattice(Lattice,L)

    ! Calculation of Initial Magnetization and Energy
    Call CalcMoment_Energy(Lattice,L,J_Ising,M,E)
    M_Per_Spin = M / Real(N) ; E_Per_Spin = E / Real(N)
    Write(*,*) "Total Initial Magnetization          : ",M
    Write(*,*) "Total Initial Energy                 : ",E
    Write(*,*) "Total Initial Magnetization per Spin : ",M_Per_Spin
    Write(*,*) "Total Initial Energy per Spin        : ",E_Per_Spin
    Do Trials = 1,Trials_Num
        Call RandomizeLattice(Lattice,L)
        ! Calculation of Initial Magnetization and Energy
        Call CalcMoment_Energy(Lattice,L,J_Ising,M,E)
        ! Ising Iteration Begins
        Do time = 1,Iter
            ! Choose a Random Lattice Site
            Call RANDOM_NUMBER(r); i = int(float(L)*r)+1
            Call RANDOM_NUMBER(r); j = int(float(L)*r)+1
            Call RANDOM_NUMBER(r); k = int(float(L)*r)+1
            Call EnergyOfSite(Lattice,L,i,j,k,J_Ising,Ei)
            !Write(*,*) "Ei = ",Ei
            Call FlipSite(Lattice,L,i,j,k)
            Call EnergyOfSite(Lattice,L,i,j,k,J_Ising,Ef)
            !Write(*,*) "Ef = ",Ef
            dE = Ef - Ei
            If (dE < 0) Then
                ! Spin Flip Acceppted and Instantaneous Energy and Magnetic Moment Calculated
                E = E + dE
                M = M + 2.0*Lattice(i,j,k)
            Else
                ! Compute Boltzman Factor
                U = Exp(-dE/T)
                Call RANDOM_NUMBER(r)
                If (r < U) Then
                    ! Spin Flip Acceppted(with a probability of Boltzman Factor) and Instantaneous Energy and Magnetic Moment Calculated
                    E = E + dE
                    M = M + 2.0*Lattice(i,j,k)
                Else
                    ! Spin Flip Not Accepted. Revert it back
                    Call FlipSite(Lattice,L,i,j,k)
                End If 
            End If
            !Write(10,*) time,E
            !Write(11,*) time,M
        End Do       
        Energy_Avg = Energy_Avg + E
        Moment_Avg = Moment_Avg + M
        M_Per_Spin = M / Real(N) ; E_Per_Spin = E / Real(N)
        Write(10,*) Trials,E_Per_Spin
        Write(11,*) Trials,M_Per_Spin
    End Do
    Energy_Avg = E / Real(Trials_Num) ; Moment_Avg = M / Real(Trials_Num)
    Write(*,*) "Magnetization Average   = ",M 
    Write(*,*) "Energy Average          = ",E
    Deallocate(Lattice)
    Close(10) ;  Close(11)  ! Close Data Files
End Program Ising

Subroutine CalcMoment_Energy(Latt,L,J_Ising,M,E)
    Implicit None
    Integer     ::  L
    Real        ::  Latt(L,L,L),E,M,J_Ising,EOS ! Energy of a site
    Integer     ::  i,j,k ! Do Loop Indices
    E  = 0.0
    M  = 0.0
    ! Calculation of Magnetization and Energy
    Do k=1,L
        Do j=1,L
            Do i=1,L

                ! Summing for Total Magnetization
                M = M + Latt(i,j,k)       
                ! Summing for Total Energy CAUTION: Every Interaction Energy is counted twice
                Call EnergyOfSite(Latt,L,i,j,k,J_Ising,EOS)
                E = E + EOS
            End Do           
        End Do  
    End Do
    E = E * 0.50 ! Multipied by 0.5 as every interaction energy is counted twice
End Subroutine CalcMoment_Energy

Subroutine RandomizeLattice(Latt,L)
    Implicit None
    Integer         :: L,i,j,k 
    Real            :: Latt(L,L,L)
    Real            :: r 
    ! Initialize the Lattice with random spin values
    Do k=1,L
        Do j=1,L
            Do i=1,L
                Call RANDOM_NUMBER(r)
                if (r > 0.50) Then
                    Latt(i,j,k) =  1.0
                Else 
                    Latt(i,j,k) = -1.0 
                End If
            End Do 
        End Do  
    End Do
End Subroutine RandomizeLattice

Subroutine EnergyOfSite(Latt,L,i,j,k,J_Ising,E_O_S) 
    Implicit None
    Integer     :: L,i,j,k 
    Real        :: Latt(L,L,L),E_O_S,J_Ising
    Integer     :: X_Plus,X_Minus,Y_Plus,Y_Minus,Z_Plus,Z_Minus !  Six dirctions
     ! Neighbours
    X_Plus = i+1; X_Minus = i-1
    Y_Plus = j+1; Y_Minus = j-1
    Z_Plus = k+1; Z_Minus = k-1
    ! Boundary Conditions on a 3D Torus
    If (X_PLus > L) X_Plus = 1; If (X_Minus < 1) X_Minus = L
    If (Y_PLus > L) Y_Plus = 1; If (Y_Minus < 1) Y_Minus = L
    If (Z_PLus > L) Z_Plus = 1; If (Z_Minus < 1) Z_Minus = L
    ! Energy of the site calculated as:
    E_O_S = - J_Ising*Latt(i,j,k)*(     Latt(X_Plus,j,k)+Latt(X_Minus,j,k)  &
                                      + Latt(i,Y_Plus,k)+Latt(i,Y_Minus,k)  &
                                      + Latt(i,j,Z_Plus)+Latt(i,j,Z_Minus))
End Subroutine EnergyOfSite 

Subroutine FlipSite(Latt,L,i,j,k)
    Implicit None
    Integer         :: L,i,j,k 
    Real            :: Latt(L,L,L)
    Latt(i,j,k) = -1.0 * Latt(i,j,k)
End Subroutine FlipSite
