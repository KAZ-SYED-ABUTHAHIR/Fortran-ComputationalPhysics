! Program to simulate a 3D Ising Lattice
Program Ising_Test
    Use ISING 
    Implicit None
    Integer             :: time                      ! Lattice and Iteration Loop indices 
    Integer             :: L,Iter,N                  ! Lattice Length,Number of Iterations and Size of the Lattice 
    Real*8 ,Allocatable :: Lattice(:,:,:)            ! 3D Ising Lattice to be simulated
    Real*8              :: T , J_Ising = 1.0d0       ! Temperature and Units
    Real*8              :: E,M,M_Per_Spin,E_Per_Spin ! Energy and Magnetic Moment
    Real*8              :: Energy_Avg,Moment_Avg,E_SqrAvg,M_SqrAvg
    Real*8              :: E_Variance,M_Variance
    Real*8              :: B = 0.0d0
    Integer,Parameter   :: EQUILIB_STEPS = 1000000
    Integer             :: T_Num_Steps,T_Index
    Real*8,Parameter    :: T_INIT = 8.0d0,T_FINISH = 1.0d0,T_STEP = -0.01d0
    Real*8              :: Cv,Chi                  ! Specific Heat Capacity and Magnetic Susceptibility
    Integer,Parameter   :: SKIP = 100              ! CAUTION: SKIP is always 1.This is Stupid Syed.OK Probem Identified. Think Why? Used Iter instead of time in the Mod expression
    Real*8              :: MSqrSqrAvg,BindersCumulant
    E = 0.0d0
    M = 0.0d0
    Energy_Avg = 0.0d0 ; Moment_Avg = 0.0d0
    M_Per_Spin = 0.0d0 ; E_Per_Spin = 0.0d0
    E_SqrAvg   = 0.0d0 ; M_SqrAvg   = 0.0d0  
    MSqrSqrAvg = 0.0d0 ; BindersCumulant = 0.0d0


    T_Num_Steps = Int(Abs(T_INIT-T_FINISH) / Abs(T_STEP))

    Open(UNIT = 10,FILE = "Phase_Varibles_L10_FineGrain_Skipped100.dat")
    L = 10
 

    Iter  = 10000000

    Allocate(Lattice(L,L,L))
    N = L * L * L

    !Lattice = 1.0d0
    Call RandomizeLattice(Lattice,L)

    ! Calculation of Initial Magnetization and Energy
    Call CalcMoment_Energy(Lattice,L,J_Ising,B,M,E)
    M_Per_Spin = M / Real(N) ; E_Per_Spin = E / Real(N)
    Write(*,*) "Lattice Size L                       : ",L
    Write(*,*) "Total Initial Magnetization          : ",M
    Write(*,*) "Total Initial Energy                 : ",E
    Write(*,*) "Total Initial Magnetization per Spin : ",M_Per_Spin
    Write(*,*) "Total Initial Energy per Spin        : ",E_Per_Spin

    Do T_Index = 0,T_Num_Steps,1
        E = 0.0d0
        M = 0.0d0 
        Energy_Avg = 0.0d0 ; Moment_Avg = 0.0d0
        M_Per_Spin = 0.0d0 ; E_Per_Spin = 0.0d0
        E_SqrAvg   = 0.0d0 ; M_SqrAvg   = 0.0d0  
        MSqrSqrAvg = 0.0d0 ; BindersCumulant = 0.0d0
        T = T_INIT + T_Index * T_STEP 

        Call RandomizeLattice(Lattice,L)
        Call Equilibrate(Lattice,L,J_Ising,T,B,M,E,EQUILIB_STEPS)

        ! Ising Iteration Begins , MCS - Mante Carlo Step
        Do time = 1,Iter,1
            Call MCS(Lattice,L,J_Ising,T,B,M,E)
            If (Mod(time,SKIP) .EQ. 0) Then 
                !Write(*,*) Energy_Avg,Moment_Avg
                Energy_Avg = Energy_Avg + E ; Moment_Avg = Moment_Avg + Abs(M)
                E_SqrAvg = E_SqrAvg + E * E ; M_SqrAvg = M_SqrAvg + (M * M)
                MSqrSqrAvg = MSqrSqrAvg +  (M * M * M * M)
            End If 
        End Do   

        Energy_Avg = Energy_Avg / Real(Iter/SKIP) ; Moment_Avg = Moment_Avg / Real(Iter/SKIP)
        E_SqrAvg   = E_SqrAvg   / Real(Iter/SKIP) ; M_SqrAvg   = M_SqrAvg   / Real(Iter/SKIP)
        MSqrSqrAvg = MSqrSqrAvg / Real(Iter/SKIP)

        E_Per_Spin = Energy_Avg / Real(N)    ; M_Per_Spin = Moment_Avg / Real(N)
        E_Variance = E_SqrAvg - (Energy_Avg * Energy_Avg)
        M_Variance = M_SqrAvg - (Moment_Avg * Moment_Avg)
        Cv = E_Variance / (T*T) ; Chi = M_Variance / T 
        BindersCumulant = 1.0d0 - (MSqrSqrAvg / (3.0d0 * M_SqrAvg * M_SqrAvg))

        Write(10,*) T,Energy_Avg,Moment_Avg,E_Per_Spin,M_Per_Spin,Cv,Chi,BindersCumulant 
        Write(*,*) T,Energy_Avg,Moment_Avg,E_Per_Spin,M_Per_Spin,Cv,Chi,BindersCumulant 
    End Do

    Deallocate(Lattice)
    Close(10) ;  !Close(11)  ! Close Data Files
End Program Ising_Test

