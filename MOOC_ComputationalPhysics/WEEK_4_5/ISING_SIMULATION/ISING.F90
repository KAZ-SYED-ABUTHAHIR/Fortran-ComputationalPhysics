Module ISING
    Implicit None
    Contains 
    
    Subroutine Equilibrate(Latt,L,J_Ising,T,B,M,E,N)
        Implicit None
        Integer     ::  L,N ! Lattice Size and Number of Equilibration Steps
        Real*8      ::  Latt(L,L,L),E,M,B,J_Ising,T
        Integer     ::  i ! Do Loop Index
        Do i = 1,N
            Call MCS(Latt,L,J_Ising,T,B,M,E)
        End Do
    End Subroutine Equilibrate

    Subroutine MCS(Latt,L,J_Ising,T,B,M,E)
        Implicit None
        Integer     ::  L
        Real*8      ::  Latt(L,L,L),E,M,B,J_Ising,T,U 
        Real*8      ::  r,Ei,Ef,dE 
        Integer     ::  i,j,k ! Site Indices
        ! Choose a Random Lattice Site
        Call RANDOM_NUMBER(r); i = int(float(L)*r)+1
        Call RANDOM_NUMBER(r); j = int(float(L)*r)+1
        Call RANDOM_NUMBER(r); k = int(float(L)*r)+1

        Call EnergyOfSite(Latt,L,i,j,k,J_Ising,B,Ei)
        Call FlipSite(Latt,L,i,j,k)
        Call EnergyOfSite(Latt,L,i,j,k,J_Ising,B,Ef)
        dE = Ef - Ei

        If (dE < 0) Then
            ! Spin Flip Acceppted and Instantaneous Energy and Magnetic Moment Calculated
            E = E + dE
            M = M + 2.0d0*Latt(i,j,k)
        Else
            ! Compute Boltzman Factor
            U = Exp(-dE/T)
            Call RANDOM_NUMBER(r)
            If (r < U) Then
                ! Spin Flip Acceppted(with a probability of Boltzman Factor) and Instantaneous Energy and Magnetic Moment Calculated
                E = E + dE
                M = M + 2.0d0*Latt(i,j,k)
            Else
                ! Spin Flip Not Accepted. Revert it back
                Call FlipSite(Latt,L,i,j,k)
            End If 
        End If
    End Subroutine MCS

    Subroutine CalcMoment_Energy(Latt,L,J_Ising,B,M,E)
        Implicit None
        Integer     ::  L
        Real*8      ::  Latt(L,L,L),E,M,B,J_Ising,EOS ! Energy of a site
        Integer     ::  i,j,k ! Do Loop Indices
        E  = 0.0d0
        M  = 0.0d0
        ! Calculation of Magnetization and Energy
        Do k=1,L
            Do j=1,L
                Do i=1,L
    
                    ! Summing for Total Magnetization
                    M = M + Latt(i,j,k)       
                    ! Summing for Total Energy CAUTION: Every Interaction Energy is counted twice
                    Call EnergyOfSite(Latt,L,i,j,k,J_Ising,B,EOS)
                    E = E + EOS
                End Do           
            End Do  
        End Do
        E = E * 0.5d0 ! Multipied by 0.5 as every interaction energy is counted twice
    End Subroutine CalcMoment_Energy
    
    Subroutine RandomizeLattice(Latt,L)
        Implicit None
        Integer         :: L,i,j,k 
        Real*8          :: Latt(L,L,L)
        Real*8          :: r 
        ! Initialize the Lattice with random spin values
        Do k=1,L
            Do j=1,L
                Do i=1,L
                    Call RANDOM_NUMBER(r)
                    if (r > 0.50) Then
                        Latt(i,j,k) =  1.0d0
                    Else 
                        Latt(i,j,k) = -1.0d0 
                    End If
                End Do 
            End Do  
        End Do
    End Subroutine RandomizeLattice
    
    Subroutine EnergyOfSite(Latt,L,i,j,k,J_Ising,B,E_O_S) 
        Implicit None
        Integer     :: L,i,j,k 
        Real*8        :: Latt(L,L,L),E_O_S,J_Ising,B
        Integer     :: X_Plus,X_Minus,Y_Plus,Y_Minus,Z_Plus,Z_Minus !  Six dirctions
        B = 0 ! CAUTION : Dummy Assignment to Avoid Warning
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
        Real*8          :: Latt(L,L,L)
        Latt(i,j,k) = -1.0d0 * Latt(i,j,k)
    End Subroutine FlipSite  
End Module ISING