! Program to Implement Molecular Dynamics For Lennard-Jones Interacting  Particles
! Code Not Optimized :-(

!___________________________________________________PARAMETERS MODULE____________________________________________________________ 
Module SimulationParameters
    Implicit None

    Integer,Parameter       ::  NumParticles = 2400,NumEquilibrations = 50000 , NumIterations = 250000,NumCalcAvg = 10
    Integer,Parameter       ::  NumRefreshNeighbours = 80
    Real*8,Parameter        ::  LX = 20.0d0,LY = 20.0d0,LZ = 20.0d0
    Real*8,Parameter        ::  Mass = 1.0d0,KbMulT = 0.25d0,Sigma = 1.0d0,Epsil = 4.0d0,delT = 0.0025d0
    Real*8,Parameter        ::  Rc = 2.50d0,Rs = 2.0d0+Rc ! Rs - Skin radius for neighbours table 
    Real*8,Parameter        ::  FourPI = 8.0d0*asin(1.0d0),SkinSphereVolume = (FourPI/3.0d0) * Rs**3
    Integer,Parameter       ::  MaxNumOfNeighbours = 8*Ceiling((NumParticles/(LX*LY*LZ)) * SkinSphereVolume)  ! CAUTION: Set the Multiplier to a high value with increased particle numbers to avoid memory segmaentation error.
    Real*8,Parameter        ::  dr = 0.050d0*Sigma ! dr - g(r) bin size
    Real*8,Parameter        ::  GOfRMax = 7.0d0, VolumeOfGOfRSphere = (FourPI * GOfRMax ** 3)/3.0d0
    Integer,Parameter       ::  NumBins = CEILING(GOfRMax/dr),GOfRCalcInterval = 100  
    Real*8,Parameter        ::  minR = 0.70d0 ! Start value for calculating g(r)
   
    Real*8,Parameter        ::  SigmaPow6 = Sigma ** 6,SigmaPow12 = Sigma ** 12
    Real*8,Parameter        ::  Fc = Epsil * ( (12.0d0*SigmaPow12/(Rc**13)) - (6.0d0*SigmaPow6/(Rc**7)) )
    Real*8,Parameter        ::  UFc = Fc*Rc + Epsil * ( ((Sigma/Rc)**12)-((Sigma/Rc)**6) )  ! Modiifed Potential

    Logical,Parameter       ::  Thermostat = .TRUE.  ! Make this Equal to .TRUE. or .FALSE to Switch ON or Off the Thermostat

    Integer                 ::  GOFRCount 
    Real*8                  ::  GOfR(NumBins),AverageDensity,SumAverageDensity

    Real*8                  ::  LXBy2 = LX/2.0d0,LYBy2 = LY/2.0d0,LZBy2 = LZ/2.0d0
    Real*8                  ::  Positions(3*NumParticles),Velocities(3*NumParticles)
    Real*8                  ::  Forces(3*NumParticles),NewForces(3*NumParticles)
    Integer                 ::  NumOfNeighbours(NumParticles),NeighboursTable(NumParticles,MaxNumOfNeighbours) 
    
    Real*8                  ::  AvrVelX,AvrVelY,AvrVelZ
    Real*8                  ::  x1,y1,z1,x2,y2,z2,dx,dy,dz,r,Lj,gX,gY,gZ,x,y,z
    Real*8                  ::  LjForce
    Real*8                  ::  NewPE,PE,KE
    Real*8                  ::  Invr,Ir2,Ir6
    Real*8                  ::  LLX,LLY,LLZ,TheoryKE,ScaleFactorForThermoStat
    Real*8                  ::  RndNum      ! Random Num Variable for Setting Random Veocities at the beginning
    Integer                 ::  i,j,k,l,m   ! Loop Indices that may come handy
    Integer                 ::  Time 
    Integer                 ::  AvgNumOfNeighbours,MostNumOfNeighbours ! Average Number of Neighbours within a Sphere of Radius  Rs  
End Module SimulationParameters

!___________________________________________________________________________________________________________________________
!                                                        MAIN PROGRAM
!___________________________________________________________________________________________________________________________

Program LennardJones
    Use SimulationParameters
    Implicit None

    Real*8              :: Density , Rad 
    Forces = 0.0d0 ; NumOfNeighbours = 0 ; AvgNumOfNeighbours = 0 ; GOFRCount = 0 ; GOfR = 0.0d0 
    SumAverageDensity = 0.0d0 


    Call InitPosRnd ! Initialize Particle Positions
    !Call InitPosLattice ! Initialize Particle Positions
    Open(UNIT=44,FILE = "InitialConfig_JMOL.xyz")   ! JMOl Movie File
    Write(44,*) NumParticles  ! NUM_OSC - Number of Oscillators in the Circular Chain
    Write(44,*) ""
    Do i = 1,NumParticles,1
        ! Writing Initial Positions to the file InitialConfig_JMOL.xyz
        Write(44,*) "N",Positions(3*i-2),Positions(3*i-1),Positions(3*i  )
    End Do 
    Close(44)

    Call InitVel        ! Initialize Particle Velocities 

    Write(*,*) "Maximum Number of Possible Neighbours: ",MaxNumOfNeighbours
    Forces = 0.0d0

    Call RefreshNeibours
    Call CalcForce

    Forces = NewForces

    Open(UNIT=33,FILE = "MovieLJ.xyz")   ! JMOl Movie File
    Open(UNIT=77,FILE = 'Neighbours.csv',STATUS = 'Unknown',FORM = 'Formatted')
    Open(UNIT=88,FILE = 'Dist.csv',STATUS = 'Unknown',FORM = 'Formatted')
    Open(UNIT=99,FILE = 'Energy.csv',STATUS = 'Unknown',FORM = 'Formatted')

    If(Thermostat) Then
        Write(*,*) "Thermostat: ON"
    Else
        Write(*,*) "Thermostat: OFF"
    End If  

    ! Equilibration Loop
    Write(*,*)  "Equilibrating..."
    Do Time = 1,NumEquilibrations,1
        If (Mod(Time,200) .EQ. 0) Then
            Write(*,*) Time 
        EndIf 
        Call UpdatePos 
        Forces = NewForces
        Call CalcForce
        Call UpdateVel        
    End Do

    ! Run Simulation
    Do Time = 1,NumIterations,1
        If (Mod(Time,100) .EQ. 0) Then
            Write(*,*) Time 
        EndIf 

        Call UpdatePos 
        Forces = NewForces
        Call CalcForce
        Call UpdateVel


        ! Time is the iteration index
        If (Mod(Time,800) == 0) Then ! Write the Position Data To The JMol Animation File Once In 800 Steps
            Write(*,*) "Writing Position to File MovieLJ.xyz ..."
            Write(33,*) NumParticles  ! NumParticles - Number of Particles in the LJ Box
            Write(33,*) ""
            Do i = 1,NumParticles,1
                ! Write the data to the file 
                Write(33,*) "N",Positions(3*i-2),Positions(3*i-1),Positions(3*i  )
            End Do 
        End If
        
        ! If (Mod(Time,NumCalcAvg) .EQ. 0) Call CalcTDAvgs

        If (Mod(Time,GOfRCalcInterval) .EQ. 0) Call CalcGofR

        Write(99,*) dFloat(Time),",",NewPE/dFloat(NumParticles),",",KE/dFloat(NumParticles), &
                                                        ",",(NewPE+KE)/dFloat(NumParticles)

        If(Mod(Time,NumRefreshNeighbours) .EQ. 0) Write(77,*) dFloat(Time),",",MostNumOfNeighbours,",",AvgNumOfNeighbours,&
                                                                           ",",AverageDensity
    End Do
    
    Open(UNIT=55,FILE = 'GOfR.csv',STATUS = 'Unknown',FORM = 'Formatted')

     Density =  SumAverageDensity / (dFloat(GOfRCount)* VolumeOfGOfRSphere)
     !Density = dFloat(NumParticles)/(LX*LY*LZ)
    GOfR = GOfR / (dFloat(GOfRCount)*Density)

    Do i = 1,NumBins,1
        Rad = dr*dFloat(i)
        GOfR(i) = GOfR(i) / (FourPI*Rad*Rad*dr)
        Write(55,*) dFloat(i)*dr,",",GOfR(i)
    End Do 

    Close(33);Close(55);Close(77);Close(88);Close(99)
End Program LennardJones

!___________________________________________________________________________________________________________________________
!                                                         SUBROUTINES
!____________________________________________________________________________________________________________________________
!________________________________________________________INIT POSITIONS____________________________________________________
Subroutine InitPosLattice
    Use SimulationParameters
    Implicit None
    
    ! Putting Particles in a Cubic Lattice
    Real*8                  ::  LatticeSpacing,xPos,yPos,zPos
    Integer                 ::  xNum,yNum,zNum,cX,cY,cZ,ParticleCounter

    Write(*,*) "Initializing Partricle Positions to Simple Cubic Lattice..."

    LatticeSpacing = (LX*LY*LZ)**(1.0d0/3.0d0) / CEILING(dFloat(NumParticles)**(1.0d0/3.0d0))
    
    xNum = FLOOR(LX/LatticeSpacing); yNum = FLOOR(LY/LatticeSpacing); zNum = FLOOR(LY/LatticeSpacing)
    
    xPos = LatticeSpacing  ; yPos = LatticeSpacing  ; zPos = LatticeSpacing  ; 

    ParticleCounter = 1

    Do cx = 1,xNum,1
        Do cy = 1,yNum,1
            Do cz = 1,zNum,1
                Positions(3*ParticleCounter-2) = dFloat(cX)*LatticeSpacing - (LatticeSpacing/2.0d0)
                Positions(3*ParticleCounter-1) = dFloat(cY)*LatticeSpacing - (LatticeSpacing/2.0d0)
                Positions(3*ParticleCounter  ) = dFloat(cZ)*LatticeSpacing - (LatticeSpacing/2.0d0)
                ParticleCounter =  ParticleCounter + 1
                If (ParticleCounter > NumParticles) EXIT            
            End Do  
            If (ParticleCounter > NumParticles) EXIT      
        End Do
        If (ParticleCounter > NumParticles) EXIT
    End Do
    Write(*,*) ParticleCounter-1,"Particles Positions Initialized"

    Open(UNIT=22,FILE='InitPositions.csv',STATUS='Unknown')

    Do i = 1,NumParticles,1
        Write(22,*) Positions(3*i-2),',',Positions(3*i-1),',',Positions(3*i  )
    End Do
    Close(22)
End Subroutine InitPosLattice 

!__________________________________________________Initialize Positions Randomly__________________________________________________
Subroutine InitPosRnd
    ! Putting Particles in a Cubic Lattice
    Use SimulationParameters
    Implicit None
    Integer                 :: ParticleCounter,IterationCount
    Real*8                  :: SigmaBy2
    Logical                 :: OverLapped = .TRUE.
    Real*8,Parameter        :: Clearance = Sigma - (Sigma/20.0d0)

    SigmaBy2 = Sigma/2.0d0 
    IterationCount = 0

    Write(*,*) "Initializing Partricle Positions Randomly..."

    ! Particle Positions are randomly initilized such that they don't overlap with their images with a tolerance of Sigma
    Do ParticleCounter = 1,NumParticles,1
        Call RANDOM_NUMBER(RndNum) ; Positions(3*ParticleCounter-2) = RndNum*(LX-Sigma) + SigmaBy2
        Call RANDOM_NUMBER(RndNum) ; Positions(3*ParticleCounter-1) = RndNum*(LY-Sigma) + SigmaBy2
        Call RANDOM_NUMBER(RndNum) ; Positions(3*ParticleCounter  ) = RndNum*(LX-Sigma) + SigmaBy2
    End Do

    ! Now comes the computationally costly part to check whether two particle postions overlap within a tolerance of Sigma
    ! and to reassign them with new random positions until no such pair does exist !!! Somehow windy logic and has to be 
    ! made clear in future
    Do While(OverLapped)
        IterationCount = IterationCount + 1
        Do i = 1,NumParticles-1,1
            Do j = i+1,NumParticles,1
                dx = Positions(3*j-2)-Positions(3*i-2);dy = Positions(3*j-1)-Positions(3*i-1);dz = Positions(3*j  )-Positions(3*i)
                r = dSqrt(dx*dx + dy*dy + dz*dz)
                if ( r .LT. Clearance ) Then
                    ! Reassign random values. Need to ressign only one of the everlapping positions.
                    Call RANDOM_NUMBER(RndNum) ; Positions(3*i-2) = RndNum*(LX-Sigma) + SigmaBy2
                    Call RANDOM_NUMBER(RndNum) ; Positions(3*i-1) = RndNum*(LY-Sigma) + SigmaBy2
                    Call RANDOM_NUMBER(RndNum) ; Positions(3*i  ) = RndNum*(LX-Sigma) + SigmaBy2
                    OverLapped = .FALSE.
                End If  
            End Do
        End Do
        if(.NOT.(OverLapped)) Then
            Overlapped = .TRUE.
        Else
            Overlapped = .FALSE.
        End If  
    End Do

   
    Write(*,*) ParticleCounter-1,"Particles Positions Initialized"

    Open(UNIT=22,FILE='InitPositionsRnd.csv',STATUS='Unknown')

    Do i = 1,NumParticles,1
        Write(22,*) Positions(3*i-2),',',Positions(3*i-1),',',Positions(3*i  )
    End Do
    Close(22)
End Subroutine InitPosRnd

!_____________________________________________________INIT VELECITIES________________________________________________________________
Subroutine InitVel
    ! Initializing Particle Velocities Randomly
    Use SimulationParameters
    Implicit None
    
    Real*8                  ::  VelConstant,AvgVx,AvgVy,AvgVz
    
    VelConstant = dSqrt(12.0d0)*KbMulT

    Write(*,*) "Initializing Partricle Velocities Randomly..."

    ! Initialize the Velocities randomly
    Do i = 1,NumParticles,1
        Call RANDOM_NUMBER(RndNum) ; Velocities(3*i-2) = VelConstant * (RndNum - 0.5d0)
        Call RANDOM_NUMBER(RndNum) ; Velocities(3*i-1) = VelConstant * (RndNum - 0.5d0)
        Call RANDOM_NUMBER(RndNum) ; Velocities(3*i  ) = VelConstant * (RndNum - 0.5d0)
    End Do

    AvgVx = 0.0d0; AvgVy = 0.0d0; AvgVz = 0.0d0

    Do i = 1,NumParticles,1
        AvgVx = AvgVx + Velocities(3*i-2); AvgVy = AvgVy + Velocities(3*i-1); AvgVz = AvgVz + Velocities(3*i  )
    End Do

    AvgVx = AvgVx/dFloat(NumParticles);AvgVy = AvgVy/dFloat(NumParticles);AvgVz = AvgVz/dFloat(NumParticles)

    Write(*,*) "Average Velocities Before CoM Stabilization: "
    Write(*,*) "Average VelX :",AvgVx,"Average VelY :",AvgVy,"Average VelZ :",AvgVz

    ! Circus to prevent the Centre of Mass from having a random initial jitter 
    Do i = 1,NumParticles,1
        Velocities(3*i-2) = AvgVx - Velocities(3*i-2)
        Velocities(3*i-1) = AvgVy - Velocities(3*i-1)
        Velocities(3*i  ) = AvgVz - Velocities(3*i  )
    End Do

    AvgVx = 0.0d0; AvgVy = 0.0d0; AvgVz = 0.0d0
    Do i = 1,NumParticles,1
        AvgVx = AvgVx + Velocities(3*i-2); AvgVy = AvgVy + Velocities(3*i-1); AvgVz = AvgVz + Velocities(3*i  )
    End Do
    Write(*,*) "Average Velocities After CoM Stabilization: "
    Write(*,*) "Average VelX :",AvgVx,"Average VelY :",AvgVy,"Average VelZ :",AvgVz

    Open(UNIT=11,FILE='InitVelocities.csv',STATUS='Unknown')

    Do i = 1,NumParticles,1
        Write(11,*) Velocities(3*i-2),',',Velocities(3*i-1),',',Velocities(3*i  )
    End Do
    Close(11)
    Write(*,*) NumParticles,"Particles Velocities Initialized"
End Subroutine InitVel

!____________________________________________________UPDATE POSITIONS________________________________________________________
Subroutine UpdatePos
    Use SimulationParameters
    Implicit None
   
    Real*8                  ::  DtSqrByTwo

    DtSqrByTwo = 0.50d0*delT*delT

    ! Intiliazing Positions
    Do i = 1,NumParticles,1
            Positions(3*i-2) = Positions(3*i-2) + Velocities(3*i-2)*delT + Forces(3*i-2)*(DtSqrByTwo/Mass)
            Positions(3*i-1) = Positions(3*i-1) + Velocities(3*i-1)*delT + Forces(3*i-1)*(DtSqrByTwo/Mass)
            Positions(3*i-0) = Positions(3*i  ) + Velocities(3*i  )*delT + Forces(3*i  )*(DtSqrByTwo/Mass)
            ! Periodic Boundary Conditions
            Positions(3*i-2) = Modulo(Positions(3*i-2),LX)
            Positions(3*i-1) = Modulo(Positions(3*i-1),LY)
            Positions(3*i  ) = Modulo(Positions(3*i  ),LZ)
    End Do
End Subroutine UpdatePos

!____________________________________________________UPDATE VELOCITIES_____________________________________________________________
Subroutine UpdateVel
    Use SimulationParameters
    Implicit None
    
    ! Updating Velocity
    KE = 0.0d0
    AvrVelX = 0.0d0 ; AvrVelY = 0.0d0 ; AvrVelZ = 0.0d0 ; 

    Do i = 1,NumParticles,1
        Velocities(3*i-2) = Velocities(3*i-2) + delT*0.50d0*(Forces(3*i-2)+NewForces(3*i-2))/Mass
        Velocities(3*i-1) = Velocities(3*i-1) + delT*0.50d0*(Forces(3*i-1)+NewForces(3*i-1))/Mass
        Velocities(3*i  ) = Velocities(3*i  ) + delT*0.50d0*(Forces(3*i  )+NewForces(3*i  ))/Mass
        

        AvrVelX = AvrVelX + Velocities(3*i-2)
        AvrVelY = AvrVelY + Velocities(3*i-1)
        AvrVelZ = AvrVelZ + Velocities(3*i  )

        KE = KE + Velocities(3*i-2)*Velocities(3*i-2) + Velocities(3*i-1)*Velocities(3*i-1) + Velocities(3*i  )*Velocities(3*i  )
    End Do 

    KE = 0.50d0*Mass*KE

    ! Thermostat Implementation
    If (Thermostat .AND. (Mod(Time,100) .EQ. 0)) Then      
            TheoryKE = 1.50d0*KbMulT*dFloat(NumParticles)
            ScaleFactorForThermoStat = dSqrt(TheoryKE/KE)
            Velocities = Velocities * ScaleFactorForThermoStat

            KE = 0.0d0
            Do i = 1,NumParticles,1
                KE = KE + Velocities(3*i-2)*Velocities(3*i-2) + Velocities(3*i-1)*Velocities(3*i-1) + &
                                                                Velocities(3*i  )*Velocities(3*i  )
            End Do 
            KE = 0.50d0*Mass*KE
    End If
End Subroutine UpdateVel

!__________________________________________________CALCULATE FORCE_____________________________________________________________
Subroutine CalcForce
    Use SimulationParameters
    Implicit None
    ! Calculation of Forces
    NewForces = 0.0d0 ; NewPE = 0.0d0
    
    If(Mod(Time,NumRefreshNeighbours) .EQ. 0) Call RefreshNeibours

    Do i = 1,NumParticles-1,1
        x1 = Positions(3*i-2) ; y1 = Positions(3*i-1) ; z1 = Positions(3*i  )
        Do j = 1,NumOfNeighbours(i),1
            k = NeighboursTable(i,j) ! Particle k is j th Neighbour of Particle i
            x2 = Positions(3*k-2) ; y2 = Positions(3*k-1) ; z2 = Positions(3*k  )

            x = x1-x2 ; y = y1-y2 ; z = z1-z2 ; 

                ! Minimum Image Convention
                If (Abs(x) .GE. LXBy2) x = (LX-Abs(x))*(-1.0d0*x/Abs(x))
                If (Abs(y) .GE. LYBy2) y = (LY-Abs(y))*(-1.0d0*y/Abs(y))
                If (Abs(z) .GE. LYBy2) z = (LZ-Abs(z))*(-1.0d0*z/Abs(z))

                r = dSqrt(x*x+y*y+z*z)

                If (r .LT. Rc) Then
                    Lj = Epsil * ((Sigma/r)**12 - (Sigma/r)**6) - UFc + Fc*r 
                    NewPE = NewPE + Lj 
                    LjForce = Epsil * ( (12.0d0*SigmaPow12/(r**13)) - (6.0d0*SigmaPow6/(r**7)) ) - Fc 
              
                    NewForces(3*i-2) = NewForces(3*i-2) + LjForce * (x/r)
                    NewForces(3*i-1) = NewForces(3*i-1) + LjForce * (y/r)
                    NewForces(3*i  ) = NewForces(3*i  ) + LjForce * (z/r)

                    NewForces(3*k-2) = NewForces(3*k-2) - LjForce * (x/r)
                    NewForces(3*k-1) = NewForces(3*k-1) - LjForce * (y/r)
                    NewForces(3*k  ) = NewForces(3*k  ) - LjForce * (z/r)
                End If
        End Do
    End Do
 End Subroutine CalcForce

!_________________________________________AUX ROUTINE: REFRESH NEIGBOURS____________________________________________________
Subroutine RefreshNeibours
    Use SimulationParameters
    Implicit None

    Integer         ::  ii,jj

    NumOfNeighbours = 0 ; NeighboursTable = 0

    Do ii = 1,NumParticles-1,1
        x1 = Positions(3*ii-2) ; y1 = Positions(3*ii-1) ; z1 = Positions(3*ii  ) 
        Do jj = ii+1,NumParticles,1
            x2 = Positions(3*jj-2) ; y2 = Positions(3*jj-1) ; z2 = Positions(3*jj  )
            x = x1-x2 ; y = y1-y2 ; z = z1-z2 ; 

            If (Abs(x) .GE. LXBy2) x = (LX-Abs(x))*(-1.0d0*x/Abs(x))
            If (Abs(y) .GE. LYBy2) y = (LY-Abs(y))*(-1.0d0*y/Abs(y))
            If (Abs(z) .GE. LYBy2) z = (LZ-Abs(z))*(-1.0d0*z/Abs(z))

            r = dSqrt(x*x+y*y+z*z)

            If (r .LT. Rs) Then
               NumOfNeighbours(ii) = NumOfNeighbours(ii) + 1
               if( NumOfNeighbours(ii) .GT. MaxNumOfNeighbours) Write(*,*) NumOfNeighbours(ii),"Out Of Bound",MaxNumOfNeighbours
               NeighboursTable(ii,NumOfNeighbours(ii)) = jj 
            End If 
        End Do ! j Loop
    End do ! i Loop   

    ! MostNumOfNeighbours = MAXVAL(NumOfNeighbours)
    ! AverageDensity = dFloat(MostNumOfNeighbours+1)/SkinSphereVolume
    ! AvgNumOfNeighbours = CEILING( dFloat(SUM(NumOfNeighbours)) / dFloat(NumParticles) )
    ! NumOfNeighbours(NumParticles),NeighboursTable(NumParticles,MaxNumOfNeighbours) 
End Subroutine RefreshNeibours
!___________________________________________CALCULATE g(r)________________________________________________________________
Subroutine CalcGofR
    Use SimulationParameters
    Implicit None

    Integer     ::  NumShell ! Number of Particles within Shell bounded by r and r+dr
    Integer     ::  ParticleNum, BinIndex,ParticleCount

    ParticleCount = 1
    NumShell = 0
    
    Call RANDOM_NUMBER(RndNum) ; ParticleNum  = INT(RndNum * NumParticles) ! Particle Randomly Selected 
    GOfRCount = GOfRCount + 1
    x1 = Positions(3*ParticleNum -2) ; y1 = Positions(3*ParticleNum -1) ; z1 = Positions(3* ParticleNum   )
        
    Do i = 1,NumParticles,1
            If (i .NE. ParticleNum) Then
                x2 = Positions(3*i -2) ; y2 = Positions(3*i -1) ; z2 = Positions(3*i   )
                x = x1-x2 ; y = y1-y2 ; z = z1-z2 ; 

                If (Abs(x) .GE. LXBy2) x = (LX-Abs(x))*(-1.0d0*x/Abs(x))
                If (Abs(y) .GE. LYBy2) y = (LY-Abs(y))*(-1.0d0*y/Abs(y))
                If (Abs(z) .GE. LYBy2) z = (LZ-Abs(z))*(-1.0d0*z/Abs(z))

                r = dSqrt(x*x+y*y+z*z)

                if (r .LE. GOfRMax) ParticleCount = ParticleCount + 1

                If ((r .GE. minR) .AND. (r .LE. GOfRMax)) Then
                    BinIndex =  CEILING(r/dr)                   
                    GOfR(BinIndex) = GOfR(BinIndex) + 1.0d0 
                End If  
            End If 
    End Do 
    AverageDensity = dFloat(ParticleCount) 
    SumAverageDensity = SumAverageDensity + AverageDensity
End Subroutine CalcGofR

!_________________________________________CALCULATE THERMODYNAMIC AVERGES__________________________________________________
Subroutine CalcTDAvgs
    Use SimulationParameters
    Implicit None
    

End Subroutine CalcTDAvgs