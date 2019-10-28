! Program to Simulate a Circular Chain of Coupled Oscillators using RK4-Coupled Mehtod
Program CircOsc
    Use DE_SUBS
    Implicit None
    Interface
            Function dydt(t,yi,v,yi_Minus_1,yi_Plus_1) Result(retval)
                Implicit None
                Real*8 :: t,yi,v,yi_Minus_1,yi_Plus_1
                Real*8 :: retval                
            End Function dydt 

            Function dvdt(t,yi,v,yi_Minus_1,yi_Plus_1) Result(retval)
                Implicit None
                Real*8 :: t,yi,v,yi_Minus_1,yi_Plus_1
                Real*8 :: retval                
            End Function dvdt 
    End Interface

    Integer,Parameter   :: NUM_OSC = 20

    Real*8              :: t,y(NUM_OSC),v(NUM_OSC),temp_y(NUM_OSC),temp_v(NUM_OSC),tempyi,tempvi,x(NUM_OSC),z(NUM_OSC),theta = 0.0d0
    Real*8,Parameter    :: dt = 0.025d0 ! Time Step 
    Integer             :: Iters = 20000,i,j,k,l,m,n
    Real*8,Parameter    :: R = 4.0d0 ! Radius of the Circular Chain of Oscillators
    Real*8,Parameter    :: TwoPI = 4.0d0*asin(1.0d0),ThetaStep = TwoPI/Real(NUM_OSC)

    !Open(UNIT = 11,FILE = "ntyv.dat")
    !Open(UNIT = 22,FILE = "tyv_particle_1.dat")
    Open(UNIT = 33,FILE = "MovieOsc.xyz")

    ! Init Positions and Velocities,All velocities are zero and all dispacements are zero
    ! Except for first and twenty sixth particles
    Do i = 1,NUM_OSC,1
        theta =  ThetaStep*Real(i-1)
        x(i) = R*cos(theta) ; z(i) = R*sin(theta)
        y(i) = 1.0d0 * sin(2.0d0*theta) ; temp_y(i) = y(i) ! Standing Wave
        ! y(i) = 4.0d0 * Exp(-15.0d0*(theta-1.570d0)**2) ; temp_y(i) = y(i) ! Gaussian Pulse
        ! y(i) = 4.0d0 * (sin(15d0*theta)+sin(2d0*theta)) ; temp_y(i) = y(i) 
        !y(i) = 0.0d0 ; !temp_y(i) = 0.0d0
        !v(i) = 0.0d0 ; !temp_v(i) = 0.0d0    
    End Do
    !y(1) = 0.80d0 ; y(INT(NUM_OSC/2)+1) = 0.80d0
    temp_y = y ; temp_v = v 

    t = 0.0d0

    !Do j = 1,NUM_OSC,1
        !Write(11,*) j,t,y(j),v(j)
    !End Do

    !Write(22,*) t,y(1),v(1)

    Do k = 1,Iters,1

        ! Update Inner points
        Do l = 2,NUM_OSC-1,1
            tempyi = temp_y(l) ; tempvi = temp_v(l)
            Call RungeKuttaIV_CoupledOsc(t,tempyi,tempvi,dt,dydt,dvdt,temp_y(l-1),temp_y(l+1))
            y(l) = tempyi ; v(l) = tempvi
        End Do
        ! Update Boundary Points in accordance with Periodic Boundary Conditions
        tempyi = temp_y(1) ; tempvi = temp_v(1) 
        Call RungeKuttaIV_CoupledOsc(t,tempyi,tempvi,dt,dydt,dvdt,temp_y(NUM_OSC),temp_y(2))
        y(1) = tempyi ; v(1) = tempvi

        tempyi = temp_y(NUM_OSC) ; tempvi = temp_v(NUM_OSC) 
        Call RungeKuttaIV_CoupledOsc(t,tempyi,tempvi,dt,dydt,dvdt,temp_y(NUM_OSC-1),temp_y(1))
        y(NUM_OSC) = tempyi ; v(NUM_OSC) = tempvi

        ! Update t 
        t = (k * dt)

        ! Update the temporary placeholder array
        temp_y = y ; temp_v = v 

        ! k is the iteration index
        If (Mod(k,10) == 0) Then ! Write the Position Data To The JMol Animation File Once In 10 Steps
            Write(33,*) NUM_OSC  ! NUM_OSC - Number of Oscillators in the Circular Chain
            Write(33,*) ""
            Do n = 1,NUM_OSC,1
                ! Write the data to the file 
                Write(33,*) "N",x(n),y(n),z(n)
            End Do 
        End If
        !Write(22,*) t,y(1),v(1)
    End Do
    !Close(11)  
    !Close(22)
    Close(33)
End Program CircOsc

Function dvdt(t,yi,v,yi_Minus_1,yi_Plus_1) Result(retval)
    Implicit None
    Real*8 :: t,yi,v,yi_Minus_1,yi_Plus_1
    Real*8 :: retval     
    Real*8,Parameter :: k_by_m = 1.0d0
    retval = k_by_m*(yi_Minus_1+yi_Plus_1-(2.0d0*yi))       
End Function dvdt 

Function dydt(t,yi,v,yi_Minus_1,yi_Plus_1) Result(retval)
    Implicit None
    Real*8 :: t,yi,v,yi_Minus_1,yi_Plus_1
    Real*8 :: retval    
    retval = v            
End Function dydt 