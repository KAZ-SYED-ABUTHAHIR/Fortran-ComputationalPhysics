! Program to solve the Non Linear Oscillator Problem using RK4 Propagator
Program NLO
    Use DE_SUBS 
    Implicit None
    Integer             :: i 
    Integer,Parameter   :: Iters = 5001
    Real*8,Parameter    :: del_T = 0.01d0
    Real*8,Parameter    :: x0 = 0.0d0,v0 = 1.999
    Real*8              :: x,v,t
    
    Interface
            Function dvdt(a,b,c) Result(retval)
                Implicit None
                Real*8 :: a,b,c
                Real*8 :: retval                
            End Function dvdt

            Function dxdt(a,b,c) Result(retval)
                Implicit None
                Real*8 :: a,b,c
                Real*8 :: retval                
            End Function dxdt
    End Interface

    x = x0 ; v = v0 ; t = 0.0d0

    open(UNIT = 11,file = "SHO_X.dat")
    open(UNIT = 22,file = "SHO_V.dat")


    Do i = 1,Iters,1
        write(11,*) t,x
        write(22,*) t,v
        ! Propagate t,x,v
        Call RungeKuttaIV_Coupled(t,x,v,del_T,dxdt,dvdt)
    End Do

    Close(11)
End Program NLO

Function dvdt(t,x,v) Result(R)
    Implicit None
    Real*8 :: t,x,v
    Real*8 :: R     
    
    R = -sin(x)    !
End Function dvdt

Function dxdt(t,x,v) Result(R)
    Implicit None
    Real*8 :: t,x,v
    Real*8 :: R     
    
    R = v       !
End Function dxdt
