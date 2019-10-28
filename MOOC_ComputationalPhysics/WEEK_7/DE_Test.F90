Program Main
    !Program to test DE Subroutines Module
    Use DE_SUBS
    Implicit None
    Real*8          :: x,y
    Real*8,Parameter:: dx = 0.01d0
    Integer         :: i,N
    Real*8,Parameter:: XL = 0.0d0,XU = 1.55d0

    Interface
            Function f(a,b) Result(retval)
                Implicit None
                Real*8 :: a,b
                Real*8 :: retval                
            End Function f
    End Interface

    N = INT((XU - XL) / dx) + 1
   
    x = XL
    y = 0.0d0

    Open(UNIT = 11,file = "RungeKutta_IV.dat")
    Open(UNIT = 22,file = "TanX.dat")

    Do i = 1,N,1
        Write(11,*) x,y
        Write(22,*) x,tan(x)
        Call RungeKuttaIV(x,y,dx,f)
    End Do
    Close(11) ; CLose(22)
End Program Main

Function f(x,y) Result(R)
    Implicit None
    Real*8,Intent(IN)        :: x,y
    Real*8                   :: R

    R = (y * y) + 1.0d0
End Function f 