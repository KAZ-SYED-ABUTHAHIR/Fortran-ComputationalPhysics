! Module Containing Suroutines to Solve Differential Equations By Different Methods
Module DE_SUBS
    Implicit None
    Contains 
    Subroutine SimpleEuler(x,y,h,func) ! Simple Euler Propagator for DE of the From y' = f(x,y)
        Interface
            Function func(x,y) Result(retval)
                Implicit None
                Real*8 :: x,y 
                Real*8 :: retval                
            End Function func
        End Interface
        
        Real*8  :: x,y
        Real*8  :: h

        y = y + h*func(x,y)  ! Propagate y to the next value
        x = x + h            ! Propagate x to the next value
    End Subroutine SimpleEuler

    Subroutine ModifiedEuler(x,y,h,func) ! Modified Euler Propagator for DE of the From y' = f(x,y)
        Interface
            Function func(x,y) Result(retval)
                Implicit None
                Real*8 :: x,y
                Real*8 :: retval                
            End Function func
        End Interface
        
        Real*8  :: x,y
        Real*8  :: h
        
        Real*8  :: temp_y

        temp_y = y + 0.50d0*h*func(x,y)  ! Propagate y to the next temp value, Half Step
        y = y + h * func(x,temp_y)       ! Now propagate y to the next value by recalcuating f(x,y) using temp y vaue
        x = x + h                        ! Propagate x to the next value
    End Subroutine ModifiedEuler

    Subroutine ImprovedEuler(x,y,h,func) ! Improved Euler Propagator for DE of the From y' = f(x,y)
        Interface
            Function func(x,y) Result(retval)
                Implicit None
                Real*8 :: x,y
                Real*8 :: retval                
            End Function func
        End Interface
        
        Real*8  :: x,y
        Real*8  :: h       
        Real*8  :: f0,fe,end_y
        ! Predictor
        f0 = func(x,y)              
        end_y = y + h*f0  
        ! Corrector
        fe = func(x+h,end_y)
        y = y + h * (f0 + fe) * 0.5d0   
        x = x + h                        
    End Subroutine ImprovedEuler

    Subroutine RungeKuttaIV(x,y,h,func) ! Fourth Order Runge-Kutta Propagator for DE of the From y' = f(x,y) 
        Interface
            Function func(x,y) Result(retval)
                Implicit None
                Real*8 :: x,y
                Real*8 :: retval                
            End Function func
        End Interface
        
        Real*8  :: x,y
        Real*8  :: h       
        Real*8  :: K1,K2,K3,K4
        K1 = func(x,y)
        K2 = func(x+0.50d0*h,y+0.50d0*h*K1)
        K3 = func(x+0.50d0*h,y+0.50d0*h*K2)
        K4 = func(x+h,y+h*K3)
        y = y + h*(K1+2.0d0*K2+2.0d0*K3+K4)/6.0d0
        x = x + h                        
    End Subroutine RungeKuttaIV

    Subroutine RungeKuttaIV_IInd_Order(t,x,v,dt,dxdt,dvdt) ! Fourth Order Runge-Kutta Propagator for a set of coupled Ist order DEs
                                                           ! Represending a IInd Order Equation x'' = f(x,x') ; v = x'
        Interface
            Function dxdt(t,x,v) Result(retval)
                Implicit None
                Real*8 :: t,x,v
                Real*8 :: retval                
            End Function dxdt 

            Function dvdt(t,x,v) Result(retval)
                Implicit None
                Real*8 :: t,x,v
                Real*8 :: retval                
            End Function dvdt 
        End Interface
        
        Real*8  :: t,x,v
        Real*8  :: dt       
        Real*8  :: KX1,KX2,KX3,KX4,KV1,KV2,KV3,KV4
        
        ! Don't Propagate Velocity and Position Simultaneously !  ??? Don't Do this ! Not Symplectic 
             
        KV1 = dvdt(t,x,v)
        KV2 = dvdt(t+0.50d0*dt,x,v+0.50d0*dt*KV1)
        KV3 = dvdt(t+0.50d0*dt,x,v+0.50d0*dt*KV2)
        KV4 = dvdt(t+dt,x,v+dt*KV3)
        v = v + dt * (KV1+2.0d0*KV2+2.0d0*KV3+KV4)/6.0d0 
        
        KX1 = dxdt(t,x,v)
        KX2 = dxdt(t+0.50d0*dt,x+0.50d0*dt*KX1,v)
        KX3 = dxdt(t+0.50d0*dt,x+0.50d0*dt*KX2,v)
        KX4 = dxdt(t+dt,x+dt*KX3,v)
        x = x + dt * (KX1+2.0d0*KX2+2.0d0*KX3+KX4)/6.0d0
        
        t = t + dt          
               
    End Subroutine RungeKuttaIV_IInd_Order

    Subroutine RungeKuttaIV_CoupledOsc(t,yi,v,dt,dydt,dvdt,yi_Minus_1,yi_Plus_1) ! Fourth Order Runge-Kutta Propagator for a set of coupled Ist order DEs
                                                                                 ! Represending a IInd Order Equation x'' = f(x,x') ; v = x'
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

        Real*8  :: t,yi,v,yi_Minus_1,yi_Plus_1
        Real*8  :: dt       
        Real*8  :: KY1,KY2,KY3,KY4,KV1,KV2,KV3,KV4

        ! Don't Propagate Velocity and Position Simultaneously !  ??? Don't Do this ! Not Symplectic !


        KV1 = dvdt(t,yi,v,yi_Minus_1,yi_Plus_1)
        KV2 = dvdt(t+0.50d0*dt,yi,v+0.50d0*dt*KV1,yi_Minus_1,yi_Plus_1)
        KV3 = dvdt(t+0.50d0*dt,yi,v+0.50d0*dt*KV2,yi_Minus_1,yi_Plus_1)
        KV4 = dvdt(t+dt,yi,v+dt*KV3,yi_Minus_1,yi_Plus_1)
        v = v + dt * (KV1+2.0d0*KV2+2.0d0*KV3+KV4)/6.0d0 

        KY1 = dydt(t,yi,v,yi_Minus_1,yi_Plus_1)
        KY2 = dydt(t+0.50d0*dt,yi+0.50d0*dt*KY1,v,yi_Minus_1,yi_Plus_1)
        KY3 = dydt(t+0.50d0*dt,yi+0.50d0*dt*KY2,v,yi_Minus_1,yi_Plus_1)
        KY4 = dydt(t+dt,yi+dt*KY3,v,yi_Minus_1,yi_Plus_1)
        yi = yi + dt * (KY1+2.0d0*KY2+2.0d0*KY3+KY4)/6.0d0

    End Subroutine RungeKuttaIV_CoupledOsc

    Subroutine JacobiIterator(nop,x,y,tolerance,fun_yi,coeff_yiminus1,coeff_yiplus1,coeff_xi,isFinished)
        Interface
            Function fun_yi(xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi) Result(R)
                Implicit None
                Real*8          :: xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi
                Real*8          :: R
            End Function fun_yi
        End Interface
        
        Integer         :: i,nop ! Number of Grid Points
        Real*8          :: x(nop),y(nop),tolerance,yold(nop),coeff_yiminus1,coeff_yiplus1,coeff_xi
        Logical         :: isFinished

        isFinished = .True.
        
        ! Copy the original Grid to the temporary Grid
        
        yold = y

        Do i = 2,nop-1,1
             y(i) = fun_yi(x(i),yold(i-1),yold(i+1),coeff_yiminus1,coeff_yiplus1,coeff_xi)
             if(Abs(y(i)-yold(i)) > tolerance) Then
                isFinished = .False.
             End If
        End Do
    End Subroutine JacobiIterator

    Subroutine GaussSeidelIterator(nop,x,y,tolerance,fun_yi,coeff_yiminus1,coeff_yiplus1,coeff_xi,isFinished)
        Interface
            Function fun_yi(xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi) Result(R)
                Implicit None
                Real*8          :: xi,yiminus1,yiplus1,coeff_yiminus1,coeff_yiplus1,coeff_xi
                Real*8          :: R
            End Function fun_yi
        End Interface
        
        Integer         :: i,nop ! Number of Grid Points
        Real*8          :: x(nop),y(nop),tolerance,yold(nop),coeff_yiminus1,coeff_yiplus1,coeff_xi
        Logical         :: isFinished

        isFinished = .True.
        
        ! Copy the original Grid to the temporary Grid  
        yold = y

        Do i = 2,nop-1,1
             y(i) = fun_yi(x(i),y(i-1),y(i+1),coeff_yiminus1,coeff_yiplus1,coeff_xi)
             if(Abs(y(i)-yold(i)) > tolerance) Then
                isFinished = .False.
             End If
        End Do
    End Subroutine GaussSeidelIterator
End Module DE_SUBS