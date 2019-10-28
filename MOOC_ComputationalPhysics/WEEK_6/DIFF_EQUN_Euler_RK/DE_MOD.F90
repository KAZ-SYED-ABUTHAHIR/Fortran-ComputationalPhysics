Module DE_SUBS
    Implicit None
    Contains 
    Subroutine SimpleEuler(x,y,h,func) ! Simple Euler Propagator 
        Interface
            Function func(a,b) Result(retval)
                Implicit None
                Real*8 :: a,b
                Real*8 :: retval                
            End Function func
        End Interface
        
        Real*8  :: x,y
        Real*8  :: h

        y = y + h*func(x,y)  ! Propagate y to the next value
        x = x + h            ! Propagate x to the next value
    End Subroutine SimpleEuler

    Subroutine ModifiedEuler(x,y,h,func) ! Modified Euler Propagator 
        Interface
            Function func(a,b) Result(retval)
                Implicit None
                Real*8 :: a,b
                Real*8 :: retval                
            End Function func
        End Interface
        
        Real*8  :: x,y
        Real*8  :: h
        
        Real*8  :: temp_y

        temp_y = y + 0.50d0*h*func(x,y)  ! Propagate y to the next temp value
        y = y + h * func(x,temp_y)       ! Now propagate y to the next value by recalcuating f(x,y) using temp y vaue
        x = x + h                        ! Propagate x to the next value
    End Subroutine ModifiedEuler

    Subroutine ImprovedEuler(x,y,h,func) ! Improved Euler Propagator 
        Interface
            Function func(a,b) Result(retval)
                Implicit None
                Real*8 :: a,b
                Real*8 :: retval                
            End Function func
        End Interface
        
        Real*8  :: x,y
        Real*8  :: h       
        Real*8  :: f0,fe,temp_y

        f0 = func(x,y)              
        temp_y = y + h*f0  
        fe = func(x+h,temp_y)

        y = y + h * (f0 + fe) * 0.5d0   
        x = x + h                        
    End Subroutine ImprovedEuler

    Subroutine RungeKuttaIV_Coupled(t,x,v,dt,dxdt,dvdt) ! Fourth Order Runge-Kutta Propagator for a set of coupled Ist order DEs
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
        
        ! Don't Propagate Velocity and Position Simultaneously !  ??? Don't Do this ! Not Symplectic You Fool 
        
        
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
               
    End Subroutine RungeKuttaIV_Coupled

    Subroutine RungeKuttaIV(x,y,h,func) ! Fourth Order Runge-Kutta Propagator 
        Interface
            Function func(a,b) Result(retval)
                Implicit None
                Real*8 :: a,b
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
End Module DE_SUBS