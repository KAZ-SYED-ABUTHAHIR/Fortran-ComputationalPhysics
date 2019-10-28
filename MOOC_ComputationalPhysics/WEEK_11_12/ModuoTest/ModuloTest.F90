Program Modulox
    Implicit None
    Real*8                      :: llx 
    Real*8                      :: x,y,z
    llx = 30
    x = MODULO(27.05d0,llx) ; y = MODULO(30.05d0,llx) ; z = MODULO(-0.03d0,llx)
    Write(*,*) "x = MODULO(27.05d0,30.0d0) -----> X = ",x
    Write(*,*) "y = MODULO(30.05d0,30.0d0) -----> Y = ",y
    Write(*,*) "z = MODULO(-0.03d0,30.0d0) -----> Z = ",z
End Program Modulox