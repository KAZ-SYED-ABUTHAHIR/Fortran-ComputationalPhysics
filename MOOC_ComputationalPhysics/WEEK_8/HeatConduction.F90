! Program to Implement Algorithm to Iterate Towards Steady State Temperature Distribuition of a Rectangualr 
! Grid Under Dirichlet Boundary Conditions

Program HeatDirichlet
    Implicit None
    Integer,Parameter   :: L = 34 ! Grid Size (Sides)
    Real*8              :: Grid(L,L),Grid_Old(L,L),Tolerance = 0.00010d0
    Integer             :: i,j      ! Loop Indices
    Logical             :: Finished = .False.

    ! Init Array Except Boundaries
    Do i = 2,L-1,1
        Do j = 2,L-1,1
            Grid_Old(j,i) = 0.0d0 ; Grid(j,i) = 0.0d0
        End Do
    End Do

    ! Init Boundaries in Accordance with Dirichlet Boundary Conditions
    Do i = 1,L,1
        Grid_Old(1,i) = 3.70d0 ; Grid(1,i) = Grid_Old(1,i)
        Grid_Old(L,i) = 0.40d0 ; Grid(L,i) = Grid_Old(L,i)
        Grid_Old(i,1) = 3.70d0 - (3.30d0/dFloat(L-1))*dFloat(i-1) ; Grid(i,1) = Grid_Old(i,1)
        Grid_Old(i,L) = 3.70d0 - (3.30d0/dFloat(L-1))*dFloat(i-1) ; Grid(i,L) = Grid_Old(i,L)
    End Do

    Open(Unit = 11,File = "ijT_Init.dat")
    Do i = 1,L,1
        Do j = 1,L,1
            Write(11,*) i,j,Grid(i,j)
        End Do
    End Do
    Close(11)

    MainLoop: Do While(.Not.(Finished))
        ! Update Grid Points
        Do i = 2,L-1,1
            Do j = 2,L-1,1
                Grid(j,i) = (Grid_Old(j-1,i)+Grid_Old(j+1,i)+Grid_Old(j,i-1)+Grid_Old(j,i+1))*0.250d0
                !Write(*,*) j,i,Grid(j,i)
            End Do
        End Do
       
        ! Let Finished =.True. Temporarily
        Finished = .True.
        ! Make Finished = .False. Even if a single difference is greater than the Tolerance level
        Do i = 2,L-1,1
            Do j = 2,L-1,1
                If (Abs(Grid(j,i)-Grid_Old(j,i)) > Tolerance) Then
                    Finished = .False.
                End If
            End Do
        End Do
        ! Save the Current Grid in Grid_Old
        Grid_Old = Grid
    End Do MainLoop 

    ! Algorithm Finished. Write the Data to the File
    Open(Unit = 11,File = "ijT.dat")
    Do i = 1,L,1
        Do j = 1,L,1
            Write(11,*) i,j,Grid(i,j)
        End Do
    End Do
    Close(11)
End Program HeatDirichlet