!Program to implement Crude Monte Carlo Integration Method
!To Integrate Exp(X): L_Lim -> U_Lim
!Author	:	KAZ
!Date	:	13.08.2019
!Context:	Computational Physics MOOC Week - 3 Exercise

Program MonteCarloCrude
	Implicit None
	Integer					:: I,N,Hits 						! Counter,Sampling Size and Counter for Num of Hits
	Real*8					:: MC,Actual,Mean,MeanSqr,R,F,SD    ! R - Random Number Generated, Function Value
	Real*8					:: L_Lim,U_Lim,Length,X,Y,Height	! X - Rnd Variable
	Real*8					:: EvaluateF						! Integrand Function

	Mean 	= 0.0D0
	MeanSqr = 0.0D0
	Hits	= 0

	Write(*,*) "Enter the Number of Random Points to be used:"
	Read(*,*)	N

	Write(*,*) "Enter the Lower Limit of the Integral   :"
	Read(*,*) L_Lim
	Write(*,*) "Enter the Upper Limit of the Integral   :"
	Read(*,*) U_Lim

	Length = U_Lim - L_Lim
	Height = Max(EvaluateF(U_Lim),EvaluateF(L_Lim))	
	Actual = EvaluateF(U_Lim) - EvaluateF(L_Lim)  		! CAUTION: Only Works for F = EXP(X)

	Do I = 1,N
		! Generate a Random Point within the rectangular area bounding the integrand curve
		Call Random_Number(R)
		X =  L_Lim + Length * R 
		Call Random_Number(R)
		Y =  Height * R 

		! Evaluate Integrand at that Point's X value. If the Y is less than the Integrand value a hit is made
		F =  EvaluateF(X)
		If (Y < F) Then
			Hits = Hits + 1
		End If 

		Mean = Mean + F
		MeanSqr = MeanSqr + F*F
	End Do

	Mean 	= Mean/Real(N)
	MeanSqr = MeanSqr/Real(N)
	SD 		= Sqrt(MeanSqr - Mean*Mean)

	MC =  Length * Height * (Real(Hits)/Real(N))

	Write(*,100) N
	100 Format ("Number of Random Points                        : ",I12)
	Write(*,101) Mean
	101 Format ("Average of Random Function Values              : ",F12.6)
	Write(*,102) SD
	102 Format ("Standard Deviation of Random Function Values   : ",F12.6)

	Write(*,110) MC
	110 Format ("Integral Value Calculated as                   : ",F12.6)
	Write(*,120) Actual	
	120 Format ("Actual Value of the Integral is                : ",F12.6)

End Program MonteCarloCrude

!Function to Evaluate the Integrand 
Function EvaluateF(X) Result(R)
	Implicit None
	Real*8	 			 ::	 R
	REAL*8, Intent(In)   ::  X
	R = Exp(X)
End Function EvaluateF