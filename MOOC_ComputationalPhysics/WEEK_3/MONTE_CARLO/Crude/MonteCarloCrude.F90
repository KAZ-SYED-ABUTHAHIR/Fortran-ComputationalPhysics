!Program to implement Crude Monte Carlo Integration Method
!To Integrate Exp(X): L_Lim -> U_Lim
!Author	:	KAZ
!Date	:	11.08.2019
!Context:	Computational Physics MOOC Week - 3 Exercise

Program MonteCarloCrude
	Implicit None
	Integer					:: I,N 
	Real*8					:: MC,Actual,Mean,MeanSqr,R,F,SD
	Real*8					:: L_Lim,U_Lim,Length,X
	Real*8					:: EvaluateF

	Mean = 0.0D0
	MeanSqr = 0.0D0
	MC = 0.0D0

	Write(*,*) "Enter the Number of Random Points to be used:"
	Read(*,*)	N

	Write(*,*) "Enter the Lower Limit of the Integral   :"
	Read(*,*) L_Lim
	Write(*,*) "Enter the Upper Limit of the Integral   :"
	Read(*,*) U_Lim

	Length = U_Lim - L_Lim
	Actual = EvaluateF(U_Lim) - EvaluateF(L_Lim)

	Do I = 1,N
		Call Random_Number(R)
		X = L_Lim + Length * R 
		F = EvaluateF(X)
		MC = MC + F
		MeanSqr = MeanSqr + F*F
	End Do

	MC  = MC/Real(N)

	Mean = MC
	MeanSqr = MeanSqr/Real(N)

	SD = Sqrt(MeanSqr - Mean*Mean)

	MC = MC * Length

	Write(*,100) N
	100 Format ("Number of Random Points                        : ",I10)
	Write(*,101) Mean
	101 Format ("Average of Random Function Values              : ",F9.6)
	Write(*,102) SD
	102 Format ("Standard Deviation of Random Function Values   : ",F9.6)

	Write(*,110) MC
	110 Format ("Integral Value Calculated as                   : ",F9.6)
	Write(*,120) Actual	
	120 Format ("Actual Value of the Integral is                : ",F9.6)

End Program MonteCarloCrude

!Function to Evaluate the Integrand 
Function EvaluateF(X) Result(R)
	Implicit None
	Real*8	 			 ::	 R
	REAL*8, Intent(In)   ::  X
	R = Exp(X)
End Function EvaluateF