MODULE MOMENT_ROUTINE

	IMPLICIT NONE
CONTAINS 
	FUNCTION MOMENT(LEN,ARRAY,K) RESULT(RET_MOMENT)
		IMPLICIT NONE
		INTEGER				::	LEN,K,I
		REAL,DIMENSION(LEN)	::	ARRAY
		REAL 				:: 	SUM_OF_POWERS
		REAL 				::	RET_MOMENT

		SUM_OF_POWERS = 0.0D0

		DO I = 1,LEN
			SUM_OF_POWERS = SUM_OF_POWERS + ARRAY(I)**K
		END DO

		RET_MOMENT = SUM_OF_POWERS / REAL(LEN)
	END FUNCTION MOMENT

END MODULE MOMENT_ROUTINE