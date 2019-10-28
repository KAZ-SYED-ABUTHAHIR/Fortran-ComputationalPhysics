PROGRAM RAN_SEED
    !USE RANDOM_SEED_MOD
    IMPLICIT NONE
    INTEGER :: I,N = 10
    REAL    :: R
    !CALL SETUP_RANDOM()
    DO I=1,N
        CALL RANDOM_NUMBER(R)
        WRITE(*,*) "[",I,"] :",R
    END DO
    PRINT *, "Hello RANDOM World!"
END PROGRAM RAN_SEED
