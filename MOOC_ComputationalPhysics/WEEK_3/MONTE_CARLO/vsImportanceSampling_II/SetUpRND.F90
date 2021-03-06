MODULE RANDOM_SEED_MOD
    IMPLICIT NONE
    CONTAINS
        SUBROUTINE SETUP_RANDOM()
            !VARIABLES FOR PORTABLE SEED SETTING
            INTEGER                             ::  I_SEED,I
            INTEGER, DIMENSION(:),ALLOCATABLE   ::  A_SEED
            INTEGER, DIMENSION(1:8)             ::  DT_SEED
            !SETUP RANDOMSEED PORTABLY
            CALL RANDOM_SEED(SIZE = I_SEED)
            !WRITE(*,*) I_SEED
            ALLOCATE(A_SEED(1:I_SEED))
            CALL RANDOM_SEED(GET = A_SEED)
            CALL DATE_AND_TIME(VALUES = DT_SEED)
            DO I=1,8
                A_SEED(I) = DT_SEED(I)
            END DO
            CALL RANDOM_SEED(PUT = A_SEED)
            DEALLOCATE(A_SEED)
            !---------------DONE SETTING UP RANDOM SEED----------------!
        END SUBROUTINE SETUP_RANDOM
    END MODULE RANDOM_SEED_MOD