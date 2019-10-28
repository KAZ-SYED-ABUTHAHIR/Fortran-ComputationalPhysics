PROGRAM OutputFile

    IMPLICIT NONE
    INTEGER:: u = 7, ios ! Unit Number of File to Written, IOSTAT == 0 is OK
    REAL:: r = 3.14159265358979
    REAL, DIMENSION(5):: v = (/ 1.1,1.4,1.5,1.9,1.4 /)
    
    ! & ==> Line Continuation
    ! IOSTAT > 0 , Previous read Failed, IOSTAT == 0, Previous read Succeeded and IOSTAT < 0, EOF Condition
    OPEN(UNIT = u,IOSTAT = ios,FILE = 'output.txt',STATUS = 'new',ACTION = 'write')
    
    IF (ios .EQ. 0) THEN
        WRITE (u,'(5f5.1)') v
        WRITE (u,'(f5.3)') r
        CLOSE(u)
    ELSE
        PRINT '(a30)' , "ERROR: File Not Open!"
    END IF
    
END PROGRAM OutputFile