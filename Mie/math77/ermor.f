      SUBROUTINE ERMOR(MSG,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1985-09-20 ERMOR  Lawson  Initial code.
C
C     --------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     MSG      Message to be printed as part of the diagnostic.
C
C     FLAG     A single character,which when set to '.' will
C              call the subroutine ERFIN and will just RETURN
C              when set to any other character.
C
C     --------------------------------------------------------------
C
      COMMON/M77ERR/IDELTA,IALPHA
      INTEGER IDELTA,IALPHA
      SAVE /M77ERR/
      CHARACTER*(*) MSG
      CHARACTER*1 FLAG
C
      IF (IALPHA.GE.-1) THEN
        WRITE (*,*) MSG
        IF (FLAG .EQ. '.') CALL ERFIN
      END IF
C
      RETURN
      END
