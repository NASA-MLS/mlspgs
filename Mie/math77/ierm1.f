      subroutine IERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1990-01-18 CLL Added Integer stmt for VALUE.  Typed all variables.
C>> 1985-08-02 IERM1  Lawson  Initial code.
C
      integer INDIC, LEVEL, VALUE
      character*(*) SUBNAM,MSG,LABEL
      character*1 FLAG
      call ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
      call IERV1(LABEL,VALUE,FLAG)
C
      return
      end
