; this processes another tex string
  PRO proc_tex_stra,lu,buff,stra,strb,sps_name
  buff = ''
  buff1 = ''
; disregard comment specifiers
  READF,lu,buff
  WHILE STRPOS(buff,'%') GT -1 DO READF,lu,buff
; read complete string and join it together
  WHILE STRPOS(buff,'\\') EQ -1 DO BEGIN
    READF,lu,buff1
    WHILE STRPOS(buff1,'%') GT -1 DO READF,lu,buff1
    buff = buff + buff1
  ENDWHILE
  IF STRPOS(buff,'\dotfill') GT -1 THEN BEGIN
    amp_pos = STRSPLIT(buff,'\\dotfill',LENGTH = len,/REGEX)
    sps_name = STRTRIM(STRMID(buff,amp_pos(1),len(1)),2)
  ENDIF ELSE BEGIN
    sps_name = ''
; get rid of trailing \\
    buff = STRJOIN(STRSPLIT(buff,'\\\\',/EXTRACT,/REGEX),'',/SINGLE)
; extract band names
    amp_pos = STRSPLIT(buff,'&',LENGTH = len)
    stra = STRTRIM(STRMID(buff,amp_pos(11),len(11)),2)
    strb = STRTRIM(STRMID(buff,amp_pos(12),len(12)),2)
; extract non character part of the string
    buff = STRMID(buff,0,amp_pos(11)-1)
; now replace every occurance of & with a space
    buff = STRJOIN(STRSPLIT(buff,'&',/EXTRACT),' ',/SINGLE)
; now replace every occurance of $\times 10^{ with e
    buff = STRJOIN(STRSPLIT(buff,'\$\\times *10\^\{',/EXTRACT,/REGEX),'e', $
           /SINGLE)
; get rid of trailing }$
    buff = STRJOIN(STRSPLIT(buff,'\}\$',/EXTRACT,/REGEX),'',/SINGLE)
  ENDELSE
  RETURN
  END
