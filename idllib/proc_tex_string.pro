; this processes a tex string
  PRO proc_tex_string,lu,sps_name,buff
  buff = ''
  buff1 = ''
; disregard comment specifiers
  READF,lu,buff
  WHILE STRPOS(buff,'%') GT -1 DO READF,lu,buff
  WHILE STRPOS(buff,'\\') EQ -1 DO BEGIN
    READF,lu,buff1
    buff = buff + buff1
  ENDWHILE
; extract molecule name
  first_amp = STRPOS(buff,'&')
  sps_name = STRTRIM(STRMID(buff,0,first_amp),2)
  buff = STRMID(buff,first_amp,STRLEN(buff)-first_amp)
; now replace every occurance of & with a space
  buff = STRJOIN(STRSPLIT(buff,'&',/EXTRACT),' ',/SINGLE)
; now replace every occurance of $\times *10^{ with e
  buff = STRJOIN(STRSPLIT(buff,'\$\\times 10\^\{',/EXTRACT,/REGEX),'e', $
         /SINGLE)
; get rid of trailing }$
  buff = STRJOIN(STRSPLIT(buff,'\}\$',/EXTRACT,/REGEX),'',/SINGLE)
; get rid of trailing \\
  buff = STRJOIN(STRSPLIT(buff,'\\\\',/EXTRACT,/REGEX),'',/SINGLE)
  RETURN
  END
