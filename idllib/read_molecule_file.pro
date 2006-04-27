;$Id$
; this reads a molecule list file
  PRO read_molecule_file,filename,table_of_molecules,spectags,spectls,hitran, $
                         hitran_iso,hitran_frac,mol_mr,dflt_wc,dflt_nc
  OPENR,lu,filename,/GET_LUN
  buff = ''
  WHILE STRPOS(buff,'format =') EQ -1 DO READF,lu,buff
; extract the format statement
  fmt = STRMID(buff,STRPOS(buff,'('),STRPOS(buff,')')-STRPOS(buff,'(')+1)
  READF,lu,buff
  str1 = ''
  str2 = ''
  str3 = ''
  str4 = ''
  i1 = 0
  r1 = 0.0
  r2 = 0.0
  r3 = 0.0
  r4 = 0.0
  READF,lu,str1,str2,str3,str4,i1,r1,r2,r3,r4, FORMAT = fmt
  table_of_molecules = [STRTRIM(str1,2)]
  spectags = [str2]
  spectls = [str3]
  hitran = [str4]
  hitran_iso = [i1]
  hitran_frac = [r1]
  mol_mr = [r2]
  dflt_wc = [r3]
  dflt_nc = [r4]
  WHILE NOT EOF(lu) DO BEGIN
    READF,lu,str1,str2,str3,str4,i1,r1,r2,r3,r4, FORMAT = fmt
    table_of_molecules = [table_of_molecules,STRTRIM(str1,2)]
    spectags = [spectags,str2]
    spectls = [spectls,str3]
    hitran = [hitran,str4]
    hitran_iso = [hitran_iso,i1]
    hitran_frac = [hitran_frac,r1]
    mol_mr = [mol_mr,r2]
    dflt_wc = [dflt_wc,r3]
    dflt_nc = [dflt_nc,r4]
  ENDWHILE
  FREE_LUN,lu
  RETURN
  END

;$Log$