; program to read spectroscopy database
  PRO read_spect_dbase,mol_file,line_file,species,sps_data
; inputs
; mol_file = input file name for molecular data
; line_file = input file name for line data
; species = string array of species wanted or blank for all of them
; sps_data = a structure array of molecular data
; find out what user wants
  buff = ''
  buff1 = ''
  buff2 = ''
  spsv = ''
  abunv = 0.0
  massv = 0.0
  freqv = 0.0d+00
  gsev  = 0.0
  istrv = 0.0
  wcv   = 0.0
  ncv   = 0.0
  psv   = 0.0
  nsv   = 0.0
  int1v = 0.0
  n1v   = 0.0
  int2v = 0.0
  n2v   = 0.0
  qv = FLTARR(3)
  contv = FLTARR(6)
  OPENR,lu1,mol_file,/GET_LUN
  OPENR,lu2,line_file,/GET_LUN
; first line is a spacer
  READF,lu2,buff1
  IF SCALARIZE(species) EQ '' THEN BEGIN
; user wants everything
    proc_tex_string,lu1,spsv,buff
    READS,buff,abunv,massv,qv,contv
    sps  = [spsv]
    abun = [abunv]
    mass = [massv]
    parf = [qv]
    conf = [contv]
    no_sps = 1
; now load lines for first molecule
    proc_tex_stra,lu2,buff,stra,strb,sps1
    no_linesv = 0
    IF qv(0) GT 0.00001 THEN BEGIN
      IF sps1 EQ spsv THEN BEGIN
        proc_tex_stra,lu2,buff,stra,strb,sps1
        READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
        freq = [freqv]
        gse  = [gsev]
        istr = [istrv]
        wc   = [wcv]
        nc   = [ncv]
        ps   = [psv]
        ns   = [nsv]
        int1 = [int1v]
        n1   = [n1v]
        int2 = [int2v]
        n2   = [n2v]
        eos_band = [stra]
        uars_band = [strb]
        no_linesv = 1
        WHILE NOT EOF(lu2) AND sps1 EQ '' DO BEGIN
          proc_tex_stra,lu2,buff,stra,strb,sps1
          IF sps1 EQ '' THEN BEGIN
            READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
            freq = [freq,freqv]
            gse  = [gse,gsev]
            istr = [istr,istrv]
            wc   = [wc,wcv]
            nc   = [nc,ncv]
            ps   = [ps,psv]
            ns   = [ns,nsv]
            int1 = [int1,int1v]
            n1   = [n1,n1v]
            int2 = [int2,int2v]
            n2   = [n2,n2v]
            eos_band  = [eos_band,stra]
            uars_band = [uars_band,strb]
            no_linesv = no_linesv + 1
          ENDIF
        ENDWHILE
      ENDIF ELSE BEGIN
        PRINT,'species match error has occured'
        PRINT,sps1,' versus ',spsv
        STOP
      ENDELSE
    ENDIF
    no_lines = [no_linesv]
    WHILE NOT EOF(lu1) DO BEGIN
      proc_tex_string,lu1,spsv,buff
      READS,buff,abunv,massv,qv,contv
      sps  = [sps,spsv]
      abun = [abun,abunv]
      mass = [mass,massv]
      parf = [parf,qv]
      conf = [conf,contv]
      no_sps = no_sps + 1
      no_linesv = 0
      IF qv(0) GT 0.00001 THEN BEGIN
        IF sps1 EQ spsv THEN BEGIN
          proc_tex_stra,lu2,buff,stra,strb,sps1
          READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
          freq = [freq,freqv]
          gse  = [gse,gsev]
          istr = [istr,istrv]
          wc   = [wc,wcv]
          nc   = [nc,ncv]
          ps   = [ps,psv]
          ns   = [ns,nsv]
          int1 = [int1,int1v]
          n1   = [n1,n1v]
          int2 = [int2,int2v]
          n2   = [n2,n2v]
          eos_band  = [eos_band,stra]
          uars_band = [uars_band,strb]
          no_linesv = 1
          WHILE NOT EOF(lu2) AND sps1 EQ '' DO BEGIN
            proc_tex_stra,lu2,buff,stra,strb,sps1
            IF sps1 EQ '' THEN BEGIN
              READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
              freq = [freq,freqv]
              gse  = [gse,gsev]
              istr = [istr,istrv]
              wc   = [wc,wcv]
              nc   = [nc,ncv]
              ps   = [ps,psv]
              ns   = [ns,nsv]
              int1 = [int1,int1v]
              n1   = [n1,n1v]
              int2 = [int2,int2v]
              n2   = [n2,n2v]
              eos_band  = [eos_band,stra]
              uars_band = [uars_band,strb]
              no_linesv = no_linesv + 1
            ENDIF
          ENDWHILE
        ENDIF ELSE BEGIN
          PRINT,'species match error has occured'
          PRINT,sps1,' versus ',spsv
          STOP
        ENDELSE
      ENDIF
      no_lines = [no_lines,no_linesv]
    ENDWHILE
  ENDIF ELSE BEGIN
; selected molecule search mode
    no_sps = N_ELEMENTS(species)
    proc_tex_string,lu1,spsv,buff
    WHILE NOT EOF(lu1) AND spsv NE species(0) DO $
        proc_tex_string,lu1,spsv,buff
    IF spsv EQ species(0) THEN BEGIN
; species is found
      READS,buff,abunv,massv,qv,contv
      sps  = [spsv]
      abun = [abunv]
      mass = [massv]
      parf = [qv]
      conf = [contv]
      sps1 = ''
      WHILE spsv NE sps1 DO BEGIN
        READF,lu2,buff
        WHILE STRPOS(buff,'\dotfill') EQ -1 DO READF,lu2,buff
        amp_pos = STRSPLIT(buff,'\\dotfill',LENGTH = len,/REGEX)
        sps1 = STRTRIM(STRMID(buff,amp_pos(1),len(1)),2)
      ENDWHILE
      proc_tex_stra,lu2,buff,stra,strb,sps1
      READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
      freq = [freqv]
      gse  = [gsev]
      istr = [istrv]
      wc   = [wcv]
      nc   = [ncv]
      ps   = [psv]
      ns   = [nsv]
      int1 = [int1v]
      n1   = [n1v]
      int2 = [int2v]
      n2   = [n2v]
      eos_band  = [stra]
      uars_band = [strb]
      no_linesv = 1
      WHILE NOT EOF(lu2) AND sps1 EQ '' DO BEGIN
        proc_tex_stra,lu2,buff,stra,strb,sps1
        IF sps1 EQ '' THEN BEGIN
          READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
          freq = [freq,freqv]
          gse  = [gse,gsev]
          istr = [istr,istrv]
          wc   = [wc,wcv]
          nc   = [nc,ncv]
          ps   = [ps,psv]
          ns   = [ns,nsv]
          int1 = [int1,int1v]
          n1   = [n1,n1v]
          int2 = [int2,int2v]
          n2   = [n2,n2v]
          eos_band  = [eos_band,stra]
          uars_band = [uars_band,strb]
          no_linesv = no_linesv + 1
        ENDIF
      ENDWHILE
      no_lines = [no_linesv]
    ENDIF ELSE BEGIN
      PRINT,'species ',species(0),' not found will try the next one'
      no_sps = no_sps - 1
    ENDELSE
    FOR i = 1 , no_sps - 1 DO BEGIN
      POINT_LUN,lu1,0
      POINT_LUN,lu2,0
      proc_tex_string,lu1,spsv,buff
      WHILE NOT EOF(lu1) AND spsv NE species(i) DO $
              proc_tex_string,lu1,spsv,buff
      IF spsv EQ species(i) THEN BEGIN
; species is found
        READS,buff,abunv,massv,qv,contv
        sps  = [sps,spsv]
        abun = [abun,abunv]
        mass = [mass,massv]
        parf = [parf,qv]
        conf = [conf,contv]
        sps1 = ''
        WHILE spsv NE sps1 DO BEGIN
          READF,lu2,buff
          WHILE STRPOS(buff,'\dotfill') EQ -1 DO READF,lu2,buff
          amp_pos = STRSPLIT(buff,'\\dotfill',LENGTH = len,/REGEX)
          sps1 = STRTRIM(STRMID(buff,amp_pos(1),len(1)),2)
        ENDWHILE
        proc_tex_stra,lu2,buff,stra,strb,sps1
        READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
        freq = [freq,freqv]
        gse  = [gse,gsev]
        istr = [istr,istrv]
        wc   = [wc,wcv]
        nc   = [nc,ncv]
        ps   = [ps,psv]
        ns   = [ns,nsv]
        int1 = [int1,int1v]
        n1   = [n1,n1v]
        int2 = [int2,int2v]
        n2   = [n2,n2v]
        eos_band  = [eos_band,stra]
        uars_band = [uars_band,strb]
        no_linesv = 1
        WHILE NOT EOF(lu2) AND sps1 EQ '' DO BEGIN
          proc_tex_stra,lu2,buff,stra,strb,sps1
          IF sps1 EQ '' THEN BEGIN
            READS,buff,freqv,gsev,istrv,wcv,ncv,psv,nsv,int1v,n1v,int2v,n2v
            freq = [freq,freqv]
            gse  = [gse,gsev]
            istr = [istr,istrv]
            wc   = [wc,wcv]
            nc   = [nc,ncv]
            ps   = [ps,psv]
            ns   = [ns,nsv]
            int1 = [int1,int1v]
            n1   = [n1,n1v]
            int2 = [int2,int2v]
            n2   = [n2,n2v]
            eos_band  = [eos_band,stra]
            uars_band = [uars_band,strb]
            no_linesv = no_linesv + 1
          ENDIF
        ENDWHILE
        no_lines = [no_lines,no_linesv]
      ENDIF ELSE BEGIN
        PRINT,'species ',species(i),' not found will try the next one'
        no_sps = no_sps - 1
      ENDELSE
    ENDFOR
  ENDELSE
  FREE_LUN,lu1
  FREE_LUN,lu2
; package into a nice structure
  IF no_sps GT 0 THEN BEGIN
    max_no_lines = MAX(no_lines)
    data = {name:'',abun:0.0, mass:0.0, q:FLTARR(3), cont:FLTARR(6), $
            no_lines:0, frq:DBLARR(max_no_lines), $
            gse:FLTARR(max_no_lines), ist:FLTARR(max_no_lines), $
            wc:FLTARR(max_no_lines), nc:FLTARR(max_no_lines), $
            ps:FLTARR(max_no_lines), ns:FLTARR(max_no_lines), $
            int1:FLTARR(max_no_lines), n1:FLTARR(max_no_lines), $
            int2:FLTARR(max_no_lines), n2:FLTARR(max_no_lines), $
            eos_bands:STRARR(max_no_lines), uars_bands:STRARR(max_no_lines)}
    sps_data = REPLICATE(data,no_sps)
    sps_data(*).name = sps
    sps_data(*).abun = abun
    sps_data(*).mass = mass
    sps_data(*).q(0:2) = REFORM(parf,3,no_sps)
    sps_data(*).cont(0:5) = REFORM(conf,6,no_sps)
    sps_data(*).no_lines = no_lines
    ptr1 = 0
    FOR i = 0 , no_sps - 1 DO BEGIN
      IF no_lines(i) GT 0 THEN BEGIN
        ptr2 = ptr1 + no_lines(i) - 1
        sps_data(i).frq        = freq(ptr1:ptr2)
        sps_data(i).gse        = gse(ptr1:ptr2)
        sps_data(i).ist        = istr(ptr1:ptr2)
        sps_data(i).wc         = wc(ptr1:ptr2)
        sps_data(i).nc         = nc(ptr1:ptr2)
        sps_data(i).ps         = ps(ptr1:ptr2)
        sps_data(i).ns         = ns(ptr1:ptr2)
        sps_data(i).int1       = int1(ptr1:ptr2)
        sps_data(i).n1         = n1(ptr1:ptr2)
        sps_data(i).int2       = int2(ptr1:ptr2)
        sps_data(i).n2         = n2(ptr1:ptr2)
        sps_data(i).eos_bands  = eos_band(ptr1:ptr2)
        sps_data(i).uars_bands = uars_band(ptr1:ptr2)
        ptr1 = ptr2 + 1
      ENDIF
   ENDFOR
  ENDIF
  RETURN
  END
