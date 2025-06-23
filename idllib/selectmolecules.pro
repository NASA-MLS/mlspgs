; $Id$

FUNCTION selectMolecules, thresh, data = data, molName=molName, lineName=lineName, moleculeFile=moleculeFile, fastRead=fastRead


;; for standalone function with default database lists invoke as:
;; moleculeList=selectmolecules(thresh)

;; if the data structure is already available then invoke as:
;; moleculeList=selectmolecules(thresh, data=data)

IF N_ELEMENTS(data) EQ 0 THEN $
  BEGIN
    if n_elements(molName) eq 0 then $
      molName = getenv('HOME') + '/mlspgs/tables/mol_data_table.tex'
    if n_elements(lineName) eq 0 then $
      lineName = getenv('HOME') + '/mlspgs/tables/line_data_table.tex'
    
    if keyword_set ( fastRead ) then begin
        lineName = getenv('HOME') +  '/mlspgs/tables/line_data_table.sav'
        restore, lineName
        data = sps_data
        lineID = 'Restored from saveset'
        readID = lineID
        molID = lineID
    endif else begin
        ;; Call Bills code to read the files
        Read_Spect_Dbase, molName, lineName, '', data,  $
          molID=molID, lineID=lineID, readID=readID
    endelse
ENDIF

if n_elements(moleculeFile) eq 0 then $
  moleculeFile = getenv('HOME') + '/mlspgs/tables/molecule_data_base.txt'

read_molecule_file,moleculeFile,table_of_molecules,spectags,spectls, $
  hitran,hitran_iso,hitran_frac,mol_mr,dflt_wc,dflt_nc

  idx2keep = [0]
  mol2keep = ['']
  alwaysKeepMatch = ['N$_{2}$', 'O$_{2}$', 'EXTINCTION','EXTINCTIONV2']
  alwaysKeepPos = ['H$_{2}$O-r', 'O$_{3}$-r', 'CLOUD\_']

  t = 240.0
  del_s = 500.0
  w_c = 3.0
  n_c = 0.7
  qtp = ALOG10([300.0,225.0,150.0])
 
  FOR i = 0, N_ELEMENTS(data) - 1 DO $
    BEGIN
      mol_ind = WHERE(table_of_molecules EQ data(i).name,no)

      IF TOTAL(STRMATCH(alwaysKeepMatch, data(i).name)) NE 0 $
        OR TOTAL(STRMATCH(alwaysKeepPos, (STRSPLIT(data(i).name,'\_',/REGEX,/EXTRACT))[0]+'\_')) NE 0 $
                          OR TOTAL(STRMATCH(alwaysKeepPos, (STRSPLIT(data(i).name,'-r',/REGEX,/EXTRACT))[0]+'-r')) NE 0 THEN $

        BEGIN
          mol2keep = [mol2keep, data(i).name]
          idx2keep = [idx2keep, i]
      END ELSE IF (no EQ 1 AND data(i).no_lines GT 0.0) THEN $
        BEGIN
          qlg = ALOG10(data(i).q)
          lt_grid = ALOG10(t)
          qrat = (qlg(1)-qlg(0)) * (lt_grid-qtp(0)) * (lt_grid GT qtp(1)) $
            / (qtp(1)-qtp(0)) + (qlg(1)-qlg(0) + (qlg(2)-qlg(1)) $
                                 * (lt_grid-qtp(1))/(qtp(2)-qtp(1))) * (lt_grid LE qtp(1))
;; a line center strength
          inds = INDGEN(data(i).no_lines)
          
          zip = WHERE(data(i).wc(inds) EQ 0.0,nozip)
          IF nozip GT 0 THEN BEGIN
              data(i).wc(zip) = w_c
              MESSAGE,/CONTINUE,'Overwriting zero WC data with default for:'+ data(i).name
          ENDIF
          
          zip = WHERE(data(i).nc(inds) EQ 0.0,nozip)
          IF nozip GT 0 THEN BEGIN
              data(i).nc(zip) = n_c 
              MESSAGE,/CONTINUE,'Overwriting zero NC data with default for:'+ data(i).name
          ENDIF
          
          beta = 2.30549e09 * (t/300.0)^data(i).nc(inds) $
            * 10.0^(data(i).ist(inds) - qrat $
                    + data(i).gse(inds)*(1.0/300.0-1.0/t)/1.600386) $
            * (1.0 - EXP(-data(i).frq(inds)/(20836.7*t))) $
            / ((1.0 - EXP(-data(i).frq(inds)/(20836.7*300.0))) * t $
               * data(i).wc(inds))

          bt = t * scalarize(mol_mr(mol_ind)) * beta * del_s
          
          IF MAX(bt) GE thresh THEN $ 
            BEGIN
              mol2keep = [mol2keep, data(i).name]
              idx2keep = [idx2keep, i]
              msg =  STRING('Accepted:', i, data(i).name, 'Max.bt =', MAX(bt), FORMAT='(A15,1X,I5,1X,A30,2X,A,1X,G10.4)')
              MESSAGE,/CONTINUE, msg
          ENDIF ELSE $
            BEGIN
              msg =  STRING('Rejected>>>:', i, data(i).name, 'Max.bt =', MAX(bt), FORMAT='(A15,1X,I5,1X,A30,2X,A,1X,G10.4)')
              MESSAGE,/CONTINUE, msg
          ENDELSE
          
      ENDIF ELSE $
        BEGIN
          msg =  STRING('Skipping>>>>>>:', i, data(i).name, 'No. of lines =', data(i).no_lines , FORMAT='(A15,1X,I5,1X,A30,2X,A,1X,I5)')          
          MESSAGE,/CONTINUE, msg
      ENDELSE
  ENDFOR
  idx2keep = idx2keep[1:*]
  mol2keep = mol2keep[1:*]
;; only return the species above the thresold
  RETURN, data[idx2keep]
END

;$Log$
;Revision 1.1  2006/04/27 23:03:20  lambert
;called by createspectroscopyl2cf.pro
;
