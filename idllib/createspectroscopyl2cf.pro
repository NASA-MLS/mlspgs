; $Id$

pro CreateSpectroscopyL2CF, $
                            molName=molName, $
                            lineName=lineName, $
                            crossRefName=crossRefName, $
                            instrument=instrument, $
                            fastRead=fastRead, $ ;; use fastRead for testing purposes only (otherwise cvs id info will be incorrect)
                            outName=outName, $
                            defsOnlyName=defsOnlyName, $
                            threshold=threshold, $
                            l2cfForPFAFwm=l2cfForPFAFwm
; 1. make full spectroscopy
; CreateSpectroscopyL2CF,molName='$HOME/mlspgs/tables/mol_data_table.tex',lineName=/users/bill/mlspgs/tables/line_data_table.tex',crossrefname='$HOME/mlspgs/tables/sps_cross_ref_table.txt',outname='$HOME/tables/MLS-Aura_L2Cal-Spectroscopy_v3-0-0_0000d000.l2cf'
; 2 make file with all PFA molecules stripped out (assuming default
; locations for standard files.
; CreateSpectroscopyL2CF,outname='$HOME/mlspgs/tables/MLS-Aura_L2Cal-Spectroscopy-PFA_v3-0-0_0000d000.l2cf',/l2cfForPFAFwm


if n_elements(instrument) eq 0 then instrument='emls'

instrumentPrefixes=['emlsSignals=[ ', 'umlsSignals=[ ', 'xptl1Signals=[ ']
noInstruments = n_elements(instrumentPrefixes)

if n_elements(outName) eq 0 then outName = 'spectroscopy.l2cf'

if n_elements(defsOnlyName) eq 0 then defsOnlyName = 'spectroscopy-defs.l2cf'

if n_elements(molName) eq 0 then $
  molName = getenv('HOME') + '/mlspgs/tables/mol_data_table.tex'
if n_elements(lineName) eq 0 then $
  lineName = getenv('HOME') + '/mlspgs/tables/line_data_table.tex'
if n_elements(crossRefName) eq 0 then $
  crossRefName = getenv('HOME') + '/mlspgs/tables/sps_cross_ref_table.txt'

if keyword_set ( fastRead ) then begin
    IF STRLEN(BYTE(fastread)) GT 1 THEN lineName = fastread ELSE $
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

;; set default threshold for inclusion of molecule (kelvin)
if N_ELEMENTS(threshold) EQ 0 THEN threshold=0.0

if N_ELEMENTS(l2cfForPFAFwm) EQ 0 THEN l2cfForPFAFwm = 0

l2cfSpectroscopyTableType = (['Full Spectroscopy', 'Reduced Spectroscopy for PFA Forward Models'])[l2cfForPFAFwm]


if instrument eq 'emls' then begin
    
;; define list of mandatory lbl molecules for lower and upper sidebands in
;; each band for use in the PFA forward models
    
    SIDSMandatoryLBLMoleculesBands = {  $
                                       B1F:{lsb:['O2', 'O_18_O', 'O_17_O', 'O2_V1'], usb:'NONE'},$
                                       B22D:{lsb:'O2', usb:'NONE'},$
                                       B32W:{lsb:'NONE', usb:'NONE'},$
                                       B21F:{lsb:['O2', 'O_18_O', 'O_17_O', 'O2_V1'], usb:'NONE'},$
                                       B26D:{lsb:'O2', usb:'NONE'},$
                                       B34W:{lsb:'NONE', usb:'NONE'},$
                                       B2F:{lsb:'H2O', usb:'NONE'},$
                                       B3F:{lsb:'H2O', usb:'NONE'},$
                                       B4F:{lsb:'NONE', usb:'NONE'},$
                                       B5F:{lsb:'NONE', usb:'NONE'},$
                                       B6F:{lsb:'NONE', usb:'NONE'},$
                                       B23D:{lsb:'H2O', usb:'NONE'},$
                                       B27M:{lsb:'NONE', usb:'NONE'},$
                                       B7F:{lsb:'O3', usb:['O3', 'O3_V2']},$
                                       B8F:{lsb:'O_18_O', usb:'NONE'},$
                                       B9F:{lsb:['CO', 'O3'], usb:'O3'},$
                                       B24D:{lsb:'O3', usb:'NONE'},$
                                       B25D:{lsb:'NONE', usb:'NONE'},$
                                       B33W:{lsb:'NONE', usb:'O3'},$
                                       B10F:{lsb:'NONE', usb:'NONE'},$
                                       B11F:{lsb:'NONE', usb:'O3'},$
                                       B12F:{lsb:'O3', usb:['N2O', 'O3', 'O3_V2']},$
                                       B13F:{lsb:['HCL_35', 'O3'], usb:'NONE'},$
                                       B14F:{lsb:['HCL_35', 'HCL_37', 'O3'], usb:'NONE'},$
                                       B28M:{lsb:'NONE', usb:'NONE'},$
                                       B29M:{lsb:'NONE', usb:'NONE'},$
                                       B30M:{lsb:'O3', usb:'NONE'},$
                                       B31M:{lsb:'NONE', usb:'NONE'},$         
                                       B15F:{lsb:'OH', usb:'NONE'},$
                                       B16F:{lsb:['OH', 'O3'], usb:'NONE'},$
                                       B17F:{lsb:['O2', 'O3_V2'], usb:'O3'},$                   
                                       B18F:{lsb:'OH', usb:'NONE'},$
                                       B19F:{lsb:['OH', 'O3'], usb:'NONE'},$
                                       B20F:{lsb:['O2', 'O3_V2'], usb:'O3'} }

    MoleculesAlwaysPresent = [ 'H2O', 'N2', 'O2', 'O3', 'EXTINCTION', 'EXTINCTIONV2' ]
    MoleculesR = ['H2O_', 'O3_']

;; create list of all mandatory LBL molecules regardless of radiometer/sideband
    AllLBLMolecules = MoleculesAlwaysPresent
    FOR ib = 0,N_TAGS(SIDSMandatoryLBLMoleculesBands)-1 DO $ 
      BEGIN
        AllLBLMolecules = [AllLBLMolecules, SIDSMandatoryLBLMoleculesBands.(ib).lsb, SIDSMandatoryLBLMoleculesBands.(ib).usb]
    ENDFOR;;ib

    ;; tidy up the list
    q = WHERE(AllLBLMolecules NE 'NONE', count)
    IF q[0] NE -1 AND count NE N_ELEMENTS(AllLBLMolecules) THEN AllLBLMolecules = AllLBLMolecules[q]                
    AllLBLMolecules = AllLBLMolecules[UNIQ(AllLBLMolecules, SORT(AllLBLMolecules))]
   
endif;;instrument EQ emls


;; scan for molecules above threshold and set data to contain just
;; those molecules
selectedData = selectMolecules(threshold, data=data, fastRead=fastRead )
data = selectedData


myID='$Id$'

; Read the cross reference information
openr, unit, crossRefName, /get_lun
line = ''
firstLine = 1
repeat begin
    readf, unit, line, format='(a)'
    if firstLine then crossRefID = line
    firstLine = 0
endrep until strmid(line,0,4) eq '----'
while not eof ( unit ) do begin
    readf, unit, line, format='(a)'
    words = strsplit ( line, ' ', /extract )
    if n_elements(words) gt 1 then begin
        if n_elements(texNames) eq 0 then begin
            texNames = words(0)
            l2Names = words(1)
            noChildren = fix ( words(3) )
        endif else begin
            texNames = [ texNames, words(0) ]
            l2Names = [ l2Names, words(1) ]
            noChildren = [ noChildren, fix ( words(3) ) ]
        endelse
    endif
  endwhile
free_lun, unit
; Now work out the parents
noCrossRefs = n_elements(texNames)
parents = intarr ( noCrossRefs ) - 1
for index = 0, noCrossRefs - 1 do begin
    if noChildren(index) ne 0 then begin
        if noChildren(index) < 0 then i0 = index else i0 = index + 1
        i1 = index + abs(noChildren(index))
        parents(i0:i1) = index
    endif
endfor  

; set dummy unity partition functions to prevent warning when taking logs 
change = where(data.name eq 'N$_{2}$' or data.name eq 'EXTINCTION' or data.name eq 'EXTINCTIONV2' $
               or strmid(data.name,0,10) eq 'H$_{2}$O-r' or strmid(data.name,0,9) eq 'O$_{3}$-r' $
               or strmid(data.name,0,7) eq 'CLOUD\_' )
data(change).q=1.0

; Open the output file
Openw, unit, outName, /get_lun

;; Open the defs only file
openw, defsUnit, defsOnlyName, /get_lun

;; Combine the units
units = [ unit, defsUnit ]

; Write in the appropriate header information
printf, unit, '; --------------------------------------------- Spectroscopy ---- '
printf, unit, '; ' + l2cfSpectroscopyTableType
printf, unit, '; This file is automatically generated by the program'
printf, unit, '; mlspgs/idllib/createspectroscopyl2cf.pro do NOT edit it by hand.'
printf, unit, STRING('; Created with extra line selection threshold on input database = ',threshold,'K', FORMAT='(A,G10.4,1X,A)')
printf, unit, ''


printf, unit, '; Source file information:'
printf, unit, '; '+strmid(molID,7,strlen(molID)-9)
printf, unit, '; '+strmid(lineID,7,strlen(lineID)-9)
printf, unit, '; '+strmid(crossRefID,5,strlen(crossRefID)-7)
printf, unit, '; '+strmid(readID,5,strlen(readID)-7)
printf, unit, '; '+strmid(myID,5,strlen(myID)-7)
printf, unit, ''

printf, unit, 'begin spectroscopy'

noEMLSBands = 34
noUMLSBands = 6
noXPTL1Bands = 1                ; XPTL1 doesn't list any.
noSidebands = 3

noMols = n_elements(data)
niceNames = strarr(noMols)
parentNames = niceNames
emlsMolecules = intarr ( noSidebands, noEMLSBands, noMols )
umlsMolecules = intarr ( noSidebands, noUMLSBands, noMols )
xptl1Molecules = intarr ( noSidebands, noXPTL1Bands, noMols )

for mol = 0, noMols - 1 do begin

    ;; First try to un TeXify the molecule names
    ;; Also, after a space, add the name of the parent
    ;; Did it after the space to avoid the case statement getting overly large
    index = ( where ( data(mol).name eq texNames ) ) ( 0 )
    if index eq -1 then begin
        if strmid ( data(mol).name, 0, 10 ) eq 'H$_{2}$O-r' then begin
            niceNames(mol) = 'H2O_R' + strupcase ( strmid ( data(mol).name, 10, 2 ) )
            parentNames(mol) = 'H2O'
        endif else if strmid ( data(mol).name, 0, 9 ) eq 'O$_{3}$-r' then begin
            niceNames(mol) = 'O3_R' + strupcase ( strmid ( data(mol).name, 9, 2 ) )
            parentNames(mol) = 'O3'
        endif else begin
            MyMessage, /error, "Unrecognized molecule: " + data(mol).name
        endelse
    endif else begin
        niceNames(mol) = l2Names(index)
        if parents(index) eq -1 then begin
            parentNames(mol) = niceNames ( mol )
        endif else begin
            parentNames(mol) = l2Names ( parents (index) )
        endelse
    endelse
    thisMolName = niceNames(mol)

    ;; Look out for molecules that are radiometer specific
    ;; Different algorithm for skipping lines on EOS vs. SMLS
    if instrument eq 'emls' then begin
        if ( min(data(mol).q) le 0.0 and data(mol).cont(0) eq 0.0 ) then begin
            print, 'Skipping ' + thisMolName
            niceNames ( mol ) = ''
            continue
        endif
    endif else begin
        if ( data(mol).no_lines eq 0 and data(mol).cont(0) eq 0.0 ) then begin
            print, 'Skipping ' + thisMolName + ' (not EOS)'
            niceNames ( mol ) = ''
            continue
        endif

        if strmid ( thisMolName, 0, 5 ) eq 'H2O_R' or $
          strmid ( thisMolName, 0, 4 ) eq 'O3_R' then begin
            print, 'Skipping radiometer based ' + thisMolName + ' as not emls'
            niceNames ( mol ) = ''
            continue
        endif
    endelse

    if n_elements(allNames) eq 0 then allNames = thisMolName $
    else allNames = allNames + ', ' + thisMolName
    
    ;; Now print out a header line
    printf,unit, ''
    printf,unit, '  ;; ----------------------------------- '+thisMolName
    if thisMolName ne 'EXTINCTION' and thisMolName ne 'EXTINCTIONV2' then begin
        isotopeLine = '!define(isotoperatio' + thisMolName + ',{' + $
          string ( data(mol).abun, format='(f10.8)' ) + '} )'
        printf, unit, '  ' + isotopeLine
        printf, defsUnit, isotopeLine
    endif
    printf,unit, ''

    ;; Lines
    ;; For l2cf PFA output we omit the line information
    IF (l2cfForPFAFwm EQ 0) OR ((l2cfForPFAFwm EQ 1) AND ((WHERE(thisMolName EQ AllLBLMolecules))[0] NE -1)) THEN $
    for line = 0, data(mol).no_Lines - 1 do begin
        text= '  '+thisMolName + '_L' + string ( line+1, format='(I0)' ) + ': line, '
        AddWordToLine,text,unit,4, $
          'v0= '+strtrim(string( data(mol).frq(line), format='(f12.4)' ),2) + ' MHz, '
        AddWordToLine,text,unit,4, $
          'el= '+strtrim(string( data(mol).gse(line), format='(g13.8)' ),2) + ', '
        AddWordToLine,text,unit,4, $
          'str= '+strtrim(string( data(mol).ist(line), format='(f8.4)' ),2) + ', '
        AddWordToLine,text,unit,4, $
          'w= '+strtrim(string( data(mol).wc(line), format='(f8.4)' ),2) + ', '
        AddWordToLine,text,unit,4, $
          'n= '+strtrim(string(data(mol).nc(line), format='(f6.3)'),2) + ', '
        AddWordToLine,text,unit,4, $
          'ps= '+strtrim(string(data(mol).ps(line), format='(g15.5)'),2) + ', '
        AddWordToLine,text,unit,4, $
          'ns= '+strtrim(string(data(mol).ns(line), format='(f6.3)'),2) + ', '
        AddWordToLine,text,unit,4, $
          'delta= '+strtrim(string(data(mol).int1(line), format='(g15.5)'),2) + ', '
        AddWordToLine,text,unit,4, $
          'n1= '+strtrim(string(data(mol).n1(line), format='(f6.3)'),2) + ', '
        AddWordToLine,text,unit,4, $
          'gamma= '+strtrim(string(data(mol).int2(line), format='(g15.5)'),2) + ', '
        AddWordToLine,text,unit,4, $
          'n2= '+strtrim(string(data(mol).n2(line), format='(f6.3)'),2) + ', '

        ;; Do the quantum numbers
        noQNs = data(mol).qnfmt(line) mod 10
        qns = data(mol).qnfmt(line)
        for q = 0, noQNs - 1 do qns = [ qns, fix ( strmid ( data(mol).qnu(line), q*2, 2 ) ) ]
        for q = 0, noQNs - 1 do qns = [ qns, fix ( strmid ( data(mol).qnl(line), q*2, 2 ) ) ]
        AddWordToLine,text,unit,4, 'qn=[ '
        for q = 0, n_elements(qns)-2 do $
          AddWordToLine,text,unit,4, string(qns(q),format='(i0)')+', '
        AddWordToLine,text,unit,3, string(qns(2*noQNS),format='(i0)') + ' ]'
        
        for i=0, noInstruments - 1 do begin
            case i of
                0 : begin 
                    prefix = 'emlsSignals=[ '
                    bands = data(mol).eos_bands(line)
                end
                1 : begin
                    prefix = 'umlsSignals=[ '
                    bands = data(mol).uars_bands(line)
                end
                2 : begin
                    prefix = 'xptl1Signals=[ '
                    bands = ''  ; data(mol).xptl1_bands(line)
                end
            endcase
            bands = strtrim(strsplit(bands,',',/extract),2)
            if bands(0) ne '' then begin
                text = text + ', '
                AddWordToLine,text,unit,4,prefix
                for j = 0, n_elements(bands)-1 do begin
                    thisBand = strupcase ( bands(j) )
                    bandNo = fix(strmid(thisBand,1,strlen(thisBand)))
                    case 1 of
                        strpos(thisBand,'U') ne -1 : sideband = 1
                        strpos(thisBand,'L') ne -1 : sideband = -1
                        else : sideband = 0
                    endcase
                    case 1 of 
                        bandNo le 21 : thisBand=thisBand+'F'
                        bandNo le 26 : thisBand=thisBand+'D'
                        bandNo le 31 : thisBand=thisBand+'M'
                        bandNo le 34 : thisBand=thisBand+'W'
                    endcase
                    thisBand="'"+thisBand+"'"
                    if j lt n_elements(bands)-1 $
                      then thisBand=thisBand+', ' $
                    else thisBand=thisBand+' ]'
                    AddWordToLine,text,unit,4,thisBand
                    ;; Keep track of the molecule by molecule usage
                    case i of
                        0 : emlsMolecules ( sideband+1, bandNo-1, mol ) = 1
                        1 : umlsMolecules ( sideband+1, bandNo-1, mol ) = 1
                        2 : xptl1Molecules ( sideband+1, bandNo-1, mol ) = 1
                    endcase
                endfor
            endif               ; Any bands
        endfor                  ; Loop over instruments

        ;; Special case for O2 polarized 118 line
        if thisMolName eq 'O2' and n_elements(qns) eq 5 and $
          min ( qns eq [ 102, 1, 1, 1, 0 ] ) then begin
            text = text + ', '
            AddWordToLine, text, unit, 4, "emlsSignalsPol=[ 'B22LD', " 
            AddWordToLine, text, unit, 4, "'B26LD' ]"
        endif
        
        if strtrim(text,2) ne '' then printf,unit,text
    endfor                      ; Loop over lines



    ;; Now print out the molecule information
    printf, unit, ''
    text = '  spectra, molecule= '+thisMolName+', Qlog=[ '

    for i = 0, 2 do begin
        text = text + strtrim(string(alog10(data(mol).q(i)), $
                                     format='(f10.4)'),2)
        if i ne 2 then text = text + ', '
    endfor

    ;; Lines
    ;; For l2cf PFA output we omit the line information
    if (data(mol).no_lines gt 0) AND ((l2cfForPFAFwm EQ 0) OR ((l2cfForPFAFwm EQ 1) AND ((WHERE(thisMolName EQ AllLBLMolecules))[0] NE -1))) then $
      begin
        text = text + ' ], $'
        printf, unit,text
        text = '    lines=[ '
        for line = 0, data(mol).no_lines - 2 do begin
            AddWordToLine, text, unit, 4, $
              thisMolName+'_L'+string(line+1,format='(I0)')+', '
        endfor
        AddWordToLine, text, unit, 4, $
          thisMolName+'_L'+string(data(mol).no_lines,format='(I0)')+' ]'
    endif else begin
        text = text + ' ]'
    endelse
    
    ;; Mass
    text = text + ', '
    AddWordToLine, text, unit, 4, $
      'mass= ' + strtrim ( string ( data(mol).mass, format='(f9.5)' ), 2 )
    
    ;; Continuum
    noNonZero = max(where(data(mol).cont ne 0.0))+1
    if noNonZero ne 0 then begin
        text = text + ', '
        AddWordToLine, text, unit, 4, $
          'continuum=[ '
        for i=0,noNonZero-2 do begin
            AddWordToLine, text, unit, 4, $
              strtrim(string(data(mol).cont(i)),2)+', '
        endfor
        AddWordToLine, text, unit, 4, $
          strtrim(string(data(mol).cont(noNonZero-1)),2)+' ]'
    endif
    
    ;; Finish off the line
    if strtrim(text,2) ne '' then printf,unit,text
endfor;;loop over molecules

printf, unit, ''
printf, unit, 'end spectroscopy'
printf, unit, ''

;; ------------------------------------- End of the spectroscopy part

;; Now write out the molecules per band stuff
case strlowcase ( instrument ) of 
    'umls'  : array=umlsMolecules
    'emls'  : array=emlsMolecules
    'xptl1' : array=xptl1Molecules
    else : MyMessage,/error,'Unknown instrument'
endcase

;; Collapse the sidebands
array(1,*,*) = array(0,*,*) or array(1,*,*) or array(2,*,*) 
printf, unit, '; ----------------------------------------- Molecules per band'
printf, defsUnit, '; ----------------------------------------- Molecules per band'
printf, unit, ''
printf, defsUnit, ''
;; Make H2O, N2, O2, O3 and extinction used everywhere.
array(*,*,where(niceNames eq 'H2O')) = 1
array(*,*,where(niceNames eq 'N2')) = 1
array(*,*,where(niceNames eq 'O2')) = 1
array(*,*,where(niceNames eq 'O3')) = 1
array(*,*,where(niceNames eq 'EXTINCTION')) = 1
array(*,*,where(niceNames eq 'EXTINCTIONV2')) = 1
sidebandStrings = ['L','','U']
noBands = (size ( array ) ) (2)

;; Get the instrument configuration
database = ReadValidSignals()

;; Make the H2O_<radiometer> species used where appropriate
if instrument eq 'emls' then begin
    for radiometer = 0, database.noRadiometers-1 do begin
        relevantBands = where ( database.bands.radiometerIndex eq radiometer )
        array(*,relevantBands, $
              where(niceNames eq 'H2O_'+database.radiometers(radiometer).prefix ) ) = 1
    endfor
    ;; Make the O3_<radiometer> species used where appropriate
    for radiometer = 0, database.noRadiometers-1 do begin
        relevantBands = where ( database.bands.radiometerIndex eq radiometer )
        o3rs = where(niceNames eq 'O3_'+database.radiometers(radiometer).prefix )
        if o3rs(0) ne -1 then array(*,relevantBands, o3rs ) = 1
    endfor
endif

;; Now, add some extra virtual bands which are the radiometers and one
;; at the end which is the whole instrument (folded only)
newArray = intarr ( noSidebands, noBands+database.noRadiometers+1, noMols )
newArray(*,0:noBands-1,*) = array
for radiometer = 0, database.noRadiometers-1 do begin
    relevantBands = where ( database.bands.radiometerIndex eq radiometer )
    collapsed = total ( array (*,relevantBands,*), 2 )
    newArray(*,noBands+radiometer,*) = collapsed gt 0
endfor
;; Do the whole instrument
newArray(1,noBands+database.noRadiometers,*) = niceNames ne ''

;; Use this new array
array=newArray

if instrument eq 'emls' then begin

;; now create the lsb,usb lists for the radiometers 

    FOR radiometer = 0, database.noRadiometers-1 DO $
      BEGIN
        
        relevantBands = where ( database.bands.radiometerIndex eq radiometer )
        RadiometerName = (database.bands(relevantBands).radiometer)[0]
        prefixList = database.bands(relevantBands).prefix

        FOR sideband = 0,2,2 DO $ ;; 'L' and 'U' 
          BEGIN
            LBLMolecules = ['NONE']

            FOR ip = 0,N_ELEMENTS(PrefixList)-1 DO $ ;; loop over bands in this radiometer
              BEGIN
                prefix = PrefixList(ip)
                ib = WHERE(prefix EQ TAG_NAMES(SIDSMandatoryLBLMoleculesBands))
                IF ib(0) EQ -1 THEN MESSAGE, 'prefix=' + prefix + ' not found in SIDSMandatoryLBLMoleculesBands'                        
                
                SIDSMandatoryLBLMolecules = (sideband EQ 0) ? SIDSMandatoryLBLMoleculesBands.(ib).lsb : SIDSMandatoryLBLMoleculesBands.(ib).usb
                
                LBLMolecules = [LBLMolecules,SIDSMandatoryLBLMolecules]                
            ENDFOR;;ip

            ;; tidy up the list
            q = WHERE(LBLMolecules NE 'NONE', count)
            IF q[0] NE -1 AND count NE N_ELEMENTS(LBLMolecules) THEN LBLMolecules = LBLMolecules[q]                
            LBLMolecules = LBLMolecules[UNIQ(LBLMolecules, SORT(LBLMolecules))]

            ;; make/append the lsb/usb structures
            w = (sideband EQ 0) ?  Create_Struct( 'lsb', LBLMolecules) : Create_Struct( w, 'usb', LBLMolecules)

        ENDFOR;;sideband

        ;; add lsb/usb structure for this radiometer
        SIDSMandatoryLBLMoleculesBands = Create_Struct( SIDSMandatoryLBLMoleculesBands, RadiometerName, w)

    END;;radiometer

endif;;instrument EQ emls

    for band = 0, noBands + database.noRadiometers do begin
        for sideband = 0, 2 do begin
            ;; Which molecules does it use
            usedMols = where ( reform ( array(sideband, band, *) ) )
            ;; What parents do they have, get a unique list
            if usedMols(0) ne -1 then begin
                usedParents = parentNames ( usedMols )
                usedParents = usedParents ( sort ( usedParents ) )
                usedParents = usedParents ( uniq ( usedParents ) )
                case 1 of
                    band lt noBands : outName = string ( band+1, format='(i0)' )
                    band ge noBands and band lt noBands + database.noRadiometers : $
                      outName = database.radiometers(band-noBands).prefix
                    else : outName = 'Instrument'
                endcase
                
                ;; First do the comprehensive lists for full forward models.
                line = '!define(moleculesFor' + outName + $
                  sidebandStrings(sideband)+',{[ '
                for p = 0, n_elements(usedParents)-1 do begin
                    if p ne 0 then line = line +', '
                    thisParent = usedParents(p)
                    if thisParent ne 'EXTINCTION' and thisParent ne 'EXTINCTIONV2' then begin
                        children = where ( parentNames eq thisParent and $
                                           reform ( array(sideband, band, *) ) )
                        AddWordToLine, line, units, 2, '[ '+thisParent
                        ;; Format is [ parent, <children> ]
                        ;; Note, having the only child be the same as the parent is fine
                        ;; (e.g. [ H2O, H2O ])
                        for c = 0, n_elements(children) - 1 do begin
                            line = line + ', '
                            AddWordToLine, line, units, 4, niceNames(children(c))
                        endfor
                        line = line + ' ]'
                    endif else begin
                        ;; For extinction just use extinction alone, no isotope
                        AddWordToLine, line, units, 2, thisParent
                    endelse
                endfor
                line = line + ' ]})'
                printf, unit, line
                printf, defsUnit, line

                if instrument eq 'emls' then begin

;; Now do the lbl molecule definitions per sideband for the PFA
;; forward models
                    
                    IF sideband NE 1 THEN $ ;; output only for individual 'L' or 'U' sidebands
                      BEGIN
                        
                        case 1 of
                            band lt noBands : BEGIN
                                prefix = database.bands(band).prefix  
                                RadiometerName = database.bands(band).radiometer
                            END
                            band ge noBands and band lt noBands + database.noRadiometers : BEGIN                                
                                prefix = database.radiometers(band-noBands).prefix 
                                RadiometerName = prefix
                            END
                            else : MESSAGE,'Incorrect band index:' + STRING(band)
                        endcase 

                        ib = WHERE(prefix EQ TAG_NAMES(SIDSMandatoryLBLMoleculesBands))
                        IF ib(0) EQ -1 THEN MESSAGE, 'prefix=' + prefix + ' not found in SIDSMandatoryLBLMoleculesBands'                        
                        
                        SIDSMandatoryLBLMolecules = (sideband EQ 0) ? SIDSMandatoryLBLMoleculesBands.(ib).lsb : SIDSMandatoryLBLMoleculesBands.(ib).usb
                        
                        IF TOTAL(SIDSMandatoryLBLMolecules EQ 'NONE') EQ 0 THEN $
                          BEGIN                            
                            RTVLMandatoryLBLMolecules = [SIDSMandatoryLBLMolecules, MoleculesAlwaysPresent, MoleculesR + RadiometerName]
                        ENDIF ELSE BEGIN 
                            RTVLMandatoryLBLMolecules = [MoleculesAlwaysPresent, MoleculesR + RadiometerName]
                        ENDELSE
                        
                        ;; sort and uniq
                        RTVLMandatoryLBLMolecules = RTVLMandatoryLBLMolecules[UNIQ(RTVLMandatoryLBLMolecules, SORT(RTVLMandatoryLBLMolecules))]

;; output the RTVL/l2pc case
                        line = '!define(RTVLMandatoryLBLMoleculesFor' + outName + $
                          sidebandStrings(sideband) + ',{ '
                        
                        doneAny = 0
                        for p = 0, n_elements(RTVLMandatoryLBLMolecules) - 1 do begin
                            thisMolecule = RTVLMandatoryLBLMolecules(p)
                            if doneAny then line = line + ', ' else doneAny = 1
                            AddWordToLine, line, units, 4, thisMolecule
                        endfor
                        
                        line = line + ' })'
                        printf, unit, line
                        printf, defsUnit, line
                        
                        
;; output the SIDS case
                        line = '!define(SIDSMandatoryLBLMoleculesFor' + outName + $
                          sidebandStrings(sideband) + ',{ '
                        
                        doneAny = 0
                        for p = 0, n_elements(SIDSMandatoryLBLMolecules) - 1 do begin
                            thisMolecule = SIDSMandatoryLBLMolecules(p)
                            if doneAny then line = line + ', ' else doneAny = 1
                            AddWordToLine, line, units, 4, thisMolecule
                        endfor
                        
                        line = line + ' })'
                        printf, unit, line
                        printf, defsUnit, line
                        
                    ENDIF;;sideband

                endif ;;instrument EQ emls
                
                
;; Now do just the parents for l2pc line forward models
;; First do the comprehensive lists for full forward models.
                line = '!define(moleculeFamiliesFor' + outName + $
                  sidebandStrings(sideband)+',{[ '
                for p = 0, n_elements(usedParents)-1 do begin
                    if p ne 0 then line = line + ', '
                    thisParent = usedParents(p)
                    children = where ( parentNames eq thisParent and $
                                       reform ( array(sideband, band, *) ) )
                    AddWordToLine, line, units, 2, thisParent
                endfor
                line = line +' ]})'
                printf, unit, line
                printf, defsUnit, line
            endif
            
        endfor                  ; Loop over sidebands
        printf, unit, ''
        printf, defsUnit, ''
    endfor                      ; Loop over bands

;; Print out a few more macro definitions
isotopicMolecules = niceNames ( where ( nicenames ne 'EXTINCTION' and $
                                        nicenames ne 'EXTINCTIONV2' and $
                                        nicenames ne '' ) )

units = [ unit, defsUnit ]
for i = 0, 1 do begin
    printf,units(i),'!define(molecules,{!moleculesFor$1})'
    printf,units(i),'!define(moleculeFamilies,{!moleculeFamiliesFor$1})'
;    printf,units(i),'!define(molecules2X,{!molecules2XFor$1})'
;    printf,units(i),'!define(moleculeFamilies2X,{!moleculeFamilies2XFor$1})'
    printf,units(i),'!define(isotopicMolecules,{' + strjoin(isotopicMolecules,',')+'})'
endfor

Free_lun, unit
Free_lun, defsUnit

end
