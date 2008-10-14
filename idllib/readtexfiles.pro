function ReadTeXFiles, $
  moleculeFile=moleculeFile, $
  lineFile=lineFile, $
  crossRefFile=crossRefFile

maxLines = 4000
;; This is a new IDL reader for Bill's tex files.
;; First read the cross reference information


;; ------------------------------------------------------------------
;; Read the cross reference information
openr, unit, crossRefFile, /get_lun
line = ''
firstLine = 1
repeat begin
  readf, unit, line, format='(a)'
  if firstLine then crossRefFileID = line
  firstLine = 0
endrep until strmid(line,0,4) eq '----'
while not eof ( unit ) do begin
  readf, unit, line, format='(a)'
  words = strsplit ( line, ' ', /extract )
  if n_elements(words) gt 1 then begin
    if n_elements(texNames) eq 0 then begin
      texNames = words[0]
      l2Names = words[1]
      hdfNames = words[2]
      noChildren = fix ( words[3] )
    endif else begin
      texNames = [ texNames, words[0] ]
      l2Names = [ l2Names, words[1] ]
      hdfNames = [ hdfNames, words[2] ]
      noChildren = [ noChildren, fix ( words[3] ) ]
    endelse
  endif
endwhile
free_lun, unit
;; Now work out the parents
noCrossRefs = n_elements(texNames)
parents = intarr ( noCrossRefs ) - 1
for index = 0, noCrossRefs - 1 do begin
  if noChildren(index) ne 0 then begin
    if noChildren(index) < 0 then i0 = index else i0 = index + 1
    i1 = index + abs(noChildren(index))
    parents(i0:i1) = index
  endif
endfor  

;; Now bundle all this up
crossRefs = { $
  noCrossRefs: noCrossRefs, $
  texNames: texNames, $
  l2Names: l2Names, $
  hdfNames: hdfNames, $
  noChildren: noChildren, $
  parents: parents }

;; ------------------------------------------------------------------
;; Read the molecule file information
openr, unit, moleculeFile, /get_lun
line = ''
readf, unit, line, format='(a)' 
molFileID = line
while not eof ( unit ) do begin
  readf, unit, line, format='(a)'
  line = strtrim ( line, 2 )
  if strmid ( line, 0, 1 ) eq '%' then continue
  words = strtrim ( strsplit ( line, '&', /extract ), 2 )
  if n_elements ( words ) ne 12 then begin
    print, 'Wrong number of words for ' + words[0]
    print, line
    stop
  end
  ;; Do some special stuff to parse any texisms in the 'cont' stuff
  contWords = words [ 6 : 11 ]
  keyWord = '$\times 10^{'
  keyLen = strlen ( keyWord )
  for i = 0, n_elements ( contWords ) - 1 do begin
    pos = strpos ( contWords[i], keyWord )
    if pos ne -1 then begin
      contWords[i] = strmid ( contWords[i], 0, pos ) + $
        'e' + strmid ( contWords[i], pos+keyLen, strlen ( contWords[i] ) - pos - keyLen - 2 )
    endif
  endfor
  thisMolecule = { $
    texName: words[0], $
    l2Name: '', $
    hdfName: '', $
    abundance: float ( words[1] ), $
    mass: float ( words[2] ), $
    q: float ( words[3:5] ), $
    cont: float ( contWords ), $
    noLines: 0l }
  ;; Note the last word includes the final ' \\' but the float function
  ;; ignores it
  if n_elements ( molecules ) eq 0 then begin
    molecules = thisMolecule
  endif else begin
    molecules = [ molecules, thisMolecule ]
  endelse
endwhile
free_lun, unit

;; Find the l2Names for each of these molecules
for i = 0, n_elements ( molecules ) - 1 do begin
  texName = molecules[i].texName
  radiometerBased = 0
  words = strsplit ( texName, '-', /extract )
  if n_elements ( words ) eq 2 then begin
    if strmid ( words[1], 0, 1 ) ne 'r' then begin
      print, 'Error: Expected r after - in cross reference file'
      stop
    endif
    texName = words[0]
    radiometerBased = 1
  endif
  nameInds = where ( texName eq crossRefs.texNames )
  if nameInds[0] eq -1 or n_elements ( nameInds ) ne 1 then begin
    print, 'Unable to find: ' + molecules[i].texName
    stop
  endif
  molecules[i].l2Name = crossRefs.l2Names [ nameInds[0] ]
  molecules[i].hdfName = crossRefs.hdfNames [ nameInds[0] ]
  if radiometerBased then begin
    molecules[i].l2Name += '_' + words[1]
    molecules[i].hdfName += '-' + words[1]
  endif
endfor

;; ------------------------------------------------------------------
;; Read the line file information
openr, unit, lineFile, /get_lun
line = ''
readf, unit, line, format='(a)'
lineFileID = line
;; Skip the next line too
readf, unit, line

;; All new molecules begin with this
molPrefix = '\multicolumn{17}{c}{\dotfill '

;; Loop over molecules
readf, unit, line, format='(a)' ;; Do the first read 
while not eof ( unit ) do begin
  ;; Got 'line' from the last time round (or above for first time round)
  if strmid ( line, 0, strlen ( molPrefix ) ) ne molPrefix then begin
    print, 'Error: expected line to begin with'
    print, molPrefix
    print, '  but got:'
    print, line
    stop
  endif
  suffix = strmid ( line, strlen ( molPrefix ), strlen ( line ) )
  words = strsplit ( suffix, ' ', /extract )
  texName = words[0]
  ;; Loop over lines
  noLines = 0
  doneThisMolecule = 0
  repeat begin
    if not eof ( unit ) then begin
      readf, unit, line
      doneThisMolecule = strmid ( line, 0, strlen ( molPrefix ) ) eq molPrefix
    endif else begin
      doneThisMolecule = 1
    endelse

    ;; Either finish this molecule or add a new line to it
    if doneThisMolecule then begin
      ;; We're done with this molecule, assemble the information
      ;; First identify the l2Name for this
      nameInds = where ( texName eq crossRefs.texNames )
      if nameInds[0] eq -1 or n_elements ( nameInds ) ne 1 then begin
        print, 'Unable to find: ' + texName
        stop
      endif
      l2Name = crossRefs.l2Names [ nameInds[0] ]
      if n_elements ( lines ) eq 0 then begin
        lines = create_struct ( l2Name, theseLines[0:noLines-1] )
      endif else begin
        lines = create_struct ( lines, l2Name, theseLines[0:noLines-1] )
      endelse
    endif else begin
      ;; Otherwise, parse this line
      words = strsplit ( line, '&', /extract )
      if n_elements ( words ) ne 17 then begin
        print, 'Error: wrong number of words for'
        print, line
        stop
      end
      thisLine = { $
        TOneLine, $
        freq:  double ( words [ 0  ] ), $
        gse:   float ( words [ 1  ] ), $
        istr:  float ( words [ 2  ] ), $
        wc:    float ( words [ 3  ] ), $
        nc:    float ( words [ 4  ] ), $
        ps:    float ( words [ 5  ] ), $
        ns:    float ( words [ 6  ] ), $
        int1:  float ( words [ 7  ] ), $
        n1:    float ( words [ 8  ] ), $
        int2:  float ( words [ 9  ] ), $
        n2:    float ( words [ 10 ] ), $
        qnfmt: long  ( words [ 11 ] ), $
        qnu:           words [ 12 ], $
        qnl:           words [ 13 ], $
        uarsBands: strtrim ( words [ 14 ], 2 ), $
        eosBands:  strtrim ( words [ 15 ], 2 ), $
        ref:           words [ 16 ] }
      if n_elements ( theseLines ) eq 0 then theseLines = replicate ( thisLine, maxLines )
      theseLines [ noLines ] = thisLine
      noLines = noLines + 1
      if noLines eq maxLines then begin
        print,'Error: Too many lines, increase maxlines in readtexfiles.pro'
        stop
      endif
    endelse
  endrep until doneThisMolecule
endwhile
free_lun, unit

;; Now we want to reorder the lines in molecule order
tagNames = tag_names ( lines )
for i = 0, n_elements ( molecules ) - 1 do begin
  name = molecules[i].l2name
  tagName = name
  if strpos ( tagName, '-' ) ne -1 then begin
    tagName = strjoin ( strsplit ( tagName, '-', /extract ), '__' )
  endif
  inds = where ( tagNames eq strupcase ( name ) )
  if inds[0] eq -1 or n_elements ( inds ) ne 1 then begin
    thisEntry = 0
  endif else begin
    thisEntry = lines.(inds[0])
  endelse
  if n_elements ( newLines ) eq 0 then begin
    newLines = create_struct ( tagName, thisEntry )
  endif else begin
    newLines = create_struct ( newLines, tagName, thisEntry )
  endelse
  molecules[i].noLines = n_tags ( thisEntry ) eq 0 ? 0 : n_elements ( thisEntry )
endfor

data = { $
  crossRefs: crossRefs, $
  molecules: molecules, $
  lines: newLines, $
  crossRefFileID: crossRefFileID, $
  moleculeFileID: molFileID, $
  lineFileID: lineFileID, $
  readID: '$Id$' }

return, data

end
