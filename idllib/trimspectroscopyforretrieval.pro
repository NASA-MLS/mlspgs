pro TrimSpectroscopyForRetrieval, data

noMolecules = n_elements(data)

for mol = 0, noMolecules - 1 do begin
  noLines = data(mol).no_lines
  for line = 0, noLines - 1 do begin
    bands = data(mol).eos_bands(line)
    frq = data(mol).frq(line)
    words = strtrim ( strsplit ( bands, ',', /extract ), 2)
    noBands = n_elements(words)
    scoreL = 0
    scoreU = 0
    for band = 0, noBands - 1 do begin
      scoreL += max ( words(band) eq [ 'B2L', 'B3L', 'B4L', 'B5L', 'B6L', 'B23L', 'B27L' ] )
      scoreU += max ( words(band) eq [ 'B2U', 'B3U', 'B4U', 'B5U', 'B6U', 'B23U', 'B27U' ] )
    endfor
    ;; Only bother with cases where all of R2 (U/L) listed
    if scoreL ne 7 and scoreU ne 7 then continue
    ;; Skip the 183GHz h2o line
    if data(mol).name eq 'H$_{2}$O' and $
      frq gt 183310.116d0 and frq lt 183310.118d0 then continue
    if scoreL eq 7 and scoreU eq 7 then begin
      print,'The retrieval spectroscopy generator hit a case it was not trained for'
      stop
    endif
    newBands = ''
    if scoreL eq 7 then begin
      if frq lt 180500.0d0 $
        then newBands = 'B5L, B6L, B27L' $
      else newBands = 'B2L, B3L, B4L, B23L'
    endif
    if scoreU eq 7 then begin
      if newBands ne '' then newBands = newBands + ', '
      if frq gt 203000.0d0 $
        then newBands = newBands + 'B5U, B6U, B27U' $
      else newBands = newBands + 'B2U, B3U, B4U, B23U'
    endif
    data(mol).eos_bands(line) = newBands
  endfor
endfor

end
