function ReadL2CFArray, filename, prefix=prefix, suffix=suffix, extraSuffix=extraSuffix

openr, unit, filename, /get_lun
found = 0
line = ''
while not found and not eof ( unit ) do begin
  readf, unit, line, format='(a)'
  line = strtrim ( line, 2 )
  start = strmid ( line, 0, strlen ( prefix ) )
  if start eq prefix then found = 1
endwhile
if not found then return, empty

done = 0

line = strmid ( line, strlen ( prefix ), strlen ( line ) )
snippet = ''

repeat begin
  lastChar = strmid ( line, strlen ( line ) - 1, 1 )
  if lastChar eq '$' then begin
    snippet = snippet + strmid ( line, 0, strlen ( line ) - 1 )
    readf, unit, line, format='(a)'
    line = strtrim ( line, 2 )
  endif else begin
    snippet = snippet + line
    done = 1
  endelse
endrep until done or eof ( unit )
if not done then return, empty

if n_elements ( suffix ) ne 0 then begin
  trimPos = strpos ( snippet, suffix, /reverse_search )
  snippet = strmid ( snippet, 0, trimPos )
endif

;; Strip out any extra suffix stuff
if n_elements ( extraSuffix ) ne 0 then begin
  for i = 0, n_elements ( extraSuffix ) - 1 do begin
    trimPos = strpos ( snippet, extraSuffix[i], /reverse_search )
    if trimPos ne -1 then snippet = strmid ( snippet, 0, trimPos )
  endfor
endif

dummy = execute ( 'result=' + snippet )
return, result

end

  

  
