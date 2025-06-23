; $Id$

pro AddWordToLine, line, units, indent, word

if strlen(line) + strlen(word) gt 72 then begin
  for i = 0, n_elements ( units ) - 1 do $
    printf, units(i), line+'$'
  line = strjoin(replicate(' ',indent))
endif
line = line + word
end
