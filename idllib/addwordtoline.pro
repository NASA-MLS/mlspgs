; $Id$

pro AddWordToLine, line, unit, indent, word

if strlen(line) + strlen(word) gt 72 then begin
  printf, unit, line+'$'
  line = strjoin(replicate(' ',indent))
endif
line = line + word
end
