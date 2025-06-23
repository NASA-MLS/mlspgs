s/integer, parameter :: \([A-Za-z_][A-Za-z_0-9]*\)/case(\1)/
s/'//g
s/^! [A-Z].*$/		print*, '&'/
s/=/!/
p
