c... Use Planck function to brightness temperature

	subroutine planck(temp,freq,tb)
	real temp,tb,freq
	real h
	real k

	h = 6.6256
	k = 1.3805

        tb=h*freq*1.e-2/(exp(h*freq*1.e-2/k/temp)-1)/k

	return 
	end
