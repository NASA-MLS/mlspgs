!... Use Planck function to brightness temperature

	subroutine planck(temp,freq,tb)
        use MLSCommon, only: r8
	real(r8) :: temp,tb,freq
	real(r8) :: h
	real(r8) :: k

	h = 6.6256
	k = 1.3805

        tb=h*freq*1.e-2_r8/(exp(h*freq*1.e-2_r8/k/temp)-1.)/k

	return 
	end

! $Log: planck.f90,v      
