module mod_init

	use mod_restart
	use mod_inOut

	implicit none

	contains

		subroutine init_phase_avg
		!--------------------------------------------------------------------
		! This subroutine initializes the phase averaging parameters
		! performs the appropriate unit conversions if requested.
		! This subroutine is called in the main program before time loop.
		!--------------------------------------------------------------------
			implicit none
			integer :: i
			character(300)  :: path

			call myMPI()

			call read_phase_average_parameters()

			! To convert from inner to outer units in case travelling
			! wave parameters are given in inner units
			if(inner_units==1) then
				cf_amp   = Re_tau/ren
				cf_omega = Re_tau**2/ren  
				cf_kappa = Re_tau
			else if(inner_units==0) then
				cf_amp   = 1.0d0
				cf_omega = 1.0d0  
				cf_kappa = 1.0d0
			end if

			amp   = Twave%amp  *cf_amp          
			kappa = Twave%kappa*cf_kappa
			omega = Twave%omega*cf_omega
			wave_speed = omega/kappa

			if(phase_averaging==1)  then
				call allocate_arrays()
				path = trim(results_path)//'phase_averaged_fields/'
				call system('mkdir -p '//trim(path))

				path = trim(results_path)//'phase_averaged_fields_recovery/'
				call system('mkdir -p '//trim(path))

				if(restart==1 .and. omega/=0) call restart_phase_avg()

			end if

			if(kappa/=0) then

				lambda = 2*Pi/kappa          ! lambda is the wavelength of the travelling wave
				xmin   = 0                   ! Starting point of the cycle: zero
				xmax   = lambda              ! Ending point of the cycle: lambda
				delX   = (xmax-xmin)/(xbins) ! Width of the bins

				! Calculate the x-coordinates of the bins in a cycle
				! Useful for plotting the phase-averaged value
				if(.not. allocated(xvals))  allocate(xvals(xbins))
				xvals=0
				do i=1,xbins 
					xvals(i) = xmin + (i-0.5)*delX
				end do

			end if
 

		end subroutine init_phase_avg


		subroutine allocate_arrays()

			implicit none
			integer :: a,b,c,d

	        a=zstart(1);b=min(zend(1),n1m)
	        c=zstart(2);d=min(zend(2),n2m)

	        allocate(counts(a:b,c:d,xbins))
	        counts=0

	        allocate(U(a:b,c:d,xbins))
	        U=0.0d0
	        allocate(V(a:b,c:d,xbins))
	        V=0.0d0
	        allocate(W(a:b,c:d,xbins))
	        W=0.0d0
	        allocate(P(a:b,c:d,xbins))
	        P=0.0d0
	        
	        
	        allocate(UU(a:b,c:d,xbins))
	        UU=0.0d0
	        allocate(VV(a:b,c:d,xbins))
	        VV=0.0d0
	        allocate(WW(a:b,c:d,xbins))
	        WW=0.0d0
	        allocate(PP(a:b,c:d,xbins))
	        PP=0.0d0


	        allocate(UV(a:b,c:d,xbins))
	        UV=0.0d0
	        allocate(UW(a:b,c:d,xbins))
	        UW=0.0d0
	        allocate(UP(a:b,c:d,xbins))
	        UP=0.0d0
	        allocate(VW(a:b,c:d,xbins))
	        VW=0.0d0
	        allocate(VP(a:b,c:d,xbins))
	        VP=0.0d0
	        
	        
	        allocate(UUV(a:b,c:d,xbins))
	        UUV=0.0d0
	        allocate(VVV(a:b,c:d,xbins))
	        VVV=0.0d0
	        allocate(WWV(a:b,c:d,xbins))
	        WWV=0.0d0
	        allocate(UVV(a:b,c:d,xbins))
	        UVV=0.0d0
	  
	        
	        allocate(dUdx(a:b,c:d,xbins))
	        dUdx=0.0d0
	        allocate(dVdx(a:b,c:d,xbins))
	        dVdx=0.0d0
	        allocate(dWdx(a:b,c:d,xbins))
	        dWdx=0.0d0
	        allocate(dPdx(a:b,c:d,xbins))
	        dPdx=0.0d0


	        allocate(PdUdx(a:b,c:d,xbins))
	        PdUdx=0.0d0
	        allocate(PdUdy(a:b,c:d,xbins))
	        PdUdy=0.0d0
	        allocate(PdVdx(a:b,c:d,xbins))
	        PdVdx=0.0d0
	        allocate(PdVdy(a:b,c:d,xbins))
	        PdVdy=0.0d0
	        allocate(PdWdz(a:b,c:d,xbins))
	        PdWdz=0.0d0


	        allocate(dUdxdUdx(a:b,c:d,xbins))
	        dUdxdUdx=0.0d0
	        allocate(dUdydUdy(a:b,c:d,xbins))
	        dUdydUdy=0.0d0
	        allocate(dUdzdUdz(a:b,c:d,xbins))
	        dUdzdUdz=0.0d0


	        allocate(dVdxdVdx(a:b,c:d,xbins))
	        dVdxdVdx=0.0d0
	        allocate(dVdydVdy(a:b,c:d,xbins))
	        dVdydVdy=0.0d0
	        allocate(dVdzdVdz(a:b,c:d,xbins))
	        dVdzdVdz=0.d0


	        allocate(dWdxdWdx(a:b,c:d,xbins))
	        dWdxdWdx=0.0d0
	        allocate(dWdydWdy(a:b,c:d,xbins))
	        dWdydWdy=0.0d0
	        allocate(dWdzdWdz(a:b,c:d,xbins))
	        dWdzdWdz=0.0d0


	        allocate(dPdxdPdx(a:b,c:d,xbins))
	        dPdxdPdx=0.0d0
	        allocate(dPdydPdy(a:b,c:d,xbins))
	        dPdydPdy=0.0d0
	        allocate(dPdzdPdz(a:b,c:d,xbins))
	        dPdzdPdz=0.0d0


	        allocate(dUdxdVdx(a:b,c:d,xbins))
	        dUdxdVdx=0.0d0
	        allocate(dUdydVdy(a:b,c:d,xbins))
	        dUdydVdy=0.0d0
	        allocate(dUdzdVdz(a:b,c:d,xbins))
	        dUdzdVdz=0.0d0

	        
	        allocate(dUUdx(a:b,c:d,xbins))
	        dUUdx=0.0d0
	        allocate(dVVdx(a:b,c:d,xbins))
	        dVVdx=0.0d0
	        allocate(dWWdx(a:b,c:d,xbins))
	        dWWdx=0.0d0
	        allocate(dUVdx(a:b,c:d,xbins))
	        dUVdx=0.0d0        

	                
	        allocate(omgX(a:b,c:d,xbins))
	        omgX=0.0d0
	        allocate(omgY(a:b,c:d,xbins))
	        omgY=0.0d0
	        allocate(omgZ(a:b,c:d,xbins))
	        omgZ=0.0d0

	      
	        allocate(omgXomgX(a:b,c:d,xbins))
	        omgXomgX=0.0d0
	        allocate(omgYomgY(a:b,c:d,xbins))
	        omgYomgY=0.0d0
	        allocate(omgZomgZ(a:b,c:d,xbins))
	        omgZomgZ=0.0d0

	        
	        allocate(omgXomgY(a:b,c:d,xbins))
	        omgXomgY=0.0d0
	        allocate(omgXomgZ(a:b,c:d,xbins))
	        omgXomgZ=0.0d0
	        allocate(omgYomgZ(a:b,c:d,xbins))
	        omgYomgZ=0.0d0

	        
	        allocate(omgXU(a:b,c:d,xbins))
	        omgXU=0.0d0
	        allocate(omgXV(a:b,c:d,xbins))
	        omgXV=0.0d0
	        allocate(omgYU(a:b,c:d,xbins))
	        omgYU=0.0d0
	        allocate(omgYV(a:b,c:d,xbins))
	        omgYV=0.0d0
	        allocate(omgZU(a:b,c:d,xbins))
	        omgZU=0.0d0
	        allocate(omgZV(a:b,c:d,xbins))
	        omgZV=0.0d0

	               
	        allocate(omgXomgXV(a:b,c:d,xbins))
	        omgXomgXV=0.0d0
	        allocate(omgYomgYV(a:b,c:d,xbins))
	        omgYomgYV=0.0d0
	        allocate(omgZomgZV(a:b,c:d,xbins))
	        omgZomgZV=0.0d0
	        allocate(omgXomgYV(a:b,c:d,xbins))
	        omgXomgYV=0.0d0

	        
	        allocate(domgXdx(a:b,c:d,xbins))
	        domgXdx=0.0d0
	        allocate(domgYdx(a:b,c:d,xbins))
	        domgYdx=0.0d0
	        allocate(domgZdx(a:b,c:d,xbins))
	        domgZdx=0.0d0

	                
	        allocate(omgXdUdx(a:b,c:d,xbins))
	        omgXdUdx=0.0d0
	        allocate(omgXdUdy(a:b,c:d,xbins))
	        omgXdUdy=0.0d0
	        allocate(omgXdUdz(a:b,c:d,xbins))
	        omgXdUdz=0.0d0  

	        
	        allocate(omgXdVdx(a:b,c:d,xbins))
	        omgXdVdx=0.0d0
	        allocate(omgXdVdy(a:b,c:d,xbins))
	        omgXdVdy=0.0d0
	        allocate(omgXdVdz(a:b,c:d,xbins))
	        omgXdVdz=0.0d0   
	        

	        allocate(omgXdWdx(a:b,c:d,xbins))
	        omgXdWdx=0.0d0
	                

	        allocate(omgYdUdx(a:b,c:d,xbins))
	        omgYdUdx=0.0d0
	        allocate(omgYdUdy(a:b,c:d,xbins))
	        omgYdUdy=0.0d0
	        allocate(omgYdUdz(a:b,c:d,xbins))
	        omgYdUdz=0.0d0   
	        

	        allocate(omgYdVdx(a:b,c:d,xbins))
	        omgYdVdx=0.0d0
	        allocate(omgYdVdy(a:b,c:d,xbins))
	        omgYdVdy=0.0d0
	        allocate(omgYdVdz(a:b,c:d,xbins))
	        omgYdVdz=0.0d0 
	        

	        allocate(omgYdWdy(a:b,c:d,xbins))
	        omgYdWdy=0.0d0 
	        

	        allocate(omgZdUdz(a:b,c:d,xbins))
	        omgZdUdz=0.0d0   
	        allocate(omgZdVdz(a:b,c:d,xbins))
	        omgZdVdz=0.0d0

	        
	        allocate(omgZdWdx(a:b,c:d,xbins))
	        omgZdWdx=0.0d0
	        allocate(omgZdWdy(a:b,c:d,xbins))
	        omgZdWdy=0.0d0
	        allocate(omgZdWdz(a:b,c:d,xbins))
	        omgZdWdz=0.0d0

	        
	        allocate(omgXomgXdUdx(a:b,c:d,xbins))
	        omgXomgXdUdx=0.0d0
	        allocate(omgXomgXdVdx(a:b,c:d,xbins))
	        omgXomgXdVdx=0.0d0

	        
	        allocate(omgXomgYdUdx(a:b,c:d,xbins))
	        omgXomgYdUdx=0.0d0
	        allocate(omgXomgYdUdy(a:b,c:d,xbins))
	        omgXomgYdUdy=0.0d0
	        allocate(omgXomgYdVdx(a:b,c:d,xbins))
	        omgXomgYdVdx=0.0d0
	        allocate(omgXomgYdVdy(a:b,c:d,xbins))
	        omgXomgYdVdy=0.0d0

	                
	        allocate(domgXdxdomgXdx(a:b,c:d,xbins))
	        domgXdxdomgXdx=0.0d0
	        allocate(domgXdydomgXdy(a:b,c:d,xbins))
	        domgXdydomgXdy=0.0d0
	        allocate(domgXdzdomgXdz(a:b,c:d,xbins))
	        domgXdzdomgXdz=0.0d0

	        
	        allocate(omgXomgZdUdz(a:b,c:d,xbins))
	        omgXomgZdUdz=0.0d0
	        allocate(omgXomgZdVdz(a:b,c:d,xbins))
	        omgXomgZdVdz=0.0d0
	        allocate(omgXomgZdWdx(a:b,c:d,xbins))
	        omgXomgZdWdx=0.0d0

	        
	        allocate(omgYomgYdUdy(a:b,c:d,xbins))
	        omgYomgYdUdy=0.0d0
	        allocate(omgYomgYdVdy(a:b,c:d,xbins))
	        omgYomgYdVdy=0.0d0
	        

	        allocate(omgYomgZdUdz(a:b,c:d,xbins))
	        omgYomgZdUdz=0.0d0
	        allocate(omgYomgZdVdz(a:b,c:d,xbins))
	        omgYomgZdVdz=0.0d0
	        allocate(omgYomgZdWdy(a:b,c:d,xbins))
	        omgYomgZdWdy=0.0d0

	        
	        allocate(omgZomgZdWdz(a:b,c:d,xbins))
	        omgZomgZdWdz=0.0d0
	        

	        allocate(domgYdxdomgYdx(a:b,c:d,xbins))
	        domgYdxdomgYdx=0.0d0
	        allocate(domgYdydomgYdy(a:b,c:d,xbins))
	        domgYdydomgYdy=0.0d0
	        allocate(domgYdzdomgYdz(a:b,c:d,xbins))
	        domgYdzdomgYdz=0.0d0
	        

	        allocate(domgZdxdomgZdx(a:b,c:d,xbins))
	        domgZdxdomgZdx=0.0d0
	        allocate(domgZdydomgZdy(a:b,c:d,xbins))
	        domgZdydomgZdy=0.0d0
	        allocate(domgZdzdomgZdz(a:b,c:d,xbins))
	        domgZdzdomgZdz=0.0d0

	        
	        allocate(domgXdxdomgYdx(a:b,c:d,xbins))
	        domgXdxdomgYdx=0.0d0
	        allocate(domgXdydomgYdy(a:b,c:d,xbins))
	        domgXdydomgYdy=0.0d0
	        allocate(domgXdzdomgYdz(a:b,c:d,xbins))
	        domgXdzdomgYdz=0.0d0
	        

		end subroutine allocate_arrays



end module mod_init
