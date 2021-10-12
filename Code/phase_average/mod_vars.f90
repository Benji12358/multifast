module mod_vars

	use mod_params

	implicit none

	real*8,  dimension(:),  allocatable :: xvals

	integer, dimension(:,:,:),allocatable :: counts

	real*8,  dimension(:,:,:),allocatable :: U
	real*8,  dimension(:,:,:),allocatable :: V
	real*8,  dimension(:,:,:),allocatable :: W
	real*8,  dimension(:,:,:),allocatable :: P

	real*8,  dimension(:,:,:),allocatable :: UU
	real*8,  dimension(:,:,:),allocatable :: VV
	real*8,  dimension(:,:,:),allocatable :: WW
	real*8,  dimension(:,:,:),allocatable :: PP

	real*8,  dimension(:,:,:),allocatable :: UV
	real*8,  dimension(:,:,:),allocatable :: UW
	real*8,  dimension(:,:,:),allocatable :: UP
	real*8,  dimension(:,:,:),allocatable :: VW
	real*8,  dimension(:,:,:),allocatable :: VP

	real*8,  dimension(:,:,:),allocatable :: UUV, VVV, WWV, UVV

	real*8,  dimension(:,:,:),allocatable :: dUdx, dVdx, dWdx, dPdx

	real*8,  dimension(:,:,:),allocatable :: PdUdx, PdUdy
	real*8,  dimension(:,:,:),allocatable :: PdVdx, PdVdy
	real*8,  dimension(:,:,:),allocatable :: PdWdz

	real*8,  dimension(:,:,:),allocatable :: dUdxdUdx, dVdxdVdx, dWdxdWdx, dPdxdPdx
	real*8,  dimension(:,:,:),allocatable :: dUdydUdy, dVdydVdy, dWdydWdy, dPdydPdy
	real*8,  dimension(:,:,:),allocatable :: dUdzdUdz, dVdzdVdz, dWdzdWdz, dPdzdPdz
	real*8,  dimension(:,:,:),allocatable :: dUdxdVdx, dUdydVdy, dUdzdVdz

	real*8,  dimension(:,:,:),allocatable :: dUUdx, dVVdx, dWWdx, dUVdx  

	real*8,  dimension(:,:,:),allocatable :: omgX 
	real*8,  dimension(:,:,:),allocatable :: omgY
	real*8,  dimension(:,:,:),allocatable :: omgZ 

	real*8,  dimension(:,:,:),allocatable :: omgXomgX 
	real*8,  dimension(:,:,:),allocatable :: omgYomgY
	real*8,  dimension(:,:,:),allocatable :: omgZomgZ

	real*8,  dimension(:,:,:),allocatable :: omgXomgY 
	real*8,  dimension(:,:,:),allocatable :: omgXomgZ
	real*8,  dimension(:,:,:),allocatable :: omgYomgZ

	real*8,  dimension(:,:,:),allocatable :: omgXU 
	real*8,  dimension(:,:,:),allocatable :: omgXV
	real*8,  dimension(:,:,:),allocatable :: omgYU
	real*8,  dimension(:,:,:),allocatable :: omgYV 
	real*8,  dimension(:,:,:),allocatable :: omgZU
	real*8,  dimension(:,:,:),allocatable :: omgZV

	real*8,  dimension(:,:,:),allocatable :: omgXomgXV 
	real*8,  dimension(:,:,:),allocatable :: omgYomgYV
	real*8,  dimension(:,:,:),allocatable :: omgZomgZV
	real*8,  dimension(:,:,:),allocatable :: omgXomgYV 

	real*8,  dimension(:,:,:),allocatable :: domgXdx 
	real*8,  dimension(:,:,:),allocatable :: domgYdx
	real*8,  dimension(:,:,:),allocatable :: domgZdx

	real*8,  dimension(:,:,:),allocatable :: omgXdUdx 
	real*8,  dimension(:,:,:),allocatable :: omgXdUdy
	real*8,  dimension(:,:,:),allocatable :: omgXdUdz
	real*8,  dimension(:,:,:),allocatable :: omgXdVdX
	real*8,  dimension(:,:,:),allocatable :: omgXdVdY
	real*8,  dimension(:,:,:),allocatable :: omgXdVdz
	real*8,  dimension(:,:,:),allocatable :: omgXdWdx

	real*8,  dimension(:,:,:),allocatable :: omgYdUdx
	real*8,  dimension(:,:,:),allocatable :: omgYdUdy
	real*8,  dimension(:,:,:),allocatable :: omgYdUdz
	real*8,  dimension(:,:,:),allocatable :: omgYdVdx 
	real*8,  dimension(:,:,:),allocatable :: omgYdVdy
	real*8,  dimension(:,:,:),allocatable :: omgYdVdz
	real*8,  dimension(:,:,:),allocatable :: omgYdWdy

	real*8,  dimension(:,:,:),allocatable :: omgZdUdz
	real*8,  dimension(:,:,:),allocatable :: omgZdVdz
	real*8,  dimension(:,:,:),allocatable :: omgZdWdx 
	real*8,  dimension(:,:,:),allocatable :: omgZdWdy
	real*8,  dimension(:,:,:),allocatable :: omgZdWdz

	real*8,  dimension(:,:,:),allocatable :: omgXomgXdUdx
	real*8,  dimension(:,:,:),allocatable :: omgXomgXdVdx

	real*8,  dimension(:,:,:),allocatable :: omgXomgYdUdx  
	real*8,  dimension(:,:,:),allocatable :: omgXomgYdUdy
	real*8,  dimension(:,:,:),allocatable :: omgXomgYdVdx
	real*8,  dimension(:,:,:),allocatable :: omgXomgYdVdy

	real*8,  dimension(:,:,:),allocatable :: omgXomgZdUdz
	real*8,  dimension(:,:,:),allocatable :: omgXomgZdVdz
	real*8,  dimension(:,:,:),allocatable :: omgXomgZdWdx

	real*8,  dimension(:,:,:),allocatable :: omgYomgYdUdy
	real*8,  dimension(:,:,:),allocatable :: omgYomgYdVdy

	real*8,  dimension(:,:,:),allocatable :: omgYomgZdUdz
	real*8,  dimension(:,:,:),allocatable :: omgYomgZdVdz  
	real*8,  dimension(:,:,:),allocatable :: omgYomgZdWdy

	real*8,  dimension(:,:,:),allocatable :: omgZomgZdWdz
	   
	real*8,  dimension(:,:,:),allocatable :: domgXdxdomgXdx 
	real*8,  dimension(:,:,:),allocatable :: domgXdydomgXdy
	real*8,  dimension(:,:,:),allocatable :: domgXdzdomgXdz

	real*8,  dimension(:,:,:),allocatable :: domgYdxdomgYdx 
	real*8,  dimension(:,:,:),allocatable :: domgYdydomgydY
	real*8,  dimension(:,:,:),allocatable :: domgYdzdomgYdz

	real*8,  dimension(:,:,:),allocatable :: domgZdxdomgZdx 
	real*8,  dimension(:,:,:),allocatable :: domgZdydomgZdy
	real*8,  dimension(:,:,:),allocatable :: domgZdzdomgZdz

	real*8,  dimension(:,:,:),allocatable :: domgXdxdomgYdx 
	real*8,  dimension(:,:,:),allocatable :: domgXdydomgYdy
	real*8,  dimension(:,:,:),allocatable :: domgXdzdomgYdz 

	
end module mod_vars