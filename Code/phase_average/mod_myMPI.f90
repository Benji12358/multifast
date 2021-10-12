module mod_myMPI

	use mod_vars

	implicit none

	integer :: ierr,mpi_er,comm, myid, numprocs

		contains

		subroutine myMPI  
			implicit none

			comm = mpi_comm_world
			! To get id of each participating process
			call mpi_comm_rank(comm,myid,ierr)
			call mpi_comm_size(comm,numprocs,ierr)

		end subroutine myMPI


end module mod_myMPI