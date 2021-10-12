module mod_restart

  use mod_inOut
  use mod_hdf5IO

  implicit none

    contains

      subroutine restart_phase_avg()

        implicit none
        integer :: i
        character(300)  :: path
        character(300)  :: filename
        logical :: fexist,check=.false.
        type(datafile)  :: file

        call myMPI()

    
               
        if(recover==1) then
          path = trim(external_fields_path)//'/phase_averaged_fields_recovery/'
            else
          path = trim(external_fields_path)//'/phase_averaged_fields/'
        endif


        filename = trim(path)//'W.h5'


        if(myid==0) then
          write(*,*)
          write(*,*)
          write(*,*) 'Searching for Phase averaged file...'
          write(*,*) filename
        end if
          
        inquire(file=filename,exist=fexist)

        if(fexist .and. myid==0) then
          write(*,*) 'Found restart file sucessfully'
          write(*,*)
        end if

        if (.not.fexist) then
          if(myid==0) then
            print '(A)', '########################################'
            print '(A)', '# Phase averaged file not found...     #'
            print '(A)', '# Exiting Solver...                    #'
            print '(A)', '########################################'
            call mpi_abort(comm,1,ierr)
          end if
        end if    

        if(myid==0) then
          write(*,*) 'Reading travelling wave parameters from restart file...'
        end if

        call read_file(filename,file)

        ! Read travelling wave parameters: amp_old, omega_old, xbins_old, tbins_old 
        call h5gopen_f(file%id,'setup',file%current_group, ierr)
        call read_hdf(file,amp_old,'amp')
        call read_hdf(file,omega_old,'omega')
        call read_hdf(file,kappa_old,'kappa')
        call read_hdf(file,xbins_old,'xbins')
        call h5gclose_f(file%current_group, ierr)

        if(myid==0) then
          write(*,*) 'amp_old =',   amp_old
          write(*,*) 'omega_old =', omega_old
          write(*,*) 'kappa_old =', kappa_old
          write(*,*) 'xbins_old =', xbins_old
        end if

        if(myid==0) then
          write(*,*) 'Successfully read travelling wave parameters from restart file'
        end if 

        if(amp_old/=amp) then
          call myMPI()
          if(myid==0) then
            write(*,*) '**************FATAL-ERROR***************'
            write(*,*) 'Change in the amplitude is not allowed with restart ON'
            write(*,*) 'Change the value of amp in travelling_wave.d to:', amp_old
            write(*,*) 'Or turn OFF restart to start phase averaging from the start.'
            write(*,*) '****************************************'
          end if
          call mpi_barrier(comm,ierr)
          call mpi_abort(comm,1,ierr)
        end if  

        if(omega_old/=omega) then
          call myMPI()
          if(myid==0) then
            write(*,*) '**************FATAL-ERROR***************'
            write(*,*) 'Change in the omega is not allowed with restart ON'
            write(*,*) 'Change the value of omega in travelling_wave.d to:', omega_old
            write(*,*) 'Or turn OFF restart to start phase averaging from the start.'
            write(*,*) '****************************************'
          end if
          call mpi_barrier(comm,ierr)
          call mpi_abort(comm,1,ierr)
        end if  

        if(kappa_old/=kappa) then
          call myMPI()
          if(myid==0) then
            write(*,*) '**************FATAL-ERROR***************'
            write(*,*) 'Change in the kappa is not allowed with restart ON'
            write(*,*) 'Change the value of omega in travelling_wave.d to:', kappa_old
            write(*,*) 'Or turn OFF restart to start phase averaging from the start.'
            write(*,*) '****************************************'
          end if
          call mpi_barrier(comm,ierr)
          call mpi_abort(comm,1,ierr)
        end if     

        if(xbins_old/=xbins) then
          call myMPI()
          if(myid==0) then
            write(*,*) '**************FATAL-ERROR***************'
            write(*,*) 'Change in number of bins is not allowed restart ON'
            write(*,*) 'Change the value of xbins in phase_averaging.d to:', xbins_old
            write(*,*) 'Or turn OFF restart to start phase averaging from the start.'
            write(*,*) '****************************************'
          end if
          call mpi_barrier(comm,ierr)
          call mpi_abort(comm,1,ierr)
        end if


        call avg_input_Xi(U,counts,'U',.false.,recover)
        call avg_input_Xi(V,counts,'V',.false.,recover)
        call avg_input_Xi(W,counts,'W',.true. ,recover)
        call avg_input_Xi(P,counts,'P',.false.,recover) 

        call avg_input_Xi(UU,counts,'UU',.false.,recover)
        call avg_input_Xi(VV,counts,'VV',.false.,recover)
        call avg_input_Xi(WW,counts,'WW',.false.,recover)
        call avg_input_Xi(PP,counts,'PP',.false.,recover)

        call avg_input_Xi(UV, counts,'UV', .false.,recover)
        call avg_input_Xi(UW, counts,'UW', .false.,recover)
        call avg_input_Xi(UP, counts,'UP', .false.,recover)
        call avg_input_Xi(VW, counts,'VW', .false.,recover)
        call avg_input_Xi(VP, counts,'VP', .false.,recover)

        call avg_input_Xi(UUV, counts, 'UUV', .false.,recover)
        call avg_input_Xi(VVV, counts, 'VVV', .false.,recover)
        call avg_input_Xi(WWV, counts, 'WWV', .false.,recover)
        call avg_input_Xi(UVV, counts, 'UVV', .false.,recover)

        call avg_input_Xi(dUdx, counts, 'dUdx', .false.,recover) 
        call avg_input_Xi(dVdx, counts, 'dVdx', .false.,recover)
        call avg_input_Xi(dWdx, counts, 'dWdx', .false.,recover)
        call avg_input_Xi(dPdx, counts, 'dPdx', .false.,recover)

        call avg_input_Xi(PdUdx, counts, 'PdUdx', .false.,recover)
        call avg_input_Xi(PdUdy, counts, 'PdUdy', .false.,recover)
        call avg_input_Xi(PdVdx, counts, 'PdVdx', .false.,recover)
        call avg_input_Xi(PdVdy, counts, 'PdVdy', .false.,recover)
        call avg_input_Xi(PdWdz, counts, 'PdWdz', .false.,recover)

        call avg_input_Xi(dUdxdUdx, counts, 'dUdxdUdx', .false.,recover)
        call avg_input_Xi(dUdydUdy, counts, 'dUdydUdy', .false.,recover)
        call avg_input_Xi(dUdzdUdz, counts, 'dUdzdUdz', .false.,recover) 

        call avg_input_Xi(dVdxdVdx, counts, 'dVdxdVdx', .false.,recover)
        call avg_input_Xi(dVdydVdy, counts, 'dVdydVdy', .false.,recover)
        call avg_input_Xi(dVdzdVdz, counts, 'dVdzdVdz', .false.,recover)

        call avg_input_Xi(dWdxdWdx, counts, 'dWdxdWdx', .false.,recover)
        call avg_input_Xi(dWdydWdy, counts, 'dWdydWdy', .false.,recover)
        call avg_input_Xi(dWdzdWdz, counts, 'dWdzdWdz', .false.,recover)

        call avg_input_Xi(dPdxdPdx, counts, 'dPdxdPdx', .false.,recover)
        call avg_input_Xi(dPdydPdy, counts, 'dPdydPdy', .false.,recover)
        call avg_input_Xi(dPdzdPdz, counts, 'dPdzdPdz', .false.,recover)

        call avg_input_Xi(dUdxdVdx, counts, 'dUdxdVdx', .false.,recover)
        call avg_input_Xi(dUdydVdy, counts, 'dUdydVdy', .false.,recover)
        call avg_input_Xi(dUdzdVdz, counts, 'dUdzdVdz', .false.,recover)

        call avg_input_Xi(dUUdx, counts, 'dUUdx', .false.,recover)
        call avg_input_Xi(dVVdx, counts, 'dVVdx', .false.,recover)
        call avg_input_Xi(dWWdx, counts, 'dWWdx', .false.,recover)
        call avg_input_Xi(dUVdx, counts, 'dUVdx', .false.,recover)

        call avg_input_Xi(omgX, counts, 'omgX', .false.,recover)
        call avg_input_Xi(omgY, counts, 'omgY', .false.,recover)
        call avg_input_Xi(omgZ, counts, 'omgZ', .false.,recover)

        call avg_input_Xi(omgXomgX, counts, 'omgXomgX', .false.,recover)
        call avg_input_Xi(omgYomgY, counts, 'omgYomgY', .false.,recover)
        call avg_input_Xi(omgZomgZ, counts, 'omgZomgZ', .false.,recover)

        call avg_input_Xi(omgXomgY, counts, 'omgXomgY', .false.,recover)
        call avg_input_Xi(omgXomgZ, counts, 'omgXomgZ', .false.,recover)
        call avg_input_Xi(omgYomgZ, counts, 'omgYomgZ', .false.,recover)

        call avg_input_Xi(omgXU, counts, 'omgXU', .false.,recover)
        call avg_input_Xi(omgXV, counts, 'omgXV', .false.,recover)
        call avg_input_Xi(omgYU, counts, 'omgYU', .false.,recover)
        call avg_input_Xi(omgYV, counts, 'omgYV', .false.,recover)
        call avg_input_Xi(omgZU, counts, 'omgZU', .false.,recover)
        call avg_input_Xi(omgZV, counts, 'omgZV', .false.,recover)

        call avg_input_Xi(omgXomgXV, counts, 'omgXomgXV', .false.,recover)
        call avg_input_Xi(omgYomgYV, counts, 'omgYomgYV', .false.,recover)
        call avg_input_Xi(omgZomgZV, counts, 'omgZomgZV', .false.,recover)
        call avg_input_Xi(omgXomgYV, counts, 'omgXomgYV', .false.,recover)
                
        call avg_input_Xi(domgXdx, counts, 'domgXdx', .false.,recover)
        call avg_input_Xi(domgYdx, counts, 'domgYdx', .false.,recover)
        call avg_input_Xi(domgZdx, counts, 'domgZdx', .false.,recover)

        call avg_input_Xi(omgXdUdx, counts, 'omgXdUdx', .false.,recover)
        call avg_input_Xi(omgXdUdy, counts, 'omgXdUdy', .false.,recover)
        call avg_input_Xi(omgXdUdz, counts, 'omgXdUdz', .false.,recover)
        call avg_input_Xi(omgXdVdx, counts, 'omgXdVdx', .false.,recover)
        call avg_input_Xi(omgXdVdy, counts, 'omgXdVdy', .false.,recover)
        call avg_input_Xi(omgXdVdz, counts, 'omgXdVdz', .false.,recover) 
        call avg_input_Xi(omgXdWdx, counts, 'omgXdWdx', .false.,recover)     

        call avg_input_Xi(omgYdUdx, counts, 'omgYdUdx', .false.,recover)
        call avg_input_Xi(omgYdUdy, counts, 'omgYdUdy', .false.,recover)
        call avg_input_Xi(omgYdUdz, counts, 'omgYdUdz', .false.,recover)
        call avg_input_Xi(omgYdVdx, counts, 'omgYdVdx', .false.,recover)
        call avg_input_Xi(omgYdVdy, counts, 'omgYdVdy', .false.,recover)
        call avg_input_Xi(omgYdVdz, counts, 'omgYdVdz', .false.,recover)
        call avg_input_Xi(omgYdWdy, counts, 'omgYdWdy', .false.,recover)

        call avg_input_Xi(omgZdUdz, counts, 'omgZdUdz', .false.,recover)
        call avg_input_Xi(omgZdVdz, counts, 'omgZdVdz', .false.,recover)
        call avg_input_Xi(omgZdWdx, counts, 'omgZdWdx', .false.,recover)
        call avg_input_Xi(omgZdWdy, counts, 'omgZdWdy', .false.,recover)
        call avg_input_Xi(omgZdWdz, counts, 'omgZdWdz', .false.,recover)

        call avg_input_Xi(omgXomgXdUdx, counts, 'omgXomgXdUdx', .false.,recover)
        call avg_input_Xi(omgXomgXdVdx, counts, 'omgXomgXdVdx', .false.,recover)

        call avg_input_Xi(omgXomgYdUdx, counts, 'omgXomgYdUdx', .false.,recover)
        call avg_input_Xi(omgXomgYdUdy, counts, 'omgXomgYdUdy', .false.,recover)
        call avg_input_Xi(omgXomgYdVdx, counts, 'omgXomgYdVdx', .false.,recover)
        call avg_input_Xi(omgXomgYdVdy, counts, 'omgXomgYdVdy', .false.,recover)

        call avg_input_Xi(omgXomgZdUdz, counts, 'omgXomgZdUdz', .false.,recover)
        call avg_input_Xi(omgXomgZdVdz, counts, 'omgXomgZdVdz', .false.,recover)
        call avg_input_Xi(omgXomgZdWdx, counts, 'omgXomgZdWdx', .false.,recover)

        call avg_input_Xi(omgYomgYdUdy, counts, 'omgYomgYdUdy', .false.,recover)
        call avg_input_Xi(omgYomgYdVdy, counts, 'omgYomgYdVdy', .false.,recover)

        call avg_input_Xi(omgYomgZdUdz, counts, 'omgYomgZdUdz', .false.,recover)
        call avg_input_Xi(omgYomgZdVdz, counts, 'omgYomgZdVdz', .false.,recover)
        call avg_input_Xi(omgYomgZdWdy, counts, 'omgYomgZdWdy', .false.,recover)

        call avg_input_Xi(omgZomgZdWdz, counts, 'omgZomgZdWdz', .false.,recover)
                
        call avg_input_Xi(domgXdxdomgXdx, counts, 'domgXdxdomgXdx', .false.,recover)
        call avg_input_Xi(domgXdydomgXdy, counts, 'domgXdydomgXdy', .false.,recover)
        call avg_input_Xi(domgXdzdomgXdz, counts, 'domgXdzdomgXdz', .false.,recover)

        call avg_input_Xi(domgYdxdomgYdx, counts, 'domgYdxdomgYdx', .false.,recover)
        call avg_input_Xi(domgYdydomgYdy, counts, 'domgYdydomgYdy', .false.,recover)
        call avg_input_Xi(domgYdzdomgYdz, counts, 'domgYdzdomgYdz', .false.,recover)

        call avg_input_Xi(domgZdxdomgZdx, counts, 'domgZdxdomgZdx', .false.,recover)
        call avg_input_Xi(domgZdydomgZdy, counts, 'domgZdydomgZdy', .false.,recover)
        call avg_input_Xi(domgZdzdomgZdz, counts, 'domgZdzdomgZdz', .false.,recover)

        call avg_input_Xi(domgXdxdomgYdx, counts, 'domgXdxdomgYdx', .false.,recover)
        call avg_input_Xi(domgXdydomgYdy, counts, 'domgXdydomgYdy', .false.,recover)
        call avg_input_Xi(domgXdzdomgYdz, counts, 'domgXdzdomgYdz', .false.,recover)
  
        end subroutine restart_phase_avg


end module mod_restart
