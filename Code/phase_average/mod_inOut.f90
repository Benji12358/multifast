module mod_inOut

    use hdf5
  use mod_hdf5IO  
  use mod_params
  use mod_vars
  use mod_myMPI


  implicit none

    contains

      subroutine read_phase_average_parameters
      !---------------------------------------------------------------
      ! This subroutine reads the phase averaging parameters
      ! from file 'phase_average.d'
      ! This subroutine is called in the init_phase_avg subroutine.
      !---------------------------------------------------------------
        implicit none

        open(15,file=trim(COMMON_settings_path)//'phase_averaging.d')
          read(15,*) phase_averaging
          read(15,*) xbins
          read(15,*) d_freq
          read(15,*) restart
          read(15,*) recover
        close(15)

      end subroutine read_phase_average_parameters

      subroutine dumpOUTPUT(recovery)

        implicit none
        integer, intent(in) :: recovery 


        call avg_output_Xi(U,counts,'U',.false.,recovery)
        call avg_output_Xi(V,counts,'V',.false.,recovery)
        call avg_output_Xi(W,counts,'W',.true. ,recovery)
        call avg_output_Xi(P,counts,'P',.false.,recovery) 

        call avg_output_Xi(UU,counts,'UU',.false.,recovery)
        call avg_output_Xi(VV,counts,'VV',.false.,recovery)
        call avg_output_Xi(WW,counts,'WW',.false.,recovery)
        call avg_output_Xi(PP,counts,'PP',.false.,recovery)

        call avg_output_Xi(UV, counts,'UV', .false.,recovery)
        call avg_output_Xi(UW, counts,'UW', .false.,recovery)
        call avg_output_Xi(UP, counts,'UP', .false.,recovery)
        call avg_output_Xi(VW, counts,'VW', .false.,recovery)
        call avg_output_Xi(VP, counts,'VP', .false.,recovery)

        call avg_output_Xi(UUV, counts, 'UUV', .false.,recovery)
        call avg_output_Xi(VVV, counts, 'VVV', .false.,recovery)
        call avg_output_Xi(WWV, counts, 'WWV', .false.,recovery)
        call avg_output_Xi(UVV, counts, 'UVV', .false.,recovery)

        call avg_output_Xi(dUdx, counts, 'dUdx', .false.,recovery) 
        call avg_output_Xi(dVdx, counts, 'dVdx', .false.,recovery)
        call avg_output_Xi(dWdx, counts, 'dWdx', .false.,recovery)
        call avg_output_Xi(dPdx, counts, 'dPdx', .false.,recovery)

        call avg_output_Xi(PdUdx, counts, 'PdUdx', .false.,recovery)
        call avg_output_Xi(PdUdy, counts, 'PdUdy', .false.,recovery)
        call avg_output_Xi(PdVdx, counts, 'PdVdx', .false.,recovery)
        call avg_output_Xi(PdVdy, counts, 'PdVdy', .false.,recovery)
        call avg_output_Xi(PdWdz, counts, 'PdWdz', .false.,recovery)

        call avg_output_Xi(dUdxdUdx, counts, 'dUdxdUdx', .false.,recovery)
        call avg_output_Xi(dUdydUdy, counts, 'dUdydUdy', .false.,recovery)
        call avg_output_Xi(dUdzdUdz, counts, 'dUdzdUdz', .false.,recovery) 

        call avg_output_Xi(dVdxdVdx, counts, 'dVdxdVdx', .false.,recovery)
        call avg_output_Xi(dVdydVdy, counts, 'dVdydVdy', .false.,recovery)
        call avg_output_Xi(dVdzdVdz, counts, 'dVdzdVdz', .false.,recovery)

        call avg_output_Xi(dWdxdWdx, counts, 'dWdxdWdx', .false.,recovery)
        call avg_output_Xi(dWdydWdy, counts, 'dWdydWdy', .false.,recovery)
        call avg_output_Xi(dWdzdWdz, counts, 'dWdzdWdz', .false.,recovery)

        call avg_output_Xi(dPdxdPdx, counts, 'dPdxdPdx', .false.,recovery)
        call avg_output_Xi(dPdydPdy, counts, 'dPdydPdy', .false.,recovery)
        call avg_output_Xi(dPdzdPdz, counts, 'dPdzdPdz', .false.,recovery)

        call avg_output_Xi(dUdxdVdx, counts, 'dUdxdVdx', .false.,recovery)
        call avg_output_Xi(dUdydVdy, counts, 'dUdydVdy', .false.,recovery)
        call avg_output_Xi(dUdzdVdz, counts, 'dUdzdVdz', .false.,recovery)

        call avg_output_Xi(dUUdx, counts, 'dUUdx', .false.,recovery)
        call avg_output_Xi(dVVdx, counts, 'dVVdx', .false.,recovery)
        call avg_output_Xi(dWWdx, counts, 'dWWdx', .false.,recovery)
        call avg_output_Xi(dUVdx, counts, 'dUVdx', .false.,recovery)

        call avg_output_Xi(omgX, counts, 'omgX', .false.,recovery)
        call avg_output_Xi(omgY, counts, 'omgY', .false.,recovery)
        call avg_output_Xi(omgZ, counts, 'omgZ', .false.,recovery)

        call avg_output_Xi(omgXomgX, counts, 'omgXomgX', .false.,recovery)
        call avg_output_Xi(omgYomgY, counts, 'omgYomgY', .false.,recovery)
        call avg_output_Xi(omgZomgZ, counts, 'omgZomgZ', .false.,recovery)

        call avg_output_Xi(omgXomgY, counts, 'omgXomgY', .false.,recovery)
        call avg_output_Xi(omgXomgZ, counts, 'omgXomgZ', .false.,recovery)
        call avg_output_Xi(omgYomgZ, counts, 'omgYomgZ', .false.,recovery)

        call avg_output_Xi(omgXU, counts, 'omgXU', .false.,recovery)
        call avg_output_Xi(omgXV, counts, 'omgXV', .false.,recovery)
        call avg_output_Xi(omgYU, counts, 'omgYU', .false.,recovery)
        call avg_output_Xi(omgYV, counts, 'omgYV', .false.,recovery)
        call avg_output_Xi(omgZU, counts, 'omgZU', .false.,recovery)
        call avg_output_Xi(omgZV, counts, 'omgZV', .false.,recovery)

        call avg_output_Xi(omgXomgXV, counts, 'omgXomgXV', .false.,recovery)
        call avg_output_Xi(omgYomgYV, counts, 'omgYomgYV', .false.,recovery)
        call avg_output_Xi(omgZomgZV, counts, 'omgZomgZV', .false.,recovery)
        call avg_output_Xi(omgXomgYV, counts, 'omgXomgYV', .false.,recovery)
                
        call avg_output_Xi(domgXdx, counts, 'domgXdx', .false.,recovery)
        call avg_output_Xi(domgYdx, counts, 'domgYdx', .false.,recovery)
        call avg_output_Xi(domgZdx, counts, 'domgZdx', .false.,recovery)

        call avg_output_Xi(omgXdUdx, counts, 'omgXdUdx', .false.,recovery)
        call avg_output_Xi(omgXdUdy, counts, 'omgXdUdy', .false.,recovery)
        call avg_output_Xi(omgXdUdz, counts, 'omgXdUdz', .false.,recovery)
        call avg_output_Xi(omgXdVdx, counts, 'omgXdVdx', .false.,recovery)
        call avg_output_Xi(omgXdVdy, counts, 'omgXdVdy', .false.,recovery)
        call avg_output_Xi(omgXdVdz, counts, 'omgXdVdz', .false.,recovery) 
        call avg_output_Xi(omgXdWdx, counts, 'omgXdWdx', .false.,recovery)     

        call avg_output_Xi(omgYdUdx, counts, 'omgYdUdx', .false.,recovery)
        call avg_output_Xi(omgYdUdy, counts, 'omgYdUdy', .false.,recovery)
        call avg_output_Xi(omgYdUdz, counts, 'omgYdUdz', .false.,recovery)
        call avg_output_Xi(omgYdVdx, counts, 'omgYdVdx', .false.,recovery)
        call avg_output_Xi(omgYdVdy, counts, 'omgYdVdy', .false.,recovery)
        call avg_output_Xi(omgYdVdz, counts, 'omgYdVdz', .false.,recovery)
        call avg_output_Xi(omgYdWdy, counts, 'omgYdWdy', .false.,recovery)

        call avg_output_Xi(omgZdUdz, counts, 'omgZdUdz', .false.,recovery)
        call avg_output_Xi(omgZdVdz, counts, 'omgZdVdz', .false.,recovery)
        call avg_output_Xi(omgZdWdx, counts, 'omgZdWdx', .false.,recovery)
        call avg_output_Xi(omgZdWdy, counts, 'omgZdWdy', .false.,recovery)
        call avg_output_Xi(omgZdWdz, counts, 'omgZdWdz', .false.,recovery)

        call avg_output_Xi(omgXomgXdUdx, counts, 'omgXomgXdUdx', .false.,recovery)
        call avg_output_Xi(omgXomgXdVdx, counts, 'omgXomgXdVdx', .false.,recovery)

        call avg_output_Xi(omgXomgYdUdx, counts, 'omgXomgYdUdx', .false.,recovery)
        call avg_output_Xi(omgXomgYdUdy, counts, 'omgXomgYdUdy', .false.,recovery)
        call avg_output_Xi(omgXomgYdVdx, counts, 'omgXomgYdVdx', .false.,recovery)
        call avg_output_Xi(omgXomgYdVdy, counts, 'omgXomgYdVdy', .false.,recovery)

        call avg_output_Xi(omgXomgZdUdz, counts, 'omgXomgZdUdz', .false.,recovery)
        call avg_output_Xi(omgXomgZdVdz, counts, 'omgXomgZdVdz', .false.,recovery)
        call avg_output_Xi(omgXomgZdWdx, counts, 'omgXomgZdWdx', .false.,recovery)

        call avg_output_Xi(omgYomgYdUdy, counts, 'omgYomgYdUdy', .false.,recovery)
        call avg_output_Xi(omgYomgYdVdy, counts, 'omgYomgYdVdy', .false.,recovery)

        call avg_output_Xi(omgYomgZdUdz, counts, 'omgYomgZdUdz', .false.,recovery)
        call avg_output_Xi(omgYomgZdVdz, counts, 'omgYomgZdVdz', .false.,recovery)
        call avg_output_Xi(omgYomgZdWdy, counts, 'omgYomgZdWdy', .false.,recovery)

        call avg_output_Xi(omgZomgZdWdz, counts, 'omgZomgZdWdz', .false.,recovery)
                
        call avg_output_Xi(domgXdxdomgXdx, counts, 'domgXdxdomgXdx', .false.,recovery)
        call avg_output_Xi(domgXdydomgXdy, counts, 'domgXdydomgXdy', .false.,recovery)
        call avg_output_Xi(domgXdzdomgXdz, counts, 'domgXdzdomgXdz', .false.,recovery)

        call avg_output_Xi(domgYdxdomgYdx, counts, 'domgYdxdomgYdx', .false.,recovery)
        call avg_output_Xi(domgYdydomgYdy, counts, 'domgYdydomgYdy', .false.,recovery)
        call avg_output_Xi(domgYdzdomgYdz, counts, 'domgYdzdomgYdz', .false.,recovery)

        call avg_output_Xi(domgZdxdomgZdx, counts, 'domgZdxdomgZdx', .false.,recovery)
        call avg_output_Xi(domgZdydomgZdy, counts, 'domgZdydomgZdy', .false.,recovery)
        call avg_output_Xi(domgZdzdomgZdz, counts, 'domgZdzdomgZdz', .false.,recovery)

        call avg_output_Xi(domgXdxdomgYdx, counts, 'domgXdxdomgYdx', .false.,recovery)
        call avg_output_Xi(domgXdydomgYdy, counts, 'domgXdydomgYdy', .false.,recovery)
        call avg_output_Xi(domgXdzdomgYdz, counts, 'domgXdzdomgYdz', .false.,recovery)


      end subroutine dumpOUTPUT


      subroutine avg_output_Xi(varIN,cntIN,varName,flag,recovery)
        implicit none

        real*8,      intent(in) :: varIN(:,:,:)
        integer,    intent(in) :: cntIN(:,:,:)
        character(*),  intent(in) :: varName
        logical,    intent(in) :: flag
        integer,    intent(in) :: recovery

        character(300)  :: path
        character(300)  :: filename
        type(datafile)  :: file
        integer(hsize_t):: dimsf(3)

        character(300)  :: fnam
        character(8)    :: suffix


        if(recovery==1) then
          path = trim(results_path)//'phase_averaged_fields_recovery/'
            else
          path = trim(results_path)//'phase_averaged_fields/'
        endif

        filename= trim(path)//trim(varName)//'.h5'
        call create_file(filename,file)   

        dimsf(1:3)=(/n1m,n2m,xbins/)


        ! Write simulation parameters: amp, kappa, omega, 
        call h5gcreate_f(file%id,'setup',file%current_group, ierr)
        call write_hdf(file,amp,'amp')
        call write_hdf(file,kappa,'kappa')
        call write_hdf(file,omega,'omega')
        call write_hdf(file,xbins,'xbins')
        call h5gclose_f(file%current_group, ierr)

        call h5gcreate_f(file%id,'grid',file%current_group,ierr)
        call write_hdf(file,xvals,'xvals')
        call write_hdf(file,Yc,'Yc')
        call write_hdf(file,Zc,'Zc')
        call h5gclose_f(file%current_group,ierr)


        call h5gcreate_f(file%id,'fields',file%current_group, ierr)
        call write_hdf(file,varIN,dimsf,trim(varName),zstart(1),min(zend(1),n1m),zstart(2),min(zend(2),n2m),1,xbins)
        if(flag) call write_hdf(file,cntIN,dimsf,'counts',zstart(1),min(zend(1),n1m),zstart(2),min(zend(2),n2m),1,xbins)
        call h5gclose_f(file%current_group, ierr)

        ! Close HDF5 file
        call close_file(file)

      end subroutine avg_output_Xi


      subroutine avg_input_Xi(varIN,cntIN,varName,flag,recovery)
        implicit none

        real*8,     intent(out)  :: varIN(:,:,:)
        integer,    intent(out)  :: cntIN(:,:,:)
        character(*),  intent(in)  :: varName
        logical,    intent(in)  :: flag
        integer,    intent(in)  :: recovery

        character(300)  :: path
        character(300)  :: filename

        type(datafile)  :: file
        integer(hsize_t):: dimsf(3)
        character(8)    :: suffix


        if(recovery==1) then
          path = trim(external_fields_path)//'/phase_averaged_fields_recovery/'
            else
          path = trim(external_fields_path)//'/phase_averaged_fields/'
        endif

        filename= trim(path)//trim(varName)//'.h5'
        call read_file(filename,file)   

        dimsf(1:3)=(/n1m,n2m,xbins/)


        call h5gopen_f(file%id,'fields',file%current_group, ierr)
        call read_hdf(file,varIN,dimsf,trim(varName),zstart(1),min(zend(1),n1m),zstart(2),min(zend(2),n2m),1,xbins)
        if(flag) call read_hdf(file,cntIN,dimsf,'counts',zstart(1),min(zend(1),n1m),zstart(2),min(zend(2),n2m),1,xbins)
        call h5gclose_f(file%current_group, ierr)

        ! Close HDF5 file
        call close_file(file)

      end subroutine avg_input_Xi


end module mod_inOut
