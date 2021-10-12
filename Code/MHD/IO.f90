module MHD_settings
    implicit none

contains

    subroutine read_settings()

        use COMMON_workspace_view
        use DNS_settings
        use boundaries

        use MHD_data

        implicit none
        integer         :: n, magnet_nb

! MHD *********************************************************************************
        open(15,file=trim(COMMON_settings_path)//'MHD.d')

        read(15,*) MHD_state                ! 0: do not solve MHD ; 1: active MHD (RK3 ou AB2) ; 2: active MHD (Euler1) ; 3: passive MHD (Fext=0)
  !      read(15,*) SCA_POT_BC1, SCA_POT_BC2, SCA_POT_BC3  ! SCA_POT est le potentiel (cf LEE CHOI JFM 2001) ! WARNING DIRECTLY LINKED TO THE PRESSURE EQ
        read(15,*) Hartmann_number                          !=> Leads also to the Hartmann number by combinng with the Reynolds number
!        read(15,*) Monodirec_Magnetic_field               ! 0:OFF 1:ON
        read(15,*) Magnetic_field_unit_vector_x
        read(15,*) Magnetic_field_unit_vector_y
        read(15,*) Magnetic_field_unit_vector_z
        read(15,*)
        read(15,*) layer_type ! Apply for dir. 2 => 0 : uniform field, 1 : sinusoïdal, 2 : bipolar
        read(15,*) delta_B  ! delta_B, length (y-axis) of the layer of magnetic field ( /!\ applied on both walls /!\ )
        read(15,*) smooth
        read(15,*) num_period ! Apply for sinusoïdal B2 field => number of periods over channel length
        read(15,*)
    ! Warning the B0 'applied magnetic field' vector should be unit vector due to non dimensionalization
    ! One direction only for the moment (eiher x , y or z)
!        read(15,*) Array_activated ! 0:OFF 1:ON
!        read(15,*) nb_magnet_dipoles  ! number of magnets dipoles,  1 by default

        read(15,*) magnet_number_pairs    ! number of magnets = EVEN only (div B = 0 )
        magnet_nb = magnet_number_pairs *2
        read(15,*) (mag_array(n)%sigma,n=1,magnet_nb)                 ! -1:Magnetic Sink   1:Magnetic Source
        read(15,*) (mag_array(n)%magnet_size_x,n=1,magnet_nb)   ! WARNING *pi in subroutine generate_fields_1_magnet
        read(15,*) (mag_array(n)%magnet_size_z,n=1,magnet_nb)    ! WARNING *pi in subroutine generate_fields_1_magnet
        read(15,*) (mag_array(n)%magnet_center_x,n=1,magnet_nb)  ! WARNING *pi in subroutine generate_fields_1_magnet
        read(15,*) (mag_array(n)%magnet_center_z,n=1,magnet_nb)  ! WARNING *pi in subroutine generate_fields_1_magnet
        read(15,*) (mag_array(n)%magnet_center_y,n=1,magnet_nb)  !
        read(15,*)
        read(15,*) nexport_MHD  ! frequency of data export concerning MHD (gradP.dat, meanJ.dat ...)
        read(15,*) MHD_export_3D ! exporting fields (Elec_curr), 0:OFF 1:ON

        close(15)

        Stuart_number = Hartmann_number*Hartmann_number/ren

 !**!**!**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       write(*,*)'Stuart_number =',  Stuart_number
 !**!**!**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end subroutine read_settings


end module MHD_settings

module MHD_dao

    use mpi
    use decomp_2d
    use decomp_2d_io
    use MHD_data

    implicit none

contains

    subroutine write_fields(fields_dir)
!
        use physical_fields
        use mesh
        use DNS_settings
        use HDF5_IO
!
        implicit none
        character(*)    :: fields_dir
!
!
        character(200)    :: file_path
!
!
        file_path=trim(fields_dir)//"/Elec_Cur_J3"
        call hdf_write_3Dfield(file_path, A3_z(:,:,:), "J3", nx_global, ny_global, nz_global, zstart(1),zend(1),zstart(2),zend(2),zstart(3),zend(3))
!
        file_path=trim(fields_dir)//"/Elec_Cur_J2"
        call hdf_write_3Dfield(file_path, A2_y(:,:,:), "J2", nx_global, ny_global, nz_global, ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))
!
        file_path=trim(fields_dir)//"/Elec_Cur_J1"
        call hdf_write_3Dfield(file_path, A1_x(:,:,:), "J1", nx_global, ny_global, nz_global, xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))



!        ! Save boundary conditions
!        file_path=trim(fields_dir)//"/ScaWall"
!        if(nrank==0)  call hdf_create_emptyfile(file_path)
!        if(nrank==0)  call hdf_addgroup(file_path, "Wall10")
!        call hdf_add_2Dfield(file_path, sca_wall10(:,:), "Wall10/sca", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
!        if(nrank==0)  call hdf_addgroup(file_path, "Wall11")
!        call hdf_add_2Dfield(file_path, sca_wall11(:,:), "Wall11/sca", ny_global, nz_global, xstart(2),xend(2), xstart(3),xend(3))
!
!        if(nrank==0)  call hdf_addgroup(file_path, "Wall20")
!        call hdf_add_2Dfield(file_path, sca_wall20(:,:), "Wall20/sca", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
!        if(nrank==0)  call hdf_addgroup(file_path, "Wall21")
!        call hdf_add_2Dfield(file_path, sca_wall21(:,:), "Wall21/sca", nx_global, nz_global, ystart(1),yend(1), ystart(3),yend(3))
!
!        if(nrank==0)  call hdf_addgroup(file_path, "Wall30")
!        call hdf_add_2Dfield(file_path, sca_wall30(:,:), "Wall30/sca", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
!        if(nrank==0)  call hdf_addgroup(file_path, "Wall31")
!        call hdf_add_2Dfield(file_path, sca_wall31(:,:), "Wall31/sca", nx_global, ny_global, zstart(1),zend(1), zstart(2),zend(2))
!
    end subroutine write_fields

end module MHD_dao



module MHD_results_writer

    use mpi
    use decomp_2d
    use decomp_2d_io
    use MHD_data

    implicit none

contains

    ! Export all physical field (velocity and scalar) in HDF5 format in individual files
    ! More generic, better than 2decomp format (not generic)
    subroutine write_fields(it)
!
        use run_ctxt_data
        use COMMON_workspace_view
        use MHD_dao, only: DAO_write_fields=>write_fields
!
        implicit none
!
        integer, intent(in) :: it
!
        character*10 tmp_str
        character*80 result_dir_path
!
        write(tmp_str, "(i10)")it
        result_dir_path=trim(COMMON_results3D_path)//'field'//trim(adjustl(tmp_str))
!
        call DAO_write_fields(result_dir_path)
!
    end subroutine write_fields

end module MHD_results_writer
