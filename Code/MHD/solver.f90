module MHD_solver
!
!*********** ASSUMING THAT STREAMWISE = 3 ******************
!
    use decomp_2d
    use mesh
    use boundaries
    use schemes3D_interface
    use irregular_derivative_coefficients
    use DNS_settings
    use MHD_data
    use physical_fields
    use mpi

    implicit none

!    real*8, dimension(:,:,:), allocatable       :: previousRHS, previousRHS1, previousRHS2

    real*8  :: gam,rom


!    logical, private     :: previousRHS_are_available=.false.

contains

    subroutine solve_MHD(q1_x, q2_y, q3_z, fb1_x, fb2_x, fb3_x,ntime,last_substep)
!    use poisson_interface, solve_Poisson_ORL => solve_Poisson
!    use poisson_020_solver
        use operators_tools
        use hdf5_io
!    use operators_tools_tests
        use numerical_methods_settings
        use poisson_020_solver, ORL_solve_Poisson=>solve_Poisson
        use poisson_interface, LAMB_solve_Poisson=>solve_Poisson
        use workspace_view, only: results3D_path

      implicit none
      integer, save   :: nbcall=1
      integer, intent(in)   :: ntime
      logical, intent(in)   :: last_substep
      integer :: i,j, k, c2, num_rank
      integer :: n1s, n1e, n2s,n2e, n3s,n3e
      real*8 :: fb1_MHD_tot, intg_fb1_x,intg_fb2_x,intg_fb3_x, Mean_J, Mean_J_span
      real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1_x, fb1_x, fb2_x, fb3_x, RHS_x
      real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q2_y, RHS_y, ucrossB1_y,ucrossB3_y,Norm_y
      real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: gradphi1_y,gradphi3_y
      real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q3_z
! Tecplot Output
      character*10 tmp_str, tmp_num_rank
      character*200 absolute_path,elec_curr_path
      real*8 NormFB
      character*1 coma
      coma=char(44)
!
!       open(155, file='MHDEBUG4', position='append')
!       write(155,*) 'sum w', sum(q3_z)

       call perform_ucrossB(q1_x, q2_y, q3_z, B01c_y, B02c_y, B03c_y, ucrossB1_x, ucrossB2_y, ucrossB3_z, Mean_J, Mean_J_span)
       ! (q2*B03 - q3*B02).e1 + (q3*B01 - q1*B03).e2 +(q1*B02 - q3*B01).e3
        if (nexport_MHD/=0)then
           if ((nrank==0).and.(last_substep).and.(mod(ntime,nexport_MHD).eq.0)) then
                open(759, file="meanJ.dat", position="append")
                write(759,*) ntime, Mean_J_span, Mean_J
                close(759)
           endif
        endif

       call perform_divergence_tool(RHS_z, ucrossB3_z, ucrossB2_y, ucrossB1_x) !div(u x B)

!       write(155,*) 'div(uxB)=', sum(RHS_z)

!     open(28,file='Divergence_SourceMHD.plt')
!     write(28,10)' TITLE = " Div"'
!     write(28,10)' VARIABLES = "X","Y","Z","Div1","Div1","Div1","Div1"'
!     write(28,11) 'ZONE T="Mag"',coma,'I=',n3-1,coma,'J=',n2m,coma,'K=',n1-1,coma,'DATAPACKING=POINT'
!
!! ************ Champs 3D  **************
!
!     do i= zstart(1), min(n1m, zend(1))   !do i=1,n1m
!      do j = zstart(2), min(n2m, zend(2))
!         do k=1,n3m
!
!         NormFB = 0.d0
!!
!         write(28,*)Xc(k),Yc(j),Zc(i),RHS_z(i,j,k),RHS_z(i,j,k),RHS_z(i,j,k),NormFB
!
!        enddo
!       enddo
!      enddo
!
!      close(28)




        if (use_generic_poisson) then
            call LAMB_solve_Poisson(RHS_z, RHS_z2)
            RHS_z=RHS_z2
        else
            call ORL_solve_Poisson(RHS_z, RHS_z2)
            RHS_z=RHS_z2
        end if
!
!
! Check that Integrale_domain(dphi/dx) =-0.5*Boy*Integrale_domain(u)
!
! from here RHS_z is Phi (electrical potential)
!  from spectral resolution of Lamballais etc, only dphi/d(x,y,z)=0 can be computed at walls and open boundary => non conducting wall (uvectoB is nul at wall)
!
!       write(155,*) 'Somme Phi', sum(RHS_z)
!
!***********************************************************************************************************
!************************************* Compute Grad Phi = RHS_z *****************************************

       !  call perform_gradphi_test

         call perform_gradphi(RHS_z,gradphi1_x,gradphi2_y,gradphi3_z)

        !absolute_path='/bettik/capognam/DNS/WORKSPACE/Codes/DNS/&
        !                &MULTIFAST_MHD_v2/TMP/CANAL_carre_lam_mhd_14401070/'
!

!************************************* Compute the MHD Body Force ****************************************
        call MHD_Body_Force(gradphi1_x, gradphi2_y, gradphi3_z, ucrossB1_x, ucrossB2_y, ucrossB3_z, fb1_x, fb2_x, fb3_x, Mean_J, Mean_J_span)
!        call add_action(fb1_MHD_x, fb2_MHD_x, fb3_MHD_x,ntime)
!*************************************************************************************************************
!

        if(mod(ntime,nexport_MHD).eq.0) then
            write(tmp_str, "(i0)")ntime
            elec_curr_path=trim(results3D_path)//'field'//trim(tmp_str)
            n1e=min(n1-1, yend(1))
            n3e=min(n3-1, yend(3))
            call transpose_z_to_y(RHS_z, RHS_y)
            call hdf_write_3Dfield(trim(elec_curr_path)//'/Phi', RHS_y, 'Phi', n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            
            call transpose_x_to_y(ucrossB1_x, ucrossB1_y)
            call transpose_z_to_y(ucrossB3_z, ucrossB3_y)
            do i = ystart(1), n1e
               do j = 1, n2-1
                  do k = ystart(3), n3e
                    Norm_y(i,j,k) = dsqrt(((ucrossB1_y(i,j,k))**2 +(ucrossB2_y(i,j,k))**2 + (ucrossB3_y(i,j,k))**2))
                  enddo
               enddo
            enddo
            
            call hdf_write_3Dfield(trim(elec_curr_path)//'/UcrossB', ucrossB1_y, 'UcrossB1', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            call hdf_add_3Dfield(trim(elec_curr_path)//'/UcrossB', ucrossB2_y, 'UcrossB2', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            call hdf_add_3Dfield(trim(elec_curr_path)//'/UcrossB', ucrossB3_y, 'UcrossB3', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            call hdf_add_3Dfield(trim(elec_curr_path)//'/UcrossB', Norm_y, 'Norm', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            
            call transpose_x_to_y(gradphi1_x, gradphi1_y)
            call transpose_z_to_y(gradphi3_z, gradphi3_y)
            do i = ystart(1), n1e
                do j = 1, n2-1
                    do k = ystart(3), n3e
                        Norm_y(i,j,k) = dsqrt(((gradphi1_y(i,j,k))**2 +(gradphi2_y(i,j,k))**2 + (gradphi3_y(i,j,k))**2))
                    enddo
                enddo
            enddo
            call hdf_write_3Dfield(trim(elec_curr_path)//'/GradPhi', gradphi1_y, 'GradPhi1', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            call hdf_add_3Dfield(trim(elec_curr_path)//'/GradPhi', gradphi2_y, 'GradPhi2', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            call hdf_add_3Dfield(trim(elec_curr_path)//'/GradPhi', gradphi3_y, 'GradPhi3', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            call hdf_add_3Dfield(trim(elec_curr_path)//'/GradPhi', Norm_y, 'Norm', &
                                    &n1-1, n2-1, n3-1, ystart(1), n1e, 1, n2-1, ystart(3), n3e)
            
            !open(28,file=trim(absolute_path)//trim(results3D_path)//&
            !                &'field'//trim(tmp_str)//'/Electric_Potential.plt')
            !do num_rank=0,nproc-1
            !if(nrank==num_rank) then
            !    write(tmp_str, "(i0)")ntime
            !    write(tmp_num_rank, "(i0)")nrank
            !    open(28,file=trim(results3D_path)//'field'//trim(tmp_str)//'/Phi_'//trim(tmp_num_rank)//'.plt')
            !         write(*,*)trim(results3D_path)//'field'//trim(tmp_str)//'/Electric_Potential.plt'
            !         write(28,10)' TITLE = "Electric Potential"'
            !         write(28,10)' VARIABLES = "X","Y","Z","phi"'
            !         write(28,11) 'ZONE T = "Mag"',coma,'I=',1,coma,'J=',n2-1,coma,'K=',n1-1,coma,'DATAPACKING=POINT'

            ! ************ Champs 3D  **************
            !     n1e=min(n1-1, yend(1))
            !     n3e=min(n3-1, yend(3))
            !    call transpose_z_to_y(RHS_z, RHS_y)
            !     do i = ystart(1), n1e
            !        do j = 1, n2-1
            !            do k = ystart(3), n3e
            !               write(28,*)Xc(k),Yc(j),Zc(i),RHS_y(i,j,k)
            !        enddo
            !       enddo
            !      enddo

            !      close(28)

         !open(29,file=trim(results3D_path)//'field'//trim(tmp_str)//'/UcrossB_'//trim(tmp_num_rank)//'.plt')
         !    write(29,10)' TITLE = "UcrossB "'
         !    write(29,10)' VARIABLES = "X","Y","Z","UcrossB1","UcrossB2","UcrossB3","Norm"'
         !    write(29,11) 'ZONE T="Mag"',coma,'I=',n3-1,coma,'J=',n2m,coma,'K=',n1-1,coma,'DATAPACKING=POINT'
         !    call transpose_x_to_y(ucrossB1_x, ucrossB1_y)
         !    call transpose_z_to_y(ucrossB3_z, ucrossB3_y)
         !do i = ystart(1), n1e
         !  do j = 1, n2-1
         !     do k = ystart(3), n3e

         !       NormFB = dsqrt(((ucrossB1_y(i,j,k))**2 +(ucrossB2_y(i,j,k))**2 + (ucrossB3_y(i,j,k))**2))
         !       write(29,*)Xc(k),Yc(j),Zc(i),ucrossB1_y(i,j,k),ucrossB2_y(i,j,k),ucrossB3_y(i,j,k),NormFB
         !   enddo
         !  enddo
         ! enddo

         ! close(29)

              !open(30,file=trim(absolute_path)//trim(results3D_path)//'field'//trim(tmp_str)//'/GradPhi.plt')
         !     open(30,file=trim(results3D_path)//'field'//trim(tmp_str)//'/GradPhi_'//trim(tmp_num_rank)//'.plt')
         !        write(30,10)' TITLE = "GradPhi "'
         !        write(30,10)' VARIABLES = "X","Y","Z","GradPhi1","GradPhi2","GradPhi3","Norm"'
         !        write(30,11) 'ZONE T="Mag"',coma,'I=',1,coma,'J=',n2m,coma,'K=',n1-1,coma,'DATAPACKING=POINT'
         !        call transpose_x_to_y(gradphi1_x, gradphi1_y)
         !        call transpose_z_to_y(gradphi3_z, gradphi3_y)
        ! ************ Champs 3D  **************

         !do i = ystart(1), n1e
         !  do j = 1, n2-1
         !     do k = ystart(3), n3e

         !       NormFB = dsqrt(((gradphi1_y(i,j,k))**2 +(gradphi2_y(i,j,k))**2 + (gradphi3_y(i,j,k))**2))
         !       write(30,*)Xc(k),Yc(j),Zc(i),gradphi1_y(i,j,k),gradphi2_y(i,j,k),gradphi3_y(i,j,k),NormFB

         !       enddo
         !      enddo
         !     enddo

         !     close(30)
        !endif
        !enddo
      end if


!
!              intg_fb1_x=0.d0
!              intg_fb2_x=0.d0
!              intg_fb3_x=0.d0
!!!
!!! Warning for Parralel: ADD the MPI SUM command
!!!
!       do j = zstart(2), min(n2m, zend(2))
!         do i=zstart(1), min(n1m, zend(1))   !do i=1,n1m
!          do k=1,n3m
!              intg_fb1_x=intg_fb1_x + fb1_x(i,j,k)*dx1*dx3*cell_size_Y(j)
!              intg_fb2_x=intg_fb2_x + fb2_x(i,j,k)*dx1*dx3*cell_size_Y(j)
!              intg_fb3_x=intg_fb3_x + fb3_x(i,j,k)*dx1*dx3*cell_size_Y(j)
!            enddo
!           enddo
!       enddo

!        open(1595, file="MHD_Force_History", position="append")
!        write(1595,*)intg_fb1_x, intg_fb2_x, intg_fb3_x
!        close(1595)
        !call check_passive_solver
!      write(155,*) 'Somme Q3 x,y,z', sum(q3_z)
!      close(155)


!************************************************************************
!******************** VISU :CHAMPS INTENSITY ****************************
!************************************************************************

!     open(28,file='Magnetic_Force.plt')
!     write(28,10)' TITLE = "Magnetic Force"'
!     write(28,10)' VARIABLES = "X","Y","Z","FB1","FB2","FB3","NormFB"'
!     write(28,11) 'ZONE T="Mag"',coma,'I=',n3-1,coma,'J=',n2m,coma,'K=',n1-1,coma,'DATAPACKING=POINT'

! ************ Champs 3D  **************

!     do i= zstart(1), min(n1m, zend(1))   !do i=1,n1m
!      do j = zstart(2), min(n2m, zend(2))
!         do k=1,n3m

!         NormFB = dsqrt(((fb1_x(i,j,k))**2 +(fb2_x(i,j,k))**2 + (fb3_x(i,j,k))**2))
!
!         write(28,*)Xc(k),Yc(j),Zc(i),fb1_x(i,j,k),fb2_x(i,j,k),fb3_x(i,j,k),NormFB

!        enddo
!       enddo
!      enddo

!      close(28)

!     open(29,file='UcrossB.plt')
!     write(29,10)' TITLE = "UcrossB "'
!     write(29,10)' VARIABLES = "X","Y","Z","UcrossB1","UcrossB2","UcrossB3","Norm"'
!     write(29,11) 'ZONE T="Mag"',coma,'I=',n3-1,coma,'J=',n2m,coma,'K=',n1-1,coma,'DATAPACKING=POINT'

! ************ Champs 3D  **************

!     do i= zstart(1), min(n1m, zend(1))   !do i=1,n1m
!      do j = zstart(2), min(n2m, zend(2))
!         do k=1,n3m

!         NormFB = dsqrt(((ucrossB1_x(i,j,k))**2 +(ucrossB2_y(i,j,k))**2 + (ucrossB3_z(i,j,k))**2))
!
!         write(29,*)Xc(k),Yc(j),Zc(i),ucrossB1_x(i,j,k),ucrossB2_y(i,j,k),ucrossB3_z(i,j,k),NormFB

!        enddo
!       enddo
!      enddo

!      close(29)



   10 format(a)
   11 format(a,a,a,i4,a,a,i3,a,a,i3,a,a)
   12 format(a,a,a,i4,a,a,i3,a,a)

    end subroutine solve_MHD
!
!*************************************************************************************************
!***************************| MHD BODY FORCE CALCULATION |******************************************
!*************************************************************************************************
!
    subroutine MHD_Body_Force(gradphi1_x, gradphi2_y, gradphi3_z, ucrossB1_x, ucrossB2_y, ucrossB3_z, fb1_x, fb2_x, fb3_x,Mean_J, Mean_J_span)

    implicit  none

        integer ::  n1e, n2e, n3e,n1s, n2s, n3s
        integer :: i,j,k
        real*8 :: Mean_J, Mean_J_span
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in) :: gradphi1_x, ucrossB1_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in) :: gradphi2_y, ucrossB2_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in) :: gradphi3_z, ucrossB3_z

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: A1int_x, fb1int_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: A2int_y, fb2int_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: A3int_z, fb3int_z, fb3_z

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: fb1_x,  fb2_x,  fb3_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: fb1_y,  fb2_y,  fb3_y
        integer, save :: cpt=1

! INITIALISATION SAUVAGE !!!

!        fb1int_x = 0.d0
        fb1_y = 0.d0
        fb2_y = 0.d0
 !       fb2int_y = 0.d0
        fb3_y =0.d0
 !       fb3_z = 0.d0
!        fb3int_z = 0.d0
        A2int_y = 0.d0
        A3int_z = 0.d0
        A1int_x = 0.d0

!            open(154, file='MHD_debug', position='append')
!            write(154,*)
!            write(154,*) 'Point Aa', cpt, sum(gradphi1_x(1:n1-1, 1:n2-1, 1:n3-1)), sum(ucrossB1_x(1:n1-1, 1:n2-1, 1:n3-1))
!            write(154,*) 'Point Ab', cpt, sum(gradphi2_y(1:n1-1, 1:n2-1, 1:n3-1)), sum(ucrossB2_y(1:n1-1, 1:n2-1, 1:n3-1))
!            write(154,*) 'Point Ac', cpt, sum(gradphi3_z(1:n1-1, 1:n2-1, 1:n3-1)), sum(ucrossB3_z(1:n1-1, 1:n2-1, 1:n3-1))
!
! Position of A1_x, A2_y, A3_z is here : FACE CENTERED
!
! A1, A2 and A3 are the compents of the electrical current
!
           A2_y = -gradphi2_y + ucrossB2_y

         if(streamwise==1) then
            A1_x = -(gradphi1_x + Mean_J)  + ucrossB1_x
            A3_z = -(gradphi3_z + Mean_J_span) + ucrossB3_z

         elseif(streamwise==3) then
            A1_x = -(gradphi1_x + Mean_J_span)  + ucrossB1_x
            A3_z = -(gradphi3_z + Mean_J) + ucrossB3_z
         endif

 !          write(154,*) 'Point A', cpt, sum(A1_x), sum(A2_y), sum(A3_z)

        ! ************ CORR_JO ************ !
        ! n1e et n3e comme dans perform_RHS !
        ! ********************************* !

        n1e=ysize(1) !borne y
        n3e=ysize(3) !borne y
 !
        ! ***** CORR_JO ****** !
        ! D0s au lieu de D0ssh !
        ! ******************** !

        call D0s_3Dy(A2_y, A2int_y, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)

        n1e=zsize(1) !borne z
        n2e=zsize(2) !borne z
        call D0s_3Dz(A3_z, A3int_z, zsize(1),n1e,zsize(2),n2e,n3, NS_Q3_BC3)
        call transpose_z_to_y(A3int_z,A3_y)

        n2e=xsize(2) !borne x
        n3e=xsize(3) !borne x
        call D0s_3Dx(A1_x, A1int_x, n1,xsize(2),n2e,xsize(3),n3e, NS_Q1_BC1)
        call transpose_x_to_y(A1int_x, A1_y)


! Position of A1_y, A2int_y, A3_y is here : CELL CENTERED

!  WARNING : For the body force fb should be transposed on x direction
! ***************************************************************************
!
       n1s=ystart(1)
       n1e=min(n1-1, yend(1))
!
       n2s=ystart(2)
       n2e=min(n2-1, yend(2))
!
       n3s=ystart(3)
       n3e=min(n3-1, yend(3))

      fb1_y=0.d0
      fb2_y=0.d0
      fb3_y=0.d0

 !    open(155, file='MHD_Fb1_debug', position='append')
 !    write(155,*)
       do k=n3s, n3e
         do i=n1s, n1e
          do j=n2s,n2e
            fb1_y(i, j, k) =  Stuart_number*(  ( A2int_y(i, j, k)*B03c_y(i, j, k) ) - ( A3_y(i, j, k)*B02c_y(i, j, k) )  )
            fb2_y(i, j, k) =  Stuart_number*(  ( A3_y(i, j, k)*B01c_y(i, j, k) ) - ( A1_y(i, j, k)*B03c_y(i, j, k) )  )
            fb3_y(i, j, k) =  Stuart_number*(  ( A1_y(i, j, k)*B02c_y(i, j, k) ) - ( A2int_y(i, j, k)*B01c_y(i, j, k) )  )
 !          write(155,*) 'B01_y,B02_y,B03_y =', B01_y(1, 1, 1), B02_y(1, 1, 1), B03_y(1, 1, 1)
 !          write(155,*) 'i,j,k,fb =',i,j,k,fb1_y(i, j, k),fb2_y(i, j, k),fb3_y(i, j, k)
!          write(155,*) 'i,j,k,gradphi=',i,j,k,gradphi1_x(i, j, k),gradphi2_y(i, j, k),gradphi3_z(i, j, k)
          enddo
         enddo
       enddo
!
  !    close(155)
!
!       Position of fb1_y, fb2_y, fb3_y is here : FACE CENTERED
!       write(154,*) 'Point B', cpt, sum(fb3_y), sum(A1_y), sum(A2_y), sum(B02_y), sum(B01_y)

        ! ************************** CORR_JO **************************** !
        ! fb needs to be interpolated from CELL CENTERED to FACE CENTERED !
        ! *************************************************************** !

!
        fb2int_y = fb2_y

        call transpose_y_to_z(fb3_y, fb3int_z)
        call transpose_y_to_x(fb1_y, fb1int_x)

       n2e=xsize(2) !borne x
       n3e=xsize(3) !borne x

      ! D0ssh because CENTER TO FACE interpolation
      ! WARNING : When D0ssh is used with a Dirichlet condition (e.g. at the walls),
      !  the value of the variable at the boundary must be specified elsewhere.
      ! By default, fb = 0.d0, which is ok at the walls (cause the velocity = 0)
      call D0ssh_3Dx(fb1int_x, fb1_x, n1, xsize(2),n2e, xsize(3),n3e, NS_Q1_BC1)
    !
           n1e=ysize(1) !borne y
           n3e=ysize(3) !borne y
      call D0ssh_3Dy(fb2int_y, fb2_y, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)
    !
           n1e=zsize(1) !borne z
           n2e=zsize(2) !borne z
      call D0ssh_3Dz(fb3int_z, fb3_z, zsize(1),n1e,zsize(2),n2e,n3, NS_Q3_BC3)
!
! ****  fb should be in x decomposition *****
        call transpose_y_to_x(fb2_y, fb2_x)
        call transpose_z_to_y(fb3_z, fb3_y)
        call transpose_y_to_x(fb3_y, fb3_x)


 !      write(154,*) 'Point C', cpt, sum(fb1_x), sum(fb2_x), sum(fb3_x)

!       cpt=cpt+1
!       close(154)

!        fb1_MHD_x=0.d0
!        fb2_MHD_x=0.d0
!        fb3_MHD_x=0.d0
    end subroutine MHD_Body_Force


end module MHD_solver
