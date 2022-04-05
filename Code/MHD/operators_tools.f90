module operators_tools

    use decomp_2d
    use mesh
    use boundaries
    use schemes3D_interface
    use irregular_derivative_coefficients
    use DNS_settings
!    use SCALAR_data
    use MHD_data
!    use physical_fields
    use schemes_loader

    implicit none


!    real*8, dimension(:,:,:), allocatable       :: previousRHS, previousRHS1, previousRHS2
!    real*8  :: gam,rom


contains
!************************************************************************************
!**************************************************************************************

     subroutine perform_divergence_tool(div_z, q3_z, q2_y, q1_x, divy_mean)
        use numerical_methods_settings
        use snapshot_writer
        use mpi
        use schemes_interface

        implicit none
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(in)    :: q3_z
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(in)    :: q2_y
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(in)    :: q1_x

        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))              :: div1_x
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(out) :: div_z

        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))              :: div3_z, div2_z, div1_z, div_z2

        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))   :: div2_y, div1_y
        real*8      :: divy_sum
        real*8, optional :: divy_mean

        integer k,j,i, mpi_err, n1e,n2e,n3e
        integer,save    :: nb=1
        !!!$    integer :: nb_taches
        !!!$    common/nb_total_threads/nb_taches


        ! X orientation ---------------------------------------------------------
        div1_x=0.d0
        !do k=1,n3m
        n2e=xsize(2)!(min(n2m, xend(2))-xstart(2))+1
        n3e=xsize(3)!(min(n3m, xend(3))-xstart(3))+1

        call D1s_3Dx(q1_x, div1_x, n1,xsize(2),n2e,xsize(3),n3e, dx1, .false., POISSON_VEL_BC1)

        ! Y orientation ---------------------------------------------------------

        div2_y=0.d0
        n1e=ysize(1)!(min(n1m, yend(1))-ystart(1))+1
        n3e=ysize(3)!(min(n3m, yend(3))-ystart(3))+1

        if (use_generic_poisson) then

            call D1s_MULT_3Dy(q2_y, div2_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .false., POISSON_VEL_BC2, Yc_to_YcTr_for_D1)

        else    ! ONLY FOR CHANNEL FLOWS

            do k=ystart(3), min(n3m, yend(3))   !do k=1,n3m
                do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m
                    div2_y(i,1,k)=(q2_y(i,2,k) - q2_y(i,1,k))*Yc_to_YcTr_for_D1(1)/dx2
                    div2_y(i,n2m,k)=(q2_y(i,n2,k) - q2_y(i,n2m,k))*Yc_to_YcTr_for_D1(n2m)/dx2
                    call D1s_Tamm_MULT(q2_y(i,2:n2m,k), div2_y(i,2:n2m-1,k), n2m-1, dx2, .false., Dirichlet, Yc_to_YcTr_for_D1(2:n2m-1))
                enddo
            enddo

        end if


        ! Z orientation ---------------------------------------------------------

        div3_z=0.d0
        n1e=zsize(1)!(min(n1m, zend(1))-zstart(1))+1
        n2e=zsize(2)!(min(n2m, zend(2))-zstart(2))+1

        call D1s_ACC_3Dz(q3_z, div3_z, zsize(1),n1e,zsize(2),n2e,n3, dx3, .false., POISSON_VEL_BC3)

        call transpose_x_to_y(div1_x, div1_y)
        call transpose_y_to_z(div1_y, div1_z)

        call transpose_y_to_z(div2_y, div2_z)

        div_z=div1_z+div2_z+div3_z


        if (present(divy_mean)) then

            divy_mean=0.d0
            divy_sum=sum(div2_z)

            call MPI_ALLREDUCE (divy_sum, divy_mean, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
            divy_mean=divy_mean/((n1-1)*(n2-1)*(n3-1))

        endif

        return
    end subroutine perform_divergence_tool


    subroutine perform_gradphi(phi_z,gradphi1_x,gradphi2_y,gradphi3_z)

     use numerical_methods_settings
     use schemes_interface

     implicit  none
     integer :: n1e, n2e, n3e
     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: gradphi1_x, phi_x
     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: gradphi2_y, phi_y
     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: gradphi3_z
     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: phi_z
     integer, save :: cpt=1
     integer :: i,k

!    open(154, file='MHDEBUG2', position='append')
!    write(154,*)
!    write(154,*) 'point A sum_phiz', cpt, sum(phi_z), sum(phi_z(1:n1-1, 1:n2-1, 1:n3-1))

!     TRANSPORT_SCA_BC3 = periodic
     n1e=(min(n1m, zend(1))-zstart(1))+1 !borne z
     n2e=(min(n2m, zend(2))-zstart(2))+1 !borne z
     call D1ssh_3Dz(phi_z, gradphi3_z, zsize(1),n1e,zsize(2),n2e,n3, dx3, .true., periodic) ! Attention cr�er des conditions aux limites sp�cifiques ex : MHD_SCA_BC3

 !    TRANSPORT_SCA_BC2 = symetric
     call transpose_z_to_y(phi_z, phi_y)
!     n1e=(min(n1m, yend(1))-ystart(1))+1 !borne y
!     n3e=(min(n3m, yend(3))-ystart(3))+1 !borne y
!     call D1ssh_MULT_3Dy(phi_y, gradphi2_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., symetric,Y_to_YTr_for_D1)

    n1e=ysize(1)
    n3e=ysize(3)

    if (use_generic_poisson) then
        call D1ssh_MULT_3Dy(phi_y, gradphi2_y, ysize(1),n1e,n2,ysize(3),n3e, dx2, .true., POISSON_PR_BC2, Y_to_YTr_for_D1)
    else    ! ONLY for channel flow. In this case, dirichlet boundary cond. are applied in Y direction

        do k=ystart(3), min(n3m, yend(3))
            do i=ystart(1), min(n1m, yend(1))   !do i=1,n1m
                call D1s_Tamm_MULT(phi_y(i,:,k), gradphi2_y(i,:,k), n2, dx2, .true., Dirichlet, Y_to_YTr_for_D1)
            enddo
        enddo

    end if
 !    TRANSPORT_SCA_BC1 = periodic
     call transpose_y_to_x(phi_y, phi_x)


     n2e=(min(n2m, xend(2))-xstart(2))+1 !borne x
     n3e=(min(n3m, xend(3))-xstart(3))+1 !borne x
     call D1ssh_3Dx(phi_x, gradphi1_x, n1,xsize(2),n2e,xsize(3),n3e, dx1, .true., periodic)


!    write(154,*) 'point B_sumgradphi2_y', cpt, sum(gradphi2_y)

!    cpt=cpt+1
!    close(154)

    end subroutine perform_gradphi

!********************************************************************
!********************************************************************
    subroutine perform_ucrossB(u_x, vint_y, w_z, B1_y, B2_y, B3_y, uvectoB1_x, uvectoB2_y, uvectoB3_z,Mean_current,Mean_spanwise_current)
! Should call B1c B2c etc.
        use mpi
        implicit  none

        integer :: i, j, k, mpi_err
        integer :: n1e, n2e, n3e, n1s, n2s, n3s
        real*8 :: Mean_current, Mean_spanwise_current, Mean_current_glob, Mean_spanwise_current_glob
        real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: u_x, uint_x, uvectoB1_x, uvectoB1int_x
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: B1_y,B2_y,B3_y,v_y, u_y, w_y
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: vint_y
        real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: uvectoB2_y,uvectoB1_y, uvectoB3_y,uvectoB2int_y,uvectoB1int_y, uvectoB3int_y
        real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: w_z, uvectoB3_z,wint_z,uvectoB3int_z
        integer, save :: cpt=1
 !       real*8, dimension(:,:,:), allocatable, save :: uvectoB1_y, uvectoB3_y

 !           if (.not. allocated(uvectoB1_y)) then
 !               allocate(uvectoB1_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
 !               allocate(uvectoB3_y(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
 !               uvectoB1_y=0.d0
 !               uvectoB3_y=0.d0
 !           endif
 !    open(154, file='MHDEBUG3', position='append')
 !    write(154,*)
 !   write(154,*) 'point A', cpt, sum(B1_x), sum(B2int_y), sum(B3_z)

        ! !*** ALL vectors variables are put at the cell center & in a  y decomposition ****
!            B2_y=0.d0
            v_y=0.d0
            n1e=(min(n1m, yend(1))-ystart(1))+1
            n3e=(min(n3m, yend(3))-ystart(3))+1
!  D0s because Face to center interpolation ! + Special Treatment for v_y and Bint_y
!********* Transposition de B à mettre dans initializer, inutile de la calculer à chaque fois ********
       call D0s_3Dy(vint_y, v_y, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)

!********* Transposition de B à mettre dans initializer, inutile de la calculer à chaque fois ********
!       call D0s_3Dy(B2int_y, B2_y, ysize(1),n1e,n2,ysize(3),n3e, Dirichlet)

!            B3int_z=0.d0
            wint_z=0.d0
            n1e=(min(n1m, zend(1))-zstart(1))+1
            n2e=(min(n2m, zend(2))-zstart(2))+1
       call D0s_3Dz(w_z, wint_z, zsize(1),n1e,zsize(2),n2e,n3, NS_Q3_BC3)
       call transpose_z_to_y(wint_z,w_y) !=> w_y is cell-centered

 !********* Transposition de B à mettre dans initializer, inutile de la calculer à chaque fois ********
!       call D0s_3Dz(B3_z, B3int_z, zsize(1),n1e,zsize(2),n2e,n3, periodic)
!       call transpose_z_to_y(B3int_z,B3_y)
!
!
!            B1int_x=0.d0
            uint_x=0.d0
            n2e=(min(n2m, xend(2))-xstart(2))+1
            n3e=(min(n3m, xend(3))-xstart(3))+1

       call D0s_3Dx(u_x, uint_x, n1, xsize(2),n2e,xsize(3),n3e, NS_Q1_BC1)
       call transpose_x_to_y(uint_x, u_y)


!********* Transposition de B à mettre dans initializer, inutile de la calculer à chaque fois ********
!       call D0s_3Dx(B1_x, B1int_x, n1,xsize(2),n2e,xsize(3),n3e, periodic)
!       call transpose_x_to_y(B1int_x, B1_y)


       ! (q2*B03 - q3*B02).e1 + (q3*B01 - q1*B03).e2 +(q1*B02 - q3*B01).e3

        n1s=ystart(1)               !A VERIFIER
        n1e=min(n1-1, yend(1))      !A VERIFIER
!
        n2s=ystart(2)               !A VERIFIER
        n2e=min(n2-1, yend(2))      !A VERIFIER
!
        n3s=ystart(3)               !
        n3e=min(n3-1, yend(3))      !
!
!
        uvectoB1int_y=0.d0
        uvectoB2int_y=0.d0
        uvectoB3int_y=0.d0

!        open(151, file='MHD_CHECKSUM5', position='append')
!       write(151,*) 'CheckSum u_y', sum(u_y)
!       write(151,*) 'CheckSum v_y', sum(v_y)
!       write(151,*) 'CheckSum w_y', sum(w_y)
!       write(151,*) 'CheckSum u_y', sum(u_y)
!       write(151,*) 'CheckSum v_y', sum(v_y)
!       write(151,*) 'CheckSum w_y', sum(w_y)
!       close(151)
!
             do k=n3s, n3e
               do i=n1s, n1e
                 do j=n2s, n2e
                   uvectoB1int_y(i, j, k)= (v_y(i, j, k)*B3_y(i, j, k)) - (w_y(i, j, k)*B2_y(i, j, k))
                   uvectoB2int_y(i, j, k)= (w_y(i, j, k)*B1_y(i, j, k)) - (u_y(i, j, k)*B3_y(i, j, k))
                   uvectoB3int_y(i, j, k)= (u_y(i, j, k)*B2_y(i, j, k)) - (v_y(i, j, k)*B1_y(i, j, k))
                 enddo
                enddo
              enddo
!

!
! *************** Compute (uvectoB)_mean to correct gradphi in order to ensure Integ_domain(J)=0
!
        Mean_current=0.d0
        Mean_spanwise_current=0.d0

        do k=ystart(3), min(n3m, yend(3))
            do j=1, n2m
                do i=ystart(1), min(n1m, yend(1))

                    if (streamwise==1) then
                        Mean_current= Mean_current+(uvectoB1int_y(i,j,k)) * cell_size_Y(j) * dx1 * dx3
                        Mean_spanwise_current=Mean_spanwise_current+(uvectoB3int_y(i, j, k)) * cell_size_Y(j) * dx1 * dx3
                    elseif (streamwise==3) then
                        Mean_current= Mean_current+(uvectoB3int_y(i, j, k)) * cell_size_Y(j) * dx1 * dx3
                        Mean_spanwise_current=Mean_spanwise_current+(uvectoB1int_y(i,j,k)) * cell_size_Y(j) * dx1 * dx3
                    endif

                enddo
            enddo
        enddo

        call MPI_ALLREDUCE (Mean_current, Mean_current_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
        call MPI_ALLREDUCE (Mean_spanwise_current, Mean_spanwise_current_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        Mean_current = Mean_current_glob/(L2*L1*L3)
        Mean_spanwise_current = Mean_spanwise_current_glob/(L2*L1*L3)

!*************************************************************************

!    write(154,*) 'point B0', cpt, sum(u_y), sum(v_y),sum(w_y)
!    write(154,*) 'point B1', cpt, sum(uvectoB1int_y), sum(uvectoB2int_y),sum(uvectoB3int_y)
!    cpt=cpt+1
!    close(154)


          call transpose_y_to_z(uvectoB3int_y, uvectoB3int_z)
!
!
          call transpose_y_to_x(uvectoB1int_y, uvectoB1int_x)


       ! + interpolation vers les FACES

               n2e=xsize(2) !borne x
               n3e=xsize(3) !borne x

          ! D0ssh because CENTER TO FACE interpolation
          ! WARNING : When D0ssh is used with a Dirichlet condition (e.g. at the walls),
          !  the value of the variable at the boundary must be specified elsewhere.
          ! By default, uvectoB = 0.d0, which is ok at the walls (cause the velocity = 0)
          call D0ssh_3Dx(uvectoB1int_x, uvectoB1_x, n1, xsize(2),n2e, xsize(3),n3e, NS_Q1_BC1)
!
               n1e=ysize(1) !borne y
               n3e=ysize(3) !borne y
          call D0ssh_3Dy(uvectoB2int_y, uvectoB2_y, ysize(1),n1e,n2,ysize(3),n3e, NS_Q2_BC2)
!
               n1e=zsize(1) !borne z
               n2e=zsize(2) !borne z
          call D0ssh_3Dz(uvectoB3int_z, uvectoB3_z, zsize(1),n1e,zsize(2),n2e,n3, NS_Q3_BC3)


    end subroutine perform_ucrossB


 end module operators_tools

! module operators_tools_tests
!
!    use decomp_2d
!    use mesh
!    use boundaries
!    use schemes3D_interface
!    use irregular_derivative_coefficients
!    use DNS_settings
!!    use SCALAR_data
!    use schemes_loader
!    use operators_tools
!
!    implicit none
!
!
!
!contains
!
!    subroutine perform_gradphi_test()
!
!     implicit  none
!     integer :: i,j,k
!     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: phi_test
!     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: grad1_x, grad1_x_theo
!     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: grad2_y, grad2_y_theo
!     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: grad3_z, grad3_z_theo
!
!
!     do i= 1, n1m
!      do k= 1,n3m
!       do j = 1, n2m
!            phi_test(i,j,k) = 2.d0*YC(j)
!
!            grad1_x_theo(i,j,k) = 0.d0
!            grad2_y_theo(i,j,k) = 2.d0
!            grad3_z_theo(i,j,k) = 0.d0
!  !          dfy_th(j)=12*y(j)**3+15*y(j)**2+4*y(j)
!  !          ddfy_th(j)=36*y(j)**2+30*y(j)+4
!        end do
!       end do
!      end do
!
!    write(6,*) 'CK1_2'
!    call perform_gradphi(phi_test,grad1_x,grad2_y,grad3_z)
!    write(6,*) 'CK1_3'
!    write(6,*)'grad1_x(64,64,64)=',grad1_x(n1m/2,n2m/2,n3m/2)
!    write(6,*)'grad2_y(64,64,64)=',grad2_y(n1m/2,n2m/2,n3m/2)
!    write(6,*)'grad3_z(64,64,64)=',grad3_z(n1m/2,n2m/2,n3m/2)
!
!    write(6,*)'grad1_x_theo(64,64,64)=',grad1_x_theo(n1m/2,n2m/2,n3m/2)
!    write(6,*)'grad2_y_theo(64,64,64)=',grad2_y_theo(n1m/2,n2m/2,n3m/2)
!    write(6,*)'grad3_z_theo(64,64,64)=',grad3_z_theo(n1m/2,n2m/2,n3m/2)
!
!    stop
!    end subroutine perform_gradphi_test
!
!
!
!
!    subroutine perform_ucrossB_test()
!
!     implicit  none
!     integer :: i,j,k
!!     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: phi_test
!     real*8 :: ucrossB1_x_theo, ucrossB2_y_theo, ucrossB3_z_theo
!     real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) :: q1_x, B01_x, ucrossB1_x
!     real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)) :: q2_y, B02_y, ucrossB2_y
!     real*8, dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)) :: q3_z, B03_z, ucrossB3_z
!
!     q1_x = 0.d0
!     q2_y = 0.d0
!     q3_z = 0.d0
!
!     do i= 1, n1m
!      do k= 1, n3m
!       do j= 1, n2m
!            q1_x(i,j,k) =0.d0
!            q2_y(i,j,k) =0.d0
!            q3_z(i,j,k) =0.4d0
!            B01_x(i,j,k) =1.d0
!            B02_y(i,j,k) =0.d0
!            B03_z(i,j,k) =0.d0
!!            phi_test(i,j,k) = 2.d0*YC(j)
!
!        end do
!       end do
!      end do
!
!
!    write(6,*) 'CK1_2'
!
!    call perform_ucrossB(q1_x, q2_y, q3_z, B01_x, B02_y, B03_z, ucrossB1_x, ucrossB2_y, ucrossB3_z)
!
!    write(6,*)'ucrossB1_x=',ucrossB1_x(n1m/2,n2m/2,n3m/2)
!    write(6,*)'ucrossB2_y=',ucrossB2_y(n1m/2,n2m/2,n3m/2)
!    write(6,*)'ucrossB3_z=',ucrossB3_z(n1m/2,n2m/2,n3m/2)
!
!    ucrossB1_x_theo = (q2_y(n1m/2,n2m/2,n3m/2)*B03_z(n1m/2,n2m/2,n3m/2)) - (q3_z(n1m/2,n2m/2,n3m/2)*B02_y(n1m/2,n2m/2,n3m/2))
!    ucrossB2_y_theo = (q3_z(n1m/2,n2m/2,n3m/2)*B01_x(n1m/2,n2m/2,n3m/2)) - (q1_x(n1m/2,n2m/2,n3m/2)*B03_z(n1m/2,n2m/2,n3m/2))
!    ucrossB3_z_theo = (q1_x(n1m/2,n2m/2,n3m/2)*B02_y(n1m/2,n2m/2,n3m/2)) - (q2_y(n1m/2,n2m/2,n3m/2)*B01_x(n1m/2,n2m/2,n3m/2))
!    write(6,*)'ucrossB1_theo_x=', ucrossB1_x_theo
!    write(6,*)'ucrossB2_theo_y=',ucrossB2_y_theo
!    write(6,*)'ucrossB3_theo_z=',ucrossB3_z_theo
!
!    stop
!    end subroutine perform_ucrossB_test
!
!end module operators_tools_tests
