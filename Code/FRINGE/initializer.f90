module FRINGE_initializer
    use decomp_2d
    use mpi
    use mesh
    use FRINGE_data
    use HDF5_IO
    use DNS_settings
    use boundaries
    use physical_fields

    use SCALAR_data, only: SCA_state, delta_T, sca_x

!
    implicit none
!
contains
!
!
    subroutine get_inflow()

        implicit none
        integer :: j, k

        if (inflow_type==2) then
            call inflow_from_file
        else if (inflow_type==3) then
            call inflow_from_prev_run
        else
            call default_inflow
        endif

        ! Assume streamwise direction is 1
        ! so streamwise velocity is q1

        initial_flowrate=0.d0

        ! do k=xstart(3),min(xend(3),n3-1)
        !     do j=xstart(2),min(xend(2),n2-1)

        !         initial_flowrate = initial_flowrate + q1_inflow(j,k)*(Y(j+1)-Y(j))*dx3

        !     enddo
        ! enddo

        return

        contains

            subroutine default_inflow()
                implicit none

                select case (inflow_type)
                    case (POISEUILLE_INFLOW)

                      if (nrank==0) write(*,*)'PERFORMING default poiseuille inflow'

                      q1_inflow = 0.d0
                      q2_inflow = 0.d0
                      q3_inflow = 0.d0
                      if (streamwise==1) call perform_stream1(q1_inflow, BC2, BC3)
                      !if (streamwise==3) call perform_stream3(q3_inflow, BC1, BC2)

                      ! q1_inflow = q1_x(1,:,:)

                    case (SQUARE_INFLOW)

                      if (nrank==0) write(*,*)'PERFORMING default square inflow'

                      q1_inflow = 0.d0
                      q2_inflow = 0.d0
                      q3_inflow = 0.d0
                      ! if (streamwise==1) q1_inflow(:,:)=1.d0
                      ! if (streamwise==3) q3_inflow(:,:)=1.d0
                      if (streamwise==1) call perform_boundary_layer_1(q1_inflow, Yc(4)) ! BLASIUS BOUNDARY LAYER ON 4 NODES 
                      if (streamwise==3) call perform_boundary_layer_3(q3_inflow, Yc(4)) ! BLASIUS BOUNDARY LAYER ON 4 NODES 

                      if (SCA_state==1) then

                        do k = xstart(3), min(n3-1, xend(3))
                            do j = xstart(2), min(n2-1, xend(2))
                                sca_inflow(j,k)=(-delta_T)+(Yc(j)/L2)*2.d0*delta_T
                            end do
                        end do

                      endif

                    case (BOUNDARY_LAYER_INFLOW)

                      if (nrank==0) write(*,*)'PERFORMING default boundary layer inflow with a boundary thickness:', delta_BL

                      q1_inflow = 0.d0
                      q2_inflow = 0.d0
                      q3_inflow = 0.d0
                      if (streamwise==1) call perform_boundary_layer_1(q1_inflow, delta_BL)
                      if (streamwise==3) call perform_boundary_layer_3(q3_inflow, delta_BL)

                endselect

            end subroutine default_inflow

            ! TO CHECK
            subroutine inflow_from_file()
                use HDF5_IO
                use start_settings, only:start_it
                use COMMON_workspace_view, only: COMMON_inflow_path

                implicit none
                integer   :: i
                character(200)       :: current_inflow_path

                real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: ts11_global
                real*8, dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))      :: ts11_tmp_x
                real*8, dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))      :: ts11_tmp_y

                character*10 tmp_str

                ! current_inflow_path=trim(COMMON_inflow_path)//'outflow_'//trim(adjustl(tmp_str))
                ! Read the streamwise velocity
                current_inflow_path=trim(COMMON_inflow_path)//'W'

                if (nrank==0) write(*,*)
                if (nrank==0) write(*,*) "Reading inflow from file:", trim(current_inflow_path)//".h5"
                if (nrank==0) write(*,*) "The inflow is get from nx =", nx_start

                call hdf_read_3Dfield(current_inflow_path, ts11_global, "W", nx_global, ny_global, nz_global, xstart(1), xend(1), xstart(2),xend(2), xstart(3),xend(3))
                !call hdf_read_3Dfield(current_inflow_path, ts12_tmp, "q2_out", 1, ny_global, nz_global, 1,1, xstart(2),xend(2), xstart(3),xend(3))
                !call hdf_read_3Dfield(current_inflow_path, ts13_tmp, "q3_out", 1, ny_global, nz_global, 1,1, xstart(2),xend(2), xstart(3),xend(3))

                q1_inflow = 0.d0
                q2_inflow = 0.d0
                q3_inflow = 0.d0

                q1_inflow(:,:)=ts11_global(nx_start,:,:)

            end subroutine inflow_from_file

            ! TO CHECK
            subroutine inflow_from_prev_run()
                use HDF5_IO
                use start_settings, only:start_it
                use COMMON_workspace_view, only: COMMON_inflow_path
                use mesh

                implicit none
            
                real*8 :: u_bulk, u_bulk_glob
                integer :: j,k,mpi_err

                real*8 :: delta, uc_inflow, coef
                logical :: pair_n2
                real*8, dimension(n2) :: inflow_profile

                real*8 :: ut1, u_bulk_theo

                u_bulk=0.d0

                do j=xstart(2),min(n2m,xend(2))
                    do k=xstart(3),min(n3m,xend(3))

                        u_bulk = u_bulk + q1_x(n1m,j,k)*(Y(j+1)-Y(j))*dx3

                    enddo
                enddo

                call MPI_ALLREDUCE(u_bulk, u_bulk_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
                u_bulk_glob = u_bulk_glob/(L3*L2)

                q1_inflow = u_bulk_glob
                q2_inflow = 0.d0
                q3_inflow = 0.d0

                !if (SCA_state==1) then
                !    sca_inflow = sca_x(1,:,:)
                !endif

                ! USING BLASIUS BOUNDARY LAYER
                delta = Y(5) ! profile on 4 nodes

                pair_n2 = (mod(n2, 2)==0)

                do j=ystart(2),n2/2
                    if (Yc(j)<delta) then
                        inflow_profile(j) = 2.d0 * (Yc(j)/delta) - 2.d0 * (Yc(j)/delta)**3 + 1.d0 * (Yc(j)/delta)**4
                        inflow_profile(n2-j) = 2.d0 * (Yc(j)/delta) - 2.d0 * (Yc(j)/delta)**3 + 1.d0 * (Yc(j)/delta)**4
                    else
                        inflow_profile(j) = 1.d0
                        inflow_profile(n2-j) = 1.d0
                    endif
                enddo

                if (pair_n2) then
                    inflow_profile(n2/2+1) = 1.d0
                endif

                ! Perform the flowrate at the inlet...
                ut1=0.d0
                do k=xstart(3),min(xend(3),n3-1)
                    do j=xstart(2),min(xend(2),n2-1)
                        ut1=ut1+inflow_profile(j)*(Y(j+1)-Y(j))*dx3
                    enddo
                enddo

                call MPI_ALLREDUCE(ut1, u_bulk_theo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)
                u_bulk_theo=u_bulk_theo/(L3*L2)

                coef = u_bulk/u_bulk_theo

                inflow_profile(:) = inflow_profile(:) * coef

                do j=xstart(2),min(xend(2),n2m)
                    q1_inflow(j,:) = inflow_profile(j)
                enddo

            end subroutine inflow_from_prev_run

            subroutine perform_stream1(stream1, BC2, BC3)
              implicit none

              real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: stream1
              integer                                                 :: BC2, BC3

              real*8                                                  :: f2(n2), f3(n3), c2, c3
              integer                                                 :: j, k

              f3=1.d0

              if (mod(n2, 2)==0) then
                  c2=Yc(n2/2)
              else
                  c2=Y((n2-1)/2+1)
              end if
              

              if (BC2==NOSLIP) then

                do j=1,n2-1
                    f2(j)=(1.d0-((Yc(j)-c2)/(0.5d0*L2))**2)
                enddo
                f2(n2)=f2(n2-1)

                  if (BC3==NOSLIP) then

                      do k=1,n3

                          if (mod(n3, 2)==0) then
                              c3=(n3/2-0.5d0)*dx3
                          else
                              c3=((n3-1)/2+1-0.5d0)*dx3
                          end if

                          f3(k)=(1.d0-(((k-0.5d0)*dx3-c3)/(0.5d0*L3))**2)

                      enddo

                  end if

                  do k = xstart(3),min(xend(3),n3-1)
                      do j= xstart(2),xend(2)
                          stream1(j,k)=f2(j)*f3(k)
                      end do
                  end do

              else
                  stream1=0.d0
              end if

          end subroutine perform_stream1

        subroutine perform_stream3(stream3, BC1, BC2)
            implicit none

            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2)) :: stream3
            integer                                                 :: BC1, BC2

            real*8                                                  :: f1(n1), f2(n2), c1, c2
            integer                                                 :: i, j,temp_test

            f1=1.d0


            if (BC2==NOSLIP) then

                if (mod(n2, 2)==0) then
                   c2=Yc(n2/2)
                else
                    c2=Y((n2-1)/2+1)
                end if

                do j=1,n2-1

                    f2(j)=(1.d0-((Yc(j)-c2)/(0.5d0*L2))**2)

                enddo

               if (BC1==NOSLIP) then

                    if (mod(n1, 2)==0) then
                       c1=(n1/2-0.5d0)*dx1
                    else
                      c1=((n1-1)/2+1-0.5d0)*dx1
                    end if

                    f1(i)=(1.d0-(((i-0.5d0)*dx1-c1)/(0.5d0*L1))**2)
                    
                 else !

                    do j = zstart(2),zend(2)
                     do i= zstart(1),min(zend(1),n1-1)
                        stream3(i,j)=f2(j)*f1(i)
                     end do
                   end do

               endif !BC1

            else !
                stream3=0.d0
            end if !BC2

        end subroutine perform_stream3

          subroutine perform_boundary_layer_1(stream1, delta_BL)
            implicit none

            real*8, dimension(xstart(2):xend(2), xstart(3):xend(3)) :: stream1
            real*8, dimension(1:n2)                                 :: f1

            real*8                                                  :: delta_BL
            integer                                                 :: j
            logical                                                 :: pair_n2

            pair_n2 = (mod(n2, 2)==0)

            do j=1,n2/2
                if (Yc(j)<delta_BL) then
                    !f1(j) = 1.5d0 * (Yc(j)/delta_BL) - 0.5d0 * (Yc(j)/delta_BL)**2
                    !f1(n2-j) = 1.5d0 * (Yc(j)/delta_BL) - 0.5d0 * (Yc(j)/delta_BL)**2
                    f1(j) = 2.d0 * (Yc(j)/delta_BL) - 2.d0 * (Yc(j)/delta_BL)**3 + 1.d0 * (Yc(j)/delta_BL)**4
                    f1(n2-j) = 2.d0 * (Yc(j)/delta_BL) - 2.d0 * (Yc(j)/delta_BL)**3 + 1.d0 * (Yc(j)/delta_BL)**4
                else
                    f1(j) = 1.d0
                    f1(n2-j) = 1.d0
                endif
            enddo

            if (pair_n2) then
                f1(n2/2+1) = 1.d0
            endif

            do j = xstart(2),xend(2)
                stream1(j,:)=f1(j)
            end do

        end subroutine perform_boundary_layer_1

        subroutine perform_boundary_layer_3(stream3, delta_BL)
            implicit none

            real*8, dimension(zstart(1):zend(1), zstart(2):zend(2)) :: stream3
            real*8, dimension(1:n2)                                 :: f3

            real*8                                                  :: delta_BL
            integer                                                 :: j
            logical                                                 :: pair_n2

            pair_n2 = (mod(n2, 2)==0)

            do j=1,n2/2
                if (Yc(j)<delta_BL) then
                    !f3(j) = (3/2) * (Yc(j)/delta_BL) - (1/2) * (Yc(j)/delta_BL)**2
                    !f3(n2-j) = (3/2) * (Yc(j)/delta_BL) - (1/2) * (Yc(j)/delta_BL)**2
                    f3(j) = 2.d0 * (Yc(j)/delta_BL) - 2.d0 * (Yc(j)/delta_BL)**3 + 1.d0 * (Yc(j)/delta_BL)**4
                    f3(n2-j) = 2.d0 * (Yc(j)/delta_BL) - 2.d0 * (Yc(j)/delta_BL)**3 + 1.d0 * (Yc(j)/delta_BL)**4
                else
                    f3(j) = 1.d0
                    f3(n2-j) = 1.d0
                endif
            enddo

            if (pair_n2) then
                f3(n2/2+1) = 1.d0
            endif

            do j = zstart(2),zend(2)
                stream3(1:min(xend(1),n1-1),j)=f3(j)
            end do

        end subroutine perform_boundary_layer_3


    end subroutine get_inflow

end module FRINGE_initializer
