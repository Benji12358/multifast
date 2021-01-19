module FRINGE_initializer
    use decomp_2d
    use mpi
    use mesh
    use FRINGE_data
    use HDF5_IO
    use DNS_settings
    use boundaries

!
    implicit none
!
contains
!
!
    subroutine get_inflow()

        implicit none

        if (inflow_type==2) then
            ! call inflow_from_file
        else
            call default_inflow
        endif

        return

        contains

            subroutine default_inflow()
                implicit none
                integer j

                select case (inflow_type)
                    case (POISEUILLE_INFLOW)

                      if (nrank==0) write(*,*)'PERFORMING default poiseuille inflow'

                        call perform_stream1(q1_inflow, BC2, BC3)
                        q2_inflow(:,:)=0.d0
                        q3_inflow(:,:)=0.d0

                    case (SQUARE_INFLOW)

                        if (nrank==0) write(*,*)'PERFORMING default square inflow'

                          q1_inflow(:,:)=1.d0
                          q2_inflow(:,:)=0.d0
                          q3_inflow(:,:)=0.d0

                endselect

            end subroutine default_inflow

            ! TO CHECK
            ! subroutine inflow_from_file()
            !     use run_ctxt_data
            !     use flow_buffer_handler
            !     use COMMON_fieldstools

            !     implicit none
            !     real*8  :: cs1, cs2, cs3

            !     call get_2dfield(id_q1, q1_x(1, :,:))
            !     call get_2dfield(id_q2, q2_x(1, :,:))
            !     call get_2dfield(id_q3, q3_x(1, :,:))

            !     call perform_checksum2D_x(q1_x(1, :,:), cs1)
            !     call perform_checksum2D_x(q2_x(1, :,:), cs2)
            !     call perform_checksum2D_x(q3_x(1, :,:), cs3)

            !     if (nrank==0) then
            !         open(15, file="debugoutflow", position="append")
            !         write(15,*) ntime, cs1, cs2, cs3
            !         close(15)
            !     endif


            ! end subroutine inflow_from_file

            subroutine perform_stream1(stream1, BC2, BC3)
              implicit none

              real*8, dimension(ystart(2):yend(2), ystart(3):yend(3)) :: stream1
              integer                                                 :: BC2, BC3

              real*8                                                  :: f2(n2), f3(n3), c2, c3
              integer                                                 :: j, k

              f3=1.d0

              if (BC2==NOSLIP) then


                  if (mod(n2, 2)==0) then
                      c2=Yc(n2/2)
                  else
                      c2=Y((n2-1)/2+1)
                  end if

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

                  do k = ystart(3),yend(3)
                      do j= ystart(2),yend(2)
                          stream1(j,k)=f2(j)*f3(k)
                      end do
                  end do

              else
                  stream1=0.d0
              end if

          end subroutine perform_stream1


    end subroutine get_inflow

end module FRINGE_initializer