module subdomains_view

    use decomp_2d
    use mesh
    use mpi

    implicit none

    ! Local decomp data (direction), specific to one processor
    real*8, dimension(3)            :: xpos_start, xpos_end
    real*8, dimension(3)            :: ypos_start, ypos_end
    real*8, dimension(3)            :: zpos_start, zpos_end


    ! Global decomp data (processor, direction), shared to each processor
    real*8, dimension(:,:), allocatable  :: xpos_start_glob, xpos_end_glob ! allocated in main > allocate_data
    real*8, dimension(:,:), allocatable  :: ypos_start_glob, ypos_end_glob ! allocated in main > allocate_data
    real*8, dimension(:,:), allocatable  :: zpos_start_glob, zpos_end_glob ! allocated in main > allocate_data

    integer, dimension(:,:), allocatable            :: xstart_glob, xend_glob
    integer, dimension(:,:), allocatable            :: ystart_glob, yend_glob
    integer, dimension(:,:), allocatable            :: zstart_glob, zend_glob

      interface get_processor_number
        module procedure get_processor_number, get_processor_number_node
      end interface

contains

    ! INITIALIZATION---------------------------------------------------------
    subroutine init_subdomains_view()
        implicit none

        integer :: i, mpi_err

        xpos_start(1)=Z(xstart(1))
        xpos_start(2)=Y(xstart(2))
        xpos_start(3)=X(xstart(3))

        ypos_start(1)=Z(ystart(1))
        ypos_start(2)=Y(ystart(2))
        ypos_start(3)=X(ystart(3))

        zpos_start(1)=Z(zstart(1))
        zpos_start(2)=Y(zstart(2))
        zpos_start(3)=X(zstart(3))

        xpos_end(1)=Z(min(n1, xend(1)+1))
        xpos_end(2)=Y(min(n2, xend(2)+1))
        xpos_end(3)=X(min(n3, xend(3)+1))

        ypos_end(1)=Z(min(n1, yend(1)+1))
        ypos_end(2)=Y(min(n2, yend(2)+1))
        ypos_end(3)=X(min(n3, yend(3)+1))

        zpos_end(1)=Z(min(n1, zend(1)+1))
        zpos_end(2)=Y(min(n2, zend(2)+1))
        zpos_end(3)=X(min(n3, zend(3)+1))

        do i = 1, 3

            call MPI_ALLGATHER (xstart(i),1, MPI_INTEGER ,xstart_glob(:,i),1, MPI_INTEGER , MPI_COMM_WORLD ,mpi_err)
            call MPI_ALLGATHER (xend(i),1, MPI_INTEGER ,xend_glob(:,i),1, MPI_INTEGER , MPI_COMM_WORLD ,mpi_err)

            call MPI_ALLGATHER (ystart(i),1, MPI_INTEGER ,ystart_glob(:,i),1, MPI_INTEGER , MPI_COMM_WORLD ,mpi_err)
            call MPI_ALLGATHER (yend(i),1, MPI_INTEGER ,yend_glob(:,i),1, MPI_INTEGER , MPI_COMM_WORLD ,mpi_err)

            call MPI_ALLGATHER (zstart(i),1, MPI_INTEGER ,zstart_glob(:,i),1, MPI_INTEGER , MPI_COMM_WORLD ,mpi_err)
            call MPI_ALLGATHER (zend(i),1, MPI_INTEGER ,zend_glob(:,i),1, MPI_INTEGER , MPI_COMM_WORLD ,mpi_err)



            call MPI_ALLGATHER (xpos_start(i),1, MPI_DOUBLE_PRECISION ,xpos_start_glob(:,i),1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD ,mpi_err)
            call MPI_ALLGATHER (xpos_end(i),1, MPI_DOUBLE_PRECISION ,xpos_end_glob(:,i),1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD ,mpi_err)

            call MPI_ALLGATHER (ypos_start(i),1, MPI_DOUBLE_PRECISION ,ypos_start_glob(:,i),1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD ,mpi_err)
            call MPI_ALLGATHER (ypos_end(i),1, MPI_DOUBLE_PRECISION ,ypos_end_glob(:,i),1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD ,mpi_err)

            call MPI_ALLGATHER (zpos_start(i),1, MPI_DOUBLE_PRECISION ,zpos_start_glob(:,i),1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD ,mpi_err)
            call MPI_ALLGATHER (zpos_end(i),1, MPI_DOUBLE_PRECISION ,zpos_end_glob(:,i),1, MPI_DOUBLE_PRECISION , MPI_COMM_WORLD ,mpi_err)


        end do
            !        do i = 1, 3
            !
            !            if (nrank==0) write(*,*) ' Direction :', i
            !
            !            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
            !            write(6,*) 'st :', nrank, xpos_start(i)
            !
            !            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
            !            write(6,*) 'end:', nrank, xpos_end(i)
            !
            !            if (nrank==0) write(*,*) '---global----:'
            !
            !            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
            !            write(6,*) 'st :', nrank, xpos_start_glob(:,i)
            !
            !            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
            !            write(6,*) 'st :', nrank, xpos_end_glob(:,i)
            !
            !        end do



!
!        do i = 1, 3
!
!            if (nrank==0) write(*,*) ' Direction :', i
!
!            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!            write(6,*) 'X :', nrank,i, xstart(i), '->', xend(i)
!
!            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!            write(6,*) 'Y :', nrank,i, ystart(i), '->', yend(i)
!
!            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!            write(6,*) 'Z :', nrank,i, zstart(i), '->', zend(i)
!
!
!
!            if (nrank==0) write(*,*) '---global----:'
!
!            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!            if (nrank==0) write(6,*) 'X :', nrank,i, xstart_glob(:,i), '->', xend_glob(:,i)
!
!            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!            if (nrank==0) write(6,*) 'Y :', nrank,i, ystart_glob(:,i), '->', yend_glob(:,i)
!
!            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!            if (nrank==0) write(6,*) 'Z :', nrank,i, zstart_glob(:,i), '->', zend_glob(:,i)
!
!        end do
    !
    !
    !    a=0.2d0
    !    b=0.d0
    !    c=0.1d0
    !    do p = 0, nproc-1
    !
    !        test(1)=xpos_start_glob(p,1)+a*(xpos_end_glob(p,1)-xpos_start_glob(p,1))
    !        test(2)=xpos_start_glob(p,2)+b*(xpos_end_glob(p,2)-xpos_start_glob(p,2))
    !        test(3)=xpos_start_glob(p,3)+c*(xpos_end_glob(p,3)-xpos_start_glob(p,3))
    !
    !        !if ((nrank==0)) write(*,*)'TEST:', test
    !
    !        call get_processor_number(test, 1, j)
    !
    !        if (nrank==0) write(*,*) 'proc_ exp:', p
    !        if (nrank==0) write(*,*) 'proc     :', j
    !
    !
    !        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    !
    !    end do
    !
    !    call get_processor_number(test, 1, j)
    !    write(6,*) xpos_start_glob(j,1), xpos_end_glob(j,1), xpos_start_glob(j,2), xpos_end_glob(j,2), xpos_start_glob(j,3), xpos_end_glob(j,3)

    end subroutine init_subdomains_view

    ! GETTING A NUMBER OF A PROCESSOR ASSIGNED TO A BUBBLE POSITION-------------------------
    subroutine get_processor_number(pt, decomp, proc_num)

        implicit none

        real*8, dimension(3), intent (in)           :: pt
        integer, intent(in)                         :: decomp
        integer, intent(out)                        :: proc_num

        real*8, dimension(0:nproc-1,3)              :: start
        real*8, dimension(0:nproc-1,3)              :: end
        integer                                     :: i, j, k

        if (decomp==1) then
            start(:,:) = xpos_start_glob(:,:)
            end(:,:)=  xpos_end_glob(:,:)
        end if

        if (decomp==2) then
            start(:,:) = ypos_start_glob(:,:)
            end(:,:) =  ypos_end_glob(:,:)
        end if

        if (decomp==3) then
            start(:,:) = zpos_start_glob(:,:)
            end(:,:) =  zpos_end_glob(:,:)
        end if

        proc_num=-1

        do i = 0, nproc-1

            if (((pt(1)>=start(i,1)).and.(pt(1)<=end(i,1))).and.((pt(2)>=start(i,2)).and.(pt(2)<=end(i,2))).and.((pt(3)>=start(i,3)).and.(pt(3)<=end(i,3)))) then
                proc_num=i
                return
            end if

        end do

    end subroutine get_processor_number

    subroutine get_processor_number_node(i,j,k, decomp, proc_num)

        implicit none

        integer, intent(in)                         :: i,j,k
        integer, intent(in)                         :: decomp
        integer, intent(out)                        :: proc_num

        integer, dimension(0:nproc-1,3)             :: start
        integer, dimension(0:nproc-1,3)             :: end
        integer                                     :: p

        if (decomp==1) then
            start(:,:) = xstart_glob(:,:)
            end(:,:)=  xend_glob(:,:)
        end if

        if (decomp==2) then
            start(:,:) = ystart_glob(:,:)
            end(:,:) =  yend_glob(:,:)
        end if

        if (decomp==3) then
            start(:,:) = zstart_glob(:,:)
            end(:,:) =  zend_glob(:,:)
        end if

        proc_num=-1

        do p = 0, nproc-1

            if (((i>=start(p,1)).and.(i<=end(p,1))).and.((j>=start(p,2)).and.(j<=end(p,2))).and.((k>=start(p,3)).and.(k<=end(p,3)))) then
                proc_num=p
                return
            end if

        end do

    end subroutine get_processor_number_node


    subroutine test()
        use mpi
        implicit none
        integer     :: mpi_err, dest, valeur1=0, valeur2=0

        if (nrank==2) dest=5

        call MPI_BCAST (dest,1, MPI_INTEGER ,2, MPI_COMM_WORLD ,mpi_err)

        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

        if (nrank==2) call MPI_SEND (10,1, MPI_INTEGER ,5,100, MPI_COMM_WORLD ,mpi_err)
        if (nrank==3) call MPI_SEND (20,1, MPI_INTEGER ,5,100, MPI_COMM_WORLD ,mpi_err)
        if (nrank==dest) call MPI_RECV (valeur1,1, MPI_INTEGER ,3,100, MPI_COMM_WORLD ,MPI_STATUS_IGNORE,mpi_err)
        if (nrank==dest) call MPI_RECV (valeur2,1, MPI_INTEGER ,2,100, MPI_COMM_WORLD ,MPI_STATUS_IGNORE,mpi_err)

    end subroutine test

end module subdomains_view
