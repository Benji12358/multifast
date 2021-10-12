module decomp2D_utils
    implicit none

    integer, parameter  :: x_pencil=1, z_pencil=3

contains


    subroutine assemble_global_x_all(local,global,nx,ny,nz)

        use decomp_2d
        use MPI

        implicit none
        integer, intent(IN) :: nx,ny,nz
        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(IN) :: local
        real(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        integer :: n

        do n = 0, nproc-1
            call assemble_global_x(local, global, nx,ny,nz, n)
        end do

    end subroutine assemble_global_x_all



    subroutine assemble_global_y_all(local,global,nx,ny,nz)

        use decomp_2d
        use MPI

        implicit none
        integer, intent(IN) :: nx,ny,nz
        real(mytype), dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(IN) :: local
        real(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        integer :: n

        do n = 0, nproc-1
            call assemble_global_y(local, global, nx,ny,nz, n)
        end do

    end subroutine assemble_global_y_all



    subroutine assemble_global_z_all(local,global,nx,ny,nz)

        use decomp_2d
        use MPI

        implicit none
        integer, intent(IN) :: nx,ny,nz
        real(mytype), dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(IN) :: local
        real(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        integer :: n

        do n = 0, nproc-1
            call assemble_global_z(local, global, nx,ny,nz, n)
        end do

    end subroutine assemble_global_z_all



    subroutine assemble_global_x(local,global,nx,ny,nz, cible)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,nz,cible
        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)), intent(IN) :: local
        real(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        real(mytype), allocatable, dimension(:,:,:) :: rbuf
        integer, dimension(9) :: sbuf1, rbuf1

        integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
        integer, dimension(MPI_STATUS_SIZE) :: status

        if (nrank==cible) then
            ! master writes its own data to a global array

            i1 = xstart(1)
            i2 = xend(1)
            j1 = xstart(2)
            j2 = xend(2)
            k1 = xstart(3)
            k2 = xend(3)

            !write(*,*)'local shape', i1+1, j1+1, k1+1
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)

            do k=k1,k2
                do j=j1,j2
                    do i=i1,i2
                        ! 'local' is assumbed shape array
                        ! but it is OK as starting index for rank 0 always 1
                        global(i,j,k)=local(i,j,k)
                    end do
                end do
            end do

            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
            ! then loop through all other ranks to collect data
            do m=0,nproc-1

                if (m/=cible) then

                    CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
                    status,ierror)
                    allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
                    rbuf1(7):rbuf1(8)))
                    CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),real_type,m, &
                    m+nproc,MPI_COMM_WORLD,status,ierror)
                    do k=rbuf1(7),rbuf1(8)
                        do j=rbuf1(4),rbuf1(5)
                            do i=rbuf1(1),rbuf1(2)
                                global(i,j,k)=rbuf(i,j,k)
                            end do
                        end do
                    end do
                    deallocate(rbuf)

                end if


            end do

            !call sleep(2)
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)
            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
        else
            ! slaves send data to mater
            sbuf1(1) = xstart(1)
            sbuf1(2) = xend(1)
            sbuf1(3) = xsize(1)
            sbuf1(4) = xstart(2)
            sbuf1(5) = xend(2)
            sbuf1(6) = xsize(2)
            sbuf1(7) = xstart(3)
            sbuf1(8) = xend(3)
            sbuf1(9) = xsize(3)
            count = xsize(1)*xsize(2)*xsize(3)
            ! send partition information
            CALL MPI_SEND(sbuf1,9,MPI_INTEGER,cible,nrank,MPI_COMM_WORLD,ierror)
            ! send data array
            CALL MPI_SEND(local,count,real_type,cible, &
            nrank+nproc,MPI_COMM_WORLD,ierror)
        end if

        return
    end subroutine assemble_global_x


    subroutine assemble_global_y(local,global,nx,ny,nz, cible)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,nz,cible
        real(mytype), dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)), intent(IN) :: local
        real(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)):: local_x

        call transpose_y_to_x(local, local_x)

        call assemble_global_x(local_x, global, nx,ny,nz, cible)

        return
    end subroutine assemble_global_y


    subroutine assemble_global_z(local,global,nx,ny,nz, cible)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,nz,cible
        real(mytype), dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)), intent(IN) :: local
        real(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        real(mytype), allocatable, dimension(:,:,:) :: rbuf
        integer, dimension(9) :: sbuf1, rbuf1

        integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
        integer, dimension(MPI_STATUS_SIZE) :: status

        if (nrank==cible) then
            ! master writes its own data to a global array
            i1 = zstart(1)
            i2 = zend(1)
            j1 = zstart(2)
            j2 = zend(2)
            k1 = zstart(3)
            k2 = zend(3)

            !write(*,*)'local shape', i1+1, j1+1, k1+1
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)

            do k=k1,k2
                do j=j1,j2
                    do i=i1,i2
                        ! 'local' is assumbed shape array
                        ! but it is OK as starting index for rank 0 always 1
                        global(i,j,k)=local(i,j,k)
                    end do
                end do
            end do

            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
            ! then loop through all other ranks to collect data
            do m=0,nproc-1

                if (m/=cible) then

                    CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
                    status,ierror)
                    allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
                    rbuf1(7):rbuf1(8)))
                    CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),real_type,m, &
                    m+nproc,MPI_COMM_WORLD,status,ierror)
                    do k=rbuf1(7),rbuf1(8)
                        do j=rbuf1(4),rbuf1(5)
                            do i=rbuf1(1),rbuf1(2)
                                global(i,j,k)=rbuf(i,j,k)
                            end do
                        end do
                    end do
                    deallocate(rbuf)

                end if


            end do

            !call sleep(2)
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)
            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
        else
            ! slaves send data to mater
            sbuf1(1) = zstart(1)
            sbuf1(2) = zend(1)
            sbuf1(3) = zsize(1)
            sbuf1(4) = zstart(2)
            sbuf1(5) = zend(2)
            sbuf1(6) = zsize(2)
            sbuf1(7) = zstart(3)
            sbuf1(8) = zend(3)
            sbuf1(9) = zsize(3)
            count = zsize(1)*zsize(2)*zsize(3)
            ! send partition information
            CALL MPI_SEND(sbuf1,9,MPI_INTEGER,cible,nrank,MPI_COMM_WORLD,ierror)
            ! send data array
            CALL MPI_SEND(local,count,real_type,cible, &
            nrank+nproc,MPI_COMM_WORLD,ierror)
        end if

        return
    end subroutine assemble_global_z




    subroutine assemble_plane(local,nx,ny, ist, ien, jst, jen, cible, global)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,cible, ist, ien, jst, jen
        real(mytype), dimension(ist:ien, jst:jen), intent(IN) :: local
        real(mytype), dimension(nx,ny), intent(OUT) :: global

        real(mytype), allocatable, dimension(:,:) :: rbuf
        integer, dimension(6) :: sbuf1, rbuf1

        integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
        integer, dimension(MPI_STATUS_SIZE) :: status

        if (nrank==cible) then
            ! master writes its own data to a global array
            i1 = ist
            i2 = ien
            j1 = jst
            j2 = jen

            !write(*,*)'local shape', i1+1, j1+1, k1+1
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)

            do j=j1,j2
                do i=i1,i2
                    ! 'local' is assumbed shape array
                    ! but it is OK as starting index for rank 0 always 1
                    global(i,j)=local(i,j)
                end do
            end do

            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
            ! then loop through all other ranks to collect data
            do m=0,nproc-1

                if (m/=cible) then

                    CALL MPI_RECV(rbuf1,6,MPI_INTEGER,m,m,MPI_COMM_WORLD,status,ierror)
                    allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5)))
                    CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6),real_type,m, &
                    m+nproc,MPI_COMM_WORLD,status,ierror)
                    do j=rbuf1(4),rbuf1(5)
                        do i=rbuf1(1),rbuf1(2)
                            global(i,j)=rbuf(i,j)
                        end do
                    end do
                    deallocate(rbuf)

                end if


            end do

            !call sleep(2)
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)
            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
        else
            ! slaves send data to mater
            sbuf1(1) = ist
            sbuf1(2) = ien
            sbuf1(3) = ien-ist+1
            sbuf1(4) = jst
            sbuf1(5) = jen
            sbuf1(6) = jen-jst+1
            count = sbuf1(3)*sbuf1(6)
            ! send partition information
            CALL MPI_SEND(sbuf1,6,MPI_INTEGER,cible,nrank,MPI_COMM_WORLD,ierror)
            ! send data array
            CALL MPI_SEND(local,count,real_type,cible, &
            nrank+nproc,MPI_COMM_WORLD,ierror)
        end if

        return
    end subroutine assemble_plane

    subroutine transpose2Dy_y_to_x(fieldy, fieldx, xstart2)

        use decomp_2d
        use mpi

        implicit none

        integer, dimension(nproc), intent(IN)                                               :: xstart2
        real(mytype), dimension(ystart(1):yend(1), ystart(3):yend(3)), intent(IN)           :: fieldy
        real(mytype), dimension(xstart(1):xend(1), xstart(3):xend(3)), intent(OUT)          :: fieldx

        real(mytype), dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))    :: tmpy
        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))    :: tmpx

        real*8  :: cs1_glob, cs2_glob, cs3_glob
        real*8  :: cs1, cs2, cs3
        integer :: mpi_err, i

        fieldx=0.d0
        tmpy=0.d0
        tmpx=0.d0


!        cs1=sum(abs(fieldy))
!        call MPI_ALLREDUCE (cs1, cs1_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        do i = 1, nproc
            tmpy(:,xstart2(i),:)=fieldy
        end do

        call transpose_y_to_x(tmpy, tmpx)

!        if (xstart(2)==1) cs2=sum(abs(tmpx))
!        call MPI_ALLREDUCE (cs2, cs2_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        fieldx(:,:)=tmpx(:, xstart(2), :)

!        cs3=0.d0
!        if (xstart(2)==1) cs3=sum(abs(fieldx))

!        call MPI_ALLREDUCE (cs3, cs3_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!        write(*,*)'transpose2Dy_y_to_x', cs1_glob, cs2_glob, cs3_glob


    end subroutine transpose2Dy_y_to_x

    subroutine transpose2Dy_z_to_x(fieldz, fieldx, xstart2)

        use decomp_2d
        use mpi

        implicit none

        integer, dimension(nproc), intent(IN)                                               :: xstart2
        real(mytype), dimension(zstart(1):zend(1), zstart(3):zend(3)), intent(IN)           :: fieldz
        real(mytype), dimension(xstart(1):xend(1), xstart(3):xend(3)), intent(OUT)          :: fieldx

        real(mytype), dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))    :: tmpz
        real(mytype), dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))    :: tmpy
        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))    :: tmpx

        real*8  :: cs1_glob, cs2_glob, cs3_glob
        real*8  :: cs1, cs2, cs3
        integer :: mpi_err, i

        fieldx=0.d0
        tmpz=0.d0
        tmpy=0.d0
        tmpx=0.d0

!        cs1=0.d0
!        if (zstart(2)==1) cs1=sum(abs(fieldz))
!        call MPI_ALLREDUCE (cs1, cs1_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        if (zstart(2)==1) tmpz(:,1,:)=fieldz

        call transpose_z_to_y(tmpz, tmpy)

        do i = 1, nproc
            tmpy(:,xstart2(i),:)=tmpy(:,1,:)
        end do

        call transpose_y_to_x(tmpy, tmpx)

!        if (xstart(2)==1) cs2=sum(abs(tmpx))
!        call MPI_ALLREDUCE (cs2, cs2_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        fieldx(:,:)=tmpx(:, xstart(2), :)

!        cs3=0.d0
!        if (xstart(2)==1) cs3=sum(abs(fieldx))

!        call MPI_ALLREDUCE (cs3, cs3_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!        write(*,*)'transpose2Dy_z_to_x', cs1_glob, cs2_glob, cs3_glob


    end subroutine transpose2Dy_z_to_x

    subroutine transpose2Dy_y_to_z(fieldy, fieldz, zstart2)

        use decomp_2d
        use mpi

        implicit none

        integer, dimension(nproc), intent(IN)                                               :: zstart2
        real(mytype), dimension(ystart(1):yend(1), ystart(3):yend(3)), intent(IN)           :: fieldy
        real(mytype), dimension(zstart(1):zend(1), zstart(3):zend(3)), intent(OUT)          :: fieldz

        real(mytype), dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))    :: tmpy
        real(mytype), dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))    :: tmpz

        real*8  :: cs1_glob, cs2_glob, cs3_glob
        real*8  :: cs1, cs2, cs3
        integer :: mpi_err, i

        fieldz=0.d0
        tmpy=0.d0
        tmpz=0.d0


!        cs1=sum(abs(fieldy))
!        call MPI_ALLREDUCE (cs1, cs1_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        do i = 1, nproc
            tmpy(:,zstart2(i),:)=fieldy
        end do

        call transpose_y_to_z(tmpy, tmpz)

!        if (zstart(2)==1) cs2=sum(abs(tmpz))
!        call MPI_ALLREDUCE (cs2, cs2_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        fieldz(:,:)=tmpz(:, zstart(2), :)

!        cs3=0.d0
!        if (zstart(2)==1) cs3=sum(abs(fieldz))

!        call MPI_ALLREDUCE (cs3, cs3_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!        write(*,*)'transpose2Dy_y_to_z', cs1_glob, cs2_glob, cs3_glob


    end subroutine transpose2Dy_y_to_z


    subroutine transpose2Dz_z_to_x(fieldz, fieldx, xstart3)

        use decomp_2d
        use mpi

        implicit none

        integer, dimension(nproc), intent(IN)                                               :: xstart3
        real(mytype), dimension(zstart(1):zend(1), zstart(2):zend(2)), intent(IN)           :: fieldz
        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2)), intent(OUT)          :: fieldx

        real(mytype), dimension(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3))    :: tmpz
        real(mytype), dimension(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3))    :: tmpy
        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))    :: tmpx

        real*8  :: cs1_glob, cs2_glob, cs3_glob
        real*8  :: cs1, cs2, cs3
        integer :: mpi_err, i


!        cs1=sum(abs(fieldz))
!        call MPI_ALLREDUCE (cs1, cs1_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        fieldx=0.d0
        tmpz=0.d0
        tmpy=0.d0
        tmpx=0.d0

        do i = 1, nproc
            tmpz(:,:,xstart3(i))=fieldz
        end do

        call transpose_z_to_y(tmpz, tmpy)
        call transpose_y_to_x(tmpy, tmpx)

!        cs2=0.d0
!        if (xstart(3)==1) cs2=sum(abs(tmpx))
!        call MPI_ALLREDUCE (cs2, cs2_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)

        fieldx(:,:)=tmpx(:, :, xstart(3))

!        cs3=0.d0
!        if (xstart(3)==1) cs3=sum(abs(fieldx))

!        call MPI_ALLREDUCE (cs3, cs3_glob, 1, MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , mpi_err)
!        write(*,*)'transpose2Dz_z_to_x', cs1_glob, cs2_glob, cs3_glob

    end subroutine transpose2Dz_z_to_x

    subroutine halo_2Dx_2(field, field_h)

        use decomp_2d

        implicit none

        real(mytype), dimension(xstart(1):xend(1), xstart(3):xend(3)), intent(IN)           :: field
        real(mytype), dimension(:,:), allocatable, intent(OUT)                              :: field_h

        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))    :: tmpx
        real(mytype), dimension(:,:,:), allocatable                                         :: tmpx_h

        tmpx(:,xstart(2),:)=field

        call update_halo(tmpx,tmpx_h,level=1)

        if (.not. allocated(field_h)) allocate(field_h(xsize(1), 0:xsize(3)+1))
!        if(nrank==0) write(*,*)'halo_2Dx_2 XS', LBOUND(field_h)
!        if(nrank==0) write(*,*)'halo_2Dx_2 XE', UBOUND(field_h)
!        if(nrank==0) write(*,*)'halo_2Dx_2 SH', SHAPE(field_h)
        field_h=tmpx_h(:,1,:)

    end subroutine halo_2Dx_2

    subroutine halo_2Dx_3(field, field_h)

        use decomp_2d

        implicit none

        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2)), intent(IN)           :: field
        real(mytype), dimension(:,:), allocatable, intent(OUT)                              :: field_h

        real(mytype), dimension(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3))    :: tmpx
        real(mytype), dimension(:,:,:), allocatable                                         :: tmpx_h

        tmpx(:,:,xstart(3))=field

        call update_halo(tmpx,tmpx_h,level=1)

        if (.not. allocated(field_h)) allocate(field_h(xsize(1), 0:xsize(2)+1))
!        if(nrank==0) write(*,*)'halo_2Dx_3 XS', LBOUND(field_h)
!        if(nrank==0) write(*,*)'halo_2Dx_3 XE', UBOUND(field_h)
!        if(nrank==0) write(*,*)'halo_2Dx_3 SH', SHAPE(field_h)
        field_h=tmpx_h(:,:,1)

    end subroutine halo_2Dx_3

    subroutine assemble_plane2(local,nx,ny, ist, ien, jst, jen, cible, global, comm)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,cible, ist, ien, jst, jen, comm
        real(mytype), dimension(ist:ien, jst:jen), intent(IN) :: local
        real(mytype), dimension(:,:), intent(OUT) :: global

        real(mytype), allocatable, dimension(:,:) :: rbuf
        integer, dimension(6) :: sbuf1, rbuf1

        integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count, nb_procs, rank, mpi_err
        integer, dimension(MPI_STATUS_SIZE) :: status

        call MPI_COMM_SIZE (comm ,nb_procs, mpi_err)
        call MPI_COMM_RANK (comm ,rank, mpi_err)

        if (rank==cible) then
            ! master writes its own data to a global array
            i1 = ist
            i2 = ien
            j1 = jst
            j2 = jen

            !write(*,*)'local shape', i1+1, j1+1, k1+1
            !write(*,*)'rank', rank, 'local ', local(i1+1, j1+1, k1+1)

            do j=j1,j2
                do i=i1,i2
                    ! 'local' is assumbed shape array
                    ! but it is OK as starting index for rank 0 always 1
                    global(i,j)=local(i,j)
                end do
            end do

            !write(*,*)'rank', rank, 'global', global(i1+1, j1+1, k1+1)
            ! then loop through all other ranks to collect data
            do m=0,nb_procs-1

                if (m/=cible) then

                    CALL MPI_RECV(rbuf1,6,MPI_INTEGER,m,m,comm,status,ierror)
                    allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5)))
                    CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6),real_type,m, &
                    m+nb_procs,comm,status,ierror)
                    do j=rbuf1(4),rbuf1(5)
                        do i=rbuf1(1),rbuf1(2)
                            global(i,j)=rbuf(i,j)
                        end do
                    end do
                    deallocate(rbuf)

                end if


            end do

            !call sleep(2)
            !write(*,*)'rank', rank, 'local ', local(i1+1, j1+1, k1+1)
            !write(*,*)'rank', rank, 'global', global(i1+1, j1+1, k1+1)
        else
            ! slaves send data to mater
            sbuf1(1) = ist
            sbuf1(2) = ien
            sbuf1(3) = ien-ist+1
            sbuf1(4) = jst
            sbuf1(5) = jen
            sbuf1(6) = jen-jst+1
            count = sbuf1(3)*sbuf1(6)
            ! send partition information
            CALL MPI_SEND(sbuf1,6,MPI_INTEGER,cible,rank,comm,ierror)
            ! send data array
            CALL MPI_SEND(local,count,real_type,cible, &
            rank+nb_procs,comm,ierror)
        end if

        return
    end subroutine assemble_plane2

    subroutine assemble_1D(local, n1s, n1e, global, n1, comm, cible)
        use mpi
        use decomp_2d
        implicit none
        integer                                 :: comm, cible
        integer                                 :: n1s, n1e, n1
        real*8, dimension(n1s:n1e)              :: local
        real*8, dimension(1:n1)                 :: global

        integer :: i,k
        integer :: mpi_err, nb_procs, rank, source, source_zst1, source_zsen1, source_zsz1
        real*8, dimension(:), allocatable   :: buf


        call MPI_COMM_SIZE (comm ,nb_procs, mpi_err)
        call MPI_COMM_RANK (comm ,rank, mpi_err)

        if (rank==cible) then

            global(n1s:n1e)=local(n1s:n1e)

            do source = 0, nb_procs-1

                if (source/=cible) then

                    call MPI_RECV(source_zst1, 1, MPI_INTEGER, source, source, comm, MPI_STATUS_IGNORE, mpi_err)
                    call MPI_RECV(source_zsen1, 1, MPI_INTEGER, source, source, comm, MPI_STATUS_IGNORE, mpi_err)

                    allocate(buf(source_zst1:source_zsen1))

                    source_zsz1=source_zsen1-source_zst1+1
                    call MPI_RECV(buf, source_zsz1, MPI_REAL8, source, source, comm, MPI_STATUS_IGNORE, mpi_err)
                    global(source_zst1:source_zsen1)=buf(source_zst1:source_zsen1)

                    deallocate(buf)

                end if

            end do


        else
            call MPI_SEND(n1s, 1, MPI_INTEGER, cible, rank, comm, mpi_err)
            call MPI_SEND(n1e, 1, MPI_INTEGER, cible, rank, comm, mpi_err)
            call MPI_SEND(local, n1e-n1s+1, MPI_REAL8, cible, rank, comm, mpi_err)
        endif

            !write(*,*)"deplacements", deplacements
           !!write(*,*)"nb_elements_recus", nb_elements_recus
            !call MPI_GATHERV (mean_loc, zsize(1), MPI_REAL8, mean, nb_elements_recus,&
            !deplacements, MPI_REAL8, 0, comm, mpi_err)



    end subroutine assemble_1D




    subroutine get_nrank(i, j, k, owner, stencil)

        use decomp_2d
        use MPI

        implicit none

        integer:: i, j, k, m, stencil, ierror, owner
        integer, dimension(MPI_STATUS_SIZE) :: status

        if (stencil==1) then
            if ((xstart(1)<=i).and.(xend(1)>=i).and.(xstart(2)<=j).and.(xend(2)>=j).and.(xstart(3)<=k).and.(xend(3)>=k) ) then

                do m = 0, nproc-1
                    if (m/=nrank) then
                        call MPI_SEND(nrank,1,MPI_INTEGER,m,m,MPI_COMM_WORLD,ierror)
                    end if
                end do

                owner=nrank

            else
                CALL MPI_RECV(owner,1,MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD,status,ierror)
            end if
        end if

        if (stencil==3) then
            if ((zstart(1)<=i).and.(zend(1)>=i).and.(zstart(2)<=j).and.(zend(2)>=j).and.(zstart(3)<=k).and.(zend(3)>=k) ) then

                do m = 0, nproc-1
                    if (m/=nrank) then
                        call MPI_SEND(nrank,1,MPI_INTEGER,m,m,MPI_COMM_WORLD,ierror)
                    end if
                end do

                owner=nrank

            else
                CALL MPI_RECV(owner,1,MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD,status,ierror)
            end if
        end if

        return
    end subroutine get_nrank


end module decomp2D_utils

module decomp2D_utils_complex
    implicit none

    integer, parameter  :: x_pencil=1, z_pencil=3

contains


    subroutine c_assemble_global_x_all(local,global,nx,ny,nz, sp_decomp)

        use decomp_2d
        use MPI

        implicit none
        type(decomp_info) :: sp_decomp
        integer, intent(IN) :: nx,ny,nz

        complex(mytype), dimension( &
        sp_decomp%xst(1):sp_decomp%xen(1),  &
        sp_decomp%xst(2):sp_decomp%xen(2),  &
        sp_decomp%xst(3):sp_decomp%xen(3)), intent(IN) :: local
        complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        integer :: n

        do n = 0, nproc-1
            call assemble_global_x(local, global, nx,ny,nz, n, sp_decomp)
        end do

    end subroutine c_assemble_global_x_all



    subroutine c_assemble_global_y_all(local,global,nx,ny,nz, sp_decomp)

        use decomp_2d
        use MPI

        implicit none
        type(decomp_info) :: sp_decomp
        integer, intent(IN) :: nx,ny,nz

        complex(mytype), dimension( &
        sp_decomp%yst(1):sp_decomp%yen(1),  &
        sp_decomp%yst(2):sp_decomp%yen(2),  &
        sp_decomp%yst(3):sp_decomp%yen(3)), intent(IN) :: local

        complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        integer :: n

        do n = 0, nproc-1
            call assemble_global_y(local, global, nx,ny,nz, n, sp_decomp)
        end do

    end subroutine c_assemble_global_y_all



    subroutine c_assemble_global_z_all(local,global,nx,ny,nz, sp_decomp)

        use decomp_2d
        use MPI

        implicit none
        integer, intent(IN) :: nx,ny,nz
        type(decomp_info) :: sp_decomp

        complex(mytype), dimension( &
        sp_decomp%zst(1):sp_decomp%zen(1),  &
        sp_decomp%zst(2):sp_decomp%zen(2),  &
        sp_decomp%zst(3):sp_decomp%zen(3)), intent(IN) :: local
        complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        integer :: n

        do n = 0, nproc-1
            call assemble_global_z(local, global, nx,ny,nz, n, sp_decomp)
        end do

    end subroutine c_assemble_global_z_all



    subroutine assemble_global_x(local,global,nx,ny,nz, cible, sp_decomp)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,nz,cible
        type(decomp_info) :: sp_decomp
        complex(mytype), dimension(sp_decomp%xst(1):sp_decomp%xen(1), &
        sp_decomp%xst(2):sp_decomp%xen(2), sp_decomp%xst(3):sp_decomp%xen(3)), intent(IN) :: local

        complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        complex(mytype), allocatable, dimension(:,:,:) :: rbuf
        integer, dimension(9) :: sbuf1, rbuf1

        integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
        integer, dimension(MPI_STATUS_SIZE) :: status

        if (nrank==cible) then
            ! master writes its own data to a global array

            i1 = sp_decomp%xst(1)
            i2 = sp_decomp%xen(1)
            j1 = sp_decomp%xst(2)
            j2 = sp_decomp%xen(2)
            k1 = sp_decomp%xst(3)
            k2 = sp_decomp%xen(3)

            !write(*,*)'local shape', i1+1, j1+1, k1+1
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)

            do k=k1,k2
                do j=j1,j2
                    do i=i1,i2
                        ! 'local' is assumbed shape array
                        ! but it is OK as starting index for rank 0 always 1
                        global(i,j,k)=local(i,j,k)
                    end do
                end do
            end do

            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
            ! then loop through all other ranks to collect data
            do m=0,nproc-1

                if (m/=cible) then

                    CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
                    status,ierror)
                    allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
                    rbuf1(7):rbuf1(8)))
                    CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),complex_type,m, &
                    m+nproc,MPI_COMM_WORLD,status,ierror)
                    do k=rbuf1(7),rbuf1(8)
                        do j=rbuf1(4),rbuf1(5)
                            do i=rbuf1(1),rbuf1(2)
                                global(i,j,k)=rbuf(i,j,k)
                            end do
                        end do
                    end do
                    deallocate(rbuf)

                end if


            end do

            !call sleep(2)
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)
            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
        else
            ! slaves send data to mater
            sbuf1(1) = sp_decomp%xst(1)
            sbuf1(2) = sp_decomp%xen(1)
            sbuf1(3) = sp_decomp%xsz(1)
            sbuf1(4) = sp_decomp%xst(2)
            sbuf1(5) = sp_decomp%xen(2)
            sbuf1(6) = sp_decomp%xsz(2)
            sbuf1(7) = sp_decomp%xst(3)
            sbuf1(8) = sp_decomp%xen(3)
            sbuf1(9) = sp_decomp%xsz(3)
            count = sp_decomp%xsz(1)*sp_decomp%xsz(2)*sp_decomp%xsz(3)
            ! send partition information
            CALL MPI_SEND(sbuf1,9,MPI_INTEGER,cible,nrank,MPI_COMM_WORLD,ierror)
            ! send data array
            CALL MPI_SEND(local,count,complex_type,cible, &
            nrank+nproc,MPI_COMM_WORLD,ierror)
        end if

        return
    end subroutine assemble_global_x


    subroutine assemble_global_y(local,global,nx,ny,nz, cible, sp_decomp)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,nz,cible
        type(decomp_info) :: sp_decomp

        complex(mytype), dimension( &
        sp_decomp%yst(1):sp_decomp%yen(1),  &
        sp_decomp%yst(2):sp_decomp%yen(2),  &
        sp_decomp%yst(3):sp_decomp%yen(3)), intent(IN) :: local

        complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        complex(mytype), dimension( &
        sp_decomp%xst(1):sp_decomp%xen(1),  &
        sp_decomp%xst(2):sp_decomp%xen(2),  &
        sp_decomp%xst(3):sp_decomp%xen(3)):: local_x

        call transpose_y_to_x(local, local_x, sp_decomp)

        call assemble_global_x(local_x, global, nx,ny,nz, cible, sp_decomp)

        return
    end subroutine assemble_global_y


    subroutine assemble_global_z(local,global,nx,ny,nz, cible, sp_decomp)

        use decomp_2d
        use MPI

        implicit none

        integer, intent(IN) :: nx,ny,nz,cible
        type(decomp_info) :: sp_decomp
        complex(mytype), dimension(sp_decomp%zst(1):sp_decomp%zen(1), sp_decomp%zst(2):sp_decomp%zen(2), sp_decomp%zst(3):sp_decomp%zen(3)), intent(IN) :: local
        complex(mytype), dimension(nx,ny,nz), intent(OUT) :: global

        complex(mytype), allocatable, dimension(:,:,:) :: rbuf
        integer, dimension(9) :: sbuf1, rbuf1

        integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
        integer, dimension(MPI_STATUS_SIZE) :: status

        if (nrank==cible) then
            ! master writes its own data to a global array
            i1 = sp_decomp%zst(1)
            i2 = sp_decomp%zen(1)
            j1 = sp_decomp%zst(2)
            j2 = sp_decomp%zen(2)
            k1 = sp_decomp%zst(3)
            k2 = sp_decomp%zen(3)

            !write(*,*)'local shape', i1+1, j1+1, k1+1
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)

            do k=k1,k2
                do j=j1,j2
                    do i=i1,i2
                        ! 'local' is assumbed shape array
                        ! but it is OK as starting index for rank 0 always 1
                        global(i,j,k)=local(i,j,k)
                    end do
                end do
            end do

            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
            ! then loop through all other ranks to collect data
            do m=0,nproc-1

                if (m/=cible) then

                    CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
                    status,ierror)
                    allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
                    rbuf1(7):rbuf1(8)))
                    CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),complex_type,m, &
                    m+nproc,MPI_COMM_WORLD,status,ierror)
                    do k=rbuf1(7),rbuf1(8)
                        do j=rbuf1(4),rbuf1(5)
                            do i=rbuf1(1),rbuf1(2)
                                global(i,j,k)=rbuf(i,j,k)
                            end do
                        end do
                    end do
                    deallocate(rbuf)

                end if


            end do

            !call sleep(2)
            !write(*,*)'nrank', nrank, 'local ', local(i1+1, j1+1, k1+1)
            !write(*,*)'nrank', nrank, 'global', global(i1+1, j1+1, k1+1)
        else
            ! slaves send data to mater
            sbuf1(1) = sp_decomp%zst(1)
            sbuf1(2) = sp_decomp%zen(1)
            sbuf1(3) = sp_decomp%zsz(1)
            sbuf1(4) = sp_decomp%zst(2)
            sbuf1(5) = sp_decomp%zen(2)
            sbuf1(6) = sp_decomp%zsz(2)
            sbuf1(7) = sp_decomp%zst(3)
            sbuf1(8) = sp_decomp%zen(3)
            sbuf1(9) = sp_decomp%zsz(3)
            count = sp_decomp%zsz(1)*sp_decomp%zsz(2)*sp_decomp%zsz(3)
            ! send partition information
            CALL MPI_SEND(sbuf1,9,MPI_INTEGER,cible,nrank,MPI_COMM_WORLD,ierror)
            ! send data array
            CALL MPI_SEND(local,count,complex_type,cible, &
            nrank+nproc,MPI_COMM_WORLD,ierror)
        end if

        return
    end subroutine assemble_global_z

end module decomp2D_utils_complex
