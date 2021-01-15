module bubble_data

    use mpi
    use decomp_2d

    implicit none

    real*8                                   :: d, d_modified, divro, rad_modified, nu_viscosity, g_modified, Re_bubble

    type bubble
        real*8  :: x1
        real*8  :: x2
        real*8  :: x3
        real*8  :: v1
        real*8  :: v2
        real*8  :: v3
        real*8  :: u01
        real*8  :: u02
        real*8  :: u03
        real*8  :: u11
        real*8  :: u12
        real*8  :: u13
        real*8  :: om01
        real*8  :: om02
        real*8  :: om03
        real*8  :: om11
        real*8  :: om12
        real*8  :: om13
        real*8  :: a1
        real*8  :: a2
        real*8  :: a3
        real*8  :: b1
        real*8  :: b2
        real*8  :: b3
        real*8  :: c1
        real*8  :: c2
        real*8  :: c3
        real*8  :: d
        real*8  :: gradtau_tau0
        real*8  :: gradtau_tau1
        real*8  :: proc
    end type bubble

    integer, parameter                      :: MAX_BUBBLES=1000000
    integer, parameter                      :: MAX_BUBBLES_TMP=10000
    integer                                 :: bubble_cpt=0
    integer                                 :: bubble_tmp_cpt=0
    integer, parameter                      :: BUBBLE_SIZE=31

    type(bubble), dimension(MAX_BUBBLES)        :: bubl, bubl_old
    type(bubble), dimension(MAX_BUBBLES_TMP)    :: bubl_tmp
    real*8, dimension(:,:,:), allocatable       :: fb1, fb2, fb3
    real*8, dimension(:,:,:), allocatable       :: fb1_alpha_x, fb2_alpha_x, fb3_alpha_x
    !real*8, dimension(:,:,:), allocatable       :: fb1_alpha_y, fb2_alpha_y
    !real*8, dimension(:,:,:), allocatable       :: fb1_alpha_z, fb2_alpha_z

    type bubble_source
        integer  :: n1
        integer  :: n2
        integer  :: n3
        integer  :: s1
        integer  :: s2
        integer  :: s3
        integer  :: nb_start
        integer  :: nb_regen
    end type bubble_source

    integer, parameter                                  :: MAX_BUBBLEGEN=10
    integer                                             :: bublgen_nb=0
    type(bubble_source), dimension(MAX_BUBBLEGEN)       :: bubl_gen

    integer                                             :: sub_step

    ! Move_bubble data
    real*8, parameter                                   :: Cd=160.d0, Cl=0.5d0

contains

    subroutine add_bubble_source(n1,n2,n3, s1,s2,s3,nb_start,nb_regen)
        implicit none
        integer, intent(in)          :: n1,n2,n3, s1,s2,s3, nb_start, nb_regen

        bublgen_nb=bublgen_nb+1
        bubl_gen(bublgen_nb)%n1=n1
        bubl_gen(bublgen_nb)%n2=n2
        bubl_gen(bublgen_nb)%n3=n3
        bubl_gen(bublgen_nb)%s1=s1
        bubl_gen(bublgen_nb)%s2=s2
        bubl_gen(bublgen_nb)%s3=s3
        bubl_gen(bublgen_nb)%nb_start=nb_start
        bubl_gen(bublgen_nb)%nb_regen=nb_regen

    end subroutine add_bubble_source

    ! Adding a bubble directly into xb and vb tables
    subroutine add_particule(bubbles, n, x1,x2,x3, v1,v2,v3, u01,u02,u03, u11,u12,u13, om01,om02,om03, om11,om12,om13, a1,a2,a3, b1,b2,b3, c1,c2,c3, d, proc)

        implicit none

        real*8, intent(in)          :: x1,x2,x3, v1,v2,v3, u01,u02,u03, u11,u12,u13, om01,om02,om03, om11,om12,om13, a1,a2,a3, b1,b2,b3, c1,c2,c3, d, proc
        integer                     :: n
        type(bubble), dimension(:)  :: bubbles

        if((x1==-1.d0).and.(x2==-1.d0).and.(x3==-1.d0)) then
            return

        else if (n+1>MAX_BUBBLES) then
            write(*,*)'MAX NUMBER OF BUBBLE REACHED'
            call exit

        else
            n=n+1

            bubbles(n)%x1=x1
            bubbles(n)%x2=x2
            bubbles(n)%x3=x3

            bubbles(n)%v1=v1
            bubbles(n)%v2=v2
            bubbles(n)%v3=v3

            bubbles(n)%u01=u01
            bubbles(n)%u02=u02
            bubbles(n)%u03=u03

            bubbles(n)%u11=u11
            bubbles(n)%u12=u12
            bubbles(n)%u13=u13

            bubbles(n)%om01=om01
            bubbles(n)%om02=om02
            bubbles(n)%om03=om03

            bubbles(n)%om11=om11
            bubbles(n)%om12=om12
            bubbles(n)%om13=om13

            bubbles(n)%a1=a1
            bubbles(n)%a2=a2
            bubbles(n)%a3=a3

            bubbles(n)%b1=b1
            bubbles(n)%b2=b2
            bubbles(n)%b3=b3

            bubbles(n)%c1=c1
            bubbles(n)%c2=c2
            bubbles(n)%c3=c3

            bubbles(n)%d=d

            bubbles(n)%proc=proc
        endif


    end subroutine add_particule


    ! Removing a bubble contained into xb & vb
    subroutine remove_bubble(i)
        implicit none
        integer :: i,  n

        n=bubble_cpt

        bubl(i)=bubl(n)

        bubble_cpt=bubble_cpt-1

    end subroutine remove_bubble

    subroutine display_all_bubbles(n)
        implicit none

        integer, intent(in) :: n
        integer             :: mpi_err, i

        if (nrank==5) write(*,*) '-----------------DISPLAY BUBBLES-----------------'
        write (*,*) '******Processor :', nrank, ' ********'
        write(*,*)'nb bulles:', bubble_cpt
        do i = 1, bubble_cpt
            write(*,*) ''
            write (*,*) bubl(i)%x1, bubl(i)%x2,bubl(i)%x3
            write(*,*) ''
        end do

    end subroutine display_all_bubbles

end module bubble_data



module bubble_parallel_data

    use decomp_2d
    use mpi
    use bubble_data

    implicit none

    real*8, dimension(:,:,:), allocatable   :: send_list !(processor, particle number, vector xb(1:3)&vb(4:6))
    integer, dimension(:), allocatable      :: send_nb, receive_nb

  interface post_bubble
    module procedure post_bubble, post_bubble_array
  end interface

contains

    ! Initializes the local tables
    subroutine reset_send_tables()

        send_nb(:)=1
        send_list(:,1,:)=0.d0
        send_list(:,1,1:3)=-1.d0

    end subroutine reset_send_tables

    ! Post_bubble prepares a bubble to be send, like a post office
    subroutine post_bubble_array(proc_num, pt)

        implicit none

        integer, intent(in)                         :: proc_num
        real*8, dimension(BUBBLE_SIZE), intent(in)  :: pt
        integer                                     :: nb_bubble

        nb_bubble=send_nb(proc_num)+1
        send_nb(proc_num)=nb_bubble

        send_list(proc_num, nb_bubble, :)=pt(:)

    end subroutine post_bubble_array

    subroutine post_bubble(proc_num, bub)

        implicit none

        type(bubble)                        :: bub
        integer, intent(in)                 :: proc_num
        real*8, dimension(BUBBLE_SIZE)      :: pt
        integer :: nb_bubble

        pt(1)=bub%x1
        pt(2)=bub%x2
        pt(3)=bub%x3

        pt(4)=bub%v1
        pt(5)=bub%v2
        pt(6)=bub%v3

        pt(7)=bub%u01
        pt(8)=bub%u02
        pt(9)=bub%u03

        pt(10)=bub%u11
        pt(11)=bub%u12
        pt(12)=bub%u13

        pt(13)=bub%om01
        pt(14)=bub%om02
        pt(15)=bub%om03

        pt(16)=bub%om11
        pt(17)=bub%om12
        pt(18)=bub%om13

        pt(19)=bub%a1
        pt(20)=bub%a2
        pt(21)=bub%a3

        pt(22)=bub%b1
        pt(23)=bub%b2
        pt(24)=bub%b3

        pt(25)=bub%c1
        pt(26)=bub%c2
        pt(27)=bub%c3

        pt(28)=bub%d

        pt(29)=bub%proc
        pt(30)=bub%gradtau_tau0
        pt(31)=bub%gradtau_tau1

        nb_bubble=send_nb(proc_num)+1
        send_nb(proc_num)=nb_bubble

        send_list(proc_num, nb_bubble, :)=pt(:)

    end subroutine post_bubble


    ! Sends a posted bubble to the assigned processors
    subroutine send_bubble(bubbles, n)

        use bubble_data

        implicit none

        ! BUBBLE_SIZE   : number of variables that defines the bubble (x1, x2, x3, v1, v2, v3, etc.)
        ! target        : The current process where the bubbles are sent
        ! sbuff, rbuff  : Send & receive buffer

        integer                                 :: target
        integer                                 :: i, j, n
        integer, dimension(BUBBLE_SIZE)         :: rcount, disp
        real*8, dimension(:), allocatable       :: rbuff, sbuff
        integer                                 :: mpi_err
        integer                                 :: received_bubbles
        type(bubble), dimension(:), intent(in)  :: bubbles

        do target = 0, nproc-1

            ! Preparing target to receive data
            call MPI_GATHER (send_nb(target), 1, MPI_INTEGER , receive_nb(:), 1, MPI_INTEGER ,target, MPI_COMM_WORLD ,mpi_err)

            ! ***************************************************************
            !_________________Setting arguments of MPI_GATHERV_______________
            ! ***************************************************************

            ! Displacement and size
            disp(1)=0
            rcount(1)=BUBBLE_SIZE*receive_nb(0)
            do i = 2, nproc
                rcount(i)=receive_nb(i-1)*BUBBLE_SIZE
                disp(i)=rcount(i-1)+disp(i-1)
            end do

            ! Send buffer
            if (allocated(sbuff)) deallocate(sbuff)
            allocate(sbuff(send_nb(target)*BUBBLE_SIZE))
            received_bubbles=sum(receive_nb) !Every particles coming from processors is targetted to nrank
            do i = 1, send_nb(target)
                sbuff( BUBBLE_SIZE*(i-1)+1 : BUBBLE_SIZE*i )=send_list(target, i, 1:BUBBLE_SIZE)
            end do

            ! Allocate receive buffer
            if (allocated(rbuff)) deallocate(rbuff)
            if (nrank==target) allocate(rbuff(received_bubbles*BUBBLE_SIZE))

            ! ***************************************************************
            !_________________Calling MPI_GATHERV____________________________
            ! ***************************************************************

            call MPI_Gatherv(sbuff(:), send_nb(target)*BUBBLE_SIZE, MPI_DOUBLE_PRECISION, rbuff(:), rcount, disp, MPI_DOUBLE_PRECISION, target, MPI_COMM_WORLD, mpi_err);


            ! ***************************************************************
            !_________________Get bubbles from receive buffer________________
            ! ***************************************************************
            if (nrank==target) then
                do i = 1, received_bubbles
                    j = BUBBLE_SIZE*(i-1)+1
                    !! //////////////////AJOUT ARGUMENT DANS BUBBLE
                    call add_particule(bubbles, n, rbuff(j), rbuff(j+1), rbuff(j+2), rbuff(j+3), rbuff(j+4), rbuff(j+5), rbuff(j+6), rbuff(j+7), rbuff(j+8), rbuff(j+9), rbuff(j+10), rbuff(j+11), rbuff(j+12), rbuff(j+13), rbuff(j+14), rbuff(j+15), rbuff(j+16), rbuff(j+17), rbuff(j+18), rbuff(j+19), rbuff(j+20), rbuff(j+21), rbuff(j+22), rbuff(j+23), rbuff(j+24), rbuff(j+25), rbuff(j+26), rbuff(j+27), rbuff(j+28))
                end do
            endif

        end do

        call reset_send_tables

    end subroutine send_bubble


    subroutine display_send_tables(p)
        implicit none
        integer :: p

        integer :: i, j, k

        if (nrank==p) then

            write(6,*) '----------------SEND TABLES------------------'
            write(6,*) ''
            write(6,*) '******Table send_nb'


            write(6,*) 'Proc',  (/ (i, i=0,nproc-1) /)
            write(6,*) 'Send', send_nb(:)

            write(6,*) ''
            write(6,*) '******Table send_list'

            do i = 0, nproc-1
                write(6,*) 'Proc', i
                do j = 1, send_nb(i)
                    write(6,*) 'Bubble', j,':', send_list(i, j, :)

                end do
            end do

            write(6,*) ''
            write(6,*) '--------------------END----------------------'

        endif

    end subroutine display_send_tables

    subroutine display_receive_tables(p)
        implicit none
        integer :: p

        integer :: i, j, kmodule

        if (nrank==p) then

            write(6,*) '----------------RECEIVE TABLES------------------'
            write(6,*) ''
            write(6,*) '******Table receive_nb'
            write(6,*) 'Proc',  (/ (i, i=0,nproc-1) /)
            write(6,*) 'Rece', receive_nb(:)
            write(6,*) ''
            write(6,*) '---------------------END------------------------'

        endif

    end subroutine display_receive_tables


end module bubble_parallel_data
