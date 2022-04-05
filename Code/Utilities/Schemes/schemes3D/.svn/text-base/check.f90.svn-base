module commun
    use O2_D0s
    use O2_D0s_sh
    use O2_D1s
    use O2_D1s_sh
    use O2_D1c
    use O2_D2c
    use DRP_D0s
    use DRP_D0s_sh
    use DRP_D1s
    use DRP_D1s_sh
    use DRP_D1c
    use DRP_D2c
    use CPT_D0s
    use CPT_D0s_sh
    use CPT_D1s
    use CPT_D1s_sh
    use CPT_D1c
    use CPT_D2c
    implicit none
    integer, parameter  :: n1=32,n2=64,n3=128
    real*8, parameter   :: PI=(2.d0*dasin(1.d0))
    real*8, parameter   :: h1=2.d0*pi/(n1-1), h2=2.d0*pi/(n2-1), h3=2.d0*pi/(n3-1)
    integer :: i,j,k

    type configuration
        procedure(D0s_DRP5_3Dx), pointer, nopass    :: D0s_3Dx
        procedure(D0ssh_DRP5_3Dx), pointer, nopass  :: D0ssh_3Dx
        procedure(D0s_DRP5_3Dy), pointer, nopass    :: D0s_3Dy
        procedure(D0ssh_DRP5_3Dy), pointer, nopass    :: D0ssh_3Dy
        procedure(D0ssh_DRP5_MULT_3Dy), pointer, nopass    :: D0ssh_MULT_3Dy
        procedure(D0s_DRP5_3Dz), pointer, nopass    :: D0s_3Dz
        procedure(D0ssh_DRP5_3Dz), pointer, nopass    :: D0ssh_3Dz
        procedure(D0ssh_DRP5_MULT_3Dz), pointer, nopass    :: D0ssh_MULT_3Dz

        procedure(D1s_DRP5_3Dx), pointer, nopass    :: D1s_3Dx
        procedure(D1s_DRP5_ACC_3Dx), pointer, nopass    :: D1s_ACC_3Dx
        procedure(D1ssh_DRP5_3Dx), pointer, nopass    :: D1ssh_3Dx
        procedure(D1ssh_DRP5_ACC_3Dx), pointer, nopass    :: D1ssh_ACC_3Dx
        procedure(D1s_DRP5_3Dy), pointer, nopass    :: D1s_3Dy
        procedure(D1ssh_DRP5_3Dy), pointer, nopass    :: D1ssh_3Dy
        procedure(D1s_DRP5_MULT_3Dy), pointer, nopass    :: D1s_MULT_3Dy
        procedure(D1s_DRP5_MULTACC_3Dy), pointer, nopass    :: D1s_MULTACC_3Dy
        procedure(D1ssh_DRP5_MULT_3Dy), pointer, nopass    :: D1ssh_MULT_3Dy
        procedure(D1ssh_DRP5_MULTACC_3Dy), pointer, nopass    :: D1ssh_MULTACC_3Dy
        procedure(D1s_DRP5_3Dz), pointer, nopass    :: D1s_3Dz
        procedure(D1ssh_DRP5_3Dz), pointer, nopass    :: D1ssh_3Dz
        procedure(D1s_DRP5_ACC_3Dz), pointer, nopass    :: D1s_ACC_3Dz
        procedure(D1ssh_DRP5_ACC_3Dz), pointer, nopass    :: D1ssh_ACC_3Dz

        procedure(D1c_DRP6_3Dx), pointer, nopass    :: D1c_3Dx
        procedure(D1c_DRP6_3Dy), pointer, nopass    :: D1c_3Dy
        procedure(D1c_DRP6_MULT_3Dy), pointer, nopass    :: D1c_MULT_3Dy
        procedure(D1c_DRP6_MULTACC_3Dy), pointer, nopass    :: D1c_MULTACC_3Dy
        procedure(D1c_DRP6_3Dz), pointer, nopass    :: D1c_3Dz

        procedure(D2c_DRP6_3Dx), pointer, nopass    :: D2c_3Dx
        procedure(D2c_DRP6_ACC_3Dx), pointer, nopass    :: D2c_ACC_3Dx
        procedure(D2c_DRP6_3Dy), pointer, nopass    :: D2c_3Dy
        procedure(D2c_DRP6_MULT_3Dy), pointer, nopass    :: D2c_MULT_3Dy
        procedure(D2c_DRP6_MULTACC_3Dy), pointer, nopass    :: D2c_MULTACC_3Dy
        procedure(D2c_DRP6_3Dz), pointer, nopass    :: D2c_3Dz
        character(20)                               :: name

    end type configuration

    type(configuration) :: O2_configuration, DRP_configuration, CPT_configuration
contains

end module commun

program check
    use commun
    implicit none

    O2_configuration%name="O2"
    O2_configuration%D0s_3Dx=>D0s_O2_3Dx
    O2_configuration%D0ssh_3Dx=>D0ssh_O2_3Dx
    O2_configuration%D0s_3Dy=>D0s_O2_3Dy
    O2_configuration%D0ssh_3Dy=>D0ssh_O2_3Dy
    O2_configuration%D0ssh_MULT_3Dy=>D0ssh_O2_MULT_3Dy
    O2_configuration%D0s_3Dz=>D0s_O2_3Dz
    O2_configuration%D0ssh_3Dz=>D0ssh_O2_3Dz
    O2_configuration%D0ssh_MULT_3Dz=>D0ssh_O2_MULT_3Dz

    O2_configuration%D1s_3Dx=>D1s_O2_3Dx
    O2_configuration%D1s_ACC_3Dx=>D1s_O2_ACC_3Dx
    O2_configuration%D1ssh_3Dx=>D1ssh_O2_3Dx
    O2_configuration%D1ssh_ACC_3Dx=>D1ssh_O2_ACC_3Dx
    O2_configuration%D1s_3Dy=>D1s_O2_3Dy
    O2_configuration%D1ssh_3Dy=>D1ssh_O2_3Dy
    O2_configuration%D1s_MULT_3Dy=>D1s_O2_MULT_3Dy
    O2_configuration%D1ssh_MULT_3Dy=>D1ssh_O2_MULT_3Dy
    O2_configuration%D1s_MULTACC_3Dy=>D1s_O2_MULTACC_3Dy
    O2_configuration%D1ssh_MULTACC_3Dy=>D1ssh_O2_MULTACC_3Dy
    O2_configuration%D1s_3Dz=>D1s_O2_3Dz
    O2_configuration%D1ssh_3Dz=>D1ssh_O2_3Dz
    O2_configuration%D1s_ACC_3Dz=>D1s_O2_ACC_3Dz
    O2_configuration%D1ssh_ACC_3Dz=>D1ssh_O2_ACC_3Dz

    O2_configuration%D1c_3Dx=>D1c_O2_3Dx
    O2_configuration%D1c_3Dy=>D1c_O2_3Dy
    O2_configuration%D1c_MULT_3Dy=>D1c_O2_MULT_3Dy
    O2_configuration%D1c_MULTACC_3Dy=>D1c_O2_MULTACC_3Dy
    O2_configuration%D1c_3Dz=>D1c_O2_3Dz

    O2_configuration%D2c_3Dx=>D2c_O2_3Dx
    O2_configuration%D2c_ACC_3Dx=>D2c_O2_ACC_3Dx
    O2_configuration%D2c_3Dy=>D2c_O2_3Dy
    O2_configuration%D2c_MULT_3Dy=>D2c_O2_MULT_3Dy
    O2_configuration%D2c_MULTACC_3Dy=>D2c_O2_MULTACC_3Dy
    O2_configuration%D2c_3Dz=>D2c_O2_3Dz




    DRP_configuration%name="DRP"
    DRP_configuration%D0s_3Dx=>D0s_DRP5_3Dx
    DRP_configuration%D0ssh_3Dx=>D0ssh_DRP5_3Dx
    DRP_configuration%D0s_3Dy=>D0s_DRP5_3Dy
    DRP_configuration%D0ssh_3Dy=>D0ssh_DRP5_3Dy
    DRP_configuration%D0ssh_MULT_3Dy=>D0ssh_DRP5_MULT_3Dy
    DRP_configuration%D0s_3Dz=>D0s_DRP5_3Dz
    DRP_configuration%D0ssh_3Dz=>D0ssh_DRP5_3Dz
    DRP_configuration%D0ssh_MULT_3Dz=>D0ssh_DRP5_MULT_3Dz

    DRP_configuration%D1s_3Dx=>D1s_DRP5_3Dx
    DRP_configuration%D1s_ACC_3Dx=>D1s_DRP5_ACC_3Dx
    DRP_configuration%D1ssh_3Dx=>D1ssh_DRP5_3Dx
    DRP_configuration%D1ssh_ACC_3Dx=>D1ssh_DRP5_ACC_3Dx
    DRP_configuration%D1s_3Dy=>D1s_DRP5_3Dy
    DRP_configuration%D1ssh_3Dy=>D1ssh_DRP5_3Dy
    DRP_configuration%D1s_MULT_3Dy=>D1s_DRP5_MULT_3Dy
    DRP_configuration%D1s_MULTACC_3Dy=>D1s_DRP5_MULTACC_3Dy
    DRP_configuration%D1ssh_MULT_3Dy=>D1ssh_DRP5_MULT_3Dy
    DRP_configuration%D1ssh_MULTACC_3Dy=>D1ssh_DRP5_MULTACC_3Dy
    DRP_configuration%D1s_3Dz=>D1s_DRP5_3Dz
    DRP_configuration%D1ssh_3Dz=>D1ssh_DRP5_3Dz
    DRP_configuration%D1s_ACC_3Dz=>D1s_DRP5_ACC_3Dz
    DRP_configuration%D1ssh_ACC_3Dz=>D1ssh_DRP5_ACC_3Dz

    DRP_configuration%D1c_3Dx=>D1c_DRP6_3Dx
    DRP_configuration%D1c_3Dy=>D1c_DRP6_3Dy
    DRP_configuration%D1c_MULT_3Dy=>D1c_DRP6_MULT_3Dy
    DRP_configuration%D1c_MULTACC_3Dy=>D1c_DRP6_MULTACC_3Dy
    DRP_configuration%D1c_3Dz=>D1c_DRP6_3Dz

    DRP_configuration%D2c_3Dx=>D2c_DRP6_3Dx
    DRP_configuration%D2c_ACC_3Dx=>D2c_DRP6_ACC_3Dx
    DRP_configuration%D2c_3Dy=>D2c_DRP6_3Dy
    DRP_configuration%D2c_MULT_3Dy=>D2c_DRP6_MULT_3Dy
    DRP_configuration%D2c_MULTACC_3Dy=>D2c_DRP6_MULTACC_3Dy
    DRP_configuration%D2c_3Dz=>D2c_DRP6_3Dz




    CPT_configuration%name="CPT"
    CPT_configuration%D0s_3Dx=>D0s_CPT_3Dx
    CPT_configuration%D0ssh_3Dx=>D0ssh_CPT_3Dx
    CPT_configuration%D0s_3Dy=>D0s_CPT_3Dy
    CPT_configuration%D0ssh_3Dy=>D0ssh_CPT_3Dy
    CPT_configuration%D0ssh_MULT_3Dy=>D0ssh_CPT_MULT_3Dy
    CPT_configuration%D0s_3Dz=>D0s_CPT_3Dz
    CPT_configuration%D0ssh_3Dz=>D0ssh_CPT_3Dz
    CPT_configuration%D0ssh_MULT_3Dz=>D0ssh_CPT_MULT_3Dz

    CPT_configuration%D1s_3Dx=>D1s_CPT_3Dx
    CPT_configuration%D1s_ACC_3Dx=>D1s_CPT_ACC_3Dx
    CPT_configuration%D1ssh_3Dx=>D1ssh_CPT_3Dx
    CPT_configuration%D1ssh_ACC_3Dx=>D1ssh_CPT_ACC_3Dx
    CPT_configuration%D1s_3Dy=>D1s_CPT_3Dy
    CPT_configuration%D1ssh_3Dy=>D1ssh_CPT_3Dy
    CPT_configuration%D1s_MULT_3Dy=>D1s_CPT_MULT_3Dy
    CPT_configuration%D1ssh_MULT_3Dy=>D1ssh_CPT_MULT_3Dy
    CPT_configuration%D1s_MULTACC_3Dy=>D1s_CPT_MULTACC_3Dy
    CPT_configuration%D1ssh_MULTACC_3Dy=>D1ssh_CPT_MULTACC_3Dy
    CPT_configuration%D1s_3Dz=>D1s_CPT_3Dz
    CPT_configuration%D1ssh_3Dz=>D1ssh_CPT_3Dz
    CPT_configuration%D1s_ACC_3Dz=>D1s_CPT_ACC_3Dz
    CPT_configuration%D1ssh_ACC_3Dz=>D1ssh_CPT_ACC_3Dz

    CPT_configuration%D1c_3Dx=>D1c_CPT_3Dx
    CPT_configuration%D1c_3Dy=>D1c_CPT_3Dy
    CPT_configuration%D1c_MULT_3Dy=>D1c_CPT_MULT_3Dy
    CPT_configuration%D1c_MULTACC_3Dy=>D1c_CPT_MULTACC_3Dy
    CPT_configuration%D1c_3Dz=>D1c_CPT_3Dz

    CPT_configuration%D2c_3Dx=>D2c_CPT_3Dx
    CPT_configuration%D2c_ACC_3Dx=>D2c_CPT_ACC_3Dx
    CPT_configuration%D2c_3Dy=>D2c_CPT_3Dy
    CPT_configuration%D2c_MULT_3Dy=>D2c_CPT_MULT_3Dy
    CPT_configuration%D2c_MULTACC_3Dy=>D2c_CPT_MULTACC_3Dy
    CPT_configuration%D2c_3Dz=>D2c_CPT_3Dz

    write(*,*)'_______________________'
    call test_3Dx(DRP_configuration)
    !call test_3Dx2(DRP_configuration)
    write(*,*)'_______________________'
    call test_3Dy(DRP_configuration)
    write(*,*)'_______________________'
    call test_3Dz(DRP_configuration)

    write(*,*)'_______________________'
    call test_3Dx(O2_configuration)
    !call test_3Dx2(O2_configuration)
    write(*,*)'_______________________'
    call test_3Dy(O2_configuration)
    write(*,*)'_______________________'
    call test_3Dz(O2_configuration)

    write(*,*)'_______________________'
    call test_3Dx(CPT_configuration)
    write(*,*)'_______________________'
    call test_3Dy(CPT_configuration)
    write(*,*)'_______________________'
    call test_3Dz(CPT_configuration)

end program check


subroutine fillx(values3D, fx, n1,n2,n3)
    implicit none
    integer :: n1,n2,n3
    real*8, dimension(n1,n2,n3) :: values3D
    real*8, dimension(n1)       :: fx
    integer                     :: i,j,k

    do k = 1, n3
        do j = 1, n2
            values3D(:,j,k)=fx
        end do
    end do

end subroutine fillx


subroutine filly(values3D, fy, n1,n2,n3)
    implicit none
    integer :: n1,n2,n3
    real*8, dimension(n1,n2,n3) :: values3D
    real*8, dimension(n2)       :: fy
    integer                     :: i,j,k

    do k = 1, n3
        do i = 1, n1
            values3D(i,:,k)=fy
        end do
    end do

end subroutine filly


subroutine fillz(values3D, fz, n1,n2,n3)
    implicit none
    integer :: n1,n2,n3
    real*8, dimension(n1,n2,n3) :: values3D
    real*8, dimension(n3)       :: fz
    integer                     :: i,j,k

    do j = 1, n2
        do i = 1, n1
            values3D(i,j,:)=fz
        end do
    end do

end subroutine fillz


subroutine diff_x(values, expected, x, n1, n1s, n1e, n2, n2s, n2e, n3, n3s, n3e, file_name)
    implicit none
    integer :: n1, n1s, n1e, n2, n2s, n2e, n3, n3s, n3e, fid
    real*8, dimension(n1,n2,n3) :: values
    real*8, dimension(n1)   :: expected, x, max_error, error
    real*8      :: error_mean
    character(*)    :: file_name
    character(30)   :: warn_flag, file_flag
    integer :: i,j,k

    error=0.d0
    max_error=0.d0

    do k = n3s, n3e
        do j = n2s, n2e
            do i = n1s, n1e
                error(i)=expected(i)-values(i,j,k)
            end do

            if (sum(abs(error))>sum(abs(max_error))) max_error=error
        end do
    end do

    fid=15
    open(fid, file=trim(file_name)//".csv")

    do i = n1s, n1e
        write(fid,*) x(i), max_error(i)
    end do

    close(fid)


    error_mean=sum(abs(max_error))/(n1e-n1s)
    file_flag=trim(file_name)
    if (error_mean>0.0001d0) then
        warn_flag="<-- WARNING"
    else
        warn_flag="<-- OK"
    endif

    write(*,*)file_flag,':', error_mean, trim(warn_flag)

end subroutine diff_x


subroutine diff_y(values, expected, y, n1, n1s, n1e, n2, n2s, n2e, n3, n3s, n3e, file_name)
    implicit none
    integer :: n1, n1s, n1e, n2, n2s, n2e, n3, n3s, n3e, fid
    real*8, dimension(n1,n2,n3) :: values
    real*8, dimension(n2)   :: expected, y, max_error, error
    real*8      :: error_mean
    character(*)    :: file_name
    character(30)   :: warn_flag, file_flag
    integer :: i,j,k

    error=0.d0
    max_error=0.d0

    do k = n3s, n3e
        do i = n1s, n1e
            do j = n2s, n2e
                error(j)=expected(j)-values(i,j,k)
            end do

            if (sum(abs(error))>sum(abs(max_error))) max_error=error
        end do
    end do

    fid=15
    open(fid, file=trim(file_name)//".csv")

    do j = n2s, n2e
        write(fid,*) y(j), max_error(j)
    end do

    close(fid)

    error_mean=sum(abs(max_error))/(n2e-n2s)
    file_flag=trim(file_name)
    if (error_mean>0.0001d0) then
        warn_flag="<-- WARNING"
    else
        warn_flag="<-- OK"
    endif

    write(*,*)file_flag,':', error_mean, trim(warn_flag)

end subroutine diff_y


subroutine diff_z(values, expected, z, n1, n1s, n1e, n2, n2s, n2e, n3, n3s, n3e, file_name)
    implicit none
    integer :: n1, n1s, n1e, n2, n2s, n2e, n3, n3s, n3e, fid
    real*8, dimension(n1,n2,n3) :: values
    real*8, dimension(n3)   :: expected, z, max_error, error
    real*8      :: error_mean
    character(*)    :: file_name
    character(30)   :: warn_flag, file_flag
    integer :: i,j,k

    error=0.d0
    max_error=0.d0

    do j = n2s, n2e
        do i = n1s, n1e
            do k = n3s, n3e
                error(k)=expected(k)-values(i,j,k)
            end do

            if (sum(abs(error))>sum(abs(max_error))) max_error=error
        end do
    end do

    fid=15
    open(fid, file=trim(file_name)//".csv")

    do k = n3s, n3e
        write(fid,*) z(k), max_error(k)
    end do

    close(fid)

    error_mean=sum(abs(max_error))/(n3e-n3s)
    file_flag=trim(file_name)
    if (error_mean>0.0001d0) then
        warn_flag="<-- WARNING"
    else
        warn_flag="<-- OK"
    endif

    write(*,*)file_flag,':', error_mean, trim(warn_flag)

end subroutine diff_z


subroutine test_3Dx(config)
    use DRP_D0s
    use DRP_D0s_sh
    use commun
    implicit none
    integer :: nb
    character(10)   :: tmp_str

    type(configuration) :: config
    real*8, dimension(n1,n2,n3) ::  fx_cos_c, fx_cos_s, fx_sin_c, fx_sin_s
    real*8, dimension(n1,n2,n3) ::  d0fx_cos_c, d0fx_cos_s, d0fx_sin_c, d0fx_sin_s
    real*8, dimension(n1,n2,n3) ::  d1fx_cos_c, d1fx_cos_s, d1fx_sin_c, d1fx_sin_s
    real*8, dimension(n1,n2,n3) ::  d2fx_cos_c, d2fx_cos_s, d2fx_sin_c, d2fx_sin_s
    real*8, dimension(n1)       :: expected, error_max, sin_s, cos_s, sin_c, cos_c, f_square_s, f_square_c, mult, zero

    real*8, dimension(n1)  :: xc, xs
    zero=0.d0
    nb=0

    do i = 1, n1
        xc(i)=((i-1.d0)/(n1-1.d0))
        xs(i)=((i-0.5d0)/(n1-1.d0))
    end do


    do i = 1, n1
        cos_c(i)=dcos(2.d0*pi*xc(i))
        cos_s(i)=dcos(2.d0*pi*xs(i))
        sin_c(i)=dsin(2.d0*pi*xc(i))
        sin_s(i)=dsin(2.d0*pi*xs(i))
    end do

    call fillx(fx_cos_c, cos_c, n1,n2,n3)
    call fillx(fx_cos_s, cos_s, n1,n2,n3)
    call fillx(fx_sin_c, sin_c, n1,n2,n3)
    call fillx(fx_sin_s, sin_s, n1,n2,n3)

    nb=nb+1
    write(tmp_str, "(i10)")nb
    ! D0s   *************************************************************************************************
    ! PERIODIC
    d0fx_cos_s=1000.d0
    call config%D0s_3Dx(fx_cos_c, d0fx_cos_s, n1,n2,n2,n3,n3, periodic)
    call diff_x(d0fx_cos_s, cos_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0S_SIMPLE_PERIODIC")
    ! PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fx_cos_c=1000.d0
    call config%D0ssh_3Dx(fx_cos_s, d0fx_cos_c, n1,n2,n2,n3,n3, periodic)
    call diff_x(d0fx_cos_c, cos_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fx_cos_s=1000.d0
    call config%D0s_3Dx(fx_cos_c, d0fx_cos_s, n1,n2,n2,n3,n3, Dirichlet)
    call diff_x(d0fx_cos_s, cos_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0S_SIMPLE_DIRICHLET")
    ! DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fx_cos_c=1000.d0
    call config%D0ssh_3Dx(fx_cos_s, d0fx_cos_c, n1,n2,n2,n3,n3, Dirichlet)
    call diff_x(d0fx_cos_c, cos_c, xc, n1, 2, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_SIMPLE_DIRICHLET")


    ! D1s   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_cos_s=1000.d0
    call config%D1s_3Dx(fx_cos_c, d1fx_cos_s, n1,n2,n2,n3,n3, h1, .false., periodic)
    call diff_x(d1fx_cos_s, -1.d0*sin_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_PERIODIC")
    ! PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_cos_c=1000.d0
    call config%D1ssh_3Dx(fx_cos_s, d1fx_cos_c, n1,n2,n2,n3,n3, h1, .false., periodic)
    call diff_x(d1fx_cos_c, -1.d0*sin_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_PERIODIC")

    ! ANTISYMETRIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_sin_s=1000.d0
    call config%D1s_3Dx(fx_sin_c, d1fx_sin_s, n1,n2,n2,n3,n3, h1, .false., antisymetric)
    call diff_x(d1fx_sin_s, cos_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_ANTISYMETRIC")
    ! SYMETRIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_cos_c=1000.d0
    call config%D1ssh_3Dx(fx_cos_s, d1fx_cos_c, n1,n2,n2,n3,n3, h1, .false., symetric)
    call diff_x(d1fx_cos_c, -1.d0*sin_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_SYMETRIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_cos_s=1000.d0
    call config%D1s_3Dx(fx_cos_c, d1fx_cos_s, n1,n2,n2,n3,n3, h1, .false., Dirichlet)
    call diff_x(d1fx_cos_s, -1.d0*sin_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_DIRICHLET")
    ! DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_cos_c=1000.d0
    call config%D1ssh_3Dx(fx_cos_s, d1fx_cos_c, n1,n2,n2,n3,n3, h1, .false., Dirichlet)
    call diff_x(d1fx_cos_c, -1.d0*sin_c, xc, n1, 2, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_DIRICHLET")


    ! ACC, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d1fx_cos_s, -1.d0*sin_s, n1,n2,n3)
    call config%D1s_ACC_3Dx(fx_cos_c, d1fx_cos_s, n1,n2,n2,n3,n3, h1, .false., periodic)
    call diff_x(d1fx_cos_s, -2.d0*sin_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_ACC_PERIODIC")
    ! ACC, PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d1fx_cos_c, -1.d0*sin_c, n1,n2,n3)
    call config%D1ssh_ACC_3Dx(fx_cos_s, d1fx_cos_c, n1,n2,n2,n3,n3, h1, .false., periodic)
    call diff_x(d1fx_cos_c, -2.d0*sin_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_ACC_PERIODIC")


    ! ACC, ANTISYMETRIC
    !nb=nb+1
    !write(tmp_str, "(i10)")nb
    !call fillx(d1fx_sin_s, cos_s, n1,n2,n3)
    !call config%D1s_ACC_3Dx(fx_sin_c, d1fx_sin_s, n1,n2,n2,n3,n3, h1, .false., antisymetric)
    !call diff_x(d1fx_sin_s, 2.d0*cos_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_ACC_ANTISYMETRIC")

    ! ACC, SYMETRIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d1fx_cos_c, -1.d0*sin_c, n1,n2,n3)
    call config%D1ssh_ACC_3Dx(fx_cos_s, d1fx_cos_c, n1,n2,n2,n3,n3, h1, .false., symetric)
    call diff_x(d1fx_cos_c, -2.d0*sin_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_ACC_SYMETRIC")

    ! ACC, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d1fx_cos_s, -1.d0*sin_s, n1,n2,n3)
    call config%D1s_ACC_3Dx(fx_cos_c, d1fx_cos_s, n1,n2,n2,n3,n3, h1, .false., Dirichlet)
    call diff_x(d1fx_cos_s, -2.d0*sin_s, xs, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_ACC_DIRICHLET")
    ! ACC, DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d1fx_cos_c, -1.d0*sin_c, n1,n2,n3)
    call config%D1ssh_ACC_3Dx(fx_cos_s, d1fx_cos_c, n1,n2,n2,n3,n3, h1, .false., Dirichlet)
    call diff_x(d1fx_cos_c, -2.d0*sin_c, xc, n1, 2, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_ACC_DIRICHLET")

    ! D1c   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_sin_c=1000.d0
    call config%D1c_3Dx(fx_sin_c, d1fx_sin_c, n1,n2,n2,n3,n3, h1, .false., periodic)
    call diff_x(d1fx_sin_c, cos_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_sin_c=1000.d0
    call config%D1c_3Dx(fx_sin_c, d1fx_sin_c, n1,n2,n2,n3,n3, h1, .false., Dirichlet)
    call diff_x(d1fx_sin_c, cos_c, xc, n1, 2, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_SIMPLE_DIRICHLET")
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_sin_s=1000.d0
    call config%D1c_3Dx(fx_sin_s, d1fx_sin_s, n1,n2,n2,n3,n3, h1, .true., Dirichlet)
    call diff_x(d1fx_sin_s, cos_s, xs, n1, 2, n1-2, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1CSH_SIMPLE_DIRICHLET")

    ! D2c   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fx_cos_c=1000.d0
    call config%D2c_3Dx(fx_cos_c, d2fx_cos_c, n1,n2,n2,n3,n3, h1, .false., periodic)
    call diff_x(d2fx_cos_c, -1.d0*cos_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fx_cos_c=1000.d0
    call config%D2c_3Dx(fx_cos_c, d2fx_cos_c, n1,n2,n2,n3,n3, h1, .false., Dirichlet)
    call diff_x(d2fx_cos_c, -1.d0*cos_c, xc, n1, 2, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_SIMPLE_DIRICHLET")
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fx_cos_s=1000.d0
    call config%D2c_3Dx(fx_cos_s, d1fx_cos_s, n1,n2,n2,n3,n3, h1, .true., Dirichlet)
    call diff_x(d1fx_cos_s, -1.d0*cos_s, xs, n1, 2, n1-2, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2CSH_SIMPLE_DIRICHLET")

    ! ACC, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d2fx_cos_c, -1.d0*cos_c, n1,n2,n3)
    call config%D2c_ACC_3Dx(fx_cos_c, d2fx_cos_c, n1,n2,n2,n3,n3, h1, .false., periodic)
    call diff_x(d2fx_cos_c, -2.d0*cos_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_ACC_PERIODIC")

    ! ACC, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d2fx_cos_c, -1.d0*cos_c, n1,n2,n3)
    call config%D2c_ACC_3Dx(fx_cos_c, d2fx_cos_c, n1,n2,n2,n3,n3, h1, .false., Dirichlet)
    call diff_x(d2fx_cos_c, -2.d0*cos_c, xc, n1, 2, n1-1, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_ACC_DIRICHLET")
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillx(d2fx_cos_s, -1.d0*cos_s, n1,n2,n3)
    call config%D2c_ACC_3Dx(fx_cos_s, d2fx_cos_s, n1,n2,n2,n3,n3, h1, .true., Dirichlet)
    call diff_x(d2fx_cos_s, -2.d0*cos_s, xs, n1, 2, n1-2, n2, 1, n2, n3, 1, n3, "X/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2CSH_ACC_DIRICHLET")


end subroutine test_3Dx


subroutine test_3Dx2(config)
    use DRP_D0s
    use DRP_D0s_sh
    use commun
    implicit none
    type(configuration) :: config
    real*8, dimension(n1)       ::  dxc, dxs
    real*8, dimension(n1,n2,n3) ::  fx_c, fx_s
    real*8, dimension(n1,n2,n3) ::  d1fx_c, d1fx_s, d1fxe_c, d1fxe_s
    real*8, dimension(n1)       :: expected, error_max, sin_s, cos_s, sin_c, cos_c, f_square_s, f_square_c, mult, zero

    real*8, dimension(n1)  :: xc, xs
    zero=0.d0

    do i = 1, n1
        xc(i)=((i-1.d0)/(n1-1.d0))
        dxc(i)=1.d0
        xs(i)=((i-0.5d0)/(n1-1.d0))
        dxs(i)=1.d0
    end do


    call fillx(fx_c, xc**2, n1,n2,n3)
    call fillx(fx_s, xs**2, n1,n2,n3)

    call fillx(d1fxe_c, 2.d0*xc, n1,n2,n3)
    call fillx(d1fxe_s, 2.d0*xs, n1,n2,n3)


    ! D1s   *************************************************************************************************

    ! SYMETRIC, shift
    call config%D1ssh_3Dx(fx_s, d1fx_c, n1,n2,n2,n3,n3, 1.d0/(n1-1.d0), .false., symetric)
    call diff_x(d1fx_c, d1fxe_c, xc, n1, 1, n1-1, n2, 1, n2, n3, 1, n3, "X2/"//trim(config%name)//"_SIMPLE_SYMETRIC")

end subroutine test_3Dx2


subroutine test_3Dy(config)
    use DRP_D0s
    use DRP_D0s_sh
    use commun
    implicit none
    type(configuration) :: config
    real*8, dimension(n1,n2,n3) ::  fy_cos_c, fy_cos_s, fy_sin_c, fy_sin_s
    real*8, dimension(n1,n2,n3) ::  d0fy_cos_c, d0fy_cos_s, d0fy_sin_c, d0fy_sin_s
    real*8, dimension(n1,n2,n3) ::  d1fy_cos_c, d1fy_cos_s, d1fy_sin_c, d1fy_sin_s
    real*8, dimension(n1,n2,n3) ::  d2fy_cos_c, d2fy_cos_s, d2fy_sin_c, d2fy_sin_s
    real*8, dimension(n2)       :: expected, error_max, sin_s, cos_s, sin_c, cos_c, f_square_s, f_square_c, mult, zero

    real*8, dimension(n2)  :: yc, ys
    integer :: nb
    character(10)   :: tmp_str
    zero=0.d0
    nb=0

    do j = 1, n2
        yc(j)=((j-1.d0)/(n2-1.d0))
        ys(j)=((j-0.5d0)/(n2-1.d0))
    end do


    do j = 1, n2
        cos_c(j)=dcos(2.d0*pi*yc(j))
        cos_s(j)=dcos(2.d0*pi*ys(j))
        sin_c(j)=dsin(2.d0*pi*yc(j))
        sin_s(j)=dsin(2.d0*pi*ys(j))

        mult(j)=0.5d0*n2*(yc(j)-0.5d0)**2
    end do

    call filly(fy_cos_c, cos_c, n1,n2,n3)
    call filly(fy_cos_s, cos_s, n1,n2,n3)
    call filly(fy_sin_c, sin_c, n1,n2,n3)
    call filly(fy_sin_s, sin_s, n1,n2,n3)

    ! D0s   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fy_cos_s=1000.d0
    call config%D0s_3Dy(fy_cos_c, d0fy_cos_s, n1,n1,n2,n3,n3, periodic)
    call diff_y(d0fy_cos_s, cos_s, ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0S_SIMPLE_PERIODIC")
    ! PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fy_cos_c=1000.d0
    call config%D0ssh_3Dy(fy_cos_s, d0fy_cos_c, n1,n1,n2,n3,n3, periodic)
    call diff_y(d0fy_cos_c, cos_c, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fy_cos_s=1000.d0
    call config%D0s_3Dy(fy_cos_c, d0fy_cos_s, n1,n1,n2,n3,n3, Dirichlet)
    call diff_y(d0fy_cos_s, cos_s, ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0S_SIMPLE_DIRICHLET")
    ! DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fy_cos_c=1000.d0
    call config%D0ssh_3Dy(fy_cos_s, d0fy_cos_c, n1,n1,n2,n3,n3, Dirichlet)
    call diff_y(d0fy_cos_c, cos_c, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_SIMPLE_DIRICHLET")


    ! MULT, PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d0fy_cos_c, cos_c, n1,n2,n3)
    call config%D0ssh_MULT_3Dy(fy_cos_s, d0fy_cos_c, n1,n1,n2,n3,n3, periodic)
    call diff_y(d0fy_cos_c, cos_c**2, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_MULT_PERIODIC")
    ! MULT, DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d0fy_cos_c, cos_c, n1,n2,n3)
    call config%D0ssh_MULT_3Dy(fy_cos_s, d0fy_cos_c, n1,n1,n2,n3,n3, Dirichlet)
    call diff_y(d0fy_cos_c, cos_c**2, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_MULT_DIRICHLET")


    ! D1s   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_s=1000.d0
    call config%D1s_3Dy(fy_cos_c, d1fy_cos_s, n1,n1,n2,n3,n3, h2, .false., periodic)
    call diff_y(d1fy_cos_s, -1.d0*sin_s, ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_PERIODIC")
    ! PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_c=1000.d0
    call config%D1ssh_3Dy(fy_cos_s, d1fy_cos_c, n1,n1,n2,n3,n3, h2, .false., periodic)
    call diff_y(d1fy_cos_c, -1.d0*sin_c, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_PERIODIC")


    ! ANTISYMETRIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_s=1000.d0
    call config%D1s_3Dy(fy_sin_c, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .false., antisymetric)
    call diff_y(d1fy_sin_s, cos_s, ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_ANTISYMETRIC")
    ! SYMETRIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_c=1000.d0
    call config%D1ssh_3Dy(fy_cos_s, d1fy_cos_c, n1,n1,n2,n3,n3, h2, .false., symetric)
    call diff_y(d1fy_cos_c, -1.d0*sin_c, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_SYMETRIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_s=1000.d0
    call config%D1s_3Dy(fy_cos_c, d1fy_cos_s, n1,n1,n2,n3,n3, h2, .false., Dirichlet)
    call diff_y(d1fy_cos_s, -1.d0*sin_s, ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_DIRICHLET")
    ! DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_c=1000.d0
    call config%D1ssh_3Dy(fy_cos_s, d1fy_cos_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet)
    call diff_y(d1fy_cos_c, -1.d0*sin_c, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_DIRICHLET")


    ! MULT, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_s=1000.d0
    call config%D1s_MULT_3Dy(fy_cos_c, d1fy_cos_s, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d1fy_cos_s, (-mult*sin_s), ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_MULT_PERIODIC")
    ! MULT, PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_c=1000.d0
    call config%D1ssh_MULT_3Dy(fy_sin_s, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d1fy_sin_c, mult*cos_c, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_MULT_PERIODIC")

    ! MULT, ANTISYMETRIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_s=1000.d0
    call config%D1s_MULT_3Dy(fy_sin_c, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .false., antisymetric, mult)
    call diff_y(d1fy_sin_s, (mult*cos_s), ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_MULT_ANTISYMETRIC")
    ! MULT, SYMETRIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_c=1000.d0
    call config%D1ssh_MULT_3Dy(fy_cos_s, d1fy_cos_c, n1,n1,n2,n3,n3, h2, .false., symetric, mult)
    call diff_y(d1fy_cos_c, -mult*sin_c, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_MULT_SYMETRIC")


    ! MULT, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_s=1000.d0
    call config%D1s_MULT_3Dy(fy_sin_c, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d1fy_sin_s, (mult*cos_s), ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_MULT_DIRICHLET")
    ! MULT, DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_c=1000.d0
    call config%D1ssh_MULT_3Dy(fy_sin_s, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d1fy_sin_c, mult*cos_c, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_MULT_DIRICHLET")


    ! MULTACC, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_cos_s, -1.d0*sin_s, n1,n2,n3)
    call config%D1s_MULTACC_3Dy(fy_cos_c, d1fy_cos_s, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d1fy_cos_s, (-1.d0*sin_s-mult*sin_s), ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_MULTACC_PERIODIC")
    ! MULTACC, PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_cos_c, -1.d0*sin_c, n1,n2,n3)
    call config%D1ssh_MULTACC_3Dy(fy_cos_s, d1fy_cos_c, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d1fy_cos_c, (-1.d0*sin_c-mult*sin_c), yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_MULTACC_PERIODIC")


    ! MULTACC, ANTISYMETRIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_sin_s, cos_s, n1,n2,n3)
    call config%D1s_MULTACC_3Dy(fy_sin_c, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .false., antisymetric, mult)
    call diff_y(d1fy_sin_s, (cos_s+mult*cos_s), ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_MULTACC_ANTISYMETRIC")
    ! MULTACC, SYMETRIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_cos_c, -1.d0*sin_c, n1,n2,n3)
    call config%D1ssh_MULTACC_3Dy(fy_cos_s, d1fy_cos_c, n1,n1,n2,n3,n3, h2, .false., symetric, mult)
    call diff_y(d1fy_cos_c, (-1.d0*sin_c-mult*sin_c), yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_MULTACC_SYMETRIC")

    ! MULTACC, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_sin_s, cos_s, n1,n2,n3)
    call config%D1s_MULTACC_3Dy(fy_sin_c, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d1fy_sin_s, (cos_s+mult*cos_s), ys, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_MULTACC_DIRICHLET")
    ! MULTACC, DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_sin_c, cos_c, n1,n2,n3)
    call config%D1ssh_MULTACC_3Dy(fy_sin_s, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d1fy_sin_c, (cos_c+mult*cos_c), yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_MULTACC_DIRICHLET")

    ! D1c   *************************************************************************************************
    ! PERIODIC
    d1fy_sin_c=1000.d0
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call config%D1c_3Dy(fy_sin_c, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., periodic)
    call diff_y(d1fy_sin_c, cos_c, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_c=1000.d0
    call config%D1c_3Dy(fy_sin_c, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet)
    call diff_y(d1fy_sin_c, cos_c, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_SIMPLE_DIRICHLET")
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_s=1000.d0
    call config%D1c_3Dy(fy_sin_s, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .true., Dirichlet)
    call diff_y(d1fy_sin_s, cos_s, ys, n1, 1, n1, n2, 2, n2-2, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1CSH_SIMPLE_DIRICHLET")

    ! MULT, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_c=1000.d0
    call config%D1c_MULT_3Dy(fy_sin_c, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d1fy_sin_c, cos_c*mult, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_MULT_PERIODIC")

    ! MULT, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_c=1000.d0
    call config%D1c_MULT_3Dy(fy_sin_c, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d1fy_sin_c, cos_c*mult, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_MULT_DIRICHLET")
    ! MULT, DIRICHLET, shifted
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_sin_s=1000.d0
    call config%D1c_MULT_3Dy(fy_sin_s, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .true., Dirichlet, mult)
    call diff_y(d1fy_sin_s, cos_s*mult, ys, n1, 1, n1, n2, 2, n2-2, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1CSH_MULT_DIRICHLET")

    ! MULTACC, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_sin_c, cos_c, n1,n2,n3)
    call config%D1c_MULTACC_3Dy(fy_sin_c, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d1fy_sin_c, cos_c+cos_c*mult, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_MULTACC_PERIODIC")

    ! MULTACC, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_sin_c, cos_c, n1,n2,n3)
    call config%D1c_MULTACC_3Dy(fy_sin_c, d1fy_sin_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d1fy_sin_c, cos_c+cos_c*mult, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_MULTACC_DIRICHLET")
    ! MULTACC, DIRICHLET, shifted
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_sin_s, cos_s, n1,n2,n3)
    call config%D1c_MULTACC_3Dy(fy_sin_s, d1fy_sin_s, n1,n1,n2,n3,n3, h2, .true., Dirichlet, mult)
    call diff_y(d1fy_sin_s, cos_s+cos_s*mult, ys, n1, 1, n1, n2, 2, n2-2, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1CSH_MULTACC_DIRICHLET")

    ! D2c   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fy_cos_c=1000.d0
    call config%D2c_3Dy(fy_cos_c, d2fy_cos_c, n1,n1,n2,n3,n3, h2, .false., periodic)
    call diff_y(d2fy_cos_c, -1.d0*cos_c, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fy_cos_c=1000.d0
    call config%D2c_3Dy(fy_cos_c, d2fy_cos_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet)
    call diff_y(d2fy_cos_c, -1.d0*cos_c, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_SIMPLE_DIRICHLET")
    ! DIRICHLET, shifted
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_s=1000.d0
    call config%D2c_3Dy(fy_cos_s, d1fy_cos_s, n1,n1,n2,n3,n3, h2, .true., Dirichlet)
    call diff_y(d1fy_cos_s, -1.d0*cos_s, ys, n1, 1, n1, n2, 2, n2-2, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2CSH_SIMPLE_DIRICHLET")



    ! MULT, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fy_cos_c=1000.d0
    call config%D2c_MULT_3Dy(fy_cos_c, d2fy_cos_c, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d2fy_cos_c, -1.d0*cos_c*mult, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_MULT_PERIODIC")

    ! MULT, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fy_cos_c=1000.d0
    call config%D2c_MULT_3Dy(fy_cos_c, d2fy_cos_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d2fy_cos_c, -1.d0*cos_c*mult, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_MULT_DIRICHLET")
    ! MULT, DIRICHLET, shifted
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fy_cos_s=1000.d0
    call config%D2c_MULT_3Dy(fy_cos_s, d1fy_cos_s, n1,n1,n2,n3,n3, h2, .true., Dirichlet, mult)
    call diff_y(d1fy_cos_s, -1.d0*cos_s*mult, ys, n1, 1, n1, n2, 2, n2-2, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2CSH_MULT_DIRICHLET")



    ! MULTACC, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d2fy_cos_c, -cos_c, n1,n2,n3)
    call config%D2c_MULTACC_3Dy(fy_cos_c, d2fy_cos_c, n1,n1,n2,n3,n3, h2, .false., periodic, mult)
    call diff_y(d2fy_cos_c, -cos_c-1.d0*cos_c*mult, yc, n1, 1, n1, n2, 1, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_MULTACC_PERIODIC")

    ! MULTACC, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d2fy_cos_c, -cos_c, n1,n2,n3)
    call config%D2c_MULTACC_3Dy(fy_cos_c, d2fy_cos_c, n1,n1,n2,n3,n3, h2, .false., Dirichlet, mult)
    call diff_y(d2fy_cos_c, -cos_c-1.d0*cos_c*mult, yc, n1, 1, n1, n2, 2, n2-1, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_MULTACC_DIRICHLET")
    ! MULTACC, DIRICHLET, shifted
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call filly(d1fy_cos_s, -cos_s, n1,n2,n3)
    call config%D2c_MULTACC_3Dy(fy_cos_s, d1fy_cos_s, n1,n1,n2,n3,n3, h2, .true., Dirichlet, mult)
    call diff_y(d1fy_cos_s, -cos_s-1.d0*cos_s*mult, ys, n1, 1, n1, n2, 2, n2-2, n3, 1, n3, "Y/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2CSH_MULTACC_DIRICHLET")



end subroutine test_3Dy


subroutine test_3Dz(config)
    use DRP_D0s
    use DRP_D0s_sh
    use commun
    implicit none
    type(configuration) :: config
    real*8, dimension(n1,n2,n3) ::  fz_cos_c, fz_cos_s, fz_sin_c, fz_sin_s
    real*8, dimension(n1,n2,n3) ::  d0fz_cos_c, d0fz_cos_s, d0fz_sin_c, d0fz_sin_s
    real*8, dimension(n1,n2,n3) ::  d1fz_cos_c, d1fz_cos_s, d1fz_sin_c, d1fz_sin_s
    real*8, dimension(n1,n2,n3) ::  d2fz_cos_c, d2fz_cos_s, d2fz_sin_c, d2fz_sin_s
    real*8, dimension(n3)       :: expected, error_max, sin_s, cos_s, sin_c, cos_c, f_square_s, f_square_c, mult, zero

    real*8, dimension(n3)  :: zc, zs
    integer :: nb
    character(10)   :: tmp_str
    zero=0.d0
    nb=0

    do k = 1, n3
        zc(k)=((k-1.d0)/(n3-1.d0))
        zs(k)=((k-0.5d0)/(n3-1.d0))
    end do


    do k = 1, n3
        cos_c(k)=dcos(2.d0*pi*zc(k))
        cos_s(k)=dcos(2.d0*pi*zs(k))
        sin_c(k)=dsin(2.d0*pi*zc(k))
        sin_s(k)=dsin(2.d0*pi*zs(k))

        mult(k)=0.5d0*n2*(zc(k)-0.5d0)**2
    end do

    call fillz(fz_cos_c, cos_c, n1,n2,n3)
    call fillz(fz_cos_s, cos_s, n1,n2,n3)
    call fillz(fz_sin_c, sin_c, n1,n2,n3)
    call fillz(fz_sin_s, sin_s, n1,n2,n3)

    ! D0s   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fz_cos_s=1000.d0
    call config%D0s_3Dz(fz_cos_c, d0fz_cos_s, n1,n1,n2,n2,n3, periodic)
    call diff_z(d0fz_cos_s, cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0S_SIMPLE_PERIODIC")
    nb=nb+1
    write(tmp_str, "(i10)")nb
    ! PERIODIC, shift
    d0fz_cos_c=1000.d0
    call config%D0ssh_3Dz(fz_cos_s, d0fz_cos_c, n1,n1,n2,n2,n3, periodic)
    call diff_z(d0fz_cos_c, cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fz_cos_s=1000.d0
    call config%D0s_3Dz(fz_cos_c, d0fz_cos_s, n1,n1,n2,n2,n3, Dirichlet)
    call diff_z(d0fz_cos_s, cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0S_SIMPLE_DIRICHLET")
    ! DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d0fz_cos_c=1000.d0
    call config%D0ssh_3Dz(fz_cos_s, d0fz_cos_c, n1,n1,n2,n2,n3, Dirichlet)
    call diff_z(d0fz_cos_c, cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 2, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_SIMPLE_DIRICHLET")


    ! MULT, PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillz(d0fz_cos_c, cos_c, n1,n2,n3)
    call config%D0ssh_MULT_3Dz(fz_cos_s, d0fz_cos_c, n1,n1,n2,n2,n3, periodic)
    call diff_z(d0fz_cos_c, cos_c**2, zc, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_MULT_PERIODIC")
    ! MULT, DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillz(d0fz_cos_c, cos_c, n1,n2,n3)
    call config%D0ssh_MULT_3Dz(fz_cos_s, d0fz_cos_c, n1,n1,n2,n2,n3, Dirichlet)
    call diff_z(d0fz_cos_c, cos_c**2, zc, n1, 1, n1, n2, 1, n2, n3, 2, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D0SSH_MULT_DIRICHLET")


    ! D1s   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_sin_s=1000.d0
    call config%D1s_3Dz(fz_sin_c, d1fz_sin_s, n1,n1,n2,n2,n3, h3, .false., periodic)
    call diff_z(d1fz_sin_s, cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_PERIODIC")
    ! PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_sin_c=1000.d0
    call config%D1ssh_3Dz(fz_sin_s, d1fz_sin_c, n1,n1,n2,n2,n3, h3, .false., periodic)
    call diff_z(d1fz_sin_c, cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_sin_s=1000.d0
    call config%D1s_3Dz(fz_sin_c, d1fz_sin_s, n1,n1,n2,n2,n3, h3, .false., Dirichlet)
    call diff_z(d1fz_sin_s, cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_SIMPLE_DIRICHLET")
    ! DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_sin_c=1000.d0
    call config%D1ssh_3Dz(fz_sin_s, d1fz_sin_c, n1,n1,n2,n2,n3, h3, .false., Dirichlet)
    call diff_z(d1fz_sin_c, cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 2, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_SIMPLE_DIRICHLET")


    ! ACC, PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillz(d1fz_sin_s, cos_s, n1,n2,n3)
    call config%D1s_ACC_3Dz(fz_sin_c, d1fz_sin_s, n1,n1,n2,n2,n3, h3, .false., periodic)
    call diff_z(d1fz_sin_s, 2.d0*cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_ACC_PERIODIC")
    ! ACC, PERIODIC, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillz(d1fz_sin_c, cos_c, n1,n2,n3)
    call config%D1ssh_ACC_3Dz(fz_sin_s, d1fz_sin_c, n1,n1,n2,n2,n3, h3, .false., periodic)
    call diff_z(d1fz_sin_c, 2.d0*cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_ACC_PERIODIC")

    ! ACC, DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillz(d1fz_sin_s, cos_s, n1,n2,n3)
    call config%D1s_ACC_3Dz(fz_sin_c, d1fz_sin_s, n1,n1,n2,n2,n3, h3, .false., Dirichlet)
    call diff_z(d1fz_sin_s, 2.d0*cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1S_ACC_DIRICHLET")
    ! ACC, DIRICHLET, shift
    nb=nb+1
    write(tmp_str, "(i10)")nb
    call fillz(d1fz_sin_c, cos_c, n1,n2,n3)
    call config%D1ssh_ACC_3Dz(fz_sin_s, d1fz_sin_c, n1,n1,n2,n2,n3, h3, .false., Dirichlet)
    call diff_z(d1fz_sin_c, 2.d0*cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 2, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1SSH_ACC_DIRICHLET")

    ! D1c   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_sin_c=1000.d0
    call config%D1c_3Dz(fz_sin_c, d1fz_sin_c, n1,n1,n2,n2,n3, h3, .false., periodic)
    call diff_z(d1fz_sin_c, cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_sin_c=1000.d0
    call config%D1c_3Dz(fz_sin_c, d1fz_sin_c, n1,n1,n2,n2,n3, h3, .false., Dirichlet)
    call diff_z(d1fz_sin_c, cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 2, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1C_SIMPLE_DIRICHLET")
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_sin_s=1000.d0
    call config%D1c_3Dz(fz_sin_s, d1fz_sin_s, n1,n1,n2,n2,n3, h3, .true., Dirichlet)
    call diff_z(d1fz_sin_s, cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 2, n3-2, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D1CSH_SIMPLE_DIRICHLET")

    ! D2c   *************************************************************************************************
    ! PERIODIC
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fz_cos_c=1000.d0
    call config%D2c_3Dz(fz_cos_c, d2fz_cos_c, n1,n1,n2,n2,n3, h3, .false., periodic)
    call diff_z(d2fz_cos_c, -1.d0*cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 1, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_SIMPLE_PERIODIC")

    ! DIRICHLET
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d2fz_cos_c=1000.d0
    call config%D2c_3Dz(fz_cos_c, d2fz_cos_c, n1,n1,n2,n2,n3, h3, .false., Dirichlet)
    call diff_z(d2fz_cos_c, -1.d0*cos_c, zc, n1, 1, n1, n2, 1, n2, n3, 2, n3-1, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2C_SIMPLE_DIRICHLET")
    ! DIRICHLET, shifted
    nb=nb+1
    write(tmp_str, "(i10)")nb
    d1fz_cos_s=1000.d0
    call config%D2c_3Dz(fz_cos_s, d1fz_cos_s, n1,n1,n2,n2,n3, h3, .true., Dirichlet)
    call diff_z(d1fz_cos_s, -1.d0*cos_s, zs, n1, 1, n1, n2, 1, n2, n3, 2, n3-2, "Z/"//trim(config%name)//"/"//trim(adjustl(tmp_str))//"_D2CSH_SIMPLE_DIRICHLET")




end subroutine test_3Dz
