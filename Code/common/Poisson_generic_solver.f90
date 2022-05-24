

module derivX

use decomp_2d, only : mytype

  real(mytype) :: alcaix6,acix6,bcix6
  real(mytype) :: ailcaix6,aicix6,bicix6,cicix6
  real(mytype) :: alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx
  real(mytype) :: cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,alsa1x,as1x,bs1x
  real(mytype) :: cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx
  real(mytype) :: asmx,alsaix,asix,bsix,csix,alsa3x,as3x,bs3x,alsatx,astx,bstx
end module derivX

module derivY

use decomp_2d, only : mytype

  real(mytype) :: alcaiy6,aciy6,bciy6
  real(mytype) :: ailcaiy6,aiciy6,biciy6,ciciy6
  real(mytype) :: alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny
  real(mytype) :: cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,alsa1y,as1y,bs1y
  real(mytype) :: cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy
  real(mytype) :: asmy,alsajy,asjy,bsjy,csjy,alsa3y,as3y,bs3y,alsaty,asty,bsty
end module derivY

module derivZ

use decomp_2d, only : mytype

  real(mytype) :: alcaiz6,aciz6,bciz6
  real(mytype) :: ailcaiz6,aiciz6,biciz6,ciciz6
  real(mytype) :: alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz
  real(mytype) :: cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,alsa1z,as1z,bs1z
  real(mytype) :: cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz
  real(mytype) :: asmz,alsakz,askz,bskz,cskz,alsa3z,as3z,bs3z,alsatz,astz,bstz
end module derivZ

module poisson_matrix_inversion
    implicit none

contains


    subroutine inversion5_v1(aaa,eee,spI, ny)

        USE decomp_2d
        USE MPI

        implicit none

        ! decomposition object for spectral space
        TYPE(DECOMP_INFO) :: spI
        integer             :: ny

        real(mytype), parameter :: epsilon = 1.e-16

        complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa
        complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
        integer :: i,j,k,m,mi,jc
        integer,dimension(2) :: ja,jb
        complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
        complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

        real(mytype) :: tmp1,tmp2,tmp3,tmp4

        do i=1,2
            ja(i)=4-i
            jb(i)=5-i
        enddo
        do m=1,ny/2-2
            do i=1,2
                mi=m+i
                do k=spI%yst(3),spI%yen(3)
                    do j=spI%yst(1),spI%yen(1)
                        if (real(aaa(j,m,k,3), kind=mytype).ne.0.) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
                        if (aimag(aaa(j,m,k,3)).ne.0.)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
                        sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                        eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
                        aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
                    enddo
                enddo
                do jc=ja(i),jb(i)
                    do k=spI%yst(3),spI%yen(3)
                        do j=spI%yst(1),spI%yen(1)
                            aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
                            aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
                        enddo
                    enddo
                enddo
            enddo
        enddo


        do k=spI%yst(3),spI%yen(3)
            do j=spI%yst(1),spI%yen(1)
                if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
                    tmp1=real(aaa(j,ny/2,k,2), kind=mytype)/real(aaa(j,ny/2-1,k,3), kind=mytype)
                else
                    tmp1=0.
                endif
                if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
                    tmp2=aimag(aaa(j,ny/2,k,2))/aimag(aaa(j,ny/2-1,k,3))
                else
                    tmp2=0.
                endif
                sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                b1(j,k)=cmplx(real(aaa(j,ny/2,k,3), kind=mytype)-tmp1*real(aaa(j,ny/2-1,k,4), kind=mytype),&
                aimag(aaa(j,ny/2,k,3))-tmp2*aimag(aaa(j,ny/2-1,k,4)), kind=mytype)

                if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
                    tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
                    tmp3=real(eee(j,ny/2,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny/2-1,k), kind=mytype)
                else
                    tmp1=0.
                    tmp3=0.
                endif
                if (abs(aimag(b1(j,k))).gt.epsilon) then
                    tmp2=aimag(sr(j,k))/aimag(b1(j,k))
                    tmp4=aimag(eee(j,ny/2,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny/2-1,k))
                else
                    tmp2=0.
                    tmp4=0.
                endif
                a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                eee(j,ny/2,k)=cmplx(tmp3,tmp4, kind=mytype)

                if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
                    tmp1=1./real(aaa(j,ny/2-1,k,3), kind=mytype)
                else
                    tmp1=0.
                endif
                if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
                    tmp2=1./aimag(aaa(j,ny/2-1,k,3))
                else
                    tmp2=0.
                endif
                b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                a1(j,k)=cmplx(real(aaa(j,ny/2-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
                aimag(aaa(j,ny/2-1,k,4))*aimag(b1(j,k)), kind=mytype)
                eee(j,ny/2-1,k)=cmplx(real(eee(j,ny/2-1,k))*real(b1(j,k))-real(a1(j,k))*real(eee(j,ny/2,k)),&
                aimag(eee(j,ny/2-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny/2,k)), kind=mytype)
            enddo
        enddo

        do i=ny/2-2,1,-1
            do k=spI%yst(3),spI%yen(3)
                do j=spI%yst(1),spI%yen(1)
                    if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
                        tmp1=1./real(aaa(j,i,k,3), kind=mytype)
                    else
                        tmp1=0.
                    endif
                    if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
                        tmp2=1./aimag(aaa(j,i,k,3))
                    else
                        tmp2=0.
                    endif
                    sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                    a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
                    aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
                    b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
                    aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
                    eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
                    real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
                    real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
                    aimag(eee(j,i,k))*aimag(sr(j,k))-&
                    aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
                enddo
            enddo
        enddo

        return

    end subroutine inversion5_v1


    subroutine inversion5_v2(aaa,eee,spI, nym)

        USE decomp_2d
        USE MPI

        implicit none

        ! decomposition object for spectral space
        TYPE(DECOMP_INFO) :: spI
        integer             :: nym

        real(mytype), parameter :: epsilon = 1.e-16

        complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
        complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
        integer :: i,j,k,m,mi,jc
        integer,dimension(2) :: ja,jb
        complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
        complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

        real(mytype) :: tmp1,tmp2,tmp3,tmp4

        do i=1,2
            ja(i)=4-i
            jb(i)=5-i
        enddo
        do m=1,nym-2
            do i=1,2
                mi=m+i
                do k=spI%yst(3),spI%yen(3)
                    do j=spI%yst(1),spI%yen(1)
                        if (real(aaa(j,m,k,3), kind=mytype).ne.0.) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
                        if (aimag(aaa(j,m,k,3)).ne.0.)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
                        sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                        eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
                        aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
                    enddo
                enddo
                do jc=ja(i),jb(i)
                    do k=spI%yst(3),spI%yen(3)
                        do j=spI%yst(1),spI%yen(1)
                            aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
                            aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        do k=spI%yst(3),spI%yen(3)
            do j=spI%yst(1),spI%yen(1)
                if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
                    tmp1=real(aaa(j,nym,k,2), kind=mytype)/real(aaa(j,nym-1,k,3), kind=mytype)
                else
                    tmp1=0.
                endif
                if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
                    tmp2=aimag(aaa(j,nym,k,2))/aimag(aaa(j,nym-1,k,3))
                else
                    tmp2=0.
                endif
                sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                b1(j,k)=cmplx(real(aaa(j,nym,k,3), kind=mytype)-tmp1*real(aaa(j,nym-1,k,4), kind=mytype),&
                aimag(aaa(j,nym,k,3))-tmp2*aimag(aaa(j,nym-1,k,4)), kind=mytype)
                if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
                    tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
                    tmp3=real(eee(j,nym,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,nym-1,k), kind=mytype)
                else
                    tmp1=0.
                    tmp3=0.
                endif
                if (abs(aimag(b1(j,k))).gt.epsilon) then
                    tmp2=aimag(sr(j,k))/aimag(b1(j,k))
                    tmp4=aimag(eee(j,nym,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,nym-1,k))
                else
                    tmp2=0.
                    tmp4=0.
                endif
                a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                eee(j,nym,k)=cmplx(tmp3,tmp4, kind=mytype)

                if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
                    tmp1=1./real(aaa(j,nym-1,k,3), kind=mytype)
                else
                    tmp1=0.
                endif
                if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
                    tmp2=1./aimag(aaa(j,nym-1,k,3))
                else
                    tmp2=0.
                endif
                b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                a1(j,k)=cmplx(real(aaa(j,nym-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
                aimag(aaa(j,nym-1,k,4))*aimag(b1(j,k)), kind=mytype)
                eee(j,nym-1,k)=cmplx(real(eee(j,nym-1,k), kind=mytype)*real(b1(j,k), kind=mytype)-&
                real(a1(j,k), kind=mytype)*real(eee(j,nym,k), kind=mytype),&
                aimag(eee(j,nym-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,nym,k)), kind=mytype)
            enddo
        enddo

        do i=nym-2,1,-1
            do k=spI%yst(3),spI%yen(3)
                do j=spI%yst(1),spI%yen(1)
                    if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
                        tmp1=1./real(aaa(j,i,k,3), kind=mytype)
                    else
                        tmp1=0.
                    endif
                    if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
                        tmp2=1./aimag(aaa(j,i,k,3))
                    else
                        tmp2=0.
                    endif
                    sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
                    a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
                    aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
                    b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
                    aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
                    eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
                    real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
                    real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
                    aimag(eee(j,i,k))*aimag(sr(j,k))-&
                    aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
                enddo
            enddo
        enddo

        return

    end subroutine inversion5_v2

end module poisson_matrix_inversion

module poisson_generic_solver

    use decomp_2d
    use decomp_2d_fft
    use poisson_matrix_inversion

    use mathematical_constants
    use schemes_interface

    implicit none
    integer nx, ny, nz
    integer nxm, nym, nzm
    real*8              :: xlx, yly, zlz
    real*8  :: dx3, dx2, dx1
    integer ::  nclx, ncly, nclz
    integer :: bcx, bcy, bcz


    integer     :: strectch_mode
    real*8      :: alpha, beta
    integer     :: arrang

    integer, parameter  :: PARTIAL_STAGGERED=0, FULL_STAGGERED=1

    private        ! Make everything private unless declared public

    !  real(mytype), private, parameter :: PI = 3.14159265358979323846_mytype

    real(mytype), parameter :: epsilon = 1.e-16

    ! decomposition object for physical space
    TYPE(DECOMP_INFO), save :: ph

    ! decomposition object for spectral space
    TYPE(DECOMP_INFO), save :: sp


    complex(mytype), save, allocatable, dimension(:,:,:) :: kxyz

    ! store sine/cosine factors
    real(mytype), save, allocatable, dimension(:) :: az,bz
    real(mytype), save, allocatable, dimension(:) :: ay,by
    real(mytype), save, allocatable, dimension(:) :: ax,bx


    complex(mytype), save, allocatable, dimension(:) :: zkz,zk2,ezs
    complex(mytype), save, allocatable, dimension(:) :: yky,yk2,eys
    complex(mytype), save, allocatable, dimension(:) :: xkx,xk2,exs

    !wave numbers for stretching in a pentadiagonal matrice
    complex(mytype), save, allocatable, dimension(:,:,:,:) :: a,a2,a3
    ! work arrays,
    ! naming convention: cw (complex); rw (real);
    !                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
    real(mytype), allocatable, dimension(:,:,:) :: rw1,rw1b,rw2,rw2b,rw3
    complex(mytype), allocatable, dimension(:,:,:) :: cw1,cw1b,cw2,cw22,cw2b,cw2c

    ! underlying FFT library only needs to be initialised once
    logical, save :: fft_initialised = .false.

    public :: decomp_2d_poisson_stg, decomp_2d_poisson_init, &
    decomp_2d_poisson_finalize, generic_poisson_infos

    ! For staggered mesh where main variables are defined in the centre of
    ! control volumes while boundary conditions are defined on interfaces
    interface decomp_2d_poisson_stg
        module procedure poisson
    end interface
contains



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise Poisson solver for given boundary conditions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine decomp_2d_poisson_init(nx1, ny1, nz1, nclx1, ncly1, nclz1, xlx1, yly1, zlz1, istret, alpha1, beta1, arrang1)

        implicit none

        integer, intent(IN) :: nclx1, ncly1, nclz1
        integer, intent(IN) :: nx1, ny1, nz1
        integer, intent(IN) :: istret, arrang1
        real*8              :: xlx1, yly1, zlz1
        real*8              :: alpha1, beta1
        integer :: i


        alpha=alpha1
        beta=beta1
        strectch_mode=istret
        arrang=arrang1

        nclx=nclx1
        ncly=ncly1
        nclz=nclz1

        if (nclx==0) then
            bcx=0
        else
            bcx=1
        endif
        if (ncly==0) then
            bcy=0
        else
            bcy=1
        endif
        if (nclz==0) then
            bcz=0
        else
            bcz=1
        endif

        xlx=xlx1
        yly=yly1
        zlz=zlz1

        nx=nx1
        ny=ny1
        nz=nz1

        nxm=nx
        nym=ny
        nzm=nz

        if (bcx==1) nxm=nx-1
        if (bcy==1) nym=ny-1
        if (bcz==1) nzm=nz-1

        dx3=xlx/nxm
        dx2=yly/nym
        dx1=zlz/nzm


        allocate(zkz(nz/2+1),zk2(nz/2+1),ezs(nz/2+1))
        allocate(yky(ny),yk2(ny),eys(ny))
        allocate(xkx(nx),xk2(nx),exs(nx))

        ! pressure-grid having 1 fewer point for non-periodic directions
        allocate(ax(nxm),bx(nxm))
        allocate(ay(nym),by(nym))
        allocate(az(nzm),bz(nzm))
        call abxyz(ax,ay,az,bx,by,bz,nxm,nym,nzm,bcx,bcy,bcz)

        call decomp_info_init(nxm, nym, nzm, ph)
        call decomp_info_init(nxm, nym, nzm/2+1, sp)


        ! allocate work space
        if (bcx==0 .and. bcy==0 .and. bcz==0) then
            allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(a(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a2(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))
        else if (bcx==1 .and. bcy==0 .and. bcz==0) then
            allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
            allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
            allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
            allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(a(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a2(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))
        else if (bcx==0 .and. bcy==1 .and. bcz==0) then
            allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
            allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
            allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(kxyz(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(a(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a2(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))
        else if (bcx==1 .and. bcy==1) then
            allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
            allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
            allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
            allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
            allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
            if (bcz==1) then
                allocate(rw3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
            end if
            allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
            allocate(a(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a2(sp%yst(1):sp%yen(1),nym/2,sp%yst(3):sp%yen(3),5))
            allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))
        end if

        call waves()

        return
    end subroutine decomp_2d_poisson_init


    subroutine generic_poisson_infos()
        implicit none
        if (nrank==0) then
            write(*,*) 'POISSON GENERIC SOLVER INFOS___________________________'
            write(*,*)'nx:',nx, 'ny:',ny, 'nz:',nz
            write(*,*)'nxm:',nxm, 'nym:',nym, 'nzm:',nzm
            write(*,*)'nclx:', nclx, 'ncly:', ncly, 'nclz:', nclz
            write(*,*)'strectch_mode', strectch_mode
            write(*,*)'alpha:', alpha, 'beta:', beta
            write(*,*)'xlx', xlx, 'yly', yly, 'zlz', zlz
            write(*,*) '_______________________________________________________'
        endif
    end subroutine generic_poisson_infos


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Release memory used by Poisson solver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine decomp_2d_poisson_finalize

        implicit none

        deallocate(ax,bx,ay,by,az,bz)

        call decomp_info_finalize(ph)
        call decomp_info_finalize(sp)

        call decomp_2d_fft_finalize
        fft_initialised = .false.

        deallocate(kxyz)

        if (bcx==0 .and. bcy==0 .and. bcz==0) then
            deallocate(cw1)
        else if (bcx==1 .and. bcy==0 .and. bcz==0) then
            deallocate(cw1,cw1b,rw1,rw1b,rw2)
        else if (bcx==0 .and. bcy==1 .and. bcz==0) then
            deallocate(cw1,cw2,cw2b,rw2,rw2b)
        else if (bcx==1 .and. bcy==1) then
            deallocate(cw1,cw1b,cw2,cw2b,rw1,rw1b,rw2,rw2b)
            if (bcz==1) then
                deallocate(rw3)
            end if
        end if

        return
    end subroutine decomp_2d_poisson_finalize


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Top level wrapper
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine poisson(rhs)

        implicit none

        real(mytype), dimension(:,:,:), intent(INOUT) :: rhs
        integer :: i

        if (bcx==0 .and. bcy==0 .and. bcz==0) then
            call poisson_000_v2(rhs)
        else if (bcx==1 .and. bcy==0 .and. bcz==0) then
            call poisson_100(rhs)
        else if (bcx==0 .and. bcy==1 .and. bcz==0) then
            call poisson_010(rhs)
        else if (bcx==1 .and. bcy==1) then   ! 110 & 111
            call poisson_11x(rhs, bcz)
        else
            stop 'boundary condition not supported'
        end if

        return
    end subroutine poisson


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solving 3D Poisson equation with periodic B.C in all 3 dimensions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine poisson_000(rhs)

        ! right-hand-side of Poisson as input
        ! solution of Poisson as output
        real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

        integer, dimension(3) :: fft_start, fft_end, fft_size

        complex(mytype) :: xyzk

        complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
        complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2


        real(mytype) :: tmp1, tmp2,x ,y, z

        integer :: i,j,k

        if (.not. fft_initialised) then
            call decomp_2d_fft_init(PHYSICAL_IN_Z,nxm,nym,nzm)
            fft_initialised = .true.
        end if

        ! compute r2c transform

        write(*,*) 'proc', nrank, 'A-poisson_000:', sum(abs(rhs))/(sp%xsz(1)*sp%xsz(2)*sp%xsz(3))

        open(unit=40, file="poisson_000_k_A.csv")

        do k = 1, nzm
            write(40,*)k, rhs(2,2,k)
        end do

        close(40)

        call decomp_2d_fft_3d(rhs,cw1)




        open(unit=40, file="poisson_000_k_B.csv")

        do k = sp%xst(3), sp%xen(3)
            write(40,*)k, cw1(2,2,k)
        end do

        close(40)




        ! normalisation
        cw1 = cw1 / real(nxm, kind=mytype) /real(nym, kind=mytype) &
        / real(nzm, kind=mytype)

        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)

                    ! Solve Poisson
                    tmp1=real(kxyz(i,j,k), kind=mytype)
                    tmp2=aimag(kxyz(i,j,k))
                    ! CANNOT DO A DIVISION BY ZERO
                    if ((tmp1.lt.epsilon).or.(tmp2.lt.epsilon)) then
                        cw1(i,j,k)=0._mytype
                    !                print *,'DIV 0',i,j,k,epsilon
                    else
                        cw1(i,j,k)=cmplx( real(cw1(i,j,k), kind=mytype) / (-tmp1), &
                        aimag(cw1(i,j,k))/(-tmp2), kind=mytype)
                    end if

                    !Print result in spectal space after Poisson
                    !     if (abs(out(i,j,k)) > 1.0e-4) then
                    !        write(*,*) 'AFTER',i,j,k,out(i,j,k),xyzk
                    !     end if


                end do
            end do
        end do

        ! compute c2r transform
        call decomp_2d_fft_3d(cw1,rhs)

        write(*,*) 'proc', nrank, 'B-poisson_000:', sum(abs(rhs))/(sp%xsz(1)*sp%xsz(2)*sp%xsz(3))
        call decomp_2d_fft_3d(rhs,cw1)

        !   call decomp_2d_fft_finalize

        return
    end subroutine poisson_000

    subroutine poisson_000_v2(rhs)

        ! right-hand-side of Poisson as input
        ! solution of Poisson as output
        real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

        integer, dimension(3) :: fft_start, fft_end, fft_size

        complex(mytype) :: xyzk

        complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
        complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2


        real(mytype) :: tmp1, tmp2,x ,y, z

        integer :: i,j,k

        if (.not. fft_initialised) then
            write(*,*)'TEST decomp_2d_fft_init', nxm,nym,nzm, 64,64,128
            call decomp_2d_fft_init(PHYSICAL_IN_Z,nxm,nym,nzm)
            fft_initialised = .true.
        end if

        ! compute r2c transform
        call decomp_2d_fft_3d(rhs,cw1)

        ! normalisation
        cw1 = cw1 / real(nxm, kind=mytype) /real(nym, kind=mytype) &
        / real(nzm, kind=mytype)

        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)

                    ! Solve Poisson
                    tmp1=real(kxyz(i,j,k), kind=mytype)
                    tmp2=aimag(kxyz(i,j,k))
                    ! CANNOT DO A DIVISION BY ZERO
                    if ((tmp1.lt.epsilon).or.(tmp2.lt.epsilon)) then
                        cw1(i,j,k)=0._mytype
                    !                print *,'DIV 0',i,j,k,epsilon
                    else
                        cw1(i,j,k)=cmplx( real(cw1(i,j,k), kind=mytype) / (-tmp1), &
                        aimag(cw1(i,j,k))/(-tmp2), kind=mytype)
                    end if

                    !Print result in spectal space after Poisson
                    !     if (abs(out(i,j,k)) > 1.0e-4) then
                    !        write(*,*) 'AFTER',i,j,k,out(i,j,k),xyzk
                    !     end if

                end do
            end do
        end do

        ! compute c2r transform
        call decomp_2d_fft_3d(cw1,rhs)

        !   call decomp_2d_fft_finalize

        return
    end subroutine poisson_000_v2


    subroutine poisson_100(rhs)

        implicit none

        real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

        complex(mytype) :: xyzk
        real(mytype) :: tmp1, tmp2, tmp3, tmp4
        real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

        integer :: i,j,k, itmp

100     format(1x,a8,3I4,2F12.6)

        ! rhs is in Z-pencil but requires global operations in X


        !  call decomp_2d_fft_finalize

        return
    end subroutine poisson_100


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solving 3D Poisson equation: Neumann in Y; periodic in X & Z
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine poisson_010(rhs)

        implicit none

        real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

        complex(mytype) :: xyzk
        real(mytype) :: tmp1, tmp2, tmp3, tmp4
        real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

        integer :: i,j,k

100     format(1x,a8,3I4,2F12.6)

        ! rhs is in Z-pencil but requires global operations in Y
        call transpose_z_to_y(rhs,rw2,ph)
        do k=ph%yst(3),ph%yen(3)
            do i=ph%yst(1),ph%yen(1)
                do j=1,nym/2
                    rw2b(i,j,k)=rw2(i,2*(j-1)+1,k)
                enddo
                do j=nym/2+1,nym
                    rw2b(i,j,k)=rw2(i,2*nym-2*j+2,k)
                enddo
            enddo
        end do
        call transpose_y_to_z(rw2b,rhs,ph)

        if (.not. fft_initialised) then
            call decomp_2d_fft_init(PHYSICAL_IN_Z,nxm,nym,nzm)
            fft_initialised = .true.
        end if
        ! compute r2c transform
        call decomp_2d_fft_3d(rhs,cw1)

        ! normalisation
        cw1 = cw1 / real(nxm, kind=mytype) /real(nym, kind=mytype) &
        / real(nzm, kind=mytype)
#ifdef DEBUG
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if (abs(cw1(i,j,k)) > 1.0e-4) then
                        write(*,100) 'START',i,j,k,cw1(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        ! post-processing in spectral space

        ! POST PROCESSING IN Z
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    tmp1 = real(cw1(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1(i,j,k))
                    cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                    tmp2*bz(k)-tmp1*az(k), kind=mytype)
#ifdef DEBUG
                    if (abs(cw1(i,j,k)) > 1.0e-4) &
                    write(*,100) 'after z',i,j,k,cw1(i,j,k)
#endif
                end do
            end do
        end do

        ! POST PROCESSING IN X
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    tmp1 = real(cw1(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1(i,j,k))
                    cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                    tmp2*bx(i)-tmp1*ax(i), kind=mytype)
                    if (i.gt.(nxm/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
                    if (abs(cw1(i,j,k)) > 1.0e-4) &
                    write(*,100) 'after x',i,j,k,cw1(i,j,k)
#endif
                end do
            end do
        end do

        ! POST PROCESSING IN Y
        ! NEED TO BE IN Y PENCILS!!!!!!!!!!!!!!!
        call transpose_x_to_y(cw1,cw2,sp)

        do k = sp%yst(3), sp%yen(3)
            do i = sp%yst(1), sp%yen(1)
                cw2b(i,1,k)=cw2(i,1,k)
                do j = 2,nym
                    tmp1 = real(cw2(i,j,k), kind=mytype)
                    tmp2 = aimag(cw2(i,j,k))
                    tmp3 = real(cw2(i,nym-j+2,k), kind=mytype)
                    tmp4 = aimag(cw2(i,nym-j+2,k))
                    xx1=tmp1*by(j)/2._mytype
                    xx2=tmp1*ay(j)/2._mytype
                    xx3=tmp2*by(j)/2._mytype
                    xx4=tmp2*ay(j)/2._mytype
                    xx5=tmp3*by(j)/2._mytype
                    xx6=tmp3*ay(j)/2._mytype
                    xx7=tmp4*by(j)/2._mytype
                    xx8=tmp4*ay(j)/2._mytype
                    cw2b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                    kind=mytype)
                end do
            end do
        end do
#ifdef DEBUG
        do k = sp%yst(3), sp%yen(3)
            do j = sp%yst(2), sp%yen(2)
                do i = sp%yst(1), sp%yen(1)
                    if (abs(cw2b(i,j,k)) > 1.0e-4) then
                        write(*,100) 'after y',i,j,k,cw2b(i,j,k)
                        print *,kxyz(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        if (strectch_mode==0) then

            ! Solve Poisson
            ! doing wave number division in Y-pencil
            do k = sp%yst(3), sp%yen(3)
                do j = sp%yst(2), sp%yen(2)
                    do i = sp%yst(1), sp%yen(1)
                        !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
                        !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
                        tmp1=real(kxyz(i,j,k), kind=mytype)
                        tmp2=aimag(kxyz(i,j,k))
                        !xyzk=cmplx(tmp1,tmp2, kind=mytype)
                        !CANNOT DO A DIVISION BY ZERO
                        if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then
                            cw2b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
                        end if
                        if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                            cw2b(i,j,k)=cmplx(0._mytype, &
                            aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
                        end if
                        if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then
                            cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                            /(-tmp1), 0._mytype, kind=mytype)
                        end if
                        if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                            cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                            /(-tmp1), &
                            aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
                        end if
                    end do
                end do
            end do

        else

            call matrice_refinement()
            !       do k = sp%yst(3), sp%yen(3)
            !          do j = 1,nym/2
            !             do i = sp%yst(1), sp%yen(1)
            !                print *,i,j,k,a(i,j,k,3)
            !!                if (nrank.le.1) print *,i,j,k,a(i,j,k,3)
            !!                if (nrank.gt.1) print *,i+4,j,k,a(i,j,k,3)
            !             enddo
            !          enddo
            !       enddo


            if (strectch_mode.ne.3) then
                cw2(:,:,:)=0.;cw2c(:,:,:)=0.
                do k = sp%yst(3), sp%yen(3)
                    do j = 1,nym/2
                        do i = sp%yst(1), sp%yen(1)
                            cw2(i,j,k)=cw2b(i,2*j-1,k)
                            cw2c(i,j,k)=cw2b(i,2*j,k)
                        enddo
                    enddo
                enddo

                !   do k = sp%yst(3), sp%yen(3)
                !      do j = 1,nym/2
                !         do i = sp%yst(1), sp%yen(1)
                !            if (abs(cw2(i,j,k)) > 1.0e-4) then
                !               write(*,*) 'before IN',i,j,k,cw2(i,j,k)!*2.
                !!            end if
                !        end do
                !     end do
                !  end do

                call inversion5_v1(a,cw2,sp, nym)
                call inversion5_v1(a2,cw2c,sp, nym)

                !         cw2(1,1,1)=cw2(1,1,1)*0.5


                !   do k = sp%yst(3), sp%yen(3)
                !       do j = 1,nym/2
                !          do i = sp%yst(1), sp%yen(1)
                !             if (abs(cw2c(i,j,k)) > 1.0e-4) then
                !                write(*,*) 'after IN',i,j,k,cw2c(i,j,k)!*2.
                !             end if
                !          end do
                !       end do
                !    end do

                cw2b(:,:,:)=0.
                do k=sp%yst(3), sp%yen(3)
                    do j=1,nym-1,2
                        do i=sp%yst(1), sp%yen(1)
                            cw2b(i,j,k)=cw2(i,(j+1)/2,k)
                        enddo
                    enddo
                    do j=2,nym,2
                        do i=sp%yst(1), sp%yen(1)
                            cw2b(i,j,k)=cw2c(i,j/2,k)
                        enddo
                    enddo
                enddo
               !do k=sp%yst(3), sp%yen(3)
               !do i=sp%yst(1), sp%yen(1)
               !   if ((xkx(i)==0).and.(zkz(k)==0)) then
               !   !   cw2b(i,1,1)=0.
               !   !   cw2b(i,nym,1)=0.
               !   endif
               !enddo
               !enddo
            else
                do k = sp%yst(3), sp%yen(3)
                    do j = 1,nym
                        do i = sp%yst(1), sp%yen(1)
                            cw2(i,j,k)=cw2b(i,j,k)
                        enddo
                    enddo
                enddo
                call inversion5_v2(a3,cw2,sp, nym)
                do k = sp%yst(3), sp%yen(3)
                    do j = 1,nym
                        do i = sp%yst(1), sp%yen(1)
                            cw2b(i,j,k)=cw2(i,j,k)
                        enddo
                    enddo
                enddo
            endif

        endif

        !    print *,nrank, sp%yst(3),sp%yen(3),sp%yst(1),sp%yen(1)


#ifdef DEBUG
        do k = sp%yst(3), sp%yen(3)
            do j = sp%yst(2), sp%yen(2)
                do i = sp%yst(1), sp%yen(1)
                    if (abs(cw2b(i,j,k)) > 1.0e-4) then
                        write(*,100) 'AFTER',i,j,k,cw2b(i,j,k)
                        print *,kxyz(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        ! post-processing backward

        ! POST PROCESSING IN Y
        do k = sp%yst(3), sp%yen(3)
            do i = sp%yst(1), sp%yen(1)
                cw2(i,1,k)=cw2b(i,1,k)
                do j = 2,nym
                    tmp1 = real(cw2b(i,j,k), kind=mytype)
                    tmp2 = aimag(cw2b(i,j,k))
                    tmp3 = real(cw2b(i,nym-j+2,k), kind=mytype)
                    tmp4 = aimag(cw2b(i,nym-j+2,k))
                    xx1=tmp1*by(j)
                    xx2=tmp1*ay(j)
                    xx3=tmp2*by(j)
                    xx4=tmp2*ay(j)
                    xx5=tmp3*by(j)
                    xx6=tmp3*ay(j)
                    xx7=tmp4*by(j)
                    xx8=tmp4*ay(j)
                    cw2(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                    kind=mytype)
                end do
            end do
        end do

        ! Back to X-pencil
        call transpose_y_to_x(cw2,cw1,sp)
#ifdef DEBUG
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if (abs(cw1(i,j,k)) > 1.0e-4) then
                        write(*,100) 'AFTER Y',i,j,k,cw1(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        ! POST PROCESSING IN X
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    tmp1 = real(cw1(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1(i,j,k))
                    cw1(i,j,k) = cmplx(tmp1*bx(i)-tmp2*ax(i), &
                    tmp2*bx(i)+tmp1*ax(i), kind=mytype)
                    if (i.gt.(nxm/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
                    if (abs(cw1(i,j,k)) > 1.0e-4) &
                    write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
#endif
                end do
            end do
        end do

        ! POST PROCESSING IN Z
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    tmp1 = real(cw1(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1(i,j,k))
                    cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                    tmp2*bz(k)+tmp1*az(k), kind=mytype)
#ifdef DEBUG
                    if (abs(cw1(i,j,k)) > 1.0e-4) &
                    write(*,100) 'END',i,j,k,cw1(i,j,k)
#endif
                end do
            end do
        end do

        ! compute c2r transform, back to physical space
        call decomp_2d_fft_3d(cw1,rhs)

        ! rhs is in Z-pencil but requires global operations in Y
        call transpose_z_to_y(rhs,rw2,ph)
        do k=ph%yst(3),ph%yen(3)
            do i=ph%yst(1),ph%yen(1)
                do j=1,nym/2
                    rw2b(i,2*j-1,k)=rw2(i,j,k)
                enddo
                do j=1,nym/2
                    rw2b(i,2*j,k)=rw2(i,nym-j+1,k)
                enddo
            enddo
        end do
        call transpose_y_to_z(rw2b,rhs,ph)

        !  call decomp_2d_fft_finalize

        return
    end subroutine poisson_010
    subroutine poisson_010_v2(rhs)

        implicit none

        real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

        complex(mytype) :: xyzk
        real(mytype) :: tmp1, tmp2, tmp3, tmp4
        real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

        integer :: i,j,k

100     format(1x,a8,3I4,2F12.6)

        ! rhs is in Z-pencil but requires global operations in Y
        call transpose_z_to_y(rhs,rw2,ph)
        do k=ph%yst(3),ph%yen(3)
            do i=ph%yst(1),ph%yen(1)
                do j=1,nym/2
                    rw2b(i,j,k)=rw2(i,2*(j-1)+1,k)
                enddo
                do j=nym/2+1,nym
                    rw2b(i,j,k)=rw2(i,2*nym-2*j+2,k)
                enddo
            enddo
        end do

        call transpose_y_to_z(rw2b,rhs,ph)

        if (.not. fft_initialised) then
            call decomp_2d_fft_init(PHYSICAL_IN_Z,nxm,nym,nzm)
            fft_initialised = .true.
        end if
        ! compute r2c transform
        call decomp_2d_fft_3d(rhs,cw1)

        ! normalisation
        cw1 = cw1 / real(nxm, kind=mytype) /real(nym, kind=mytype) &
        / real(nzm, kind=mytype)


        ! post-processing in spectral space

        ! POST PROCESSING IN Y
        ! NEED TO BE IN Y PENCILS!!!!!!!!!!!!!!!
        call transpose_x_to_y(cw1,cw2,sp)

        do k = sp%yst(3), sp%yen(3)
            do i = sp%yst(1), sp%yen(1)
                cw2b(i,1,k)=cw2(i,1,k)
                do j = 2,nym
                    tmp1 = real(cw2(i,j,k), kind=mytype)
                    tmp2 = aimag(cw2(i,j,k))
                    tmp3 = real(cw2(i,nym-j+2,k), kind=mytype)
                    tmp4 = aimag(cw2(i,nym-j+2,k))
                    xx1=tmp1*by(j)/2._mytype
                    xx2=tmp1*ay(j)/2._mytype
                    xx3=tmp2*by(j)/2._mytype
                    xx4=tmp2*ay(j)/2._mytype
                    xx5=tmp3*by(j)/2._mytype
                    xx6=tmp3*ay(j)/2._mytype
                    xx7=tmp4*by(j)/2._mytype
                    xx8=tmp4*ay(j)/2._mytype
                    cw2b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                    kind=mytype)
                end do
            end do
        end do

        ! Solve Poisson
        ! doing wave number division in Y-pencil
        do k = sp%yst(3), sp%yen(3)
            do j = sp%yst(2), sp%yen(2)
                do i = sp%yst(1), sp%yen(1)
                    !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
                    !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
                    tmp1=real(kxyz(i,j,k), kind=mytype)
                    tmp2=aimag(kxyz(i,j,k))
                    !xyzk=cmplx(tmp1,tmp2, kind=mytype)
                    !CANNOT DO A DIVISION BY ZERO
                    if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then
                        cw2b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
                    end if
                    if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                        cw2b(i,j,k)=cmplx(0._mytype, &
                        aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
                    end if
                    if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then
                        cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                        /(-tmp1), 0._mytype, kind=mytype)
                    end if
                    if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                        cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                        /(-tmp1), &
                        aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
                    end if
                end do
            end do
        end do

        ! post-processing backward

        ! POST PROCESSING IN Y
        do k = sp%yst(3), sp%yen(3)
            do i = sp%yst(1), sp%yen(1)
                cw2(i,1,k)=cw2b(i,1,k)
                do j = 2,nym
                    tmp1 = real(cw2b(i,j,k), kind=mytype)
                    tmp2 = aimag(cw2b(i,j,k))
                    tmp3 = real(cw2b(i,nym-j+2,k), kind=mytype)
                    tmp4 = aimag(cw2b(i,nym-j+2,k))
                    xx1=tmp1*by(j)
                    xx2=tmp1*ay(j)
                    xx3=tmp2*by(j)
                    xx4=tmp2*ay(j)
                    xx5=tmp3*by(j)
                    xx6=tmp3*ay(j)
                    xx7=tmp4*by(j)
                    xx8=tmp4*ay(j)
                    cw2(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                    kind=mytype)
                end do
            end do
        end do

        ! Back to X-pencil
        call transpose_y_to_x(cw2,cw1,sp)

        ! compute c2r transform, back to physical space
        call decomp_2d_fft_3d(cw1,rhs)

        ! rhs is in Z-pencil but requires global operations in Y
        call transpose_z_to_y(rhs,rw2,ph)
        do k=ph%yst(3),ph%yen(3)
            do i=ph%yst(1),ph%yen(1)
                do j=1,nym/2
                    rw2b(i,2*j-1,k)=rw2(i,j,k)
                enddo
                do j=1,nym/2
                    rw2b(i,2*j,k)=rw2(i,nym-j+1,k)
                enddo
            enddo
        end do
        call transpose_y_to_z(rw2b,rhs,ph)

        !  call decomp_2d_fft_finalize

        return
    end subroutine poisson_010_v2


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solving 3D Poisson equation: Neumann in X, Y; Neumann/periodic in Z
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine poisson_11x(rhs, nclz1)

        implicit none

        integer, intent(IN) :: nclz1
        real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

        complex(mytype) :: xyzk
        real(mytype) :: tmp1, tmp2, tmp3, tmp4
        real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

        integer :: i,j,k

100     format(1x,a8,3I4,2F12.6)

        if (nclz1==1) then
            do j=1,ph%zsz(2)
                do i=1,ph%zsz(1)
                    do k=1,nzm/2
                        rw3(i,j,k)=rhs(i,j,2*(k-1)+1)
                    end do
                    do k=nzm/2+1,nzm
                        rw3(i,j,k)=rhs(i,j,2*nzm-2*k+2)
                    end do
                end do
            end do
            call transpose_z_to_y(rw3,rw2,ph)
        else if (nclz1==0) then
            call transpose_z_to_y(rhs,rw2,ph)
        end if


        do k=ph%yst(3),ph%yen(3)
            do i=ph%yst(1),ph%yen(1)
                do j=1,nym/2
                    rw2b(i,j,k)=rw2(i,2*(j-1)+1,k)
                end do
                do j=nym/2+1,nym
                    rw2b(i,j,k)=rw2(i,2*nym-2*j+2,k)
                end do
            end do
        end do

        ! the global operations in X
        call transpose_y_to_x(rw2b,rw1,ph)

        do k=ph%xst(3),ph%xen(3)
            do j=ph%xst(2),ph%xen(2)
                do i=1,nxm/2
                    rw1b(i,j,k)=rw1(2*(i-1)+1,j,k)
                end do
                do i=nxm/2+1,nxm
                    rw1b(i,j,k)=rw1(2*nxm-2*i+2,j,k)
                end do
            end do
        end do

        ! back to Z-pencil
        call transpose_x_to_y(rw1b,rw2,ph)
        call transpose_y_to_z(rw2,rhs,ph)

        if (.not. fft_initialised) then
            call decomp_2d_fft_init(PHYSICAL_IN_Z,nxm,nym,nzm)
            fft_initialised = .true.
        end if

        ! compute r2c transform

        call decomp_2d_fft_3d(rhs,cw1)



        ! normalisation
        cw1 = cw1 / real(nxm, kind=mytype) /real(nym, kind=mytype) &
        / real(nzm, kind=mytype)

        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if ((i==1).and.(j==1).and.(k==1)) then
                        write(*,*) 'mode 111 dans DIV', cw1(i,j,k)
                        cw1(i,j,k)=0.
                    end if
                end do
            end do
        end do
#ifdef DEBUG
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if (abs(cw1(i,j,k)) > 1.0e-4) then
                        write(*,100) 'START',i,j,k,cw1(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        ! post-processing in spectral space

        ! POST PROCESSING IN Z
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    tmp1 = real(cw1(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1(i,j,k))
                    cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                    tmp2*bz(k)-tmp1*az(k), kind=mytype)
#ifdef DEBUG
                    if (abs(cw1(i,j,k)) > 1.0e-4) &
                    write(*,100) 'after z',i,j,k,cw1(i,j,k)
#endif
                end do
            end do
        end do

        ! POST PROCESSING IN Y
        ! WE HAVE TO BE IN Y PENCILS
        call transpose_x_to_y(cw1,cw2,sp)
        do k = sp%yst(3), sp%yen(3)
            do i = sp%yst(1), sp%yen(1)
                cw2b(i,1,k)=cw2(i,1,k)
                do j = 2,nym
                    tmp1 = real(cw2(i,j,k), kind=mytype)
                    tmp2 = aimag(cw2(i,j,k))
                    tmp3 = real(cw2(i,nym-j+2,k), kind=mytype)
                    tmp4 = aimag(cw2(i,nym-j+2,k))
                    xx1=tmp1*by(j)/2._mytype
                    xx2=tmp1*ay(j)/2._mytype
                    xx3=tmp2*by(j)/2._mytype
                    xx4=tmp2*ay(j)/2._mytype
                    xx5=tmp3*by(j)/2._mytype
                    xx6=tmp3*ay(j)/2._mytype
                    xx7=tmp4*by(j)/2._mytype
                    xx8=tmp4*ay(j)/2._mytype
                    cw2b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                    kind=mytype)
                end do
            end do
        end do

        ! back to X-pencil
        call transpose_y_to_x(cw2b,cw1,sp)
#ifdef DEBUG
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if (abs(cw1(i,j,k)) > 1.0e-4) then
                        write(*,100) 'after y',i,j,k,cw1(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        ! POST PROCESSING IN X
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                cw1b(1,j,k)=cw1(1,j,k)
                do i = 2,nxm
                    tmp1 = real(cw1(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1(i,j,k))
                    tmp3 = real(cw1(nxm-i+2,j,k), kind=mytype)
                    tmp4 = aimag(cw1(nxm-i+2,j,k))
                    xx1=tmp1*bx(i)/2._mytype
                    xx2=tmp1*ax(i)/2._mytype
                    xx3=tmp2*bx(i)/2._mytype
                    xx4=tmp2*ax(i)/2._mytype
                    xx5=tmp3*bx(i)/2._mytype
                    xx6=tmp3*ax(i)/2._mytype
                    xx7=tmp4*bx(i)/2._mytype
                    xx8=tmp4*ax(i)/2._mytype
                    cw1b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                    kind=mytype)
                end do
            end do
        end do

#ifdef DEBUG
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if (abs(cw1b(i,j,k)) > 1.0e-4) then
                        write(*,*) 'BEFORE',i,j,k,cw1b(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        if (strectch_mode==0) then

            ! Solve Poisson
            do k = sp%xst(3),sp%xen(3)
                do j = sp%xst(2),sp%xen(2)
                    do i = sp%xst(1),sp%xen(1)
                        !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
                        !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
                        tmp1=real(kxyz(i,j,k), kind=mytype)
                        tmp2=aimag(kxyz(i,j,k))
                        !xyzk=cmplx(tmp1,tmp2, kind=mytype)
                        !CANNOT DO A DIVISION BY ZERO
                        if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then
                            cw1b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
                        end if
                        if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                            cw1b(i,j,k)=cmplx(0._mytype, &
                            aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
                        end if
                        if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then
                            cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                            /(-tmp1), 0._mytype, kind=mytype)
                        end if
                        if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                            cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                            /(-tmp1), &
                            aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
                        end if
                    end do
                end do
            end do

        else
            call matrice_refinement()
            ! the stretching is only working in Y pencils
            call transpose_x_to_y(cw1b,cw2b,sp)
            !we are now in Y pencil

            if (strectch_mode.ne.3) then
                cw2(:,:,:)=0.;cw2c(:,:,:)=0.
                do k = sp%yst(3), sp%yen(3)
                    do j = 1,nym/2
                        do i = sp%yst(1), sp%yen(1)
                            cw2(i,j,k)=cw2b(i,2*j-1,k)
                            cw2c(i,j,k)=cw2b(i,2*j,k)
                        enddo
                    enddo
                enddo
                call inversion5_v1(a,cw2,sp, nym)
                call inversion5_v1(a2,cw2c,sp, nym)

                cw2b(:,:,:)=0.
                do k=sp%yst(3), sp%yen(3)
                    do j=1,nym-1,2
                        do i=sp%yst(1), sp%yen(1)
                            cw2b(i,j,k)=cw2(i,(j+1)/2,k)
                        enddo
                    enddo
                    do j=2,nym,2
                        do i=sp%yst(1), sp%yen(1)
                            cw2b(i,j,k)=cw2c(i,j/2,k)
                        enddo
                    enddo
                enddo
            else
                cw2(:,:,:)=0.
                do k = sp%yst(3), sp%yen(3)
                    do j = sp%yst(2), sp%yen(2)
                        do i = sp%yst(1), sp%yen(1)
                            cw2(i,j,k)=cw2b(i,j,k)
                        enddo
                    enddo
                enddo

                call inversion5_v2(a3,cw2,sp, nym)

                do k = sp%yst(3), sp%yen(3)
                    do j = sp%yst(2), sp%yen(2)
                        do i = sp%yst(1), sp%yen(1)
                            cw2b(i,j,k)=cw2(i,j,k)
                        enddo
                    enddo
                enddo
            endif
            !we have to go back in X pencils
            call transpose_y_to_x(cw2b,cw1b,sp)
        endif

#ifdef DEBUG
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if (abs(cw1b(i,j,k)) > 1.0e-6) then
                        write(*,*) 'AFTER',i,j,k,cw1b(i,j,k)
                    end if
                end do
            end do
        end do
#endif
        !stop
        ! post-processing backward

        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                cw1(1,j,k)=cw1b(1,j,k)
                do i = 2,nxm
                    tmp1 = real(cw1b(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1b(i,j,k))
                    tmp3 = real(cw1b(nxm-i+2,j,k), kind=mytype)
                    tmp4 = aimag(cw1b(nxm-i+2,j,k))
                    xx1=tmp1*bx(i)
                    xx2=tmp1*ax(i)
                    xx3=tmp2*bx(i)
                    xx4=tmp2*ax(i)
                    xx5=tmp3*bx(i)
                    xx6=tmp3*ax(i)
                    xx7=tmp4*bx(i)
                    xx8=tmp4*ax(i)
                    cw1(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                    kind=mytype)
                end do
            end do
        end do
#ifdef DEBUG
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if (abs(cw1(i,j,k)) > 1.0e-4) then
                        write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
                    end if
                end do
            end do
        end do
#endif

        ! POST PROCESSING IN Y
        ! NEED to be in Y-pencil
        call transpose_x_to_y(cw1,cw2,sp)
        do k = sp%yst(3), sp%yen(3)
            do i = sp%yst(1), sp%yen(1)
                cw2b(i,1,k)=cw2(i,1,k)
                do j = 2,nym
                    tmp1 = real(cw2(i,j,k), kind=mytype)
                    tmp2 = aimag(cw2(i,j,k))
                    tmp3 = real(cw2(i,nym-j+2,k), kind=mytype)
                    tmp4 = aimag(cw2(i,nym-j+2,k))
                    xx1=tmp1*by(j)
                    xx2=tmp1*ay(j)
                    xx3=tmp2*by(j)
                    xx4=tmp2*ay(j)
                    xx5=tmp3*by(j)
                    xx6=tmp3*ay(j)
                    xx7=tmp4*by(j)
                    xx8=tmp4*ay(j)
                    cw2b(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                    kind=mytype)
                end do
            end do
        end do
#ifdef DEBUG
        do k = sp%yst(3), sp%yen(3)
            do j = sp%yst(2), sp%yen(2)
                do i = sp%yst(1), sp%yen(1)
                    if (abs(cw2b(i,j,k)) > 1.0e-4) then
                        write(*,100) 'AFTER Y',i,j,k,cw2b(i,j,k)
                    end if
                end do
            end do
        end do
#endif
        ! back to X-pencil
        call transpose_y_to_x(cw2b,cw1,sp)

        ! POST PROCESSING IN Z
        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    tmp1 = real(cw1(i,j,k), kind=mytype)
                    tmp2 = aimag(cw1(i,j,k))
                    cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                    tmp2*bz(k)+tmp1*az(k), kind=mytype)
#ifdef DEBUG
                    if (abs(cw1(i,j,k)) > 1.0e-4) &
                    write(*,100) 'END',i,j,k,cw1(i,j,k)
#endif
                end do
            end do
        end do

        ! compute c2r transform, back to physical space



        do k = sp%xst(3),sp%xen(3)
            do j = sp%xst(2),sp%xen(2)
                do i = sp%xst(1),sp%xen(1)
                    if ((i==1).and.(j==1).and.(k==1)) then
                        write(*,*) 'mode 111 dans CORRECTION', cw1(i,j,k)
                    end if
                end do
            end do
        end do

        call decomp_2d_fft_3d(cw1,rhs)

        if (nclz1==1) then
            do j=1,ph%zsz(2)
                do i=1,ph%zsz(1)
                    do k=1,nzm/2
                        rw3(i,j,2*k-1)=rhs(i,j,k)
                    end do
                    do k=1,nzm/2
                        rw3(i,j,2*k)=rhs(i,j,nzm-k+1)
                    end do
                end do
            end do
            call transpose_z_to_y(rw3,rw2,ph)
        else if (nclz1==0) then
            call transpose_z_to_y(rhs,rw2,ph)
        end if

        do k=ph%yst(3),ph%yen(3)
            do i=ph%yst(1),ph%yen(1)
                do j=1,nym/2
                    rw2b(i,2*j-1,k)=rw2(i,j,k)
                end do
                do j=1,nym/2
                    rw2b(i,2*j,k)=rw2(i,nym-j+1,k)
                end do
            enddo
        end do
        call transpose_y_to_x(rw2b,rw1,ph)
        do k=ph%xst(3),ph%xen(3)
            do j=ph%xst(2),ph%xen(2)
                do i=1,nxm/2
                    rw1b(2*i-1,j,k)=rw1(i,j,k)
                enddo
                do i=1,nxm/2
                    rw1b(2*i,j,k)=rw1(nxm-i+1,j,k)
                enddo
            enddo
        end do
        call transpose_x_to_y(rw1b,rw2,ph)
        call transpose_y_to_z(rw2,rhs,ph)

        !  call decomp_2d_fft_finalize



        return
    end subroutine poisson_11x



    ! ***********************************************************
    !
    subroutine waves ()
        !
        !***********************************************************

        use d2c_from_d1s
        use d1s_schemes

        implicit none

        integer :: i,j,k
        real(mytype) :: w,wp,w1,w1p



        real*8                                  :: wp2

        xkx(:)=0. ; xk2(:)=0. ; yky(:)=0. ; yk2(:)=0.
        zkz(:)=0. ; zk2(:)=0.

        !WAVE NUMBER IN X
        if (nclx==0) then
            do i=1,nx/2+1
                w=2.d0*PI*(i-1)/nx
                wp=keq(D1s, nx+1, i-1)
                wp=2.d0*PI*wp/nx

                xkx(i)=cmplx(nx*wp/xlx,nx*wp/xlx, kind=mytype)
                exs(i)=cmplx(nx*w/xlx,nx*w/xlx, kind=mytype)
                xk2(i)=cmplx((nx*wp/xlx)**2,(nx*wp/xlx)**2, kind=mytype)
            enddo
            !call exit
            do i=nx/2+2,nx
                xkx(i)=xkx(nx-i+2)
                exs(i)=exs(nx-i+2)
                xk2(i)=xk2(nx-i+2)
            enddo
        else
            do i=1,nx
                w=2.*PI*0.5*(i-1)/nxm
                wp=keq(D1s, 2*nx-1, i-1)
                wp=2.d0*PI*0.5*wp/nxm

                xkx(i)=cmplx(nxm*wp/xlx,nxm*wp/xlx, kind=mytype)
                exs(i)=cmplx(nxm*w/xlx,nxm*w/xlx, kind=mytype)
                xk2(i)=cmplx((nxm*wp/xlx)**2,(nxm*wp/xlx)**2, kind=mytype)
            enddo

            xkx(1)=0.
            exs(1)=0.
            xk2(1)=0.
        endif

        !WAVE NUMBER IN Y
        if (ncly==0) then
            do j=1,ny/2+1
                w=2.*PI*(j-1)/ny
!                wp=keq(D1s_ExpCtr_O0Fp5, ny+1, j-1)
                wp=keq(D1s, ny+1, j-1)
                wp=2.d0*PI*wp/ny

                if (strectch_mode==0) yky(j)=cmplx(ny*wp/yly,ny*wp/yly, kind=mytype)
                if (strectch_mode.ne.0) yky(j)=cmplx(ny*wp,ny*wp, kind=mytype)
                eys(j)=cmplx(ny*w/yly,ny*w/yly, kind=mytype)
                yk2(j)=cmplx((ny*wp/yly)**2,(ny*wp/yly)**2, kind=mytype)
            enddo
            do j=ny/2+2,ny
                yky(j)=yky(ny-j+2)
                eys(j)=eys(ny-j+2)
                yk2(j)=yk2(ny-j+2)
            enddo
        else
            do j=1,ny
                w=2.*PI*0.5*(j-1)/nym
!                wp=keq(D1s_ExpCtr_O0Fp5, 2*ny-1, j-1)
                wp=keq(D1s, 2*ny-1, j-1)
                wp=2.d0*PI*0.5*wp/nym

                if (strectch_mode==0) yky(j)=cmplx(nym*wp/yly,nym*wp/yly, kind=mytype)
                if (strectch_mode.ne.0) yky(j)=cmplx(nym*wp,nym*wp, kind=mytype)
                eys(j)=cmplx(nym*w/yly,nym*w/yly, kind=mytype)
                yk2(j)=cmplx((nym*wp/yly)**2,(nym*wp/yly)**2, kind=mytype)
            enddo
            yky(1)=0.
            eys(1)=0.
            yk2(1)=0.
        endif

        !WAVE NUMBER IN Z
        if (nclz==0) then
            do k=1,nz/2+1
                w=2.*PI*(k-1)/nz
                wp=keq(D1s, nz+1, k-1)
                wp=2.d0*PI*wp/nz

                zkz(k)=cmplx(nz*wp/zlz,nz*wp/zlz, kind=mytype)
                ezs(k)=cmplx(nz*w/zlz,nz*w/zlz, kind=mytype)
                zk2(k)=cmplx((nz*wp/zlz)**2,(nz*wp/zlz)**2, kind=mytype)
            enddo
        else
            do k=1,nz/2+1
                w=2.*PI*0.5*(k-1)/nzm
                wp=keq(D1s, 2*nz-1, k-1)
                wp=2.d0*PI*0.5*wp/nzm


                w1=2.*PI*0.5*(nzm-k+1)/nzm
                w1p=keq(D1s, 2*nz-1, (nzm-k+1))
                w1p=2.d0*PI*0.5*w1p/nzm

                zkz(k)=cmplx(nzm*wp/zlz,-nzm*w1p/zlz, kind=mytype)
                ezs(k)=cmplx(nzm*w/zlz,nzm*w1/zlz, kind=mytype)
                zk2(k)=cmplx((nzm*wp/zlz)**2,(nzm*w1p/zlz)**2, kind=mytype)
            enddo
        endif
        !
        !if (nrank==0) then
        !   do i=1,nx
        !      print *,i,ezs(i)
        !   enddo
        !endif
        !stop

        call perform_kxyz


    !          do k=1,1!nz
    !          do j=1,ny
    !          do i=1,1!!nx
    !             print *,j,a(i,j,k,3),kxyz(i,j,k)
    !          enddo
    !          enddo
    !          enddo


    contains



        real(kind=8) function keq(d1s_scheme, N, k)

            implicit none
            real(kind=8) :: dx3, PI
            integer :: i, N, k
            procedure(D0s_ExpCtr_O2Fp0) :: d1s_scheme
            real(kind=8)  x(N),f1(N), f2(N), df1_calc(N), df2_calc(N), xs(N)

            PI=dacos(-1.d0)

            dx3=(2.d0*PI/(N-1))

            x(1) = 1.d0

            do i=2, N
                x(i)=x(i-1)+dx3
            enddo

            do i=1, N-1
                xs(i)=(x(i)+x(i+1))*0.5d0
            enddo

            do i=1, N
                f1(i)=dcos(k*x(i))
                f2(i)=dsin(k*x(i))
            enddo

            call d1s_scheme(f1, df1_calc, N, dx3, .false., periodic)
            call d1s_scheme(f2, df2_calc, N, dx3, .false., periodic)

            i=N/2
            keq =  (  dcos(k*xs(i))*df2_calc(i) - dsin(k*xs(i))*df1_calc(i) )


            return
        end function

    end subroutine waves

    !**************************************************************************
    !
    subroutine matrice_refinement()
        !
        !**************************************************************************
        implicit none

        integer :: i,j,k

        complex(mytype),dimension(sp%yst(1):sp%yen(1)) :: transx
        complex(mytype),dimension(sp%yst(2):sp%yen(2)) :: transy
        complex(mytype),dimension(sp%yst(3):sp%yen(3)) :: transz
        real(mytype) :: xa0,xa1
        complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
        complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3


        ! Calcul des coefficients de transfert
        call perform_TF_coeffs(transx, transy, transz)

        cw2=cmplx(0.d0, 0.d0)

        if ((strectch_mode==1).or.(strectch_mode==2)) then


            xa0=alpha/PI+1./2./beta/PI
            if (strectch_mode==1) xa1=1./4./beta/PI
            if (strectch_mode==2) xa1=-1./4./beta/PI
            !
            !construction of the pentadiagonal matrice
            !
            do k=sp%yst(3),sp%yen(3)
                do j=1,ny/2
                    do i=sp%yst(1),sp%yen(1)
                        ! Pour chaque mode
                        !cw22 <= keq*Trx*Trz
                        cw22(i,j,k)=cmplx(real(yky(2*j-1), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind=mytype),&
                        aimag(yky(2*j-1))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
                        cw2(i,j,k)=cmplx(real(yky(2*j), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind=mytype),&
                        aimag(yky(2*j))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
                    enddo
                enddo
            enddo




            !main diagonal
            ! LES CALCULS CI-DESSOUS IMPLEMENTENT LES RELATIONS 64 de la publi
            do k=sp%yst(3),sp%yen(3)
                do j=2,ny/2-1
                    do i=sp%yst(1),sp%yen(1)
                        a(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2*j-1), kind=mytype)*real(transy(2*j-1), kind=mytype)&
                        *real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                        -(real(zk2(k), kind=mytype)*real(transy(2*j-1), kind=mytype)*real(transy(2*j-1), kind=mytype)*&
                        real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
                        -real(cw22(i,j,k), kind=mytype)*real(cw22(i,j,k), kind=mytype)*xa0*xa0-&
                        xa1*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j,k), kind=mytype)*&
                        real(cw22(i,j+1,k), kind=mytype)),&
                        -(aimag(xk2(i))*aimag(transy(2*j-1))*aimag(transy(2*j-1))*aimag(transz(k))*aimag(transz(k)))&
                        -(aimag(zk2(k))*aimag(transy(2*j-1))*aimag(transy(2*j-1))*aimag(transx(i))*aimag(transx(i)))&
                        -aimag(cw22(i,j,k))*aimag(cw22(i,j,k))*xa0*xa0-&
                        xa1*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))), kind=mytype)

                        a2(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2*j), kind=mytype)*real(transy(2*j), kind=mytype)*&
                        real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                        -(real(zk2(k), kind=mytype)*real(transy(2*j), kind=mytype)*real(transy(2*j), kind=mytype)*&
                        real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
                        -real(cw2(i,j,k), kind=mytype)*real(cw2(i,j,k), kind=mytype)*xa0*xa0-&
                        xa1*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j-1,k), kind=mytype)+real(cw2(i,j,k), kind=mytype)*&
                        real(cw2(i,j+1,k), kind=mytype)),&
                        -(aimag(xk2(i))*aimag(transy(2*j))*aimag(transy(2*j))*aimag(transz(k))*aimag(transz(k)))&
                        -(aimag(zk2(k))*aimag(transy(2*j))*aimag(transy(2*j))*aimag(transx(i))*aimag(transx(i)))&
                        -aimag(cw2(i,j,k))*aimag(cw2(i,j,k))*xa0*xa0-&
                        xa1*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j-1,k))+aimag(cw2(i,j,k))*aimag(cw2(i,j+1,k))), kind=mytype)
                    enddo
                enddo
                do i=sp%yst(1),sp%yen(1)
                    a(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*&
                    real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                    -(real(zk2(k), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*&
                    real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
                    -real(cw22(i,1,k), kind=mytype)*real(cw22(i,1,k), kind=mytype)*xa0*xa0-xa1*xa1*(real(cw22(i,1,k), kind=mytype)*&
                    real(cw22(i,2,k), kind=mytype)),&
                    -(aimag(xk2(i))*aimag(transy(1))*aimag(transy(1))*aimag(transz(k))*aimag(transz(k)))&
                    -(aimag(zk2(k))*aimag(transy(1))*aimag(transy(1))*aimag(transx(i))*aimag(transx(i)))&
                    -aimag(cw22(i,1,k))*aimag(cw22(i,1,k))*xa0*xa0-xa1*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))), kind=mytype)
                    a(i,ny/2,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(ny-2), kind=mytype)*real(transy(ny-2), kind=mytype)&
                    *real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                    -(real(zk2(k), kind=mytype)*real(transy(ny-2), kind=mytype)*real(transy(ny-2), kind=mytype)*&
                    real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
                    -real(cw22(i,ny/2,k), kind=mytype)*real(cw22(i,ny/2,k), kind=mytype)*xa0*xa0-&
                    xa1*xa1*(real(cw22(i,ny/2,k), kind=mytype)*real(cw22(i,ny/2-1,k), kind=mytype)),&
                    -(aimag(xk2(i))*aimag(transy(ny-2))*aimag(transy(ny-2))*aimag(transz(k))*aimag(transz(k)))&
                    -(aimag(zk2(k))*aimag(transy(ny-2))*aimag(transy(ny-2))*aimag(transx(i))*aimag(transx(i)))&
                    -aimag(cw22(i,ny/2,k))*aimag(cw22(i,ny/2,k))*xa0*xa0-&
                    xa1*xa1*(aimag(cw22(i,ny/2,k))*aimag(cw22(i,ny/2-1,k))), kind=mytype)
                    a2(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2), kind=mytype)*real(transy(2), kind=mytype)*&
                    real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                    -(real(zk2(k), kind=mytype)*real(transy(2), kind=mytype)*real(transy(2), kind=mytype)*&
                    real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
                    -real(cw2(i,1,k), kind=mytype)*real(cw2(i,1,k), kind=mytype)*(xa0-xa1)*(xa0+xa1)-xa1*xa1*(real(cw2(i,1,k), kind=mytype)*&
                    real(cw2(i,2,k), kind=mytype)),&
                    -(aimag(xk2(i))*aimag(transy(2))*aimag(transy(2))*aimag(transz(k))*aimag(transz(k)))&
                    -(aimag(zk2(k))*aimag(transy(2))*aimag(transy(2))*aimag(transx(i))*aimag(transx(i)))&
                    -aimag(cw2(i,1,k))*aimag(cw2(i,1,k))*(xa0-xa1)*(xa0+xa1)-xa1*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,2,k))), kind=mytype)
                    a2(i,ny/2,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(ny-1), kind=mytype)*real(transy(ny-1), kind=mytype)*&
                    real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                    -(real(zk2(k), kind=mytype)*real(transy(ny-1), kind=mytype)*real(transy(ny-1), kind=mytype)*&
                    real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
                    -real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype)*(xa0+xa1)*(xa0+xa1)-&
                    xa1*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype)),&
                    -(aimag(xk2(i))*aimag(transy(ny-1))*aimag(transy(ny-1))*aimag(transz(k))*aimag(transz(k)))&
                    -(aimag(zk2(k))*aimag(transy(ny-1))*aimag(transy(ny-1))*aimag(transx(i))*aimag(transx(i)))&
                    -aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2,k))*(xa0+xa1)*(xa0+xa1)-&
                    xa1*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2-1,k))), kind=mytype)
                enddo
            enddo





            !sup diag +1
            do k=sp%yst(3),sp%yen(3)
                do j=2,ny/2-1
                    do i=sp%yst(1),sp%yen(1)
                        a(i,j,k,4)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j+1,k), kind=mytype)+real(cw22(i,j+1,k), kind=mytype)*&
                        real(cw22(i,j+1,k), kind=mytype)),&
                        xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))+aimag(cw22(i,j+1,k))*aimag(cw22(i,j+1,k))), kind=mytype)
                        a2(i,j,k,4)=cmplx(xa0*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j+1,k), kind=mytype)+real(cw2(i,j+1,k), kind=mytype)*&
                        real(cw2(i,j+1,k), kind=mytype)),&
                        xa0*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j+1,k))+aimag(cw2(i,j+1,k))*aimag(cw2(i,j+1,k))), kind=mytype)
                    enddo
                enddo
                do i=sp%yst(1),sp%yen(1)
                    a(i,1,k,4)=2.*cmplx((xa0*xa1*(real(cw22(i,1,k), kind=mytype)*real(cw22(i,2,k), kind=mytype)+real(cw22(i,2,k), kind=mytype)*&
                    real(cw22(i,2,k), kind=mytype))),&
                    (xa0*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))+aimag(cw22(i,2,k))*aimag(cw22(i,2,k)))), kind=mytype)
                    a2(i,1,k,4)=cmplx((xa0-xa1)*xa1*(real(cw2(i,1,k), kind=mytype)*real(cw2(i,2,k), kind=mytype))&
                    +xa0*xa1*(real(cw2(i,2,k), kind=mytype)*real(cw2(i,2,k), kind=mytype)),&
                    (xa0-xa1)*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,2,k)))&
                    +xa0*xa1*(aimag(cw2(i,2,k))*aimag(cw2(i,2,k))), kind=mytype)
                    a2(i,ny/2-1,k,4)=cmplx(xa0*xa1*(real(cw2(i,ny/2-1,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype))&
                    +(xa0+xa1)*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype)),&
                    xa0*xa1*(aimag(cw2(i,ny/2-1,k))*aimag(cw2(i,ny/2,k)))&
                    +(xa0+xa1)*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2,k))), kind=mytype)
                    a2(i,ny/2,k,4)=0.
                enddo
            enddo



            !sup diag +2
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    do j=1,ny/2-2
                        a(i,j,k,5)=cmplx(-real(cw22(i,j+1,k), kind=mytype)*real(cw22(i,j+2,k), kind=mytype)*xa1*xa1,&
                        -aimag(cw22(i,j+1,k))*aimag(cw22(i,j+2,k))*xa1*xa1, kind=mytype)
                        a2(i,j,k,5)=cmplx(-real(cw2(i,j+1,k), kind=mytype)*real(cw2(i,j+2,k), kind=mytype)*xa1*xa1,&
                        -aimag(cw2(i,j+1,k))*aimag(cw2(i,j+2,k))*xa1*xa1, kind=mytype)
                    enddo
                    a(i,1,k,5)=cmplx(real(a(i,1,k,5), kind=mytype)*2.,aimag(a(i,1,k,5))*2., kind=mytype)
                    a(i,ny/2-1,k,5)=0.
                    a(i,ny/2,k,5)=0.
                    a2(i,ny/2-1,k,5)=0.
                    a2(i,ny/2,k,5)=0.
                enddo
            enddo



            !inf diag -1
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    do j=2,ny/2
                        a(i,j,k,2)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j-1,k), kind=mytype)*&
                        real(cw22(i,j-1,k), kind=mytype)),&
                        xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j-1,k))*aimag(cw22(i,j-1,k))), kind=mytype)
                        a2(i,j,k,2)=cmplx(xa0*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j-1,k), kind=mytype)+real(cw2(i,j-1,k), kind=mytype)*&
                        real(cw2(i,j-1,k), kind=mytype)),&
                        xa0*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j-1,k))+aimag(cw2(i,j-1,k))*aimag(cw2(i,j-1,k))), kind=mytype)
                    enddo
                    a(i,1,k,2)=0.
                    a2(i,1,k,2)=0.
                    a2(i,2,k,2)=cmplx(xa0*xa1*(real(cw2(i,2,k), kind=mytype)*real(cw2(i,1,k), kind=mytype))&
                    +(xa0+xa1)*xa1*(real(cw2(i,1,k), kind=mytype)*real(cw2(i,1,k), kind=mytype)),&
                    xa0*xa1*(aimag(cw2(i,2,k))*aimag(cw2(i,1,k)))&
                    +(xa0+xa1)*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,1,k))), kind=mytype)
                    a2(i,ny/2,k,2)=cmplx((xa0+xa1)*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype))&
                    +xa0*xa1*(real(cw2(i,ny/2-1,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype)),&
                    (xa0+xa1)*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2-1,k)))&
                    +xa0*xa1*(aimag(cw2(i,ny/2-1,k))*aimag(cw2(i,ny/2-1,k))), kind=mytype)
                enddo
            enddo
            !inf diag -2
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    do j=3,ny/2
                        a(i,j,k,1)=cmplx(-real(cw22(i,j-1,k), kind=mytype)*real(cw22(i,j-2,k), kind=mytype)*xa1*xa1,&
                        -aimag(cw22(i,j-1,k))*aimag(cw22(i,j-2,k))*xa1*xa1, kind=mytype)
                        a2(i,j,k,1)=cmplx(-real(cw2(i,j-1,k), kind=mytype)*real(cw2(i,j-2,k), kind=mytype)*xa1*xa1,&
                        -aimag(cw2(i,j-1,k))*aimag(cw2(i,j-2,k))*xa1*xa1, kind=mytype)
                    enddo
                    a(i,1,k,1)=0.
                    a(i,2,k,1)=0.
                    a2(i,1,k,1)=0.
                    a2(i,2,k,1)=0.
                enddo
            enddo
            !not to have a singular matrice
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    if ((real(xk2(i), kind=mytype)==0.).and.(real(zk2(k), kind=mytype)==0)) then
                        a(i,1,k,3)=cmplx(1.,1., kind=mytype)
                        a(i,1,k,4)=0.
                        a(i,1,k,5)=0.
                    endif
                enddo
            enddo

        else
            xa0=alpha/PI+1./2./beta/PI
            xa1=-1./4./beta/PI
            !
            !construction of the pentadiagonal matrice
            !
            do k=sp%yst(3),sp%yen(3)
                do j=1,nym
                    do i=sp%yst(1),sp%yen(1)
                        cw22(i,j,k)=cmplx(real(yky(j), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind=mytype),&
                        aimag(yky(j))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
                    enddo
                enddo
            enddo

            !main diagonal
            do k=sp%yst(3),sp%yen(3)
                do j=2,nym-1
                    do i=sp%yst(1),sp%yen(1)
                        a3(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(j), kind=mytype)*real(transy(j), kind=mytype)*&
                        real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                        -(real(zk2(k), kind=mytype)*real(transy(j), kind=mytype)*real(transy(j), kind=mytype)*&
                        real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
                        -real(cw22(i,j,k), kind=mytype)*real(cw22(i,j,k), kind=mytype)*xa0*xa0-&
                        xa1*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j,k), kind=mytype)*&
                        real(cw22(i,j+1,k), kind=mytype)),&
                        -(aimag(xk2(i))*aimag(transy(j))*aimag(transy(j))*aimag(transz(k))*aimag(transz(k)))&
                        -(aimag(zk2(k))*aimag(transy(j))*aimag(transy(j))*aimag(transx(i))*aimag(transx(i)))&
                        -aimag(cw22(i,j,k))*aimag(cw22(i,j,k))*xa0*xa0-&
                        xa1*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))), kind=mytype)
                    enddo
                enddo
            enddo

            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    a3(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*&
                    real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                    -(real(zk2(k), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*real(transx(i), kind=mytype)*&
                    real(transx(i), kind=mytype))&
                    -real(cw22(i,1,k), kind=mytype)*real(cw22(i,1,k), kind=mytype)*xa0*xa0-xa1*xa1*(real(cw22(i,1,k), kind=mytype)*&
                    real(cw22(i,2,k), kind=mytype)),&
                    -(aimag(xk2(i))*aimag(transy(1))*aimag(transy(1))*aimag(transz(k))*aimag(transz(k)))&
                    -(aimag(zk2(k))*aimag(transy(1))*aimag(transy(1))*aimag(transx(i))*aimag(transx(i)))&
                    -aimag(cw22(i,1,k))*aimag(cw22(i,1,k))*xa0*xa0-xa1*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))), kind=mytype)
                    a3(i,nym,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(nym), kind=mytype)*real(transy(nym), kind=mytype)*&
                    real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
                    -(real(zk2(k), kind=mytype)*real(transy(nym), kind=mytype)*real(transy(nym), kind=mytype)*real(transx(i), kind=mytype)*&
                    real(transx(i), kind=mytype))&
                    -real(cw22(i,nym,k), kind=mytype)*real(cw22(i,nym,k), kind=mytype)*xa0*xa0-&
                    xa1*xa1*(real(cw22(i,nym,k), kind=mytype)*real(cw22(i,nym-1,k), kind=mytype)),&
                    -(aimag(xk2(i))*aimag(transy(nym))*aimag(transy(nym))*aimag(transz(k))*aimag(transz(k)))&
                    -(aimag(zk2(k))*aimag(transy(nym))*aimag(transy(nym))*aimag(transx(i))*aimag(transx(i)))&
                    -aimag(cw22(i,nym,k))*aimag(cw22(i,nym,k))*xa0*xa0-&
                    xa1*xa1*(aimag(cw22(i,nym,k))*aimag(cw22(i,nym-1,k))), kind=mytype)
                enddo
            enddo




            !sup diag +1
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    do j=2,nym-1
                        a3(i,j,k,4)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j+1,k), kind=mytype)+real(cw22(i,j+1,k), kind=mytype)*&
                        real(cw22(i,j+1,k), kind=mytype)),&
                        xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))+aimag(cw22(i,j+1,k))*aimag(cw22(i,j+1,k))), kind=mytype)
                    enddo
                    a3(i,1,k,4)=cmplx((xa0*xa1*(real(cw22(i,1,k), kind=mytype)*real(cw22(i,2,k), kind=mytype)+real(cw22(i,2,k), kind=mytype)*&
                    real(cw22(i,2,k), kind=mytype))),&
                    (xa0*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))+aimag(cw22(i,2,k))*aimag(cw22(i,2,k)))), kind=mytype)
                enddo
            enddo
            !sup diag +2
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    do j=1,nym-2
                        a3(i,j,k,5)=cmplx(-real(cw22(i,j+1,k), kind=mytype)*real(cw22(i,j+2,k), kind=mytype)*xa1*xa1,&
                        -aimag(cw22(i,j+1,k))*aimag(cw22(i,j+2,k))*xa1*xa1, kind=mytype)
                    enddo
                    !a3(i,1,k,5)=a3(i,1,k,5)*2.
                    !a3(i,1,k,5)=0.
                    a3(i,nym-1,k,5)=0.
                    a3(i,nym,k,5)=0.
                enddo
            enddo


            !inf diag -1
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    do j=2,nym
                        a3(i,j,k,2)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j-1,k), kind=mytype)*&
                        real(cw22(i,j-1,k), kind=mytype)),&
                        xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j-1,k))*aimag(cw22(i,j-1,k))), kind=mytype)
                    enddo
                    a3(i,1,k,2)=0.
                enddo
            enddo
            !inf diag -2
            do k=sp%yst(3),sp%yen(3)
                do i=sp%yst(1),sp%yen(1)
                    do j=3,nym
                        a3(i,j,k,1)=cmplx(-real(cw22(i,j-1,k), kind=mytype)*real(cw22(i,j-2,k), kind=mytype)*xa1*xa1,&
                        -aimag(cw22(i,j-1,k))*aimag(cw22(i,j-2,k))*xa1*xa1, kind=mytype)
                    enddo
                    a3(i,1,k,1)=0.
                    a3(i,2,k,1)=0.
                enddo
            enddo

            !not to have a singular matrice
            !   do k=sp%yst(3),sp%yen(3)
            !   do i=sp%yst(1),sp%yen(1)
            !      if ((xkx(i)==0.).and.(zkz(k)==0)) then
            if (nrank==0) then
                a3(1,1,1,3)=cmplx(1.,1., kind=mytype)
                a3(1,1,1,4)=0.
                a3(1,1,1,5)=0.
            endif
        !      endif
        !   enddo
        !   enddo
        endif





        return
    end subroutine matrice_refinement





    subroutine abxyz(ax,ay,az,bx,by,bz,nx0,ny0,nz0,bcx,bcy,bcz)

        use mathematical_constants

        implicit none

        integer, intent(IN) :: nx0,ny0,nz0
        integer, intent(IN) :: bcx,bcy,bcz
        real(mytype), dimension(:), intent(OUT) :: ax,bx
        real(mytype), dimension(:), intent(OUT) :: ay,by
        real(mytype), dimension(:), intent(OUT) :: az,bz

        integer :: i,j,k

        if (bcx==0) then
            do i=1,nx0
                ax(i) = sin(real(i-1, kind=mytype)*PI/real(nx0, kind=mytype))
                bx(i) = cos(real(i-1, kind=mytype)*PI/real(nx0, kind=mytype))
            end do
        else if (bcx==1) then
            do i=1,nx0
                ax(i) = sin(real(i-1, kind=mytype)*PI/2.0_mytype/ &
                real(nx0, kind=mytype))
                bx(i) = cos(real(i-1, kind=mytype)*PI/2.0_mytype/ &
                real(nx0, kind=mytype))
            end do
        end if

        if (bcy==0) then
            do j=1,ny0
                ay(j) = sin(real(j-1, kind=mytype)*PI/real(ny0, kind=mytype))
                by(j) = cos(real(j-1, kind=mytype)*PI/real(ny0, kind=mytype))
            end do
        else if (bcy==1) then
            do j=1,ny0
                ay(j) = sin(real(j-1, kind=mytype)*PI/2.0_mytype/ &
                real(ny0, kind=mytype))
                by(j) = cos(real(j-1, kind=mytype)*PI/2.0_mytype/ &
                real(ny0, kind=mytype))
            end do
        end if

        if (bcz==0) then
            do k=1,nz0
                az(k) = sin(real(k-1, kind=mytype)*PI/real(nz0, kind=mytype))
                bz(k) = cos(real(k-1, kind=mytype)*PI/real(nz0, kind=mytype))
            end do
        else if (bcz==1) then
            do k=1,nz0
                az(k) = sin(real(k-1, kind=mytype)*PI/2.0_mytype/ &
                real(nz0, kind=mytype))
                bz(k) = cos(real(k-1, kind=mytype)*PI/2.0_mytype/ &
                real(nz0, kind=mytype))
            end do
        end if

        return
    end subroutine abxyz

    subroutine perform_kxyz()

        USE derivX
        USE derivY
        USE derivZ
        USE decomp_2d

        implicit none

        integer :: i,j,k
        real(mytype) :: w,wp,w1,w1p
        complex(mytype) :: xyzk
        complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
        complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3
        complex(mytype) :: tmp4,tmp5,tmp6
        real*8      :: tp1, tp2, tp3
        real*8     :: ast, bst

        tmp4=cmplx(1.d0, 0.d0)
        tmp5=cmplx(1.d0, 0.d0)
        tmp6=cmplx(1.d0, 0.d0)

        if (arrang==FULL_STAGGERED) then
            ast=1.d0
            bst=0.d0
        else
            ast=0.d0
            bst=1.d0
        end if


        if ((nclx==0).and.(nclz==0).and.((ncly==1).or.(ncly==2))) then
            do k = sp%yst(3), sp%yen(3)
                do j = sp%yst(2), sp%yen(2)
                    do i = sp%yst(1), sp%yen(1)

                        ! Relation 21, partie 1 (b,c)
                        xtt=cmplx((bicix6*2.*cos(real(exs(i))*3.*dx3/2.)+&
                        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)),&
                        (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx3/2.)+&
                        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)), kind=mytype)

                        ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
                        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)),&
                        (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
                        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)), kind=mytype)

                        ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                        ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)),&
                        (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                        ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)), kind=mytype)


                        ! Relation 21, partie 2 (a)
                        xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)),&
                        (aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)), kind=mytype)

                        ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)),&
                        (aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)), kind=mytype)

                        ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)),&
                        (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)), kind=mytype)

                        ! Relation 21, partie 2 (cote gauche)
                        xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx3)),&
                        (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx3)), kind=mytype)

                        yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)),&
                        (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)), kind=mytype)

                        zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)),&
                        (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)), kind=mytype)


                        ! Relation 44 (interpolations pour maillage uniforme)
                        !(kx*Ty*Tz)**2
                        ! NOTE: on utilise la propriété mathématique: (a+ia)/(b+ib)= (a/b+i0)
                        xt2=ast*cmplx(1.d0, 0.d0)+bst*((((ytt1+ytt)/yt1)*((ztt1+ztt)/zt1))**2)
                        !(ky*Tx*Tz)**2
                        yt2=ast*cmplx(1.d0, 0.d0)+bst*((((xtt1+xtt)/xt1)*((ztt1+ztt)/zt1))**2)
                        !(kz*Tx*Ty)**2
                        zt2=ast*cmplx(1.d0, 0.d0)+bst*((((xtt1+xtt)/xt1)*((ytt1+ytt)/yt1))**2)

                        xyzk=xt2*xk2(i) + yt2*yk2(j) + zt2*zk2(k)
                        kxyz(i,j,k)=xyzk
                    !   print *,i,j,k, kxyz(i,j,k)
                    enddo
                enddo
            enddo
        else
            if (nclz==0) then
                do k = sp%xst(3),sp%xen(3)
                    do j = sp%xst(2),sp%xen(2)
                        do i = sp%xst(1),sp%xen(1)

                            xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx3/2.)+&
                            cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)),&
                            (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx3/2.)+&
                            cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)), kind=mytype)

                            ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
                            ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)),&
                            (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
                            ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)), kind=mytype)

                            ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                            ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)),&
                            (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                            ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)), kind=mytype)

                            xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)),&
                            (aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)), kind=mytype)

                            ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)),&
                            (aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)), kind=mytype)

                            ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)),&
                            (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)), kind=mytype)

                            xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx3)),&
                            (1.+2.*ailcaix6*cos(real(exs(i))*dx3)), kind=mytype)

                            yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)),&
                            (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)), kind=mytype)

                            zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)),&
                            (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)), kind=mytype)

                            !(kx*Ty*Tz)**2
                            xt2=ast*cmplx(1.d0, 0.d0)+bst*((((ytt1+ytt)/yt1)*((ztt1+ztt)/zt1))**2)
                            !(ky*Tx*Tz)**2
                            yt2=ast*cmplx(1.d0, 0.d0)+bst*((((xtt1+xtt)/xt1)*((ztt1+ztt)/zt1))**2)
                            !(kz*Tx*Ty)**2
                            zt2=ast*cmplx(1.d0, 0.d0)+bst*((((xtt1+xtt)/xt1)*((ytt1+ytt)/yt1))**2)

                            xyzk=xt2*xk2(i) + yt2*yk2(j) + zt2*zk2(k)
                            kxyz(i,j,k)=xyzk
                        !   print *,i,j,k, kxyz(i,j,k)
                        enddo
                    enddo
                enddo
            else
                do k = sp%xst(3),sp%xen(3)
                    do j = sp%xst(2),sp%xen(2)
                        do i = sp%xst(1),sp%xen(1)
                            xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx3/2.)+&
                            cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)),&
                            (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx3/2.)+&
                            cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)), kind=mytype)

                            ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
                            ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)),&
                            (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
                            ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)), kind=mytype)
                            !
                            ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                            ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)),&
                            (biciz6*2.*cos(aimag(ezs(k))*3.*dx1/2.)+&
                            ciciz6*2.*cos(aimag(ezs(k))*5.*dx1/2.)), kind=mytype)
                            !
                            xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)),&
                            (aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)), kind=mytype)
                            ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)),&
                            (aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)), kind=mytype)
                            !
                            ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)),&
                            (aiciz6*2.*cos(aimag(ezs(k))*dx1/2.)), kind=mytype)
                            !
                            xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx3)),&
                            (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx3)), kind=mytype)
                            yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)),&
                            (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)), kind=mytype)
                            zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)),&
                            (1.+2.*ailcaiz6*cos(aimag(ezs(k))*dx1)), kind=mytype)

                            ! tmp1<=Tr(kx)
                            tmp1=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
                            aimag(ztt1+ztt)/aimag(zt1), kind=mytype)
                            ! tmp2<=Tr(ky)
                            tmp2=cmplx(real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype),&
                            real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype), kind=mytype)
                            ! tmp3<=Tr(kz)
                            tmp3=cmplx(real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype),&
                            real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype), kind=mytype)

                            ! tmp4<=(Tr(kx)Tr(ky))**2
                            tmp4=cmplx((real(tmp1, kind=mytype)*real(tmp2, kind=mytype))**2,(aimag(tmp1)*aimag(tmp2))**2, kind=mytype)

                            ! tmp5<=(Tr(kx)Tr(kz))**2
                            tmp5=cmplx((real(tmp1, kind=mytype)*real(tmp3, kind=mytype))**2,(aimag(tmp1)*aimag(tmp3))**2, kind=mytype)

                            ! tmp6<=(Tr(ky)Tr(kz))**2
                            tmp6=cmplx((real(tmp3, kind=mytype)*real(tmp2, kind=mytype))**2,(aimag(tmp3)*aimag(tmp2))**2, kind=mytype)

                            ! On peut retrouver le cas collocated avec tmp[4,5,6]=cmplx(1.d0, 0.d0)
                            tmp1=cmplx(real(tmp4, kind=mytype)*real(xk2(i), kind=mytype),aimag(tmp4)*aimag(xk2(i)), kind=mytype)
                            tmp2=cmplx(real(tmp5, kind=mytype)*real(yk2(j), kind=mytype),aimag(tmp5)*aimag(yk2(j)), kind=mytype)
                            tmp3=cmplx(real(tmp6, kind=mytype)*real(zk2(k), kind=mytype),aimag(tmp6)*aimag(zk2(k)), kind=mytype)

                            xyzk=tmp1+tmp2+tmp3
                            kxyz(i,j,k)=xyzk
                        !         print *,i,j,k,zt1,yt1
                        enddo
                    enddo
                enddo
            endif
        endif
    end subroutine perform_kxyz

    subroutine perform_TF_coeffs(transx, transy, transz)
        USE decomp_2d
        USE derivX
        USE derivY
        USE derivZ

        implicit none

        integer :: i,j,k

        complex(mytype),dimension(sp%yst(1):sp%yen(1)) :: transx
        complex(mytype),dimension(sp%yst(2):sp%yen(2)) :: transy
        complex(mytype),dimension(sp%yst(3):sp%yen(3)) :: transz

        complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
        complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3

        do i = sp%yst(1),sp%yen(1)

            ! transx <= Tr(kdx)
            xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx3/2.)+&
            cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)),&
            (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx3/2.)+&
            cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx3/2.)), kind=mytype)

            xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)),&
            (aicix6*2.*cos(real(exs(i), kind=mytype)*dx3/2.)), kind=mytype)

            xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx3)),&
            (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx3)), kind=mytype)

            transx(i)=cmplx(real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype),&
            real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype), kind=mytype)!(xtt+xtt)/xt1
        enddo

        do j = sp%yst(2),sp%yen(2)

            ! transy <= Tr(kdy)
            ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
            ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)),&
            (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dx2/2.)+&
            ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dx2/2.)), kind=mytype)

            ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)),&
            (aiciy6*2.*cos(real(eys(j), kind=mytype)*dx2/2.)), kind=mytype)

            yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)),&
            (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dx2)), kind=mytype)

            transy(j)=cmplx(real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype),&
            real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype), kind=mytype)!(ytt+ytt)/yt1
        enddo

        if (nclz==0) then
            do k = sp%yst(3),sp%yen(3)

                ! transz <= Tr(kdz)
                ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)),&
                (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)), kind=mytype)

                ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)),&
                (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)), kind=mytype)

                zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)),&
                (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)), kind=mytype)

                transz(k)=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
                aimag(ztt1+ztt)/aimag(zt1), kind=mytype)!(ztt+ztt)/zt1
            enddo

        else
            do k = sp%yst(3),sp%yen(3)

                ! transz <= Tr(kdz)
                ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dx1/2.)+&
                ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dx1/2.)),&
                (biciz6*2.*cos(aimag(ezs(k))*3.*dx1/2.)+&
                ciciz6*2.*cos(aimag(ezs(k))*5.*dx1/2.)), kind=mytype)

                ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dx1/2.)),&
                (aiciz6*2.*cos(aimag(ezs(k))*dx1/2.)), kind=mytype)

                zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dx1)),&
                (1.+2.*ailcaiz6*cos(aimag(ezs(k))*dx1)), kind=mytype)

                transz(k)=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
                aimag(ztt1+ztt)/aimag(zt1), kind=mytype)!(ztt+ztt)/zt1
            enddo
        endif

        if (arrang==FULL_STAGGERED) then
            transx=cmplx(1.d0, 1.d0)
            transy=cmplx(1.d0, 1.d0)
            transz=cmplx(1.d0, 1.d0)
        end if

    end subroutine perform_TF_coeffs

end module poisson_generic_solver



module poisson_interface
    USE decomp_2d
    implicit none

    public :: poisson_init, solve_Poisson, solve_Poisson_ibm
contains


    subroutine convert_press_000_to_001(pr_000, pr_001, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                     :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable :: pr_000, pr_001
        real(mytype), dimension(:,:,:), allocatable, save :: pr_001_z

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_000  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_001  ! decomposition object


        call decomp_info_init(n1, n2, n3, ph_000)
        call decomp_info_init(n1, n2, n3+a3, ph_001)

        if (.not. allocated(pr_001_z)) then
            allocate (pr_001_z(ph_001%zst(1):ph_001%zen(1),ph_001%zst(2):ph_001%zen(2),ph_001%zst(3):ph_001%zen(3)))
        endif

        pr_001_z=0.d0
        do i = ph_000%zst(1), ph_000%zen(1)
            do j = ph_000%zst(2), ph_000%zen(2)
                do k = ph_000%zst(3), ph_000%zen(3)
                    pr_001_z(i,j,k)=pr_000(i,j,k)
                end do
            end do
        end do

        call transpose_z_to_y(pr_001_z, pr_001, ph_001)



    end subroutine convert_press_000_to_001


    subroutine convert_press_001_to_000(pr_001, pr_000, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                     :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable :: pr_000
        real(mytype), dimension(:,:,:), allocatable :: pr_001
        real(mytype), dimension(:,:,:), allocatable, save :: pr_001_z

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_000  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_001  ! decomposition object


        call decomp_info_init(n1, n2, n3, ph_000)
        call decomp_info_init(n1, n2, n3+a3, ph_001)

        if (.not. allocated(pr_001_z)) then
            allocate (pr_001_z(ph_001%zst(1):ph_001%zen(1),ph_001%zst(2):ph_001%zen(2),ph_001%zst(3):ph_001%zen(3)))
        endif

        call transpose_y_to_z(pr_001, pr_001_z, ph_001)

        do i = ph_000%zst(1), ph_000%zen(1)
            do j = ph_000%zst(2), ph_000%zen(2)
                do k = ph_000%zst(3), ph_000%zen(3)
                    pr_000(i,j,k)=pr_001_z(i,j,k)
                end do
            end do
        end do

    end subroutine convert_press_001_to_000


    subroutine convert_press_001_to_011(pr_001, pr_011, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                             :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable         :: pr_001, pr_011
        real(mytype), dimension(:,:,:), allocatable, save   :: pr_011_y

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_001  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_011  ! decomposition object


        call decomp_info_init(n1, n2, n3+a3, ph_001)
        call decomp_info_init(n1, n2+a2, n3+a3, ph_011)

        if (.not. allocated(pr_011_y)) then
            allocate (pr_011_y(ph_011%yst(1):ph_011%yen(1),ph_011%yst(2):ph_011%yen(2),ph_011%yst(3):ph_011%yen(3)))
        endif

        pr_011_y=0.d0
        do i = ph_001%yst(1), ph_001%yen(1)
            do j = ph_001%yst(2), ph_001%yen(2)
                do k = ph_001%yst(3), ph_001%yen(3)
                    pr_011_y(i,j,k)=pr_001(i,j,k)
                end do
            end do
        end do

        call transpose_y_to_x(pr_011_y, pr_011, ph_011)



    end subroutine convert_press_001_to_011


    subroutine convert_press_011_to_001(pr_011, pr_001, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                             :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable         :: pr_011, pr_001
        real(mytype), dimension(:,:,:), allocatable, save   :: pr_011_y

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_001  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_011  ! decomposition object


        call decomp_info_init(n1, n2, n3+a3, ph_001)
        call decomp_info_init(n1, n2+a2, n3+a3, ph_011)

        if (.not. allocated(pr_011_y)) then
            allocate (pr_011_y(ph_011%yst(1):ph_011%yen(1),ph_011%yst(2):ph_011%yen(2),ph_011%yst(3):ph_011%yen(3)))
        endif


        call transpose_x_to_y(pr_011, pr_011_y, ph_011)

        do i = ph_001%yst(1), ph_001%yen(1)
            do j = ph_001%yst(2), ph_001%yen(2)
                do k = ph_001%yst(3), ph_001%yen(3)
                    pr_001(i,j,k)=pr_011_y(i,j,k)
                end do
            end do
        end do

    end subroutine convert_press_011_to_001


    subroutine convert_press_011_to_111(pr_011, pr_111, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                     :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable :: pr_011, pr_111

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_011  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_111  ! decomposition object


        call decomp_info_init(n1, n2+a2, n3+a3, ph_011)
        call decomp_info_init(n1+a1, n2+a2, n3+a3, ph_111)

        pr_111=0.d0
        do i = ph_011%xst(1), ph_011%xen(1)
            do j = ph_011%xst(2), ph_011%xen(2)
                do k = ph_011%xst(3), ph_011%xen(3)
                    pr_111(i,j,k)=pr_011(i,j,k)
                end do
            end do
        end do


    end subroutine convert_press_011_to_111


    subroutine convert_press_111_to_011(pr_111, pr_011, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                     :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable :: pr_011, pr_111

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_011  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_111  ! decomposition object


        call decomp_info_init(n1, n2+a2, n3+a3, ph_011)
        call decomp_info_init(n1+a1, n2+a2, n3+a3, ph_111)

        do i = ph_011%xst(1), ph_011%xen(1)
            do j = ph_011%xst(2), ph_011%xen(2)
                do k = ph_011%xst(3), ph_011%xen(3)
                    pr_011(i,j,k)=pr_111(i,j,k)
                end do
            end do
        end do


    end subroutine convert_press_111_to_011

    subroutine convert_000_to_111(pr_000, pr_111, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                             :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable         :: pr_000, pr_111
        real(mytype), dimension(:,:,:), allocatable, save   :: pr_011_x
        real(mytype), dimension(:,:,:), allocatable, save   :: pr_001_y

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_000  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_001  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_011  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_111  ! decomposition object


        call decomp_info_init(n1, n2, n3, ph_000)
        call decomp_info_init(n1, n2, n3+a3, ph_001)
        call decomp_info_init(n1, n2+a2, n3+a3, ph_011)
        call decomp_info_init(n1+a1, n2+a2, n3+a3, ph_111)

        if (.not. allocated(pr_001_y)) then
            allocate (pr_001_y(ph_001%yst(1):ph_001%yen(1),ph_001%yst(2):ph_001%yen(2),ph_001%yst(3):ph_001%yen(3)))
            allocate (pr_011_x(ph_011%xst(1):ph_011%xen(1),ph_011%xst(2):ph_011%xen(2),ph_011%xst(3):ph_011%xen(3)))
        endif


        call  convert_press_000_to_001(pr_000, pr_001_y, n1, n2, n3, a1, a2, a3)
        call  convert_press_001_to_011(pr_001_y, pr_011_x, n1, n2, n3, a1, a2, a3)
        call  convert_press_011_to_111(pr_011_x, pr_111, n1, n2, n3, a1, a2, a3)

    end subroutine convert_000_to_111

    subroutine convert_111_to_000(pr_111, pr_000, n1, n2, n3, a1, a2, a3)
        implicit none
        integer                                             :: n1, n2, n3, a1, a2, a3
        real(mytype), dimension(:,:,:), allocatable         :: pr_000, pr_111
        real(mytype), dimension(:,:,:), allocatable, save   :: pr_011_x
        real(mytype), dimension(:,:,:), allocatable, save   :: pr_001_y

        integer :: i,j,k

        TYPE(DECOMP_INFO) :: ph_000  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_001  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_011  ! decomposition object
        TYPE(DECOMP_INFO) :: ph_111  ! decomposition object


        call decomp_info_init(n1, n2, n3, ph_000)
        call decomp_info_init(n1, n2, n3+a3, ph_001)
        call decomp_info_init(n1, n2+a2, n3+a3, ph_011)
        call decomp_info_init(n1+a1, n2+a2, n3+a3, ph_111)

        if (.not. allocated(pr_001_y)) then
            allocate (pr_001_y(ph_001%yst(1):ph_001%yen(1),ph_001%yst(2):ph_001%yen(2),ph_001%yst(3):ph_001%yen(3)))
            allocate (pr_011_x(ph_011%xst(1):ph_011%xen(1),ph_011%xst(2):ph_011%xen(2),ph_011%xst(3):ph_011%xen(3)))
        endif



        call convert_press_111_to_011(pr_111, pr_011_x, n1, n2, n3, a1, a2, a3)
        call convert_press_011_to_001(pr_011_x, pr_001_y, n1, n2, n3, a1, a2, a3)
        call convert_press_001_to_000(pr_001_y, pr_000, n1, n2, n3, a1, a2, a3)

    end subroutine convert_111_to_000

    subroutine poisson_init(n1, n2, n3, nclx1, ncly1, nclz1, lx, ly, lz, istret, alpha, beta, arrang)
        use boundaries
        USE poisson_generic_solver
        implicit none

        integer, intent(IN) :: nclx1, ncly1, nclz1
        integer, intent(IN) :: n1, n2, n3
        integer, intent(IN) :: istret, arrang
        real*8, intent(IN)  :: lx, ly, lz
        real*8, intent(IN)  :: alpha, beta

        integer             :: nclx, ncly, nclz
        integer             :: nx, ny, nz

        select case (nclx1)
            case (UNBOUNDED, FRINGE)
                nclx=0
                nx=n1-1
            case (FREESLIP)
                nclx=1
                nx=n1
            ! Noslip or open boundary conditions
            case default
                nclx=2
                nx=n1

        end select

        select case (ncly1)
            case (UNBOUNDED, FRINGE)
                ncly=0
                ny=n2-1
            case (FREESLIP)
                ncly=1
                ny=n2
            ! Noslip or open boundary conditions
            case default
                ncly=2
                ny=n2

        end select

        select case (nclz1)
            case (UNBOUNDED, FRINGE)
                nclz=0
                nz=n3-1
            case (FREESLIP)
                nclz=1
                nz=n3
            ! Noslip or open boundary conditions
            case default
                nclz=2
                nz=n3

        end select

        call decomp_2d_poisson_init(nx,ny,nz, nclx, ncly, nclz, lx, ly, lz, istret, alpha, beta, arrang)

    end subroutine poisson_init

    subroutine solve_Poisson(RHS_z, pr_EDT_z)
        use mesh
        USE decomp_2d
        USE poisson_generic_solver
        implicit none
        real(mytype), dimension(:,:,:), intent(in)      :: RHS_z
        real(mytype), dimension(:,:,:), intent(out)     :: pr_EDT_z

        integer     :: i,j,k
        real(mytype), dimension(:,:,:), allocatable, save :: pr_EDT_x, pr_EDT_y, pp3


        TYPE(DECOMP_INFO), save :: ph_111  ! decomposition object
        TYPE(DECOMP_INFO), save :: ph_000  ! decomposition object

        pr_EDT_z=RHS_z

        if (.not. allocated(pp3)) then

            call decomp_info_init(n1, n2, n3, ph_111)
            call decomp_info_init(n1-1,n2-1,n3-1,ph_000)

            allocate (pr_EDT_x(ph_111%xst(1):ph_111%xen(1),ph_111%xst(2):ph_111%xen(2),ph_111%xst(3):ph_111%xen(3)))
            allocate (pr_EDT_y(ph_111%yst(1):ph_111%yen(1),ph_111%yst(2):ph_111%yen(2),ph_111%yst(3):ph_111%yen(3)))
            allocate (pp3(ph_000%zst(1):ph_000%zen(1),ph_000%zst(2):ph_000%zen(2),ph_000%zst(3):ph_000%zen(3)))

        end if

        call transpose_z_to_y(pr_EDT_z, pr_EDT_y, ph_111)
        call transpose_y_to_x(pr_EDT_y, pr_EDT_x, ph_111)
        call convert_111_to_000(pr_EDT_x, pp3, n1-1, n2-1, n3-1, 1, 1, 1)

        call decomp_2d_poisson_stg(pp3)

        call convert_000_to_111(pp3, pr_EDT_x, n1-1, n2-1, n3-1, 1, 1, 1)
        pr_EDT_z=0.d0
        pr_EDT_y=0.d0
        pp3=0.d0
        call transpose_x_to_y(pr_EDT_x, pr_EDT_y, ph_111)
        call transpose_y_to_z(pr_EDT_y, pr_EDT_z, ph_111)

    end subroutine solve_Poisson

    subroutine solve_Poisson_ibm(RHS_z_ibm, pr_EDT_z_ibm)
        use IBM_data, only:n1_ibm,n2_ibm,n3_ibm
        USE decomp_2d
        USE poisson_generic_solver
        implicit none
        real(mytype), dimension(:,:,:), intent(in)      :: RHS_z_ibm
        real(mytype), dimension(:,:,:), intent(out)     :: pr_EDT_z_ibm

        integer     :: i,j,k
        real(mytype), dimension(:,:,:), allocatable, save :: pr_EDT_x_ibm, pr_EDT_y_ibm, pp3_ibm


        TYPE(DECOMP_INFO), save :: ph_111_ibm  ! decomposition object
        TYPE(DECOMP_INFO), save :: ph_000_ibm  ! decomposition object

        pr_EDT_z_ibm=RHS_z_ibm

        if (.not. allocated(pp3_ibm)) then

            call decomp_info_init(n1_ibm, n2_ibm, n3_ibm, ph_111_ibm)
            call decomp_info_init(n1_ibm-1,n2_ibm-1,n3_ibm-1,ph_000_ibm)

            allocate (pr_EDT_x_ibm(ph_111_ibm%xst(1):ph_111_ibm%xen(1),ph_111_ibm%xst(2):ph_111_ibm%xen(2),ph_111_ibm%xst(3):ph_111_ibm%xen(3)))
            allocate (pr_EDT_y_ibm(ph_111_ibm%yst(1):ph_111_ibm%yen(1),ph_111_ibm%yst(2):ph_111_ibm%yen(2),ph_111_ibm%yst(3):ph_111_ibm%yen(3)))
            allocate (pp3_ibm(ph_000_ibm%zst(1):ph_000_ibm%zen(1),ph_000_ibm%zst(2):ph_000_ibm%zen(2),ph_000_ibm%zst(3):ph_000_ibm%zen(3)))

        end if

        call transpose_z_to_y(pr_EDT_z_ibm, pr_EDT_y_ibm, ph_111_ibm)
        call transpose_y_to_x(pr_EDT_y_ibm, pr_EDT_x_ibm, ph_111_ibm)
        call convert_111_to_000(pr_EDT_x_ibm, pp3_ibm, n1_ibm-1, n2_ibm-1, n3_ibm-1, 1, 1, 1)

        call decomp_2d_poisson_stg(pp3_ibm)

        call convert_000_to_111(pp3_ibm, pr_EDT_x_ibm, n1_ibm-1, n2_ibm-1, n3_ibm-1, 1, 1, 1)
        pr_EDT_z_ibm=0.d0
        pr_EDT_y_ibm=0.d0
        pp3_ibm=0.d0
        call transpose_x_to_y(pr_EDT_x_ibm, pr_EDT_y_ibm, ph_111_ibm)
        call transpose_y_to_z(pr_EDT_y_ibm, pr_EDT_z_ibm, ph_111_ibm)

    end subroutine solve_Poisson_ibm

end module poisson_interface
