




module O2_D0s
    use boundaries_types
    implicit none

contains

    subroutine D0s_O2_3Dx(f, ff, n1,n2,n2e,n3,n3e, bc_type)

        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: ff

        integer :: i,j,k

        do k=1,n3e
            do j=1,n2e
                do i=1,n1-2
                    ff(i, j, k)= 0.5d0*(f(i+1, j, k) + f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e
                    ff(n1-1, j, k)= 0.5d0*(f(1, j, k) + f(n1-1, j, k))
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    ff(n1-1,j,k)=0.5d0*(f(n1,j,k)+f(n1-1,j,k))
                enddo
            enddo

        endif

        return

    end subroutine D0s_O2_3Dx



    subroutine D0s_O2_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k
        real(kind=8) A

        A =  1.d0       /2.d0

        do k=1, n3e
            do i=1, n1e
                do j=1,n2-2
                    ff(i, j, k)= A*(f(i, j+1, k) + f(i, j, k))
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e
                    ff(i, n2-1, k)= A*(f(i, 1, k) + f(i, n2-1, k))
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1, n3e
                do i=1, n1e
                    ff(i,n2-1,k)=0.5d0*(f(i,n2,k)+f(i,n2-1,k))
                enddo
            enddo

        endif

        return

    end subroutine D0s_O2_3Dy

    subroutine D0s_O2_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        integer :: i,j,k
        real(kind=8) A

        A =  1.d0     /2.d0

        do j=1, n2e
            do i=1, n1e
                do k=1,n3-2
                    ff(i, j, k)= A*(f(i, j, k+1) + f(i, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    ff(i, j, n3-1)= A*(f(i, j, 1) + f(i, j, n3-1))
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do j=1, n2e
                do i=1, n1e
                    ff(i,j,n3-1)=0.5d0*(f(i,j,n3)+f(i,j,n3-1))
                enddo
            enddo

        endif

        return

    end subroutine D0s_O2_3Dz

end module O2_D0s

module O2_D0s_sh
    use boundaries_types
    implicit none

contains

    subroutine D0ssh_O2_3Dx(f, ff, n1,n2,n2e,n3,n3e, bc_type)
        implicit none


        integer, intent(in) :: n1,n2,n2e,n3,n3e, bc_type
        real*8, dimension(n1,n2,n3), intent(in):: f
        real*8, dimension(n1,n2,n3), intent(out):: ff

        integer :: i,j,k

        do k=1,n3e
            do j=1,n2e
                do i=2,n1-1
                    ff(i, j, k)= 0.5d0*(f(i, j, k) + f(i-1, j, k))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do k=1,n3e
                do j=1,n2e
                    ff(1, j, k)=   0.5d0*(f(1, j, k) + f(n1-1, j, k))
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then

            do k=1,n3e
                do j=1,n2e
                    ff(2,j,k)=0.5d0*(f(2,j,k)+f(1,j,k))
                enddo
            enddo

        endif

        return
    end subroutine D0ssh_O2_3Dx



    subroutine D0ssh_O2_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     ::n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k
        real(kind=8) A

        A =  1.d0     /2.d0

        do k=1, n3e
            do i=1, n1e
                do j=2,n2-1
                    ff(i, j, k)= A*(f(i, j, k) + f(i, j-1, k))
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e
                    ff(i, 1, k)=   A*(f(i, 1, k) + f(i, n2-1, k))
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then
            ! DO NOTHING
        endif

        return

    end subroutine D0ssh_O2_3Dy



    subroutine D0ssh_O2_MULT_3Dy(f, ff, n1,n1e,n2,n3,n3e, bc_type)

        implicit none

        integer, intent(in)     ::n1,n1e,n2,n3,n3e, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)      :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)     :: ff

        integer :: i,j,k
        real(kind=8) :: A

        A =  1.d0     /2.d0

        do k=1, n3e
            do i=1, n1e
                do j=2,n2-1
                    ff(i, j, k)= A*(f(i, j, k) + f(i, j-1, k))*ff(i, j, k)
                enddo
            enddo
        enddo


        if (bc_type.eq.periodic) then

            do k=1, n3e
                do i=1, n1e
                    ff(i, 1, k)= A*(f(i, 1, k) + f(i, n2-1, k))*ff(i, 1, k)
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then
            ! DO NOTHING
        endif

        return

    end subroutine D0ssh_O2_MULT_3Dy

    subroutine D0ssh_O2_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        integer :: i,j,k
        real(kind=8) A

        A =  1.d0     /2.d0

        do j=1, n2e
            do i=1, n1e
                do k=2,n3-1
                    ff(i, j, k)= A*(f(i, j, k) + f(i, j, k-1))
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    ff(i, j, 1)=   A*(f(i, j, 1) + f(i, j, n3-1))
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then
            ! DO NOTHING
        endif

        return

    end subroutine D0ssh_O2_3Dz

    subroutine D0ssh_O2_MULT_3Dz(f, ff, n1,n1e,n2,n2e,n3, bc_type)

        implicit none

        integer, intent(in)     :: n1,n1e,n2,n2e,n3, bc_type
        real(kind=8), dimension(n1,n2,n3), intent(in)        :: f
        real(kind=8), dimension(n1,n2,n3), intent(out)       :: ff

        integer :: i,j,k
        real(kind=8) A

        A =  1.d0     /2.d0

        do j=1, n2e
            do i=1, n1e
                do k=2,n3-1
                    ff(i, j, k)= (A*(f(i, j, k) + f(i, j, k-1)))*ff(i, j, k)
                enddo
            enddo
        enddo

        if (bc_type.eq.periodic) then

            do j=1, n2e
                do i=1, n1e
                    ff(i, j, 1)=   (A*(f(i, j, 1) + f(i, j, n3-1)))*ff(i, j, 1)
                enddo
            enddo

        endif


        if ((bc_type.eq.Dirichlet)) then
            ! DO NOTHING
        endif

        return

    end subroutine D0ssh_O2_MULT_3Dz

end module O2_D0s_sh
