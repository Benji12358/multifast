





module CPT_INIT
    implicit none

contains

    subroutine TS_initDIR(A, B, C, bp, m, n)

        ! AA - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n - number of equations

        implicit none
        integer,intent(in)                      :: n
        real(8),dimension(n),intent(in)         :: A, B, C
        real(8),dimension(n),intent(out)        :: bp,m
        integer i

        bp(1) = 1.d0/B(1)

        !   The first pass (setting coefficients):
        do i = 2,n
            m(i) = A(i)*bp(i-1)
            bp(i) = 1.d0/(B(i) - m(i)*C(i-1))
        enddo

    end subroutine


    subroutine TS_initPR(A, B, C, bp, cp, m, q,v,vq, n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer,intent(in)                      :: n_max
        real(kind=8), dimension(:), intent(inout)       :: A, B, C, bp, cp, m, q,v
        real*8  :: vq

        integer i

        q=0.d0
        v=0.d0
        vq=0.d0

        q(1)=2.d0
        q(n_max)=C(n_max)
        v(1)=1.d0
        v(n_max)=A(1)/q(1)

        B(n_max)=B(n_max) - q(n_max)*v(n_max)
        B(1)=B(1) - q(1)

        bp(1) = 1.d0/B(1)
        cp(1)=C(1)

        !   The first pass (setting coefficients):
        do i = 2,n_max
            m(i) = A(i)*bp(i-1)
            bp(i) = 1.d0/(B(i) - m(i)*C(i-1))
            cp(i)=C(i)
        enddo

        do i = 2,n_max
            q(i)=q(i) - m(i)*q(i-1)
        enddo

        q(n_max) = q(n_max)*bp(n_max)

        !   The second pass (back-substition)
        do i = n_max-1, 1, -1
            q(i)=(q(i) - cp(i)*q(i+1))*bp(i)
        enddo


        do i=1, n_max
            vq=vq+v(i)*q(i)
        enddo

    end subroutine

end module CPT_INIT




