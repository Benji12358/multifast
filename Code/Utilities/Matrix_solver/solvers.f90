
module LU_solvers

contains

    subroutine LU_solver_Non_Periodic(A, B, C, d, x, n_max)

        implicit none

        real(kind=8),dimension (:), intent(in)    :: A, B, C, d
        real(kind=8), dimension(:), intent(out)   :: x

        real(kind=8),dimension (n_max)       :: btj,gmkj
        integer i,n_max,je

        je=n_max-1
        ! LU factorisation -----------------------------------------
        btj(1)=B(1)
        do i=2, n_max
            gmkj(i)=C(i-1)/btj(i-1)
            btj(i)=B(i)-A(i)*gmkj(i)
        enddo

        ! Forward substitution -------------------------------------
        x(1)=d(1)/btj(1)

        do i=2, n_max
            x(i)=(d(i)-A(i)*x(i-1))/btj(i)
        enddo

        ! Backward substitution ------------------------------------

        do i=je,1,-1
            x(i)=x(i)-gmkj(i+1)*x(i+1)
        enddo

        return
    end subroutine

    subroutine LU_solver_Periodic(A,B,C, d, x, n_max)

        ! vectorized for right hand side and coefficients

        implicit none

        real(kind=8), dimension(:), intent(in)    ::  A, B, C
        real(kind=8), dimension(:), intent(in)    :: d
        real(kind=8), dimension(:), intent(out)   :: x

        real(kind=8), dimension(n_max)    ::  fei, ppi, ssi, qqi
        integer i, n_max

        ! Initialization ----------------------------------
        qqi(1) = -C(1)/B(1)
        ssi(1) = -A(1)/B(1)

        ! 1st Forward elimination sweep -------------------

        do i=2, n_max-1
            ppi(i) =1.d0/( B(i) + A(i)*qqi(i-1))
            qqi(i) = - C(i)*ppi(i)
            ssi(i) = - A(i)*ssi(i-1)*ppi(i)
        enddo


        ! 1st Backward substitution --------------------------
        ssi(n_max) = 1.d0
        do i=n_max-1, 1, -1
            ssi(i) = ssi(i) + qqi(i)*ssi(i+1)
        enddo
        ppi(n_max)=1.d0/(C(n_max)*ssi(1)+A(n_max)*ssi(n_max-1)+B(n_max))


        fei(n_max) = 0.d0

        ! 2nd Forward elimination sweep -------------------
        x(1) = d(1)/B(1)

        do i=2, n_max-1
            x(i) = (d(i) - A(i)*x(i-1))*ppi(i)
        enddo


        do i=n_max-1, 1, -1
            fei(i) = x(i) + qqi(i)*fei(i+1)
        enddo

        x(n_max)  =   (x(n_max)-C(n_max)*fei(1) - A(n_max)*fei(n_max-1))*ppi(n_max)

        ! 3. backward elimination pass ---------------------

        do i=n_max-1,1,-1
            x(i)=x(n_max)*ssi(i) + fei(i)
        enddo

        return

    end subroutine

end module LU_solvers

module Lapack_wrappers
    implicit none

contains

    subroutine create_band_storage_matrix(matrix, matrix_bs, N, kl, ku)

        implicit none
        integer                 ::  N, kl, ku, i, j, je, je_min, je_max, i_max, i_min
        integer                 :: dkl
        real(kind=8), dimension(:,:)  :: matrix
        real(kind=8), dimension(:,:)  :: matrix_bs

        matrix_bs=0.d0
        dkl=kl+1

        do i=1, N
            do j = -kl, ku
                je=j+i
                je_min=1
                je_max=N
                if (je.ge.je_min.and.je.le.je_max) then
                    i_min=max(1,je-ku)
                    i_max=min(N,je+kl)
                    if (i.ge.(i_min).and.i.le.(i_max)) then
                        matrix_bs(kl+ku+1+i-je,je) = matrix(i,j+dkl)
                    endif

                end if

            end do
        end do

    end subroutine

    subroutine Lapack_solver(A, B, C, d, x, n_max)

        implicit none

        integer, intent(in)                 :: n_max
        real(kind=8),dimension (:), intent(in)    :: d, A, B, C
        real(kind=8),dimension (:), intent(out)   :: x

        real(kind=8),dimension (n_max)    :: ipiv, du2, rhs
        integer :: info

        call DGTTRF( n_max, A(2:n_max), B, C(1:n_max-1), du2(1:n_max-2), ipiv, info )

        x=d
        call DGTTRS( 'N', n_max, 1, A(2:n_max), B, C(1:n_max-1), du2(1:n_max-2), ipiv, x, n_max, info )

    end subroutine

    subroutine Lapack_band_solver(A, kl, ku, n_max, sol)

        implicit none

        integer, intent(in)                     :: n_max, kl, ku
        real(kind=8), dimension(:, :), intent(in)     :: A
        real(kind=8),dimension (:), intent(out)       :: sol

        integer :: info, i
        real(kind=8),dimension (n_max) :: ipiv, du2, rhs
        real(kind=8), dimension(n_max, n_max)  :: AB

        call create_band_storage_matrix(A, AB(1:2*kl+ku+1,:), n_max, kl, ku)

        call DGBTRF( n_max, n_max, kl, ku, AB, n_max, ipiv, info)

        call DGBTRS( 'N', n_max, kl, ku, 1, AB, n_max, ipiv, sol, n_max, info )

    end subroutine

    subroutine Lapack_Generic_solver(A, n_max, sol)


        implicit none

        real(kind=8), dimension(:, :), intent(in) :: A
        integer, intent(in)                 :: n_max
        real(kind=8),dimension (:), intent(out)   :: sol

        integer :: info
        real(kind=8),dimension (n_max) :: ipiv, du2, rhs

        call DGETRF( n_max, n_max, A, n_max, ipiv, info )

        call DGETRS( 'N', n_max, 1, A, n_max, ipiv, sol, n_max, info )

    end subroutine

end module Lapack_wrappers

subroutine penta_solver(E, A, D, C, F, B, X, n_max)

    !   RESULTS:  matrix has 5 bands, EADCF, with D being the main diagonal,
    !   E and A are the lower diagonals, and C and F are the upper diagonals.

    !     E is defined for rows i = 3:N, but is defined as E(1) to E(N-2)
    !     A is defined for rows i = 2:N, but is defined as A(1) to A(N-1)
    !     D is defined for rows i = 1:N
    !     C is defined for rows i = 1:N-1, but the last element isn't used
    !     F is defined for rows i = 1:N-2, but the last 2 elements aren't used

    !   B is the right-hand side
    !   X is the solution vector

    implicit none
    integer i,n_max
    double precision E(n_max),A(n_max),D(n_max),C(n_max),F(n_max),B(n_max),X(n_max),XMULT
    do i = 2,n_max-1
        XMULT = A(i-1)/D(i-1)
        D(i) = D(i) - XMULT*C(i-1)
        C(i) = C(i) - XMULT*F(i-1)
        B(i) = B(i) - XMULT*B(i-1)
        XMULT = E(i-1)/D(i-1)
        A(i) = A(i) - XMULT*C(i-1)
        D(i+1) = D(i+1) - XMULT*F(i-1)
        B(i+1) = B(i+1) - XMULT*B(i-1)
    enddo
    XMULT = A(n_max-1)/D(n_max-1)
    D(n_max) = D(n_max) - XMULT*C(n_max-1)
    X(n_max) = (B(n_max) - XMULT*B(n_max-1))/D(n_max)
    X(n_max-1) = (B(n_max-1) - C(n_max-1)*X(n_max))/D(n_max-1)
    do i = n_max-2,1,-1
        X(i) = (B(i) - F(i)*X(i+2) - C(i)*X(i+1))/D(i)
    enddo
    return
end




module Thomas_solver
    implicit none
    real(kind=8), allocatable, dimension(:) :: m, bp, upper

contains

    subroutine TS_init(A, B, C, n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer,intent(in)                      :: n_max
        real(kind=8), dimension(:), intent(in)       :: A, B, C

        integer i

        if (allocated(m)) then
            deallocate(m, bp, upper)
        end if

        allocate(m(n_max), bp(n_max), upper(n_max))

        bp(1) = 1.d0/B(1)
        upper(1)=C(1)

        !   The first pass (setting coefficients):
        do i = 2,n_max
            m(i) = A(i)*bp(i-1)
            bp(i) = 1.d0/(B(i) - m(i)*C(i-1))
            upper(i)=C(i)
        enddo

    end subroutine

    subroutine TS_NPr(v, x, n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer,intent(in)                  :: n_max
        real(kind=8),dimension(:), intent(in)    :: v
        real(kind=8),dimension(:), intent(out)   :: x
        real(kind=8),dimension(n_max)                :: vp
        real*8  :: tmp

        integer i

        !   Make copies of the b and v variables so that they are
        !   unaltered by this sub
        vp(1) = v(1)

        tmp=vp(1)
        do i = 2,n_max
            vp(i) = v(i) - m(i)*tmp
            tmp=vp(i)
        enddo

        x(n_max) = vp(n_max)*bp(n_max)

        !   The second pass (back-substition)

        tmp=x(n_max)
        do i = n_max-1, 1, -1
            x(i) = (vp(i) - upper(i)*tmp)*bp(i)
            tmp=x(i)
        enddo

    end subroutine

    subroutine TS_Pr(A,B,C,d,x,n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer, intent(in)                     :: n_max
        real(kind=8), dimension(:), intent(in)       :: d
        real(kind=8), dimension(:), intent(out)      :: x
        real(kind=8), dimension(:), intent(inout)    :: A, B, C
        real(kind=8), dimension(n_max)               :: bp, u, v, vp, q, y
        real(kind=8)                                 :: vy, vq, coef_q

        real(kind=8)  :: B1, B_nmax

        integer i

        u=0.d0
        v=0.d0

        u(1)=0.9d0*B(1)
        u(n_max)=C(n_max)
        v(1)=1.d0
        v(n_max)=A(1)/u(1)
        vy=0.d0
        vq=0.d0

        B_nmax=B(n_max)
        B1=B(1)

        B(n_max)=B(n_max) - u(n_max)*v(n_max)
        B(1)=B(1) - u(1)

        call TS_init(A, B, C, n_max)

        call TS_NPr(d, y, n_max)
        call TS_NPr(u, q, n_max)

        do i=1, n_max
            vy=vy+v(i)*y(i)
            vq=vq+v(i)*q(i)
        enddo
        coef_q = vy/(1+vq)

        do i=1, n_max
            x(i)=y(i)-coef_q*q(i)
        enddo

        B(n_max)=B_nmax
        B(1)=B1

    end subroutine

    subroutine TS_init_Pr(A,B,C,u,v,vq,n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer, intent(in)                         :: n_max
        real(kind=8), dimension(:), intent(inout)   :: A, B, C
        real(kind=8), dimension(:), intent(inout)   :: u, v
        real(kind=8), intent(out)                   :: vq
        real(kind=8), dimension(n_max)              :: bp, vp, q, y
        real(kind=8)                                :: vy, coef_q

        real(kind=8)  :: B1, B_nmax

        integer i

        u=0.d0
        v=0.d0
        vq=0.d0

        u(1)=0.9d0*B(1)
        u(n_max)=C(n_max)
        v(1)=1.d0
        v(n_max)=A(1)/u(1)

        B_nmax=B(n_max)
        B1=B(1)

        B(n_max)=B(n_max) - u(n_max)*v(n_max)
        B(1)=B(1) - u(1)

        call TS_init(A, B, C, n_max)
        call TS_NPr(u, u, n_max)

        do i=1, n_max
            vq=vq+v(i)*u(i)
        enddo

    end subroutine

    subroutine TS_Pr_opt(A,B,C,d,x,u,v,vq,n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer, intent(in)                         :: n_max
        real(kind=8), dimension(:), intent(in)      :: d
        real(kind=8), dimension(:), intent(out)     :: x
        real(kind=8), dimension(:), intent(inout)   :: A, B, C
        real(kind=8), intent(out)                   :: vq
        real(kind=8), dimension(n_max)              :: bp, u, v, vp, q, y
        real(kind=8)                                :: vy, coef_q

        real(kind=8)  :: B1, B_nmax

        integer i
        vy=0.d0

        call TS_NPr(d, x, n_max)

        do i=1, n_max
            vy=vy+v(i)*x(i)
        enddo
        coef_q = vy/(1+vq)

        do i=1, n_max
            x(i)=x(i)-coef_q*u(i)
        enddo

    end subroutine

#ifdef IFORT

    subroutine TS_Pr2(A,B,C,d,x,n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer, intent(in)                     :: n_max
        real(kind=8), dimension(:), intent(in)       :: d
        real(kind=8), dimension(:), intent(out)      :: x
        real(kind=8), dimension(:), intent(inout)    :: A, B, C
        real(kind=8), dimension(n_max)               :: bp, u, v, vp, q, y
        real(kind=8)                                 :: vy, vq, coef_q
        integer     :: info

        real(kind=8)  :: B1, B_nmax

        integer i

        u(1)=0.9d0*B(1)
        u(n_max)=C(n_max)
        v(1)=1.d0
        v(n_max)=A(1)/u(1)
        vy=0.d0
        vq=0.d0

        B_nmax=B(n_max)
        B1=B(1)

        B(n_max)=B(n_max) - u(n_max)*v(n_max)
        B(1)=B(1) - u(1)

        !call TS_init(A, B, C, n_max)
        call ddttrfb( n_max, A, B, C, info)

        do i=2, n_max-1
            u(i)=0.d0
            v(i)=0.d0
        enddo

        !call TS_NPr(d, y, n_max)
        y=d
        call ddttrsb( 'N', n_max, 1, A, B, C, y, n_max, info )

        !call TS_NPr(u, q, n_max)
        q=u
        call ddttrsb( 'N', n_max, 1, A, B, C, q(1:n_max), n_max, info )

        do i=1, n_max
            vy=vy+v(i)*y(i)
            vq=vq+v(i)*q(i)
        enddo
        coef_q = vy/(1+vq)

        do i=1, n_max
            x(i)=y(i)-coef_q*q(i)
        enddo

        B(n_max)=B_nmax
        B(1)=B1

    end subroutine

    subroutine TS_Pr3(A,B,C,d,x,n_max)

        ! a - sub-diagonal (means it is the diagonal below the main diagonal)
        ! b - the main diagonal
        ! c - sup-diagonal (means it is the diagonal above the main diagonal)
        ! v - right part
        ! x - the answer
        ! n_max - number of equations

        implicit none
        integer, intent(in)                     :: n_max
        real(kind=8), dimension(:), intent(in)       :: d
        real(kind=8), dimension(:), intent(out)      :: x
        real(kind=8), dimension(:), intent(inout)    :: A, B, C
        real(kind=8), dimension(n_max)               :: bp, u, v, vp, q, y
        real(kind=8), dimension(n_max,2)               :: yq
        real(kind=8)                                 :: vy, vq, coef_q
        integer     :: info

        real(kind=8)  :: B1, B_nmax

        integer i

        u(1)=0.9d0*B(1)
        u(n_max)=C(n_max)
        v(1)=1.d0
        v(n_max)=A(1)/u(1)
        vy=0.d0
        vq=0.d0

        B_nmax=B(n_max)
        B1=B(1)

        B(n_max)=B(n_max) - u(n_max)*v(n_max)
        B(1)=B(1) - u(1)

        !call TS_init(A, B, C, n_max)
        call ddttrfb( n_max, A, B, C, info)

        do i=2, n_max-1
            u(i)=0.d0
            v(i)=0.d0
        enddo

        !call TS_NPr(d, y, n_max)
        yq(:,1)=d
        yq(:,2)=u
        call ddttrsb( 'N', n_max, 2, A, B, C, yq, n_max, info )

        do i=1, n_max
            vy=vy+v(i)*yq(i,1)
            vq=vq+v(i)*yq(i,2)
        enddo
        coef_q = vy/(1+vq)

        do i=1, n_max
            x(i)=yq(i,1)-coef_q*yq(i,2)
        enddo

        B(n_max)=B_nmax
        B(1)=B1

    end subroutine

#endif


end module Thomas_solver

