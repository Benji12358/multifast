module d1s_schemes_ibm
use boundaries_types
use schemes_settings

implicit none

contains

    subroutine D1s_ExpCtr_O2Fp0_ibm(f, d1f, n_max, h, shifted, n_start_ibm, n_end_ibm, n_objects)

        implicit none
        integer, intent(in)                 :: n_max, n_objects, n_start_ibm(:), n_end_ibm(:)
        logical, intent(in)                 :: shifted
        real(kind=8), intent(in)            :: h, f(:)
        real(kind=8), intent(out)           :: d1f(:)

        integer :: i, n, ith, ith_previous, ith_following
        real(kind=8) A

        A=1.d0/h

        if (n_objects.gt.0) then

            ! First, we apply the classic D1s scheme
            ! Then, we correct the node upstream the object

            if (shifted) then

                do i=2,n_max-1
                    d1f(i)= ( f(i) - f(i-1) ) * A
                enddo

                d1f(1)= ( f(1) - f(n_max-1) ) * A

                do n = 1,n_objects

                    ! Basically, we apply D1s = ( f(i-1) - f(i) )/h
                    ! This is applied on 4 nodes upstream the object

                    ith = n_start_ibm(n)
                    ith_previous = ith - 1
                    if (ith_previous.le.0) ith_previous=n_max-1

                    do i=1,4

                        d1f(ith)= ( f(ith_previous) - f(ith) ) * A

                        ith = ith-1
                        ith_previous = ith_previous-1
                        if (ith.le.0) ith=n_max-1

                    enddo

                    do i=n_start_ibm(n),min(n_end_ibm(n),n_max)

                        d1f(i)= 0

                    enddo

                enddo

            else

                do i=1,n_max-2
                    d1f(i)= ( f(i+1) - f(i) ) * A
                enddo

                d1f(n_max-1)= ( f(1) - f(n_max-1) ) * A

                do n = 1,n_objects

                    ! Basically, we apply D1s = ( f(i) - f(i+1) )/h
                    ! This is applied on 4 nodes upstream the object

                    ith = n_start_ibm(n)
                    ith_following = ith + 1

                    do i=1,4

                        d1f(ith)= ( f(ith) - f(ith_following) ) * A

                        ith = ith-1
                        ith_following = ith_following-1
                        if (ith.le.0) ith=n_max-1
                        if (ith_following.le.0) ith_following=n_max-1

                    enddo

                    do i=n_start_ibm(n),min(n_end_ibm(n),n_max)

                        d1f(i)= 0

                    enddo

                enddo

            end if

        else

            ! There is no object in the domain
            ! Thus, we apply the classic D1s scheme

            if (shifted) then

                do i=2,n_max-1
                    d1f(i)= ( f(i) - f(i-1) ) * A
                enddo

                d1f(1)= ( f(1) - f(n_max-1) ) * A

            else

                do i=1,n_max-2
                    d1f(i)= ( f(i+1) - f(i) ) * A
                enddo

                d1f(n_max-1)= ( f(1) - f(n_max-1) ) * A

            end if

        endif

        return

    end subroutine

end module