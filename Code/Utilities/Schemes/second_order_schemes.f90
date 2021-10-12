module second_order_schemes

contains

    subroutine compute_backward_diff_order_1(f0, f1, f2, d1f, X, nend)
        implicit none
        real*8                              :: tmpR1, alpha, dzN
        real*8                              :: b1d2c0, b1d2c1, b1d2c2
        integer, intent(in)                 :: nend
        real*8, intent(in)                  :: f0, f1, f2
        real*8, dimension(:), intent(in)    :: X
        real*8, intent(inout)               :: d1f

        !Backward Derivative
        dzN = X(nend)-X(nend-1)
        alpha = (X(nend-1)-X(nend-2))/dzN

        ! Backward first derivative second-order accurate
        tmpR1 = dzN*alpha*(alpha + 1.d0)
        b1d2c0 = ((alpha + 1.d0)**2.d0 - 1.d0)/tmpR1
        b1d2c1 = -((alpha + 1.d0)**2.d0)/tmpR1  
        b1d2c2 = 1.d0 / tmpR1

        d1f = b1d2c0 * f0 + b1d2c1 * f1 + b1d2c2 * f2

    end subroutine

    subroutine compute_backward_diff_order_2(f0, f1, f2, f3, d2f, X, nend)
        implicit none
        real*8                              :: db12, db23, db13, db24, db14, db34
        real*8                              :: b2d2c0, b2d2c1, b2d2c2, b2d2c3
        integer, intent(in)                 :: nend
        real*8, intent(in)                  :: f0, f1, f2, f3
        real*8, dimension(:), intent(in)    :: X
        real*8, intent(inout)               :: d2f

        !Backward Derivative
        db12 = X(nend)-X(nend-1)
        db23 =X(nend-1)-X(nend-2) 
        db13 = X(nend)-X(nend-2)
        db24 =X(nend-1)-X(nend-3)
        db14 = X(nend)-X(nend-3)
        db34 =X(nend-2)-X(nend-3) 

        ! Backward second derivative second-order accurate
        b2d2c0 = (2.d0*(db12 + db13 + db14))/(  db12 *db13*db14      )        
        b2d2c1 = (2.d0*(db13 + db14       ))/((-db12)*db23*db24      )
        b2d2c2 = (2.d0*(db12 + db14       ))/((-db13)*(-db23)*db34   )
        b2d2c3 = (2.d0*(db12 + db13       ))/((-db14)*(-db24)*(-db34))

        d2f = b2d2c0 * f0 + b2d2c1 * f1 + b2d2c2 * f2 + b2d2c3 * f3

    end subroutine

    subroutine compute_upward_diff_order_1(f0, f1, f2, d1f, X, nstart)
        implicit none
        real*8                              :: tmpR1, alpha, dz1
        real*8                              :: f1d2c0, f1d2c1, f1d2c2
        integer, intent(in)                 :: nstart
        real*8, intent(in)                  :: f0, f1, f2
        real*8, dimension(:), intent(in)    :: X
        real*8, intent(inout)               :: d1f

        ! Forward Derivatives    
        dz1 = X(nstart+1)-X(nstart)
        alpha = (X(nstart+2)-X(nstart+1))/dz1

        ! Forward first derivative second-order accurate
        tmpR1 = dz1*alpha*(alpha + 1.d0)
        f1d2c0 = (1.d0 - (alpha + 1.d0)**2.d0)/tmpR1
        f1d2c1 = ((alpha + 1.d0)**2.d0)/tmpR1
        f1d2c2 = -1.d0 / tmpR1

        d1f = f1d2c0 * f0 + f1d2c1 * f1 + f1d2c2 * f2

    end subroutine

    subroutine compute_upward_diff_order_2(f0, f1, f2, f3, d2f, X, nstart)
        implicit none
        real*8                              :: df21, df32, df31, df42, df41, df43
        real*8                              :: f2d2c0, f2d2c1, f2d2c2, f2d2c3
        integer, intent(in)                 :: nstart
        real*8, intent(in)                  :: f0, f1, f2, f3
        real*8, dimension(:), intent(in)    :: X
        real*8, intent(inout)               :: d2f

        !Upward Derivative
        df21 = X(nstart+1)-X(nstart)
        df32 = X(nstart+2)-X(nstart+1)   
        df31 = X(nstart+2)-X(nstart)
        df42 = X(nstart+3)-X(nstart+1)
        df41 = X(nstart+3)-X(nstart)
        df43 = X(nstart+3)-X(nstart+2)

        ! Backward second derivative second-order accurate
        f2d2c0 = (2.d0*(-df21-df31-df41))/((-df21)*(-df31)*(-df41))         
        f2d2c1 = (2.d0*(-df31-df41)     )/(  df21 *(-df32)*(-df42)) 
        f2d2c2 = (2.d0*(-df21-df41)     )/(  df31 *  df32 *(-df43))
        f2d2c3 = (2.d0*(-df21-df31)     )/(  df41 *  df42 *  df43 )

        d2f = f2d2c0 * f0 + f2d2c1 * f1 + f2d2c2 * f2 + f2d2c3 * f3

    end subroutine

end module second_order_schemes