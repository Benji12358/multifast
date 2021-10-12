module mod_phase_average

  use mod_init

  implicit none


  contains


    subroutine perform_phase_averaging(var1s_x, var2s_y, var3s_z, var4s_x, w_wall, v_wall, u_wall)

    !---------------------------------------------------------------
    ! This subroutine performs appropriate transposition of arrays
    ! which is required before calling the phase averaging 
    ! subroutine depending on the type of boundary condition.
    ! This subroutine is called in the multiphysics module.
    !---------------------------------------------------------------
        implicit none

        real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: var1s_x,var4s_x
        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), intent(in) :: var2s_y
        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)), intent(in) :: var3s_z
        
        real*8, dimension(ystart(1):yend(1),ystart(3):yend(3)), intent(in) :: w_wall, v_wall, u_wall
        real*8, dimension(ystart(1):yend(1),ystart(3):yend(3)) :: p_wall

        real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: U_x, V_x, W_x, P_x
        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: U_y, V_y, W_y, P_y
        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: U_z, V_z, W_z, P_z

        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dUdx_y, dUdy_y, dUdz_y
        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dVdx_y, dVdy_y, dVdz_y
        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dWdx_y, dWdy_y, dWdz_y

        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dUdx_z, dUdy_z, dUdz_z
        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dVdx_z, dVdy_z, dVdz_z
        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dWdx_z, dWdy_z, dWdz_z

        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dPdx_y, dPdy_y, dPdz_y
        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dPdx_z, dPdy_z, dPdz_z
        
        real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dUUdx_y, dVVdx_y, dWWdx_y, dUVdx_y
        real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dUUdx_z, dVVdx_z, dWWdx_z, dUVdx_z

          
        integer :: i,j,k
        character(len=300) :: fnam
        character(len=300) :: path
        
        ! Local variables
        integer :: a,b,c,d,e,f
        integer :: nw,bin
        real*8  :: time,Tref
        real*8  :: avg
        real*8  :: Xi

        call myMPI()

        !Interpolating velocity components to the cell centers
        !Needed to compute products of cross quantities and derivatives, etc
        call perform_velocity_at_center(var3s_z, var2s_y, var1s_x, U_z, V_y, W_x)
        P_x = var4s_x ! Pressure is already located at the cell center

        call spread_to_all_configs(W_x,'X',W_y,W_z)
        call spread_to_all_configs(V_y,'Y',V_x,V_z)
        call spread_to_all_configs(U_z,'Z',U_y,U_x)
        call spread_to_all_configs(P_x,'X',P_y,P_z)


        !------------------------------------
        
        call x_derivative(W_y,dWdx_y)
        call x_derivative(V_y,dVdx_y)
        call x_derivative(U_y,dUdx_y)
        call x_derivative(P_y,dPdx_y)

        call transpose_y_to_z(dWdx_y,dWdx_z)
        call transpose_y_to_z(dVdx_y,dVdx_z)
        call transpose_y_to_z(dUdx_y,dUdx_z)
        call transpose_y_to_z(dPdx_y,dPdx_z)
        
        !------------------------------------
        
        call compute_P_wall(V_y,P_y,p_wall)
        
        call y_derivative(W_y,dWdy_y,w_wall)
        call y_derivative(V_y,dVdy_y,v_wall)
        call y_derivative(U_y,dUdy_y,u_wall)        
        call y_derivative(P_y,dPdy_y,p_wall)

        call transpose_y_to_z(dWdy_y,dWdy_z)
        call transpose_y_to_z(dVdy_y,dVdy_z)
        call transpose_y_to_z(dUdy_y,dUdy_z)
        call transpose_y_to_z(dPdy_y,dPdy_z)

        !------------------------------------
        
        call z_derivative(W_y,dWdz_y)
        call z_derivative(V_y,dVdz_y)
        call z_derivative(U_y,dUdz_y)
        call z_derivative(P_y,dPdz_y)

        call transpose_y_to_z(dWdz_y,dWdz_z)
        call transpose_y_to_z(dVdz_y,dVdz_z)
        call transpose_y_to_z(dUdz_y,dUdz_z)
        call transpose_y_to_z(dPdz_y,dPdz_z)
        
        !------------------------------------
        
        call x_derivative(W_y*W_y,dWWdx_y)
        call x_derivative(V_y*V_y,dVVdx_y)
        call x_derivative(U_y*U_y,dUUdx_y)
        call x_derivative(U_y*V_y,dUVdx_y)

        call transpose_y_to_z(dWWdx_y,dWWdx_z)
        call transpose_y_to_z(dVVdx_y,dVVdx_z)
        call transpose_y_to_z(dUUdx_y,dUUdx_z)
        call transpose_y_to_z(dUVdx_y,dUVdx_z)
        
        
        !------------------------------------


        a=zstart(1); b=min(zend(1),n1m)
        c=zstart(2); d=min(zend(2),n2m)
        e=zstart(3); f=min(zend(3),n3m)

        !To move along streamwise direction
        do i=e,f

          Xi = Xc(i) - wave_speed*t

          ! Calculate to which cycle this 'Xi' belongs    
          if(Xi .ge. 0.0d0) nw = int(floor(Xi/lambda) + 1)
          if(Xi .lt. 0.0d0) nw = int(floor(Xi/lambda))

          ! Center the 'Xi' between 0 and lambda
          if(Xi .ge. 0.0d0) Xi = Xi-(nw-1)*lambda
          if(Xi .lt. 0.0d0) Xi = Xi-(nw+1)*lambda

          ! Calculate to which bin this 'xref' belongs
          if(Xi .ge. 0.0d0) bin = int(floor(Xi/delX) + 1)
          if(Xi .lt. 0.0d0) bin = xbins + int(floor(Xi/delX) + 1)      
             


          counts(a:b,c:d,bin) = counts(a:b,c:d,bin) + 1

          U(a:b,c:d,bin) = U(a:b,c:d,bin) + U_z(a:b,c:d,i)
          V(a:b,c:d,bin) = V(a:b,c:d,bin) + V_z(a:b,c:d,i)
          W(a:b,c:d,bin) = W(a:b,c:d,bin) + W_z(a:b,c:d,i)
          P(a:b,c:d,bin) = P(a:b,c:d,bin) + P_z(a:b,c:d,i)


          UU(a:b,c:d,bin) = UU(a:b,c:d,bin) + U_z(a:b,c:d,i)*U_z(a:b,c:d,i)
          VV(a:b,c:d,bin) = VV(a:b,c:d,bin) + V_z(a:b,c:d,i)*V_z(a:b,c:d,i)
          WW(a:b,c:d,bin) = WW(a:b,c:d,bin) + W_z(a:b,c:d,i)*W_z(a:b,c:d,i)
          PP(a:b,c:d,bin) = PP(a:b,c:d,bin) + P_z(a:b,c:d,i)*P_z(a:b,c:d,i)


          UV(a:b,c:d,bin) = UV(a:b,c:d,bin) + U_z(a:b,c:d,i)*V_z(a:b,c:d,i)
          UW(a:b,c:d,bin) = UW(a:b,c:d,bin) + U_z(a:b,c:d,i)*W_z(a:b,c:d,i) 
          UP(a:b,c:d,bin) = UP(a:b,c:d,bin) + U_z(a:b,c:d,i)*P_z(a:b,c:d,i) 
          VW(a:b,c:d,bin) = VW(a:b,c:d,bin) + V_z(a:b,c:d,i)*W_z(a:b,c:d,i) 
          VP(a:b,c:d,bin) = VP(a:b,c:d,bin) + V_z(a:b,c:d,i)*P_z(a:b,c:d,i)
          
          
          UUV(a:b,c:d,bin) = UUV(a:b,c:d,bin) + U_z(a:b,c:d,i)*U_z(a:b,c:d,i)*V_z(a:b,c:d,i)
          VVV(a:b,c:d,bin) = VVV(a:b,c:d,bin) + V_z(a:b,c:d,i)*V_z(a:b,c:d,i)*V_z(a:b,c:d,i)                   
          WWV(a:b,c:d,bin) = WWV(a:b,c:d,bin) + W_z(a:b,c:d,i)*W_z(a:b,c:d,i)*V_z(a:b,c:d,i) 
          UVV(a:b,c:d,bin) = UVV(a:b,c:d,bin) + U_z(a:b,c:d,i)*V_z(a:b,c:d,i)*V_z(a:b,c:d,i) 


          dUdx(a:b,c:d,bin) = dUdx(a:b,c:d,bin) + dUdx_z(a:b,c:d,i) 
          dVdx(a:b,c:d,bin) = dVdx(a:b,c:d,bin) + dVdx_z(a:b,c:d,i) 
          dWdx(a:b,c:d,bin) = dWdx(a:b,c:d,bin) + dWdx_z(a:b,c:d,i) 
          dPdx(a:b,c:d,bin) = dPdx(a:b,c:d,bin) + dPdx_z(a:b,c:d,i) 


          PdUdx(a:b,c:d,bin) = PdUdx(a:b,c:d,bin) + P_z(a:b,c:d,i)*dUdx_z(a:b,c:d,i) 
          PdUdy(a:b,c:d,bin) = PdUdy(a:b,c:d,bin) + P_z(a:b,c:d,i)*dUdy_z(a:b,c:d,i) 
          PdVdx(a:b,c:d,bin) = PdVdx(a:b,c:d,bin) + P_z(a:b,c:d,i)*dVdx_z(a:b,c:d,i) 
          PdVdy(a:b,c:d,bin) = PdVdy(a:b,c:d,bin) + P_z(a:b,c:d,i)*dVdy_z(a:b,c:d,i)  
          PdWdz(a:b,c:d,bin) = PdWdz(a:b,c:d,bin) + P_z(a:b,c:d,i)*dWdz_z(a:b,c:d,i) 


          dUdxdUdx(a:b,c:d,bin) = dUdxdUdx(a:b,c:d,bin) + dUdx_z(a:b,c:d,i)*dUdx_z(a:b,c:d,i) 
          dUdydUdy(a:b,c:d,bin) = dUdydUdy(a:b,c:d,bin) + dUdy_z(a:b,c:d,i)*dUdy_z(a:b,c:d,i) 
          dUdzdUdz(a:b,c:d,bin) = dUdzdUdz(a:b,c:d,bin) + dUdz_z(a:b,c:d,i)*dUdz_z(a:b,c:d,i) 


          dVdxdVdx(a:b,c:d,bin) = dVdxdVdx(a:b,c:d,bin) + dVdx_z(a:b,c:d,i)*dVdx_z(a:b,c:d,i) 
          dVdydVdy(a:b,c:d,bin) = dVdydVdy(a:b,c:d,bin) + dVdy_z(a:b,c:d,i)*dVdy_z(a:b,c:d,i)
          dVdzdVdz(a:b,c:d,bin) = dVdzdVdz(a:b,c:d,bin) + dVdz_z(a:b,c:d,i)*dVdz_z(a:b,c:d,i) 


          dWdxdWdx(a:b,c:d,bin) = dWdxdWdx(a:b,c:d,bin) + dWdx_z(a:b,c:d,i)*dWdx_z(a:b,c:d,i) 
          dWdydWdy(a:b,c:d,bin) = dWdydWdy(a:b,c:d,bin) + dWdy_z(a:b,c:d,i)*dWdy_z(a:b,c:d,i) 
          dWdzdWdz(a:b,c:d,bin) = dWdzdWdz(a:b,c:d,bin) + dWdz_z(a:b,c:d,i)*dWdz_z(a:b,c:d,i) 


          dPdxdPdx(a:b,c:d,bin) = dPdxdPdx(a:b,c:d,bin) + dPdx_z(a:b,c:d,i)*dPdx_z(a:b,c:d,i) 
          dPdydPdy(a:b,c:d,bin) = dPdydPdy(a:b,c:d,bin) + dPdy_z(a:b,c:d,i)*dPdy_z(a:b,c:d,i) 
          dPdzdPdz(a:b,c:d,bin) = dPdzdPdz(a:b,c:d,bin) + dPdz_z(a:b,c:d,i)*dPdz_z(a:b,c:d,i) 


          dUdxdVdx(a:b,c:d,bin) = dUdxdVdx(a:b,c:d,bin) + dUdx_z(a:b,c:d,i)*dVdx_z(a:b,c:d,i) 
          dUdydVdy(a:b,c:d,bin) = dUdydVdy(a:b,c:d,bin) + dUdy_z(a:b,c:d,i)*dVdy_z(a:b,c:d,i)
          dUdzdVdz(a:b,c:d,bin) = dUdzdVdz(a:b,c:d,bin) + dUdz_z(a:b,c:d,i)*dVdz_z(a:b,c:d,i) 

          
          dUUdx(a:b,c:d,bin) = dUUdx(a:b,c:d,bin) + dUUdx_z(a:b,c:d,i) 
          dVVdx(a:b,c:d,bin) = dVVdx(a:b,c:d,bin) + dVVdx_z(a:b,c:d,i)
          dWWdx(a:b,c:d,bin) = dWWdx(a:b,c:d,bin) + dWWdx_z(a:b,c:d,i)
          dUVdx(a:b,c:d,bin) = dUVdx(a:b,c:d,bin) + dUVdx_z(a:b,c:d,i)

          
        end do

          


    end subroutine perform_phase_averaging


     subroutine perform_phase_averaging_vort(var1s_x, var2s_y, var3s_z, w_wall, v_wall, u_wall)

      !---------------------------------------------------------------
      ! This subroutine performs appropriate transposition of arrays
      ! which is required before calling the phase averaging 
      ! subroutine depending on the type of boundary condition.
      ! This subroutine is called in the multiphysics module.
      !---------------------------------------------------------------
      implicit none

      real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)), intent(in) :: var1s_x
      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)), intent(in) :: var2s_y
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)), intent(in) :: var3s_z
      
      real*8, dimension(ystart(1):yend(1),ystart(3):yend(3)), intent(in) :: w_wall, v_wall, u_wall

      real*8, dimension(ystart(1):yend(1),ystart(3):yend(3))  :: OMGX_wall, OMGY_wall, OMGZ_wall

      real*8, dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: U_x, V_x, W_x
      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: U_y, V_y, W_y
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: U_z, V_z, W_z

      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dUdx_y, dUdy_y, dUdz_y
      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dVdx_y, dVdy_y, dVdz_y
      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dWdx_y, dWdy_y, dWdz_y

      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dUdx_z, dUdy_z, dUdz_z
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dVdx_z, dVdy_z, dVdz_z
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dWdx_z, dWdy_z, dWdz_z


      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: omgX_y, omgY_y, omgZ_y
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: omgX_z, omgY_z, omgZ_z

      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: domgXdx_y, domgYdx_y, domgZdx_y
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: domgXdx_z, domgYdx_z, domgZdx_z

      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: domgXdy_y, domgYdy_y, domgZdy_y
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: domgXdy_z, domgYdy_z, domgZdy_z

      real*8, dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: domgXdz_y, domgYdz_y, domgZdz_y
      real*8, dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: domgXdz_z, domgYdz_z, domgZdz_z 


      integer :: i,j,k
      character(len=300) :: fnam
      character(len=300) :: path

      ! Local variables
      integer :: a,b,c,d,flag
      integer :: nw,bin
      real*8  :: time,Tref
      real*8  :: avg
      real*8  :: Xi

      U_x=0.0d0; U_y=0.0d0; U_z=0.0d0
      V_x=0.0d0; V_y=0.0d0; V_z=0.0d0
      W_x=0.0d0; W_y=0.0d0; W_z=0.0d0

      call myMPI()

      !Interpolating velocity components to the cell centers
      !Needed to compute products of cross quantities and derivatives, etc
      call perform_velocity_at_center(var3s_z, var2s_y, var1s_x, U_z, V_y, W_x)

      call spread_to_all_configs(W_x,'X',W_y,W_z)
      call spread_to_all_configs(V_y,'Y',V_x,V_z)
      call spread_to_all_configs(U_z,'Z',U_y,U_x)

      !------------------------------------

      call x_derivative(W_y,dWdx_y)
      call x_derivative(V_y,dVdx_y)
      call x_derivative(U_y,dUdx_y)


      call transpose_y_to_z(dWdx_y,dWdx_z)
      call transpose_y_to_z(dVdx_y,dVdx_z)
      call transpose_y_to_z(dUdx_y,dUdx_z)

      !------------------------------------

      call y_derivative(W_y,dWdy_y,w_wall)
      call y_derivative(V_y,dVdy_y,v_wall)
      call y_derivative(U_y,dUdy_y,u_wall)


      call transpose_y_to_z(dWdy_y,dWdy_z)
      call transpose_y_to_z(dVdy_y,dVdy_z)
      call transpose_y_to_z(dUdy_y,dUdy_z)

      !------------------------------------

      call z_derivative(W_y,dWdz_y)
      call z_derivative(V_y,dVdz_y)
      call z_derivative(U_y,dUdz_y)

      call transpose_y_to_z(dWdz_y,dWdz_z)
      call transpose_y_to_z(dVdz_y,dVdz_z)
      call transpose_y_to_z(dUdz_y,dUdz_z)
         
      !------------------------------------

      if (streamwise==3) then
        call perform_vorticity(W_y, V_y, U_y, omgX_y, omgY_y, omgZ_y)
        omgX_Y = -1.0d0*omgX_Y
        omgY_Y = -1.0d0*omgY_Y
        omgZ_Y = -1.0d0*omgZ_Y
      end if

      call transpose_y_to_z(omgX_y,omgX_z)
      call transpose_y_to_z(omgY_y,omgY_z)
      call transpose_y_to_z(omgZ_y,omgZ_z)
      
      
      call compute_OMG_wall(U_y,V_y,W_y,w_wall,amp,kappa,omega,t,OMGX_wall,OMGY_wall,OMGZ_wall)

      !------------------------------------

      call x_derivative(omgX_y,domgXdx_y)
      call x_derivative(omgY_y,domgYdx_y)
      call x_derivative(omgZ_y,domgZdx_y)

      call transpose_y_to_z(domgXdx_y,domgXdx_z)
      call transpose_y_to_z(domgYdx_y,domgYdx_z)
      call transpose_y_to_z(domgZdx_y,domgZdx_z)

      !------------------------------------

      call y_derivative(omgX_y,domgXdy_y,OMGX_wall)
      call y_derivative(omgY_y,domgYdy_y,OMGY_wall)
      call y_derivative(omgZ_y,domgZdy_y,OMGZ_wall)

      call transpose_y_to_z(domgXdy_y,domgXdy_z)
      call transpose_y_to_z(domgYdy_y,domgYdy_z)
      call transpose_y_to_z(domgZdy_y,domgZdy_z)

      !------------------------------------

      call z_derivative(omgX_y,domgXdz_y)
      call z_derivative(omgY_y,domgYdz_y)
      call z_derivative(omgZ_y,domgZdz_y)

      call transpose_y_to_z(domgXdz_y,domgXdz_z)
      call transpose_y_to_z(domgYdz_y,domgYdz_z)
      call transpose_y_to_z(domgZdz_y,domgZdz_z)
    
      !------------------------------------

      a=zstart(1); b=min(zend(1),n1m)
      c=zstart(2); d=min(zend(2),n2m)

      ! To move along streamwise direction
      do i=zstart(3),min(zend(3),n3m)

        Xi = Xc(i) - wave_speed*t

        ! Calculate to which cycle this 'Xi' belongs    
        if(Xi .ge. 0.0d0) nw = int(floor(Xi/lambda) + 1)
        if(Xi .lt. 0.0d0) nw = int(floor(Xi/lambda))

        ! Center the 'Xi' between 0 and lambda
        if(Xi .ge. 0.0d0) Xi = Xi-(nw-1)*lambda
        if(Xi .lt. 0.0d0) Xi = Xi-(nw+1)*lambda

        ! Calculate to which bin this 'xref' belongs
        if(Xi .ge. 0.0d0) bin = int(floor(Xi/delX) + 1)
        if(Xi .lt. 0.0d0) bin = xbins + int(floor(Xi/delX) + 1)

        if(bin .gt. xbins) bin=1 ! It will not happen anyway but just in case 


        omgX (a:b,c:d,bin) = omgX (a:b,c:d,bin) + omgX_z (a:b,c:d,i) 
        omgY (a:b,c:d,bin) = omgY (a:b,c:d,bin) + omgY_z (a:b,c:d,i) 
        omgZ (a:b,c:d,bin) = omgZ (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) 


        omgXomgX (a:b,c:d,bin) = omgXomgX (a:b,c:d,bin) + omgX_z (a:b,c:d,i)**2
        omgYomgY (a:b,c:d,bin) = omgYomgY (a:b,c:d,bin) + omgY_z (a:b,c:d,i)**2
        omgZomgZ (a:b,c:d,bin) = omgZomgZ (a:b,c:d,bin) + omgZ_z (a:b,c:d,i)**2


        omgXomgY (a:b,c:d,bin) = omgXomgY (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)
        omgXomgZ (a:b,c:d,bin) = omgXomgZ (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)
        omgYomgZ (a:b,c:d,bin) = omgYomgZ (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)


        omgXU (a:b,c:d,bin) = omgXU (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * U_z (a:b,c:d,i) 
        omgXV (a:b,c:d,bin) = omgXV (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * V_z (a:b,c:d,i) 
        omgYU (a:b,c:d,bin) = omgYU (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * U_z (a:b,c:d,i) 
        omgYV (a:b,c:d,bin) = omgYV (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * V_z (a:b,c:d,i)
        omgZU (a:b,c:d,bin) = omgZU (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) * U_z (a:b,c:d,i) 
        omgZV (a:b,c:d,bin) = omgZV (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) * V_z (a:b,c:d,i)


        omgXomgXV (a:b,c:d,bin) = omgXomgXV (a:b,c:d,bin) + V_z (a:b,c:d,i) * omgX_z (a:b,c:d,i)**2
        omgYomgYV (a:b,c:d,bin) = omgYomgYV (a:b,c:d,bin) + V_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)**2
        omgZomgZV (a:b,c:d,bin) = omgZomgZV (a:b,c:d,bin) + V_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)**2
        omgXomgYV (a:b,c:d,bin) = omgXomgYV (a:b,c:d,bin) + V_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)


        domgXdx (a:b,c:d,bin) = domgXdx (a:b,c:d,bin) + domgXdx_z (a:b,c:d,i)
        domgYdx (a:b,c:d,bin) = domgYdx (a:b,c:d,bin) + domgYdx_z (a:b,c:d,i)
        domgZdx (a:b,c:d,bin) = domgZdx (a:b,c:d,bin) + domgZdx_z (a:b,c:d,i)


        omgXdUdx (a:b,c:d,bin) = omgXdUdx (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * dUdx_z (a:b,c:d,i)
        omgXdUdy (a:b,c:d,bin) = omgXdUdy (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * dUdy_z (a:b,c:d,i)
        omgXdUdz (a:b,c:d,bin) = omgXdUdz (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * dUdz_z (a:b,c:d,i)
        omgXdVdx (a:b,c:d,bin) = omgXdVdx (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * dVdx_z (a:b,c:d,i)
        omgXdVdy (a:b,c:d,bin) = omgXdVdy (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * dVdy_z (a:b,c:d,i)
        omgXdVdz (a:b,c:d,bin) = omgXdVdz (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * dVdz_z (a:b,c:d,i)
        omgXdWdx (a:b,c:d,bin) = omgXdWdx (a:b,c:d,bin) + omgX_z (a:b,c:d,i) * dWdx_z (a:b,c:d,i)


        omgYdUdx (a:b,c:d,bin) = omgYdUdx (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * dUdx_z (a:b,c:d,i)
        omgYdUdy (a:b,c:d,bin) = omgYdUdy (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * dUdy_z (a:b,c:d,i)
        omgYdUdz (a:b,c:d,bin) = omgYdUdz (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * dUdz_z (a:b,c:d,i)
        omgYdVdx (a:b,c:d,bin) = omgYdVdx (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * dVdx_z (a:b,c:d,i)
        omgYdVdy (a:b,c:d,bin) = omgYdVdy (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * dVdy_z (a:b,c:d,i)
        omgYdVdz (a:b,c:d,bin) = omgYdVdz (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * dVdz_z (a:b,c:d,i)
        omgYdWdy (a:b,c:d,bin) = omgYdWdy (a:b,c:d,bin) + omgY_z (a:b,c:d,i) * dWdy_z (a:b,c:d,i)


        omgZdUdz (a:b,c:d,bin) = omgZdUdz (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) * dUdz_z (a:b,c:d,i)
        omgZdVdz (a:b,c:d,bin) = omgZdVdz (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) * dVdz_z (a:b,c:d,i)
        omgZdWdx (a:b,c:d,bin) = omgZdWdx (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) * dWdx_z (a:b,c:d,i)
        omgZdWdy (a:b,c:d,bin) = omgZdWdy (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) * dWdy_z (a:b,c:d,i)
        omgZdWdz (a:b,c:d,bin) = omgZdWdz (a:b,c:d,bin) + omgZ_z (a:b,c:d,i) * dWdz_z (a:b,c:d,i)


        omgXomgXdUdx (a:b,c:d,bin) = omgXomgXdUdx (a:b,c:d,bin) + dUdx_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgX_z (a:b,c:d,i)
        omgXomgXdVdx (a:b,c:d,bin) = omgXomgXdVdx (a:b,c:d,bin) + dVdx_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgX_z (a:b,c:d,i)


        omgXomgYdUdx (a:b,c:d,bin) = omgXomgYdUdx (a:b,c:d,bin) + dUdx_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)
        omgXomgYdUdy (a:b,c:d,bin) = omgXomgYdUdy (a:b,c:d,bin) + dUdy_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)
        omgXomgYdVdx (a:b,c:d,bin) = omgXomgYdVdx (a:b,c:d,bin) + dVdx_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)
        omgXomgYdVdy (a:b,c:d,bin) = omgXomgYdVdy (a:b,c:d,bin) + dVdy_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)


        omgXomgZdUdz (a:b,c:d,bin) = omgXomgZdUdz (a:b,c:d,bin) + dUdz_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)
        omgXomgZdVdz (a:b,c:d,bin) = omgXomgZdVdz (a:b,c:d,bin) + dVdz_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)
        omgXomgZdWdx (a:b,c:d,bin) = omgXomgZdWdx (a:b,c:d,bin) + dWdx_z (a:b,c:d,i) * omgX_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)


        omgYomgYdUdy (a:b,c:d,bin) = omgYomgYdUdy (a:b,c:d,bin) + dUdy_z (a:b,c:d,i) * omgY_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)
        omgYomgYdVdy (a:b,c:d,bin) = omgYomgYdVdy (a:b,c:d,bin) + dVdy_z (a:b,c:d,i) * omgY_z (a:b,c:d,i) * omgY_z (a:b,c:d,i)


        omgYomgZdUdz (a:b,c:d,bin) = omgYomgZdUdz (a:b,c:d,bin) + dUdz_z (a:b,c:d,i) * omgY_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)
        omgYomgZdVdz (a:b,c:d,bin) = omgYomgZdVdz (a:b,c:d,bin) + dVdz_z (a:b,c:d,i) * omgY_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)
        omgYomgZdWdy (a:b,c:d,bin) = omgYomgZdWdy (a:b,c:d,bin) + dWdy_z (a:b,c:d,i) * omgY_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)


        omgZomgZdWdz (a:b,c:d,bin) = omgZomgZdWdz (a:b,c:d,bin) + dWdz_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i) * omgZ_z (a:b,c:d,i)


        domgXdxdomgXdx (a:b,c:d,bin) = domgXdxdomgXdx (a:b,c:d,bin) + domgXdx_z (a:b,c:d,i)**2
        domgXdydomgXdy (a:b,c:d,bin) = domgXdydomgXdy (a:b,c:d,bin) + domgXdy_z (a:b,c:d,i)**2
        domgXdzdomgXdz (a:b,c:d,bin) = domgXdzdomgXdz (a:b,c:d,bin) + domgXdz_z (a:b,c:d,i)**2


        domgYdxdomgYdx (a:b,c:d,bin) = domgYdxdomgYdx (a:b,c:d,bin) + domgYdx_z (a:b,c:d,i)**2
        domgYdydomgYdy (a:b,c:d,bin) = domgYdydomgYdy (a:b,c:d,bin) + domgYdy_z (a:b,c:d,i)**2
        domgYdzdomgYdz (a:b,c:d,bin) = domgYdzdomgYdz (a:b,c:d,bin) + domgYdz_z (a:b,c:d,i)**2


        domgZdxdomgZdx (a:b,c:d,bin) = domgZdxdomgZdx (a:b,c:d,bin) + domgZdx_z (a:b,c:d,i)**2
        domgZdydomgZdy (a:b,c:d,bin) = domgZdydomgZdy (a:b,c:d,bin) + domgZdy_z (a:b,c:d,i)**2
        domgZdzdomgZdz (a:b,c:d,bin) = domgZdzdomgZdz (a:b,c:d,bin) + domgZdz_z (a:b,c:d,i)**2


        domgXdxdomgYdx (a:b,c:d,bin) = domgXdxdomgYdx (a:b,c:d,bin) + domgXdx_z (a:b,c:d,i) * domgYdx_z (a:b,c:d,i)
        domgXdydomgYdy (a:b,c:d,bin) = domgXdydomgYdy (a:b,c:d,bin) + domgXdy_z (a:b,c:d,i) * domgYdy_z (a:b,c:d,i)
        domgXdzdomgYdz (a:b,c:d,bin) = domgXdzdomgYdz (a:b,c:d,bin) + domgXdz_z (a:b,c:d,i) * domgYdz_z (a:b,c:d,i)

      
      end do


      if(mod(ntime,d_freq)==0)        call dumpOUTPUT(0)
      if(mod(ntime,d_freq)==d_freq-1)  call dumpOUTPUT(1)

          
    end subroutine perform_phase_averaging_vort




    subroutine spread_to_all_configs(varIN,config,varOUT1,varOUT2)

          real*8, intent(in)  :: varIN(:,:,:)
          real*8, intent(out) :: varOUT1(:,:,:)
          real*8, intent(out) :: varOUT2(:,:,:)

          character(*), intent(in) :: config

          select case (config)

            case('X')

              call transpose_x_to_y(varIN,  varOUT1)
              call transpose_y_to_z(varOUT1,varOUT2)

            case('Y')

              call transpose_y_to_x(varIN,varOUT1)
              call transpose_y_to_z(varIN,varOUT2)

            case('Z')

              call transpose_z_to_y(varIN,  varOUT1)
              call transpose_y_to_x(varOUT1,varOUT2)

          end select

      end subroutine spread_to_all_configs


end module mod_phase_average
