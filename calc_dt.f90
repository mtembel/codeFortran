real*8 FUNCTION calc_dt(umax,vmax)

	USE Parametres
	USE Constantes

	implicit none

	real*8, intent(in) :: umax,vmax

	real*8 :: conv_cfl,visc_cfl,grav_cfl,diff_cfl
	
	conv_cfl = dmax1(umax*ddx,vmax*ddy) !cfl = 1
	visc_cfl = nu*(2.*ddx*ddx+2.*ddy*ddy) ! r=0.5
    diff_cfl = kappa*(2.*ddx*ddx+2.*ddy*ddy) ! r=0.5
	grav_cfl = dsqrt(dmax1(dabs(grav_x),dabs(grav_y))*ddy)
	

	calc_dt = 2.0d0/(2.*conv_cfl+visc_cfl+diff_cfl+grav_cfl)


	return

END FUNCTION calc_dt

real*8 FUNCTION max_comp(vit,nx1,ny1)

	USE Parametres
	USE Constantes

	implicit none

	real*8, dimension(-2:nx+3,-2:ny+3), intent(in) :: vit

	integer :: i,j,nx1,ny1

	max_comp = 0.0
	
	do j=1,ny1
		do i=1,nx1
			max_comp = dmax1(max_comp,dabs(vit(i,j)))
		end do
	end do

	return

END FUNCTION max_comp
