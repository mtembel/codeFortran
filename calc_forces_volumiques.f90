SUBROUTINE calc_forces_volumiques(force_u,force_v,temp)

	USE Parametres
	USE Constantes
	implicit none

real*8 ,dimension(-2:nx+3,-2:ny+3),intent(in)  :: temp  
real*8 ,dimension(-2:nx+3,-2:ny+3),intent(out) :: force_u,force_v

Integer :: i,j 

do j=1,ny
	do i=1,nx
		force_u(i,j)=rho_ref*(1.d0-alpha*(0.5*(temp(i,j)+temp(i+1,j))-temp_ref))*grav_x
		force_v(i,j)=rho_ref*(1.d0-alpha*(0.5*(temp(i,j)+temp(i,j+1))-temp_ref))*grav_y
	End do
End do

END SUBROUTINE 
