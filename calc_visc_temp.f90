SUBROUTINE calc_visc_temp(temp,visct)

	USE Parametres
	USE Constantes
	implicit none

real*8 ,dimension(-2:nx+3,-2:ny+3),intent(in)  :: temp  
real*8 ,dimension(-2:nx+3,-2:ny+3),intent(out)  :: visct  

real*8 :: d_e,d_o,d_n,d_s ! Flux diffusifs est ouest nord sud
integer:: i,j

Do j=-2,ny+3
	do i=-2,nx+3
		visct(i,j)=0
	End do
End do

Do j=1,ny
	do i=1,nx
		d_o=-1*kappa*(temp(i,j)-temp(i-1,j))/(dx*dx)
		d_e=kappa*(temp(i+1,j)-temp(i,j))/(dx*dx)
		d_n=kappa*(temp(i,j+1)-temp(i,j))/(dy*dy)
		d_s=-1*kappa*(temp(i,j)-temp(i,j-1))/(dy*dy)
		visct(i,j)=(d_e+d_o+d_n+d_s)
	End do
End do


END SUBROUTINE 
