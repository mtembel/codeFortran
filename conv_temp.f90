SUBROUTINE conv_temp (temp,u,v,ltemp)

	USE Parametres
	USE Constantes
	implicit none

real*8 ,dimension(-2:nx+3,-2:ny+3),intent(in)  :: temp  
real*8 ,dimension(-2:nx+3,-2:ny+3),intent(in) :: u,v
real*8 ,dimension(-2:nx+3,-2:ny+3),intent(out)  :: ltemp 

real*8 :: a_e,a_o,a_n,a_s ! Flux advectifs est ouest nord sud
integer:: i,j

Do j=-2,ny+3
	do i=-2,nx+3
		ltemp(i,j)=0
	End do
End do

! schema centre
do i=1,nx
     do j=1,ny
     
       a_e = -1*u(i,j)*0.5*(temp(i+1,j)+temp(i,j))/dx
       a_o = u(i-1,j)*0.5*(temp(i-1,j)+temp(i,j))/dx
       a_n = -1*v(i,j)*0.5*(temp(i,j)+temp(i,j+1))/dy
       a_s = v(i,j-1)*0.5*(temp(i,j)+temp(i,j-1))/dy
       ltemp(i,j)=a_e+a_o+a_n+a_s
       
     end do
end do



!schema upwind checker les moins

!do i=1,nx
!     do j=1,ny
     
!       a_e = -1*u(i,j)*(0.5*(1+sign(1.,u(i,j)))*temp(i,j)+ 0.5*(1-sign(1.,u(i,j)))*temp(i+1,j))/ dx
!       a_o = u(i-1,j)*(0.5*(1+sign(1.,u(i-1,j)))*temp(i-1,j)+ 0.5*(1-sign(1.,u(i-1,j)))*temp(i,j))/dx
!       a_n = -1*v(i,j)*(0.5*(1+sign(1.,v(i,j)))*temp(i,j)+ 0.5*(1-sign(1.,v(i,j)))*temp(i,j+1))/dy
!       a_s = v(i,j-1)*(0.5*(1+sign(1.,v(i,j-1)))*temp(i,j-1)+ 0.5*(1-sign(1.,v(i,j-1)))*temp(i,j))/dy
!       ltemp(i,j)=a_e+a_o+a_n+a_s
       
!     end do
!end do



END SUBROUTINE 



