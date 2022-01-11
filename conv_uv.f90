SUBROUTINE conv_u(u,v,l)

	USE Parametres
	USE Constantes

	implicit none

	real*8, dimension(-2:nx+3,-2:ny+3), intent(in ) :: u,v
	real*8, dimension(-2:nx+3,-2:ny+3), intent(out) :: l


	integer :: i,j
	real*8 :: weno5_z,uph,umh,vph,vmh
	real*8 :: fimh,fiph,gjmh,gjph

    

   do j=1,nyu
     do i=1,nxu
       
				uph = 0.5*(u(i+1,j)+u(i,j))
				umh = 0.5*(u(i-1,j)+u(i,j))

				if ( uph > zero ) then
					fiph=uph*weno5_z(u(i-2,j),u(i-1,j),u(i,j),u(i+1,j),u(i+2,j))
				else
					fiph=uph*weno5_z(u(i+3,j),u(i+2,j),u(i+1,j),u(i,j),u(i-1,j))
				end if

				if ( umh > zero ) then
					fimh=umh*weno5_z(u(i-3,j),u(i-2,j),u(i-1,j),u(i,j),u(i+1,j))
				else
					fimh=umh*weno5_z(u(i+2,j),u(i+1,j),u(i,j),u(i-1,j),u(i-2,j))
				end if

				vph = 0.5*(v(i,j  )+v(i+1,j  ))
				vmh = 0.5*(v(i,j-1)+v(i+1,j-1))

				if ( vph > zero ) then
					gjph=vph*weno5_z(u(i,j-2),u(i,j-1),u(i,j),u(i,j+1),u(i,j+2))
				else
					gjph=vph*weno5_z(u(i,j+3),u(i,j+2),u(i,j+1),u(i,j),u(i,j-1))
				end if

				if ( vmh > zero ) then
				 	gjmh=vmh*weno5_z(u(i,j-3),u(i,j-2),u(i,j-1),u(i,j),u(i,j+1))
				else
					gjmh=vmh*weno5_z(u(i,j+2),u(i,j+1),u(i,j),u(i,j-1),u(i,j-2))
				end if

				l(i,j) = - ( (fiph - fimh)*ddx + (gjph - gjmh)*ddy )

     enddo
   enddo
          
END SUBROUTINE conv_u 

SUBROUTINE conv_v(u,v,l)

	USE Parametres
	USE Constantes

	implicit none

	real*8, dimension(-2:nx+3,-2:ny+3), intent(in ) :: u,v
	real*8, dimension(-2:nx+3,-2:ny+3), intent(out) :: l


	integer :: i,j
	real*8 :: weno5_z,uph,umh,vph,vmh
	real*8 :: fimh,fiph,gjmh,gjph

    

   do j=1,nyv
     do i=1,nxv
       
				uph = 0.5*(u(i,j  )+u(i  ,j+1))
				umh = 0.5*(u(i-1,j)+u(i-1,j+1))

				if ( uph > zero ) then
					fiph=uph*weno5_z(v(i-2,j),v(i-1,j),v(i,j),v(i+1,j),v(i+2,j))
				else
					fiph=uph*weno5_z(v(i+3,j),v(i+2,j),v(i+1,j),v(i,j),v(i-1,j))
				end if

				if ( umh > zero ) then
					fimh=umh*weno5_z(v(i-3,j),v(i-2,j),v(i-1,j),v(i,j),v(i+1,j))
				else
					fimh=umh*weno5_z(v(i+2,j),v(i+1,j),v(i,j),v(i-1,j),v(i-2,j))
				end if

				vph = 0.5*(v(i,j  )+v(i,j+1))
				vmh = 0.5*(v(i,j-1)+v(i,j  ))

				if ( vph > zero ) then
					gjph=vph*weno5_z(v(i,j-2),v(i,j-1),v(i,j),v(i,j+1),v(i,j+2))
				else
					gjph=vph*weno5_z(v(i,j+3),v(i,j+2),v(i,j+1),v(i,j),v(i,j-1))
				end if

				if ( vmh > zero ) then
				 	gjmh=vmh*weno5_z(v(i,j-3),v(i,j-2),v(i,j-1),v(i,j),v(i,j+1))
				else
					gjmh=vmh*weno5_z(v(i,j+2),v(i,j+1),v(i,j),v(i,j-1),v(i,j-2))
				end if

				l(i,j) = - ( (fiph - fimh)*ddx + (gjph - gjmh)*ddy )

     enddo
   enddo
          
END SUBROUTINE conv_v    

real*8 FUNCTION weno5(a,b,c,d,e)

	USE Constantes

	implicit none

	real*8, intent(in) :: a,b,c,d,e

	real*8 :: q1,q2,q3
	real*8 :: is1,is2,is3
	real*8 :: alpha1,alpha2,alpha3

	q1 = a * d1p3  + b * d7p6m + c * d11p6
	q2 = b * d1p6m + c * d5p6  + d * d1p3
	q3 = c * d1p3  + d * d5p6  + e * d1p6m

	is1 = c13*( a - two*b + c )**2 + c3*( a - c4*b + c3*c )**2  + eps_weno
	is2 = c13*( b - two*c + d )**2 + c3*( d - b )**2            + eps_weno
	is3 = c13*( c - two*d + e )**2 + c3*( c3*c - c4*d + e )**2  + eps_weno

	is1 = is1**2
	is2 = is2**2
	is3 = is3**2

	alpha1 = is2*is3*d1p18
	alpha2 = is1*is3*d1p3
	alpha3 = is1*is2*d1p6

	weno5 = alpha1 * q1 + alpha2 * q2 + alpha3 * q3
	weno5 = weno5 / ( alpha1 + alpha2 + alpha3 )
	
	return

END FUNCTION weno5


real*8 FUNCTION weno5_z(a,b,c,d,e)

	USE Constantes

	implicit none

	real*8, intent(in) :: a,b,c,d,e

	real*8 :: q1,q2,q3
	real*8 :: is1,is2,is3
	real*8 :: alpha1,alpha2,alpha3
	real*8 :: tau5

	q1 = a * d1p3  + b * d7p6m + c * d11p6
	q2 = b * d1p6m + c * d5p6  + d * d1p3
	q3 = c * d1p3  + d * d5p6  + e * d1p6m

	is1 = c13*( a - two*b + c )**2 + c3*( a - c4*b + c3*c )**2  + eps_weno
	is2 = c13*( b - two*c + d )**2 + c3*( d - b )**2            + eps_weno
	is3 = c13*( c - two*d + e )**2 + c3*( c3*c - c4*d + e )**2  + eps_weno

	tau5 = dabs( is1 - is3 )

	is1  = is1 **2
	is2  = is2 **2
	is3  = is3 **2
	tau5 = tau5**2

	alpha1 = is2*is3*(is1+tau5)*d1p18
	alpha2 = is1*is3*(is2+tau5)*d1p3
	alpha3 = is1*is2*(is3+tau5)*d1p6

	weno5_z = alpha1 * q1 + alpha2 * q2 + alpha3 * q3
	weno5_z = weno5_z / ( alpha1 + alpha2 + alpha3 )
	
	return

END FUNCTION weno5_z
  