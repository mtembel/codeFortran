SUBROUTINE calc_visc(u,v,viscu,viscv)

	USE Parametres
	USE Constantes

	implicit none
    
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(in)  :: u,v       
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(out) :: viscu,viscv       
     
    integer:: i,j
!
!   CALCUL DES TERMES VISQUEUX POUR L'EQUATION DE U
!

    do j=1,nyu
      do i=1,nxu

         viscu(i,j)=nu*( &
         (u(i+1,j)-2.d0*u(i,j)+u(i-1,j))*ddx*ddx+ &
         (u(i,j+1)-2.d0*u(i,j)+u(i,j-1))*ddy*ddy)
      enddo
    enddo



!   CALCUL DES TERMES VISQUEUX POUR L'EQUATION DE V
!

    do j=1,nyv
      do i=1,nxv

         viscv(i,j)=nu*( &
         (v(i+1,j)-2.d0*v(i,j)+v(i-1,j))*ddx*ddx+&
         (v(i,j+1)-2.d0*v(i,j)+v(i,j-1))*ddy*ddy)
      enddo
    enddo

END SUBROUTINE calc_visc