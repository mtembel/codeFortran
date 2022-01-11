SUBROUTINE calc_vort_uv_centre(u,v,vort,u_cent,v_cent)

	USE Parametres
	USE Constantes

	implicit none
    
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(in)  :: u,v       
    real*8 ,dimension(1:nx,1:ny),intent(out) :: u_cent,v_cent,vort    

! 
!   tableau de travail pour le calcul de la vorticite
!
 
     real*8 ,dimension(-2:nx+3,-2:ny+3)  :: ome     
    
    integer:: i,j

    do j=1,nyu
      do i=1,nxu

         u_cent(i,j)=0.5*(u(i,j)+u(i-1,j))

      enddo
    enddo

     do j=1,nyv
      do i=1,nxv

         v_cent(i,j)=0.5*(v(i,j)+v(i,j-1))

      enddo
    enddo
!
!
!  
    do j=2,ny-1
       do i=2,nx-1

          ome(i,j)=(v(i,j)-v(i-1,j))*ddx-(u(i,j)-u(i,j-1))*ddy

      enddo
    enddo

    do i=2,nx-1
  
       ome(i,1   )=-2.*u(i,1 )*ddy
       ome(i,ny+1)= 2.*u(i,ny)*ddy
       
    enddo


    do j=2,ny-1
  
       ome(1   ,j)= 2.*v(1 ,j)*ddx
       ome(nx+1,j)=-2.*v(nx,j)*ddx
       
    enddo
     
    ome(1,1)=0.;ome(1,ny+1)=0.;ome(nx+1,1)=0.;ome(nx+1,ny+1)=0.;

    do j=1,ny
      do i=1,nx

 !        vort(i,j)=0.25*(ome(i,j)+ome(i+1,j)+ome(i,j+1)+ome(i+1,j+1))
        vort(i,j)=0.25*((v(i+1,j  )-v(i,j    ))/dx-(u(i  ,j+1)-u(i  ,j  ))/dy + &
                        (v(i  ,j  )-v(i-1,j  ))/dx-(u(i-1,j+1)-u(i-1,j  ))/dy + &
                        (v(i+1,j-1)-v(i,j-1  ))/dx-(u(i  ,j  )-u(i  ,j-1))/dy + &
                        (v(i  ,j-1)-v(i-1,j-1))/dx-(u(i-1,j  )-u(i-1,j-1))/dy)
      enddo
    enddo

    return

END SUBROUTINE calc_vort_uv_centre
