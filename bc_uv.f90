SUBROUTINE bc_uv(u,v)
	USE Parametres
	USE Constantes

	implicit none
    
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(inout)  :: u,v       
     
    integer:: i, j

!
!   CONDITION LIMITE SUR LA PAROI GAUCHE : ADHERENCE A LA PAROI u=v=0
!

!
!   CONDITION SUR U
!

    do j=-2,ny+3
        
       u( 0,j)=0.
       u(-1,j)=0.
       u(-2,j)=0.
    enddo

!
!   CONDITION SUR V
!
     do j=-2,ny+3
        
       v( 0,j)=-v(1,j)
       v(-1,j)=-v(2,j)
       v(-2,j)=-v(3,j)
    enddo

!
!   CONDITION LIMITE SUR LA PAROI DROITE : ADHERENCE A LA PAROI u=v=0
!

!
!   CONDITION SUR U
!

    do j=-2,ny+3
        
       u(nxu+1,j)=0.
       u(nxu+2,j)=0.
       u(nxu+3,j)=0.
    enddo

!
!   CONDITION SUR V
!
     do j=-2,ny+3
        
		v(nxv+1,j)= - v(nxv  ,j)
		v(nxv+2,j)= - v(nxv-1,j)
		v(nxv+3,j)= - v(nxv-2,j)
    enddo
 
!
!   CONDITION LIMITE SUR LA PAROI INFERIEURE : ADHERENCE A LA PAROI u=v=0
!

!
!   CONDITION SUR U
!

    do i=-2,nx+3
        
       u(i,0 )= - u(i,1)
       u(i,-1)= - u(i,2)
       u(i,-2)= - u(i,3)
    enddo

!
!   CONDITION SUR V
!
     do i=-2,nx+3
        
		v(i, 0)= 0.
		v(i,-1)= 0.
		v(i,-2)= 0.
    enddo
 
!
!   CONDITION LIMITE SUR LA PAROI SUPERIEURE : ADHERENCE A LA PAROI u=v=0
!

!
!   CONDITION SUR U
!

    do i=-2,nx+3
        
       u(i,nyu+1)=-u(i,nyu  )
       u(i,nyu+2)=-u(i,nyu-1)
       u(i,nyu+3)=-u(i,nyu-2)
!       u(i,nyu+1)= 2.- u(i,nyu  )
!       u(i,nyu+2)= u(i,nyu+1)
!       u(i,nyu+3)= u(i,nyu+1)
    enddo

!
!   CONDITION SUR V
!
     do i=-2,nx+3
        
		v(i,nyv+1)= 0.
		v(i,nyv+2)= 0.
		v(i,nyv+3)= 0.
    enddo

END SUBROUTINE bc_uv
