SUBROUTINE bc_temp(temp)
	USE Parametres
	USE Constantes

	implicit none
    
    real*8 , dimension(-2:nx+3,-2:ny+3),intent(inout)  :: temp       
     
    integer:: i,j

!
!   CONDITION LIMITE SUR LA PAROI GAUCHE : flux de temperature nul
!

    do j=-2,ny+3
        
       temp( 0,j)=temp(1,j)
       temp(-1,j)=temp(2,j)
       temp(-2,j)=temp(3,j)
    enddo


!
!   CONDITION LIMITE SUR LA PAROI DROITE :  flux de temperature nul
!

!
!   CONDITION SUR U
!

    do j=-2,ny+3
        
       temp(nx+1,j)=temp(nx,j)
       temp(nx+2,j)=temp(nx-1,j)
       temp(nx+3,j)=temp(nx-2,j)
    enddo

 
!
!   CONDITION LIMITE SUR LA PAROI INFERIEURE : temperature imposee temp_bas
!

!
!   CONDITION SUR U
!

    do i=-2,nx+3
        
       temp(i,0 )= 2.*temp_bas-temp(i,1)
       temp(i,-1)= 2.*temp_bas-temp(i,2)
       temp(i,-2)= 2.*temp_bas-temp(i,3)
    enddo

 
!
!   CONDITION LIMITE SUR LA PAROI SUPERIEURE : temperature imposee temp_haut
!

    do i=-2,nx+3
        
       temp(i,ny+1)= 2.*temp_haut-temp(i,ny  )
       temp(i,ny+2)= 2.*temp_haut-temp(i,ny-1)
       temp(i,ny+3)= 2.*temp_haut-temp(i,ny-2)
    enddo


END SUBROUTINE bc_temp
