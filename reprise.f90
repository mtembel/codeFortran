SUBROUTINE read_reprise(u,v,temp,lu,lv,ltemp,nstep_fin)

	USE Parametres
	USE Constantes

	implicit none
    
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(inout)  :: u,v ,temp      
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(inout)  :: lu,lv ,ltemp 
	integer :: nstep_fin     
    

    open(unit=30,file='fichier_reprise',form='unformatted')
    
    read(30) nstep_fin
    read(30) t
    read(30) u
	read(30) v
	read(30) temp
	read(30) lu
	read(30) lv
	read(30) ltemp

	close(30)

END SUBROUTINE read_reprise

SUBROUTINE write_reprise(u,v,temp,lu,lv,ltemp,nstep_fin)

	USE Parametres
	USE Constantes

	implicit none
    
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(inout)  :: u,v ,temp      
    real*8 ,dimension(-2:nx+3,-2:ny+3),intent(inout)  :: lu,lv ,ltemp      
	integer :: nstep_fin     
    

    open(unit=30,file='fichier_reprise',status='replace',form='unformatted')

    write(30) nstep_fin
    write(30) t
    write(30) u
	write(30) v
	write(30) temp
	write(30) lu
	write(30) lv
	write(30) ltemp

	close(30)

END SUBROUTINE write_reprise
