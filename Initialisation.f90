MODULE Parametres

	implicit none

!***************************************************************!
!	Variables Generales pour le Solveur Navier-Stokes			!
!***************************************************************!

	integer, parameter :: nx  =  30   !nombre de cellules de calcul en x
	integer, parameter :: ny  =  15    !nombre de cellules de calcul en y
	integer, parameter :: nz  =  1      !nombre de cellules de calcul en z

	integer, parameter :: nstep  =  30000 !nombre de pas de temps de la simulation
	integer, parameter :: isto   =  50 !frequence de stockage de la solution ( tous les isto pas de temps)

	integer, parameter :: reprise = 0     !indicateur de reprise de calcul : 0 on part de t=0, 1 on lit un fichier de reprise

	real*8 , parameter :: tmax = 10.       ! temps maxi de simulation

	real*8 , parameter :: cfl     =  0.5 ! condition de stabilité pour l'intégration en temps des eqts de Navier-Stokes


!***************************************************************!
!	Domaine + Conditions aux limites + Parametres Physiques		!
!***************************************************************!

	real*8 , parameter :: xmin  =   0.0 !en m
	real*8 , parameter :: xmax  =   0.04   !en m

	real*8 , parameter :: ymin  =  0.   !en m
	real*8 , parameter :: ymax  =  0.01   !en m

	real*8 , parameter :: nu       =  0.000018  !coefficient de viscosite cinematique du fluide en m^2/s
	real*8 , parameter :: rho_ref  =  1.  !densite de reference du fluide  en kg/m^3
	real*8 , parameter :: temp_ref =  300.   !température de reference du fluide  en K
	real*8 , parameter :: temp_haut=  295.   !température de reference du fluide  en K
	real*8 , parameter :: temp_bas =  2850   !température de reference du fluide  en K
  	real*8 , parameter :: kappa    =  0.000026    !coefficient de diffusivite thermique du fluide
  	real*8 , parameter :: alpha    =  0.01  !coefficient de dilatation thermique en kg/m^3/K
	
	real*8 , parameter :: grav_x  =  0.0  !
	real*8 , parameter :: grav_y  =  -9.81 ! en m/s^2

!***************************************************************!
!	Pas de discretisation spatiale dans les 2 directions		!
!***************************************************************!

	real*8 , parameter :: dx = (xmax-xmin)/nx
	real*8 , parameter :: dy = (ymax-ymin)/ny
!***************************************************************!
!   Parametres et tableaux de travail pour le solveur de pression
!***************************************************************!
    integer ,parameter ::ndim = nx*ny,mdim=3
    real*8  ,dimension(1:ndim,1:mdim) :: coef 
	real*8  ,dimension(1:ndim) :: rhs1,p_s,r_s,r2_s
	real*8  ,dimension(1:ndim) :: q_s,s_s,x1
	real*8  ,dimension(1:ndim,1:5)::l_s
	real*8  :: zeta=1.e-10 !precision du solveur de pression
	integer :: itmax=300   !nombre d'iterations max du solveur
    integer, dimension(1:mdim):: jcoef
	integer, dimension(1:5):: ldiag

!***************************************************************!
!	Variables Autres - a ne pas modifier						!
!***************************************************************!

	real*8 :: dt = 0.0, t = 0.0  !pas de temps et temps

	integer :: istep = 1         !indice de pas de temps pour la boucle en temps

	integer :: nxu = nx-1        !taille des tableaux pour chaque composante de  la vitesse (u,v)
	integer :: nyu = ny          !sur maillages décalés 
	integer :: nxv = nx
	integer :: nyv = ny-1
    
	real*8, parameter  :: ddx = 1.d0/dx
	real*8, parameter  :: ddy = 1.d0/dy


END MODULE Parametres


MODULE Constantes

	implicit none

	real*8 :: d1p3		=		1.d0/3.d0
	real*8 :: d2p3		=		2.d0/3.d0
	real*8 :: d1p4		=		1.d0/4.d0
	real*8 :: d3p4		=		3.d0/4.d0
	real*8 :: d1p6		=		1.d0/6.d0
	real*8 :: d5p6		=		5.d0/6.d0
	real*8 :: d1p6m		=	-	1.d0/6.d0
	real*8 :: d7p6		=		7.d0/6.d0
	real*8 :: d7p6m		=	-	7.d0/6.d0
	real*8 :: d11p6		=		11.d0/6.d0
	real*8 :: d13p12	=		13.d0/12.d0
	real*8 :: d1p12		=		1.d0/12.d0
	real*8 :: d1p8		=		1.d0/8.d0
	real*8 :: d1p18		=		1.d0/18.d0

	real*8, parameter :: eps_weno	=	1.e-16

	real*8 :: zero = 0.0d0
	real*8 :: one  = 1.0d0
	real*8 :: two  = 2.0d0
	real*8 :: demi = 0.5d0

	real*8 :: c13  = 13.0d0
	real*8 :: c3   = 3.0d0
	real*8 :: c4   = 4.0d0

	real*8, parameter :: pi = 3.1415926535897932384626
	real*8, parameter :: dpi = 1.d0/3.1415926535897932384626


END MODULE Constantes
