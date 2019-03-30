MODULE variables

implicit none
public
!implicit double precision (a-h,o-z)


!Version sans couplage
integer :: mbufa, mstif


integer ::  nnode,      &   ! Number of NODes par Element
            npoin,      &   ! Total number of nodal points in the structure/problem
            mdouble,    &
            melem,      &   ! Maximum of ELEMents
            mevab,      &   ! Maxium of Element VAriaBles
            mfron,      &   ! Maximum front width
            mmats,      &   ! ?? Maximum number of material sets??
            mpoin,      &   ! Maximum permissible number of nodal points in the program
        	mgaus,      &   ! Utilisée dans problemes transitoire cf. Owen p.415
            mvfix,      &   ! Maximum number of prescribed boundary nodes
            ndofn,      &   ! Number of Degrees Of Freedom per Node. =2: X et Y
            mtotv,      &   ! Maximum total number of degrees of freedom =16
            mtotg,      &   ! =mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices
            nprop,      &   ! Number of material parameters required to define characteristics of a material completely
                                !=4 for elasto-plastic problem; =2 for other applications
            ndim,       &   ! nombre de dimensions du problème: ici 2 => ndim=2
            nelem,      &   ! Total number of elements in the structure. =1 pour le moment
            nmats           ! Number of different MATerialS in the structure
    

integer ::  nvfix,      &   ! Total number of boudary points, ie nodal points at which one or more dof are restrained
            ntype,      &   ! Problem type parameter; 1-Plane stress, 2-Plane strain, 3- Axial symmetry
            ngaus,      &   ! Order of gaussian quadrature rule to be employed for numerical integration of the element stiffness matrices, etc...
                                ! cf. section 6.3.2. p.173; 2: two-point Gauss, 3: three-point rule
            nalgo,      &   ! Parameter controlling nonlinear solution algorithm - Owen p. 206
                                ! 1-Initial stiffness method, 2-Tangential, 3-Combined algo. 1st ite., 4-Combined algo. 2nd ite.
            ncrit,      &   ! Yield criterion. 1-Tresca, 2-Von Mises, 3-Mohr-Coulomb, 4-Drucker-Prager
            nincs,      &   ! Total number of increments in which the final loading is to be applied
            nstre           ! Number of indpt stress components for the application: 3- Plane stress/strain, 4-Axial symmetry


integer ::  nstep,      &   ! Number of iteration - Version Elkassas
            ntotv,      &   ! Total number of degrees of freedom ? - =npoin*ndofn=mtotv
            amodel,     &   ! ??? Elkassas
            nevab,      &   ! Number of Element VAriaBles - nevab=ndofn*nnode
            nstr1,      &   ! ??
            ngaus2,     &   ! ngaus^2
            ntotg         ! Total number of Gauss points in the structure = nelem * ngaus^2
            !kgaus           ! iterator over each gauss points
integer :: nstage           ! Number of stages

! Conditions aux limites
integer ::  iplod,		&   ! Indicate point loads
            igrav,      &   ! Indicate gravity loads
            iedge,      &   ! Indicate distributed edge loading
            nedge           ! Nombre total de cotés (edge) ayant une condition aux limites 

! FIN Conditions aux limites 

double precision ::       thick, dense

integer, dimension(24):: neror

double precision :: einitial,saturi,winitial   

!double precision :: exisp, etasp
double precision :: akll
double precision :: coslx, cosly

!integer:: lprop, ielem, kgasp

! Subroutine jacob2
double precision :: djacb

double precision :: tfact=0.0

double precision , dimension (:),allocatable     :: posgp         !posgp(4)
double precision , dimension (:) ,allocatable    :: weigp         !weigp(4)
double precision , dimension (2,2)   :: xjacm         !xjacm(2,2) - Jacobian matrix
double precision , dimension (2,2)   :: xjaci         !xjaci(2,2) - Inverse of the Jacobian matrix
double precision , dimension (9)     :: shap          !shap(9)
double precision , dimension (2,9)   :: deriv         !deriv(2,9)
double precision , dimension (2,9)   :: gpcod         !gpcod(2,9) - Cartesian coordinates of the Gauss point


integer, dimension (50)    :: amode 		!Tableau probablement inutile 09/09/12
double precision   , dimension (10)    :: tinstp 		!tinstp(10) Tableau pour valeurs temporelles? read 5
double precision   , dimension (10)    :: aiso 			!aiso(10)
double precision   , dimension (10)    :: iampl 		!iampl(10)
double precision   , dimension (10)    :: antsstp 		!antsstp(10)
double precision   , dimension (10)    :: aifreq 		!aifreq(10)
double precision   , dimension (1,8)   :: shap1         !shap1(1,8)
double precision   , dimension (9,2)   :: derivt        !derivt(9,2)
double precision   , dimension (9,2)   :: cartdt        !cartdt(9,2)
double precision   , dimension (:,:),allocatable   :: elcod         !elcod(2,9)
double precision   , dimension (2,2)   :: conm          !conm(2,2)
double precision   , dimension (2)     :: deltaz        !deltaz(2)
double precision   , dimension (2,9)   :: cartd         !cartd(2,9)
double precision   , dimension (16,1)  :: cardc         !cardc(16,1)
double precision   , dimension (2,5000)   :: gpcodg     !gpcodg(2,5000)

integer, dimension (:)  , allocatable   :: alloc_stat       !Vecteurs pour checker allocation ie assez RAM
integer, dimension (:)  , allocatable   :: matno            !matno(melem), modifier matno(nelem) - Stores material property identification number (each element may be assigned different material properties)
integer, dimension (:,:), allocatable   :: lnods            !lnods(melem,9) - ELement NODe numberS listed for each element = Table de connectivite. Remplacé par lnods(nelem,nnode) p.37 - elements nodal connections and the property numbers 
integer, dimension (:)	, allocatable   :: iffix            !iffix(mtotv) mtotv=ndofn*npoin=Total number of degrees of freedom of the structure - Stock conditions aux limites 
integer, dimension (:)	, allocatable   :: nofix            !nofix(mvfix) - Restrained node number. N° global des noeuds restraints en condition limite. Stores nodes at which one or more degrees of freedom are restrained / Temporary matrix
integer, dimension (:)	, allocatable   :: ndfro            !ndfro(melem)



double precision , dimension (:,:), allocatable   :: coord            !coord(mpoin,2), remplace par coord(npoin,2) -- Coord. of the nodes
double precision , dimension (:,:), allocatable   :: presc0           !(mvfix,ndofn) - presc(ivfix,1): prescribed value of the x (or r) component of nodal displacement
                                                                        !  - presc(ivfix,2): prescribed value of the y (or z) component of nodal displacement
double precision , dimension (:,:), allocatable   :: props            !props(mmats,nprop)
double precision , dimension (:,:), allocatable   :: rload            !rload(nelem,nevab) - Stores consistent nodal loads evaluated for each element separately
double precision , dimension (:,:), allocatable   :: eload            !eload(nelem,nevab) - Contains the loading to be applied to the structure for each iteration of the solution proces. For techniques other than the direct iteration method, 
                                                                            !this vector will contain the residual nodal forces and thus differs from the vector of applied loads.
double precision , dimension (:,:), allocatable   :: tload            !tload(nelem,nevab) - Accumulates the total loading applied to the structure at any stage of the analysis.

! Subroutine zero( )
double precision , dimension (:),   allocatable   :: tdisp            !tdisp(ntotv)
double precision , dimension (:,:), allocatable   :: treac            !treac(nvfix,ndofn)
double precision , dimension (:,:), allocatable   :: treac2            !treac(npoin,ndofn)
double precision , dimension (:),   allocatable   :: epstn            !epstn(ntotg)
double precision , dimension (:),   allocatable   :: effst            !effst(ntotg)
double precision , dimension (:,:), allocatable   :: strsg            !strsg(nstr1,ntotg)
double precision , dimension (:,:), allocatable   :: defsg            !strsg(nstr1,ntotg)
double precision , dimension (:,:), allocatable   :: strsg_out        !strsg(3*nelem,(ngaus^2)+1)

! Subroutine increm( )
integer  ::  kresl,			&   ! Resolution counter which indicates if the element stiffness matrix is to be reformulated or not
									! 1- Reformulation of the element stiffnesses accompanied by a full equation solution
									! 2- Element stiffnesses are not to be modified and only equation resolution takes place
			 miter,			&	! Maximum number of iterations. Safety mesure if solution process does not converge.
			 iiter,			&
			 iincs,         &
             istage

integer , dimension(2) :: noutp ! Controls the output of the unconverged results after the first iteration. Cf. p.211
					
double precision :: facto,	&	! Controls magnitude of the load increment. If loading is prescribed by displacements the same factoring process holds
					toler		! Controls tolerance permitted on the convergence process. Cf. 3.9.3
double precision , dimension (:),   allocatable   :: fixed            !fixed(ntotv) - Vector of prescribed displacements
double precision :: abeta, young, poiss
double precision , dimension (4)  :: dvect 
!double precision , dimension (:,:), allocatable :: dbmat 
!double precision , dimension (:,:), allocatable   :: bmatx
!double precision , dimension (:,:), allocatable    :: dmatx
double precision , dimension (4,18)  :: dbmat 
double precision , dimension (:,:),allocatable   :: bmatx  ! B Matrix B(4,2*nnode)
double precision , dimension (:,:),allocatable    :: dmatx  !D Matrix D(4,4)
double precision , dimension (4)  :: devia
double precision , dimension(:),allocatable      :: stres
double precision , dimension(4)      :: avect
double precision, dimension(:,:),allocatable :: bmatx_tr   !Transpose of Matrix B 
double precision, dimension(:),allocatable :: eload_temp   !Temporary vector of load forces to store them in the big eload

double precision :: smean  

double precision :: eeta
double precision :: yield   ! Effective stress (stress level) calculated using invariants and F 

!Subroutine invar
double precision :: phira, sint3

!Subroutine Yieldf
double precision, dimension(4)  :: veca1, veca2, veca3, veca4
double precision :: steff
double precision :: theta
double precision :: I1,varj2, varj3
double precision :: frict, snphi, tanth, tant3, sinth, costh, cost3
double precision :: abthe, cons1, cons2, cons3, plumi,gseta


! Subroutine front
double precision, dimension(:,:),allocatable  :: estif


double precision   , dimension (:,:), allocatable   :: gradb            !gradb(mpoin,4)
double precision   , dimension (:,:), allocatable   :: ngrad            !ngrad(melem,4)
double precision   , dimension (:)  , allocatable   :: delz             !delz(mpoin)
double precision   , dimension (:)	, allocatable   :: asdis            !asdis(mtotv)
double precision   , dimension (:,:)	, allocatable   :: eldis            !eldis(2,9) remplacé par eldis(2,nnode)



double precision, dimension (:,:), allocatable :: gstif
! double precision, dimension (:,:), allocatable :: gstif_backup
double precision, dimension (:,:), allocatable :: gstif_temp
double precision, dimension (:,:), allocatable :: rgstif
double precision, dimension (:), allocatable :: gload, gload_temp, asdis_temp
double precision, dimension (:), allocatable :: gload_backup
! double precision, dimension (:), allocatable :: rgload_backup

double precision, dimension (:), allocatable :: rgload
double precision, dimension (:), allocatable :: vec_residu
double precision, dimension (:), allocatable :: vec_residu_pd !residual vector for prescribed displacements
double precision, dimension (:), allocatable :: vec_residu_backup
double precision, dimension (:), allocatable :: rasdis
integer :: nvlib, nnodefix



! Subroutine residu
double precision   , dimension (:),allocatable :: stran            !stran(4)
double precision   , dimension (2,9) :: dlcod          !dlcod(2,9)
double precision   , dimension (4) :: desig            !desig(4)
double precision   , dimension (:),allocatable :: sigma            !sigma(4)
double precision   , dimension (4) :: sgtot            !sgtot(4)

! Subroutine conver
double precision   , dimension (:), allocatable  :: tofor            !tofor(mtotv)
double precision   , dimension (:), allocatable  :: stfor            !tofor(mtotv)

! Subroutine ouput
integer  ::  ncheck

integer :: drapeau


double precision, parameter :: root3=1.7320508075688772d0            !dsqrt(3)
double precision, parameter :: twopi=2*3.14159265358979323846d0      !2*pi
double precision, parameter :: radian=57.29577951308232d0           !180/pi radian en [degres]
double precision, parameter :: degres=0.0174532925199432957692369d0   !pi/180 degres en [radians]
double precision, parameter :: root13=0.57735026918962576450914878050196d0  !racine(1/3)

!Version V2 
!     double precision, parameter :: epsilon=1.0d-6
    double precision :: epsilon
    double precision :: pvalue, norme

! Variables for time calculation of the simulation
double precision :: start, finish, start_inc, finish_inc

! Subroutine front( )
!double precision, dimension (mstif,mstif) :: gstif
!double precision, dimension (ntotv,ntotv) :: gstif_temp
!double precision, dimension (nvfix,nvfix) :: rgstif
!double precision, dimension (ntotv) :: gload, gload_temp, asdis_temp

!double precision, dimension (nvfix) :: rgload
!double precision, dimension (nvfix) :: rasdis
integer, dimension (:), allocatable :: gcoord ! coordonnées globales des elements de rgload dans gload par rapport aux déplacements bloqués

integer, parameter :: SIZE2=24

integer :: drap_pres

!For lapack solve
integer,dimension(:),allocatable :: ipiv
double precision, dimension (:,:), allocatable :: rgstif_lapack
double precision, dimension (:), allocatable :: vec_residu_lapack
double precision, dimension (:), allocatable :: rasdis_lapack
double precision :: alpha_dgemm, beta_dgemm


! integer :: num_example

integer :: compt_alloc
integer :: plast_code
character(len=30):: name_disp

double precision, dimension(:),allocatable :: rstfor
integer, dimension (:), allocatable :: riffix
double precision, dimension (:), allocatable :: rfixed


! Matrices for output_python
! double precision, dimension (:,:), allocatable :: gpcod_all_global !gpcod_all_global(nelem,2*ngaus)
integer :: count_file

! Cam Clay
double precision :: CClambda, CCkappa, Nlambda, Nkappa, csl
double precision, dimension (:), allocatable :: volspegp
double precision, dimension (:), allocatable :: epstnp
double precision, dimension (:), allocatable :: yieldgp
double precision, dimension (:), allocatable :: young_gp
double precision, dimension (:), allocatable :: deviatoric_stress
double precision, dimension (:), allocatable :: effective_mean_stress
double precision, dimension (:), allocatable :: epstnp_pre
! Lime treated
double precision :: Dei, Dec, pyI, pyII, pb, beta
double precision, dimension (:), allocatable :: Deigp


END MODULE variables
 
