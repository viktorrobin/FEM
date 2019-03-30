
!***********************************************************************************************
!
!   This file contains subroutines for geometry of the problem and boundaries conditions
!   
!   By Victor ROBIN - 09.06.2013, Exeter
!**********************************************************************************************


!**********************************************************************************************


!***********************************************************************************************
!   This subroutine presets variables associated with dynamic dimensioning
!   Owen p.238
!**********************************************************************************************
 

subroutine dimen( )

    use variables

    !implicit double precision (a-h,o-z)
    implicit none
    integer :: num_example
    ! num_example=0: Example OWEN p.542
    ! num_example=1: 2 elements Q4
    ! num_example=2: 1 element Q4
    ! num_example=3: 2 elements T3
    ! num_example=4: 1 element Q8
    ! num_example=5: 4 elements Q8
    ! num_example=6: 4 elements Q8 + Distributed loading
    ! num_example=7: 4 elements Q8 + Deplacements imposes
    ! num_example=8: 1 element Q4 + Deplacements imposes - DOESN'T WORK
    ! num_example=9: 4 elements Q8 + Distributed loading - Axisymetry
    ! num_example=10: 4 elements Q8 + Point Load - AXIS
    ! num_example=11: 1 elements Q4 + Point Load - AXIS
    ! num_example=12: 2 elements T3 + Deplacements inposes
    ! num_example=13: 1 elements Q4 + Point Load - PLASTICITE PARFAITE
    ! num_example=14: 1 elements Q4 - PLASTICITE PARFAITE + Deplacements Imposes
    ! num_example=15: 2T3 + Deplacements Imposes - Calibration pour deplacements imposes
    ! num_example=16: 4Q4 - Plasticité parfaite
    ! num_example=17: 4Q4 - Plasticité parfaite - Deplacements imposés

    ! num_example=18: 1Q4 - Cam Clay
    ! num_example=19: 1Q4 - Cam Clay Lime Treated Soils

    num_example=18

    call initial_cond(num_example)


    return

end



!*******************************************************************************

!   Subroutine input: this subroutine accepts most of the input data Owen p.208

!*******************************************************************************

subroutine input( )

    use variables
    implicit none

    integer :: ipoin, idime, nb_tab_dyn, ivfix, ifpre, iprop, ielem, igaus
    integer :: i, j, idofn, ifdof, imats, istep, ngash, nloca, numat, numel, inode
    integer :: icompt, itotv

    real :: cond(mpoin,5)
    character (len=60) :: title

    include 'subroutines/interface.f90'

    ! If you want to print console output in a txt file
    !   open(unit=6,file='output_terminal.txt')         !Writing of the variables values


    ! Create file to store element stiffness matrices
    open(unit=1,form='unformatted',file='temp/estif.tmp')

    !1 Lecture et ecriture du titre du fichier d'entree
    read(5,920) title 
    write(6,920) title
    920 format(a)      

    ! Read the first data card, aand echo it immediately

    !2 Caracteristiques du maillage
    read(5,*) npoin,nelem,nvfix,ntype,nnode,nmats,ngaus,nalgo,ncrit,nincs,nstre  !2
    900 format(11i5)
    
    nevab=ndofn*nnode
    nstr1=nstre+1

    if(ntype == 3) then
        nstr1=nstre
    end if

    ntotv=npoin*ndofn
    ngaus2=ngaus*ngaus
    ntotg=nelem*ngaus2


    !  ***      DYNAMIC ALLOCATION OF MEMORY  ****
    include 'workspace/Alloc.f90'

    if (ncrit==5 .or. ncrit==6) then
        include 'workspace/Alloc_CamClay.f90'
    end if

    !! **** PARTIE COUPLAGE
    !3 Conditions initiales de saturation, indice des vides, etc...
    read(5,*) nstep,einitial,saturi,winitial,(amode(istep), istep=1,nstep) !3 

    !4 Valeur amodel??
    read(5,*) amodel !4

    !5 Valeurs liées à la discrétisation temporelle??
    read(5,*) (iampl(i),aiso(i),tinstp(i),antsstp(i),aifreq(i),i=1,nstep) !5
    !! **** FIN PARTIE COUPLAGE




        write(6,901) npoin,nelem,nvfix,ntype,nnode,nmats,ngaus,nevab,nalgo,ncrit,nincs,nstre !2
    901 format(//8h npoin =,i4,4x,8h nelem =,i4,4x,8h nvfix =,i4,4x,        &
                      8h ntype =,i4,4x,8h nnode =,i4,                       &   
                      //8h nmats =,i4,4x,8h ngaus= ,i4,                     &
                      4x,8h nevab =,i4,4x,8h nalgo =,i4,4x,                 &
                      //8h ncrit =,i4,4x,8h nincs =,i4,4x,8h nstre =,i4)


    !call check1(ndofn,nelem,ngaus,nmats,nnode,npoin,nstre,ntype,nvfix,ncrit,nalgo,nincs)   
    call check1( )

    ! On indique le type de problem:
    write(6,*) ''
    write(6,*) '--------------------------------------------'
    if (ntype==1) then
        write(6,*), "CALCUL CHARACTERISTICS: PLANE STRESS"
    end if
    
    if (ntype==2) then
        write(6,*) "CALCUL CHARACTERISTICS: PLANE STRAIN"
    end if

    if (ntype==3) then
        write(6,*) "CALCUL CHARACTERISTICS: AXISYMMETRIC PROBLEM"
    end if

    select case(ncrit)
        case(1)     !ncrit=1: Tresca yield criterion
            write(6,*) "TRESCA"
            
        case(2)     !ncrit=2: Von Mises
            write(6,*) "VON MISES"
    !       
        case(3)     !ncrit=3: Mohr-Coulomb
            write(6,*) "MOHR-COULOMB"

        case(4)     !ncrit=4: Drucker-Prager
            write(6,*) "DRUCKER-PRAGER"

        case(5)     !ncrit=5: Cam-Clay
            write(6,*) "CAM CLAY"
    end select


    write(6,*) '--------------------------------------------'

    !Read the elements nodal connections and the property numbers
        write(6,902) !3
    902 format(//8h element,3x,8hproperty,6x,12hnode numbers)
          
    do ielem=1,nelem
        read(5,*) numel,matno(numel),(lnods(numel,inode),inode=1,nnode)   !6
        write(6,903) numel,matno(numel),(lnods(numel,inode),inode=1,nnode)!4  
    end do

    903 format(1x,i5,i9,6x,8i5)
    
!     do i=1,numel
!         write(6,*) (lnods(i,j),j=1,nnode)
!     end do


    ! Set to zero all the nodal coordinates, prior to reading some of them
    do ipoin=1,npoin
        do idime=1,ndim
            coord(ipoin,idime)=0.0
        end do
    end do



    ! Read some nodal coordinates, finishing with the last node of all
    ! Numeros des noeuds + coord. dans ficher de sortie
        write(6,904) !5
    904 format(//5h node,10x,1hx,10x,1hy)

    !7 Lecture des coordonnees
    do while (ipoin/=npoin)
        read(5,*) ipoin,(coord(ipoin,idime), idime=1,ndim) !7
    end do

    905 format(i5,6f10.5)

    ! interpolate coordinates of mid-side nodes
    !call nodexy(coord,lnods,melem,mpoin,nelem,nnode)

    if (nnode>3) then
        call nodexy( )
    end if

    
    do ipoin=1,npoin
        write(6,906) ipoin,(coord(ipoin,idime), idime=1,ndim) !6
    end do

    906 format(1x,i5,3f10.4)



    !Read the fixed values
    write(6,907) !7
    ! original: write(6,*) !7
    907 format(//5h node,6x,4hcode,6x,12h fixed values)


    ! Version Owen: presc au lieu de presc0 -> pour couplage
    ! ifpre: Restraint code: 
        ! 10 Nodal displacement restrained in the x (or r) direction
        ! 01 Nodal displacement restrained in the y (or z) direction
        ! 11 Nodal displacement restrained in both coordinate direction
        ! 0  No boundary condition
        ! ifpre controls which degrees of freedom of a particualr node are to have a specified displacement value

    !*** Version Owen
    !do ivfix=1,nvfix
    !    read(5,908) nofix(ivfix),ifpre,(presc0(nofix(ivfix),idofn),idofn=1,ndofn) !8
    !    write(6,908) nofix(ivfix),ifpre,(presc0(nofix(ivfix),idofn),idofn=1,ndofn) !8
    !    nloca=(nofix(ivfix)-1)*ndofn ! Donne le numéro de ligne dans iffix de la premiere composante (selon x) du noeud nofix(ivfix)
    !    ifdof=10**(ndofn-1) !=10


    !    do idofn=1,ndofn
    !        ngash=nloca+idofn
    !        if(ifpre >= ifdof) then !si oui, ie ifpre>= 10 alors c'est bloqué selon x, si >=1, bloqué selon y
    !            iffix(ngash)=1  ! si =1, ca veut dire que ce noeud à une condition aux limites pour déplacements dans direction idofn
    !            ifpre=ifpre-ifdof
    !        end if
    !       ifdof=ifdof/10  !ifdof=1
    !    end do
    !end do
    !*** FIN Version Owen

    do itotv=1,ntotv
        iffix(itotv)=0
    end do

    !*** Version Victor pour 2D
    do ivfix=1,nvfix
        read(5,*) nofix(ivfix),ifpre,(presc0(ivfix,idofn),idofn=1,ndofn) !8
        write(6,908) nofix(ivfix),ifpre,(presc0(ivfix,idofn),idofn=1,ndofn) !8
        nloca=(nofix(ivfix)-1)*ndofn+1 ! Donne le numéro de ligne dans iffix de la premiere composante (selon x) du noeud nofix(ivfix)

        if (ifpre==10 .or. ifpre==11) then
            iffix(nloca)=1
        end if
        
        if (ifpre==01 .or. ifpre==11) then
            iffix(nloca+1)=1
        end if
    end do
    !*** FIN Version Victor

   

    908 format(1x,i4,5x,i5,5x,f10.4,1x,f10.4)

    ! Read the available selection of element properties
    write(6,923)
    923  format( //20hMaterials properties)
    16  write(6,910)
    910  format(6hNumber,6x,19h Element properties)

    
        do imats=1,nmats
            read(5,900) numat !9
            read(5,*)(props(numat,iprop),iprop=1,nprop) !10
    !     930 format(8f10.5)
            write(6,*) numat,(props(numat,iprop),iprop=1,nprop)
        end do
         

!     911 format(1x,i4,3x,8e14.6)

!       911 format(8h Young =,f10.2,4x,11h Poisson = ,f4.2,4x,14h Thickness t = ,f3.1,4x,   &
!                      10h Density =,f5.1,4x,15h Sigma_y\c\py =,f7.2,4x,4h H =,f5.2,4x,   &
!                      6h phi =,f6.2/,                                                    &
!                      /6h v0 = ,f5.3,4x,10h Lambda = ,f6.4,4x,9h Kappa = ,f6.4,4x,5h M = ,f6.4)

    ! Version plus lisible: 911 format(1x,i4,3x,8f14.2)

    !911 format(//8h Young =,i4,4x,8h Poisson ratio =,i4,4x,8h Thickness t =,i4,4x,        &
    !                  8h Density rho =,i4,4x,8h Sigma_y/Cohesion c =,i4,                       &   
    !                  //8h Hardening parameter =,i4,4x,8h Angle of friction phi (°)= ,i4)
    !911 format(//11h Materiau =,i4,4x,8h Young =,f5.2,8hPoisson=,4f5.2,4x,8hThickness,i4,4x)

    !write(6,*) props(1,:)

  

    ! ****** Version Elkassas
    ! 11 Read initial conditions
    !do j=1,npoin
    !    read(5,*) ipoin,(cond(ipoin,i),i=1,ndofn) !11 ATTENTION MODIFICATION de l'iterateur
    !end do                        

    !Put in vector of displacements
    !do ipoin=1,npoin
    !    do idofn=1,ndofn
    !        nloca=ndofn*(ipoin-1)+idofn
    !        asdis(nloca)=cond(ipoin,idofn)
    !    end do
    !end do
    ! ****** Fin Version Elkassas



    ! Set up gaussian integration constants

    !call gaussq(ngaus,posgp,weigp)
    call gaussq( )


    !call check2(coord,iffix,lnods,matno,melem,mfron,mpoin,mtotv,mvfix,ndfro,ndofn,nelem,nmats,nnode,nofix,npoin,nvfix)
    call check2( )


    ! On cherche le nombre de conditions aux limites, ie nombres de noeuds bloqués pour dimensionner ultérieurement rgstif, rgload
        ! D'abord on cherche le nombre de directions qui sont bloquées
    nnodefix=0
    do i=1,ntotv
        if (iffix(i)==1) then
            nnodefix=nnodefix+1
        end if
    end do

        ! Ensuite on regarde s'il y a des déplacements imposes
!     do i=1,nvfix
!         do j=1,ndim
!             if (presc0(i,j)/=0.d0) then
!                 nnodefix=nnodefix-1
!             end if
!         end do
!     end do

        ! On en déduit le nombre de noeuds libres = inconnues du systeme
    nvlib=ntotv-nnodefix
    allocate(gcoord(nvlib))    !Tableau 28
                

    return

end














