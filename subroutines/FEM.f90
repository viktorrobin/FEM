!**********************************************************************************************
!***********************************************************************************************
!
!   This file contains subroutines for general FE procedures
!
!   By Victor ROBIN - 09.06.2013, Exeter
!
!**********************************************************************************************
!**********************************************************************************************

!************************************************************************
!   Subroutine greduc - Owen p.52
!   Modify global load vector to account prescribed displacements
!***********************************************************************

subroutine greduc()

use variables
implicit none

integer :: ievab, irows, ieqn1, kount, icols, jevab
double precision :: factr, pivot


do ievab=1,ntotv
    if (iffix(ievab)==1) then
        do irows=1,ntotv
                ! Modify gload
                gload(irows)=gload(irows)-gstif(irows,ievab)*fixed(ievab)
        end do
    end if
!     70 continue
end do

! stop

end

!************************************************************************
!   Subroutine prescdisp - Owen p.52
!   Modify global load vector to account prescribed displacements during iteration process
!***********************************************************************

subroutine prescdisp()

use variables
implicit none

integer :: ievab, irows, ieqn1, kount, icols, jevab,i
double precision :: factr, pivot


! We copy the unmodified residual vector
do ievab=1,nvlib
    vec_residu_pd(ievab)=vec_residu(ievab)
end do


! We add the prescribed displacements
kount=1
do ievab=1,ntotv
    if (iffix(ievab) == 1 .and. fixed(ievab) /= 0.d0) then

        do i=1,nvlib
            irows=gcoord(i)
                vec_residu_pd(i)=vec_residu_pd(i)-gstif(irows,ievab)*fixed(ievab)
                kount=kount+1
        end do
    end if
end do
return
end

!************************************************************************
!   Subroutine modps - Owen p.192
!   Evaluates the elasto-plastic matrix D for either plane stress for either plane stress, plane strain, or axisymmetric situations.
!   D matrix is stored in dmatx.
!***********************************************************************

subroutine modps(lprop)

    use variables

    implicit none

    integer :: istr1, jstr1, lprop, kgaus, istre,jstre
    double precision :: const, conss

    young=props(lprop,1)    ! Young modulus
    poiss=props(lprop,2)    ! Poisson's ratio

    ! Set to zero elasto-plastic matrix
    dmatx=0.d0

    !*** ntype = 1 ---> Plane stress
    !*** ntype = 2 ---> Plane strain
    !*** ntype = 3 ---> Axial symmetry

    ! Plane stress
    if(ntype == 1) then
        !D matrix for plane stress case
        const=young/(1.d0-poiss*poiss)
        dmatx(1,1)=const
        dmatx(2,2)=const
        dmatx(1,2)=const*poiss
        dmatx(2,1)=const*poiss
        dmatx(3,3)=const*(1.d0-poiss)/2.d0
        return
    end if

    ! Plane strain
    if(ntype == 2) then
        !D matrix for plane strain case
        const=young*(1.d0-poiss)/((1.d0+poiss)*(1.d0-2.d0*poiss))
        dmatx(1,1)=const
        dmatx(2,2)=const
        dmatx(1,2)=const*poiss/(1.d0-poiss)
        dmatx(2,1)=const*poiss/(1.d0-poiss)
        dmatx(3,3)=const*(1.d0-2.d0*poiss)/(2.d0*(1.d0-poiss))
        return
    end if

    ! Axial symmetry
    if(ntype == 3) then
        !D matrix for axisymmetric case
        const=young*(1.d0-poiss)/((1.d0+poiss)*(1.d0-2.d0*poiss))
        conss=poiss/(1.d0-poiss)
        dmatx(1,1)=const
        dmatx(2,2)=const
        dmatx(3,3)=const*(1.d0-2.d0*poiss)/(2.d0*(1.d0-poiss))
        dmatx(1,2)=const*conss
        dmatx(1,4)=const*conss
        dmatx(2,1)=const*conss
        dmatx(2,4)=const*conss
        dmatx(4,1)=const*conss
        dmatx(4,2)=const*conss
        dmatx(4,4)=const
    end if


    return
end


!************************************************************************
!   Subroutine bmatps - Owen p.191
!   This subroutine evaluates the strain matrix B at any position within an element.
!   The B matrix is stored in bmatx.
!*************************************************************************

subroutine bmatps(kgasp)

    use variables

    implicit none
    integer :: ngash, mgash, inode, i, j, kgasp

    ngash=0

    do inode=1,nnode
        mgash=ngash+1
        ngash=mgash+1
        bmatx(1,mgash)=cartd(1,inode)
        bmatx(1,ngash)=0.d0
        bmatx(2,mgash)=0.d0
        bmatx(2,ngash)=cartd(2,inode)
        bmatx(3,mgash)=cartd(2,inode)
        bmatx(3,ngash)=cartd(1,inode)

        if(ntype == 3) then
            bmatx(4,mgash)=shap(inode)/gpcod(1,kgasp)
            bmatx(4,ngash)=0.d0
        end if
    end do

    return
end



!*******************************************************************************
!   Subroutine stiffp( ) - Owen p.244
!   Evaluates the stiffness matrix for each element in turn. It uses elasto-plastic matrix Dep (7.47).
!*******************************************************************************


subroutine stiffp( )

    use variables
    implicit none

    integer :: inode, idime, ievab, jevab, igaus, jgaus, iposn, lnode, istr1, istre, jstre
    double precision :: dvolu, etasp, exisp
    integer :: i, j, kgaus, idofn, ielem, kgasp, lprop

    include 'subroutines/interface.f90'

    rewind(1)
    kgaus=0

    ! Loop over each element
    do ielem=1,nelem
        lprop=matno(ielem)

        ! Evaluate the coordinates of the element nodal points
        do inode=1,nnode
            lnode=iabs(lnods(ielem,inode))
            iposn=(lnode-1)*2

            do idime=1,ndim
                iposn=iposn+1
                elcod(idime,inode)=coord(lnode,idime)
            end do
        end do

        ! Read thickness of the plate
        thick=props(lprop,3)

        ! Initialize the element stiffness matrix estif
        do ievab=1,nevab
            do jevab=1,nevab
                estif(ievab,jevab)=0.d0
            end do
        end do

        kgasp=0

        ! Enter loops for for area numerical integration
        do igaus=1,ngaus
            exisp=posgp(igaus)
            do jgaus=1,ngaus
                etasp=posgp(jgaus)
                kgasp=kgasp+1
                kgaus=kgaus+1

                ! Evaluate the D-Matrix
                call modps(lprop)

                ! Evaluate the shape functions, elemental volume, etc...
                call sfr2(etasp,exisp)

                ! Evaluate cartesian coordinates, Jacobian Matrix/Inverse/Determinant,cartesian derivatives of element shape functions
                call jacob2(ielem,kgasp)

                ! End of the integral to calculate stiffness matrix using Gauss-Legendre integration
                dvolu=djacb*weigp(igaus)*weigp(jgaus)

                ! If axisymmetry problem (Triaxial)
                if (ntype==3) then
                    dvolu=dvolu*twopi*gpcod(1,kgasp)
                end if

                ! If plate case
                if (thick /= 0.d0) then
                    dvolu=dvolu*thick
                end if

                ! Evaluate the B matrices
                call bmatps(kgasp)

                dbmat=matmul(dmatx,bmatx)

                ! Calculate the element stiffness
                do ievab=1,nevab
                    do jevab=ievab,nevab
                        do istre=1,nstre
                            estif(ievab,jevab)=estif(ievab,jevab)+bmatx(istre,ievab)*dbmat(istre,jevab)*dvolu
                        end do
                    end do
                end do

            end do
        end do

        ! Construct the lower triangle of the stiffness matrix
        do ievab=1,nevab
            do jevab=1,nevab
                estif(jevab,ievab)=estif(ievab,jevab)
            end do
        end do


        if (nnode==3) then
            estif=estif/4.d0
        end if

        ! Store the stiffness matrix, stress matrix, and sampling point coordinates for each element on disc file
        write(1) estif

    end do

    return
end

!***********************************************************************
!   Subroutine dbe - Owen p.194
!   this subroutine multiplies the D-matx by the B-matx
!   Possible replacement by matmull?
!*******************************************************************************

subroutine dbe( )

use variables

implicit none

integer :: istre, ievab, jstre, i, j


dbmat=matmul(dmatx,bmatx)
return
do istre=1,nstre
    do ievab=1,nevab
        dbmat(istre,ievab)=0.d0
        do jstre=1,nstre
            dbmat(istre,ievab)=dbmat(istre,ievab)+dmatx(istre,jstre)*bmatx(jstre,ievab)
        end do
    end do
end do

return
end



!*******************************************************************************
!   Subroutine assem( ) - Based on front( ), Owen p.194
!   This subroutine assembles the contributions of each element to form the global stiffness matrix and global load vector
!       and to solve the resulting set of simultaneous equations by Gaussian direct elimination.
!   The main feature of the frontal solution technique is that it assembles the equations and eliminates the variables at
!       the same time.
!*******************************************************************************

!subroutine assemb(gstif,gstif_temp,rgstif,gload, gload_temp, asdis_temp,rgload,rasdis,gcoord,a,y,x)
subroutine assemb( )

    use variables
    implicit none


    integer :: i, j, ibufa, istif, ifron, idofn, nikno, ievab, nfron
    integer :: kevab, inode, ipoin, kelva, kexis, klast, locno, nbufa, nlast
    integer :: nposi, idest, jdest, ngash, ngish, jevab, jfron, nloca
    integer :: kboun, itotv, ielva, lfron, mgash, ngush, nfunc

    integer:: nrowe,nrows, nodej,jnode,nodei,ncole,ncols,jdofn
    integer :: fehler, compteur, loca, loca2, ielem

    !  ***      DYNAMIC ALLOCATION OF MEMORY  ****
    if (iincs==1) then
        include 'workspace/Alloc2.f90'
    end if
    ! Reposition at the beginning of the file storing element stiffness matrices
    rewind(1)

    ! Set to zero global vectors and matrices
    gload=0.d0
    asdis=0.d0
    gstif=0.d0

    ! Set to zero reduced vector and matrices
    rasdis=0.d0
    rgload=0.d0
    gcoord=0.d0
    rgstif=0.d0

    !Assembly global stiffness matrix
    do ielem=1,nelem
        read(1) estif
        do inode=1,nnode
            nodei=lnods(ielem,inode)
            do idofn=1,ndofn

                !Assembly gload
                nrows=(nodei-1)*ndofn+idofn     !Row in the global stiffness matrix and load vector
                nrowe=(inode-1)*ndofn+idofn     !Row in the element stiffness matrix and load vector
                gload(nrows)=gload(nrows)+eload(ielem,nrowe)

                !Assembly gstif
                if (kresl==2) cycle
                do jnode=1,nnode
                    nodej=lnods(ielem,jnode)
                    do jdofn=1,ndofn
                        ncols=(nodej-1)*ndofn+jdofn
                        ncole=(jnode-1)*ndofn+jdofn
                        gstif(nrows,ncols)=gstif(nrows,ncols)+estif(nrowe,ncole)
                    end do
                end do
            end do
        end do
    end do

!     Display global stiffness matrix
!     call disp_gstif()

    return
end


!*******************************************************************************
!   Subroutine reduc( ) - Victor
!   this subroutine reduces stiffness matrix and load vector using boundary conditions
!*******************************************************************************

subroutine reduc( )

    use variables
    implicit none

    integer :: i, j, loca, loca2, compteur

    compteur=0
    do i=1,ntotv
        if (iffix(i)==0) then
!         if (iffix(i) == 0.d0 .or. fixed(i)/=0.d0) then
            compteur=compteur+1
            gcoord(compteur)=i
        end if
    end do


    ! Reduce global stiffness matrix gstif to rgstif
    do i=1,nvlib
        loca=gcoord(i)
        do j=1,nvlib
            loca2=gcoord(j)
            rgstif(i,j)=gstif(loca,loca2)
        end do
    end do


    ! Reduce vector rload to rgload
    do i=1,nvlib
        loca=gcoord(i)
        rgload(i)=gload(loca)
    end do

    !Reduce fixity vector
    do i=1,nvlib
        loca=gcoord(i)
        riffix(i)=iffix(loca)
    end do

    !Reduce prescribed displacements vectors
    do i=1,nvlib
        loca=gcoord(i)
        rfixed(i)=fixed(loca)
    end do


    return
end


! **************************************************
!   Subroutine increm - Owen p.210
!   This subroutine increments the applied load or any prescribed displacements according to the load factors specified as input.
!   The final role is to insert appropriate values in the fixity array to control any prescribed displacements.
! **************************************************

subroutine increm( )

    use variables

    implicit none

    integer :: ievab, itotv, nloca, ngash, idofn, ivfix, ielem
    double precision :: dnincs
    ! Read load factor for this increment, toler (useless), maximum number of iteration, and parameters for output (useless)
!     read(5,950) facto,toler,miter,noutp(1),noutp(2)
!     950 format(2f10.5,3i5)
    noutp(1)=0
    noutp(2)=0

    if (iincs==nincs) then
        noutp(1)=3
        noutp(2)=3
    end if


        dnincs=nincs
        facto=1/dnincs

    ! Update total load factor
    tfact=tfact+facto

    write(6,960) iincs, 100*tfact, epsilon!,miter!,noutp(1),noutp(2)
    960 format(/11hIncrement= ,i6,5x,14hload factor = ,f6.2,' %',5x,24h convergence tolerance =,e12.4)!,5x,24hMax. No. of iterations =,i5)!, &
                        !//27h initial output paramater =,i5,5x,24hfinal output paramater =,i5)


    do ielem=1,nelem
        do ievab=1,nevab
            eload(ielem,ievab)=rload(ielem,ievab)*facto
            tload(ielem,ievab)=rload(ielem,ievab)*facto
        end do
    end do


    ! Interpret fixity data in vector form
    do itotv=1,ntotv
        fixed(itotv)=0.d0
    end do

    ! Flag to know if prescribed displacements
    drap_pres=0
    do ivfix=1,nvfix
        nloca=(nofix(ivfix)-1)*ndofn
        do idofn=1,ndofn
            ngash=nloca+idofn
            fixed(ngash)=presc0(ivfix,idofn)*facto
!             print *,fixed(ngash)
            if (presc0(ivfix,idofn)/=0.d0) then
                drap_pres=1
            end if

        end do
    end do

    return
end



!*******************************************************************************
!   Subroutine residu( ) - Owen p.249
!   This subroutine evaluates the nodal forces which are statically equivalent to the stress field satisfying
!       elasto-plastic conditions. Comparison of these equivalent nodal forces with the applied loads gives the residual
!       forces.
!*******************************************************************************

subroutine residu( )

    use variables
    implicit none

    integer :: ievab, igaus, jgaus, idofn, inode, nposn, lnode, i, j
    integer :: istr1, istep, istre, mgash, mstep, kgaus, ielem, kgasp, lprop

    double precision :: hards, uniax, dvolu, agash, pi, etasp, exisp
    double precision :: astep, bgash, bring, curys, dlamd, escur, espre, preys, reduc, rfact, bgashCC
    double precision :: part1, part2, part3, part4

    double precision   , dimension (4) :: strsg_pre            !desig(4)

    include 'subroutines/interface.f90'

    plast_code=0
    ! Zero array in which the equivalent nodal forces, calculated at the end of the subroutine, will be stored.
    do ielem=1,nelem
        do ievab=1,nevab
            eload(ielem,ievab)=0.d0
        end do
    end do


    !Zero gauss point counter over all elements
    kgaus=0

    ! Loop over each element
    do ielem=1,nelem

        !Identify the element material property number
        lprop=matno(ielem)

        !Identify the initial uniaxial yield stress sigma_Y° (or c for Mohr-Coulomb material or Drucker-Prager criteria)
        uniax=props(lprop,5) !Yield Stress
        hards=props(lprop,6)
        frict=props(lprop,7)

        if (ncrit==5 .or. ncrit==6) then
            CCkappa=props(lprop,10)
            CClambda=props(lprop,9)
            csl=props(lprop,11)
            Nlambda=props(lprop,8)
            Nkappa=Nlambda-(CClambda - CCkappa)*dlog(uniax)
            pb=props(lprop,16)
            pyI=props(lprop,5)
            pyII=props(lprop,12)
            Dei=props(lprop,14)
            Dec=props(lprop,15)
            beta=props(lprop,13)
        end if

        ! Mohr-Coulomb: evaluates the equivalent yield stress as c.cos(phi)
        if (ncrit == 3) then
            uniax=props(lprop,5)*dcos(frict*degres)
        end if

        ! Drucker-Prager: evaluates the equivalent yield stress as k' according to 7.18
        if (ncrit ==4) then
            uniax=6.d0*props(lprop,5)*dcos(frict*degres)/(root3*(3.d0-dsin(frict*degres)))
        end if



        ! Store nodal coordinates in elcod and the nodal displacements due to the application of the residual forces in array eldis
        do inode=1,nnode
            lnode=iabs(lnods(ielem,inode))
            nposn=(lnode-1)*ndofn
            do idofn=1,ndofn
                nposn=nposn+1
                elcod(idofn,inode)=coord(lnode,idofn)
                eldis(idofn,inode)=asdis(nposn)
            end do
        end do

!         name_disp= "eldis"
!         call disp_mat(eldis)

!         name_disp='elcod'
!         call disp_mat(elcod)

        !Evaluate the elastic D matrix
        call modps(lprop)

!         name_disp='D matrix'
!         call disp_mat(dmatx)

        !Identitfy the element thickness
        thick=props(lprop,3)

        !Zero local gauss point counter
        kgasp=0

        !Enter loop for numerical integration
        do igaus=1,ngaus
            do jgaus=1,ngaus
                !write(6,*) igaus,jgaus
                exisp=posgp(igaus)
                etasp=posgp(jgaus)
                kgaus=kgaus+1
                kgasp=kgasp+1

                !Evaluate shape functions Ni and their derivatives dNi/d_ksi, etc...
                call sfr2(etasp,exisp)

                !Evaluate the Gauss point coordinates, the determinant of the Jacobian matrix, and the cartesian derivatives of the shape functions
                call jacob2(ielem,kgasp)

                !Calculate the elemental volume for numerical integration
                dvolu=djacb*weigp(igaus)*weigp(jgaus)

                !Axisymmetric problem
                if (ntype == 3) then
                    dvolu=dvolu*twopi*gpcod(1,kgasp)
                end if

                !if thickness (plane stress?)
                if (thick /= 0.d0) then
                    dvolu=dvolu*thick
                end if

                !Evaluates B matrix for the Gauss Point
                call bmatps(kgasp)

!                 name_disp='bmatx'
!                 call disp_mat(bmatx)


                !Computes stress increment stres(istr1) assuming ELASTIC behaviour
                    ! 1) strain components thanks to B matrix: {Epsilon}=[B].{Un}
                    ! 2) stress components thanks to Hooke: {sigma}=[D].{Epsilon}
                    ! Renvoie matrices stres(sigma) and stran(epsilon)
                call linear(lprop,kgasp,kgaus)



!                 name_disp='Tenseur des deformations'
!                 call disp_vec(stran)

!                 name_disp='Tenseur des contraintes elastic'
!                 call disp_vec(stres)

                ! Compute yield stress for (r-1)th iteration
                if (ncrit < 5) then
                    ! Compute yield stress for (r-1)th iteration as: sigma_y +H'*eps_p^(r-1)
                    preys = uniax+epstn(kgaus)*hards
                else if (ncrit == 5) then
                    ! Compute yield stress for (r-1)th iteration for Cam Clay
                    ! For first increment, we assume elastic behaviour so no need to calculate specific volume (epstnp=0)
!                    print *,"epstnp=",epstnp(kgaus),"vol spe=",volspegp(kgaus)
                   preys = uniax*dexp( epstnp(kgaus)*(volspegp(kgaus)/(CClambda - CCkappa)) )
!                    print *, 'preys=',preys,'ieldgp=',yieldgp(kgaus)
!                    print *,"uniax=",uniax,"preys=",preys,"epstnp(kgaus)=",epstnp(kgaus),"v=",volspegp(kgaus)
                else if (ncrit == 6) then
                    preys = yieldgp(kgaus)
                end if


                strsg_pre=0.d0

                !Store desigma_e^r et sigma_e^r
                !strsg: stock état de contrainte précédent sigma_e^(r-1)? vaut 0.0 à la premiere iteration (initialisé dans subroutine ZERO)
                do istr1=1,nstr1
                    desig(istr1)=stres(istr1)   !dsigma_e^r = increment de contrainte a iteration r
                    strsg_pre(istr1)=strsg(istr1,kgaus) !pour pouvoir evaluer le vrai desig a la fin du processus de projection
                    sigma(istr1)=strsg(istr1,kgaus)+stres(istr1)  !sigma_e^r = etat de contrainte actuel à iteration r
                end do

                !** xxxx  /!\ PLASTICITY AREA /!\ xxxx **!

                !Evaluate the effective stress and store as YIELD
                call invar(lprop,sigma)

                !Check if Gauss point had yielded on the the previous iteration = first operation of step d
                espre=effst(kgaus)-preys

                    !if >0, Gauss point had yielded on the previous iteration: YES p.250
                    !if <0, Gauss point had not previously yielded: NO p.250

                if (espre >= 0.d0) then
                    !Gauss point had yielded previously

                    escur=yield-effst(kgaus)

                    if (escur <= 0.d0) then
                        !Gauss point unloading elastically: R=0
!                         rfact=0.d0
                        go to 60

                    else if (escur > 0.d0) then
                        !Gauss point had yielded previously and stress is still increasing
                        rfact=1.d0
                    end if

                else if (espre < 0.d0) then
                    !If Gauss point had not yielded previously (was elastic)
                    !Check to see if it has yielded during this iteration

                    escur=yield-preys

                    if (escur <= 0.d0) then
                        !Gauss point still elastic
                        go to 60

                    else if (escur > 0.d0) then
                        !Gauss point has yielded during this iteration
                        rfact = escur/(yield - effst(kgaus))
                    end if

                end if

!                 print *,"preys=",preys,"effst(kgaus)=",effst(kgaus),"espre=",espre,'yield=',yield,'escur=',escur,'fract=',rfact

                !Evaluate number of steps into which the excess stress (plastic part) R*dsigma_e^r is to be divided
                mstep = (escur*8.d0/uniax) + 1.d0
                astep = mstep

                !Calculate portion of stress which satisfies yield criterion (sgtot) and the one which has to be reduced on F (stres)
                do istr1=1,nstr1
!                     sgtot(istr1) = strsg(istr1,kgaus) + (1.d0 - rfact)*stres(istr1)
!                     stres(istr1) = rfact*stres(istr1)/astep
                    sgtot(istr1) = strsg(istr1,kgaus) + (1.d0 - rfact)*desig(istr1)
                    stres(istr1) = rfact*desig(istr1)/astep !Stress to be applied at each loop
                end do


                !Loop over each stress reduction step
                do istep=1,mstep
                    !Compute vectors a and d_D
                    call invar(lprop,sgtot)
                    call yieldf(lprop,kgaus)
                    call flowpl(lprop,kgaus)

                    !Compute dLambda and store as dlamd
                    agash=0.d0


                            do istr1=1,nstr1
                                agash = agash + avect(istr1)*stres(istr1)
                            end do

                            dlamd = agash*abeta

                            if (dlamd < 0.d0) then
                                dlamd = 0.d0
                            end if





                    !Compute the value of the corrected stress
                    !However, it is not yet on the yield function!
                    bgash=0.d0
                    bgashCC=0.d0

                    do istr1=1,nstr1
                        bgash = bgash + avect(istr1)*sgtot(istr1)
                        bgashCC = bgashCC + avect(istr1)
                        sgtot(istr1) = sgtot(istr1) + stres(istr1) - dlamd*dvect(istr1)
                    end do


                    !Compute effective plastic strain to scale new stress and bring it back on F
                    if (ncrit<5) then
                        epstn(kgaus) = epstn(kgaus) + dlamd*bgash/yield
                    else if (ncrit==5 .or. ncrit==6) then
                        epstn(kgaus) = epstn(kgaus) + dlamd*bgash/yield
                        epstnp(kgaus) = epstnp(kgaus) - dlamd*bgashCC !'-' because Cam Clay
!                         print *,"dlamd=",dlamd,"bgashCC=",bgashCC,"ds=",stres(1),"v=",volspegp(kgaus),'epp=',epstnp(kgaus)
                    end if

                    !Scale the obtained stress so it lies on F
                    !Additional refinement!
                    call invar(lprop,sgtot)

!                     curys = uniax + epstn(kgaus)*hards

                    if (ncrit < 5) then
                        ! Compute yield stress for (r-1)th iteration as: sigma_y +H'*eps_p^(r-1)
                        curys = uniax+epstn(kgaus)*hards
                    else if (ncrit == 5) then
                        ! Compute yield stress for (r-1)th iteration for Cam Clay
    !                    print *,"epstnp=",epstnp(kgaus),"vol spe=",volspegp(kgaus)
                       curys = uniax*dexp( epstnp(kgaus)*(volspegp(kgaus)/(CClambda - CCkappa)) )
    !                    print *,"preys=",preys
                    else if (ncrit == 6) then
                        part1 = (dlamd*bgashCC)*volspegp(kgaus)
                        part2 = exp(beta*yieldgp(kgaus))*beta*(exp(pyI*beta) + exp(pyII*beta))*(Dec-Dei)
                        part3 = (exp(yieldgp(kgaus)*beta)+exp(pyII*beta))*(exp(yieldgp(kgaus)*beta)+exp(pyII*beta))
                        part4 = part1/( -part2/part3 + (CClambda-CCkappa)/yieldgp(kgaus) )
                        curys = yieldgp(kgaus) - part4
!                         print *,'curys=',curys,"beta=",beta,"pyII=",pyII,'pyI=',pyI
                    end if

!                     print *,'pb=',pb,"Dei=",Dei,"Dec=",Dec
!                     print *,'curys=',curys
!                     print *,'epstnp=',epstnp(kgaus),'epstn=',epstn(kgaus),'error=',(epstnp(kgaus)-epstn(kgaus))*100.d0/epstn(kgaus)
!                     print *,'epstn=',epstn(kgaus)
!                     print *,"agash=",agash,'abeta=',abeta
!                     print *,"stran=",stran
!                     print *,'dvect=',dvect
!                     print *,'avect=',avect
!                     print *,"sgtot=",sgtot
!                     print *,'I1=',sgtot(1)+sgtot(2)+sgtot(4)
!                     print *,'stres=',stres
!                     print *,'I1/3=',I1/3.d0,"smean=",smean,'varj2=',varj2,'dlamd=',dlamd,"yield=",yield,'v=',volspegp(kgaus)
!                     pause

                    bring = 1.d0


                    if (yield > curys) then
                        !the corrected stress is still outside F so we scale it and bring it on F
                        bring = curys/yield
                    end if


                    !On incremente le vecteur contrainte corrigé au tenseur des contraintes total
                    !Apres m steps, on aura incremente toute la partie plastique du vecteur contrainte calcule au debut
                    do istr1=1,nstr1
!                         desig(istr1) = bring*sgtot(istr1) - strsg(istr1,kgaus)
                        strsg(istr1,kgaus) = bring*sgtot(istr1)
                    end do

                    !Store effective stress sigma^r in effst
                    effst(kgaus) = bring*yield

                    if (ncrit==5 .or. ncrit==6) then
                        yieldgp(kgaus)=curys
                    end if
                end do


                !We calculate the actual stress increment after the projection of the yield function
!                 desig=0.d0

                do istr1=1,nstr1
                    desig(istr1) = strsg(istr1,kgaus)-strsg_pre(istr1)
                end do

                if (ncrit==5 .or. ncrit==6) then
                    ! For Cam Clay, update the new specific volume
                    volspegp(kgaus)=volspegp(kgaus)*( 1.d0 + (stran(1)+stran(2)+stran(4)) )
                    call invar(lprop,strsg)
                    deviatoric_stress(kgaus)=dsqrt(3.d0*varj2)
                    effective_mean_stress(kgaus)=-smean
                end if

                ! Variable to know if plasticity has occured during this iteration
                plast_code=1



                go to 190

                !** xxxx  /!\ END PLASTICITY AREA /!\ xxxx **!

                60 continue !Exit point in case Gauss point is not yielding and is elastic

                !Dans le cas ou il n'y a pas eu de plasticité, on rajoute l'increment de contrainte au tenseur des contraintes totales
                do istr1=1,nstr1
                    !strsg(istr1,kgaus)=strsg(istr1,kgaus)+desig(istr1)
                    strsg(istr1,kgaus)=strsg(istr1,kgaus)+desig(istr1)
                end do

                effst(kgaus) = yield


                ! For Cam Clay, Calculate new specific void ratio
                if (ncrit==5) then
                    call invar(lprop,strsg)

                    if (volspegp(kgaus)==0.d0) then
                        volspegp(kgaus)=Nkappa - CCkappa*dlog(dabs(smean))
                    else
!                         print *,(stran(1)+stran(2)+stran(4))
!                         print *,volspegp(kgaus)*( 1.d0 + (stran(1)+stran(2)+stran(4)) )
                        volspegp(kgaus)=volspegp(kgaus)*( 1.d0 + (stran(1)+stran(2)+stran(4)) )
                    end if
!                     yieldgp(kgaus)=yield
                        yieldgp(kgaus)=uniax*dexp( epstnp(kgaus)*(volspegp(kgaus)/(CClambda - CCkappa)) )
!                     print *,yield
!                     volspegp(kgaus)=volspegp(kgaus)*( 1.d0 + (stran(1)+stran(2)+stran(4)) )
!                     print *,"volspe=",volspegp(kgaus),"p=",dabs(smean)
!                     print *,'epstnp=',epstnp(kgaus)
!                     print *,'stran=',stran
!                     print *,(stran(1)+stran(2)+stran(4))
                        call invar(lprop,strsg)
                        deviatoric_stress(kgaus)=dsqrt(3.d0*varj2)
                        effective_mean_stress(kgaus)=-smean

                end if

                ! TO MODIFY
                if (ncrit==6) then
                    call invar(lprop,strsg)
                        volspegp(kgaus)=Nkappa - CCkappa*dlog(dabs(smean))


                end if


                190 continue

                !Initialisation de bmatx_tr
                do i=1,size(bmatx_tr,1)
                    do j=1,size(bmatx_tr,2)
                        bmatx_tr(i,j)=0.d0
                    end do
                end do


                !Calcul de la transposé de B au point de Gauss pour pouvoir recalculer les forces nodales à incrementer
                bmatx_tr=transpose(bmatx)


                !eload_temp=matmul(bmatx_tr,strsg(:,kgaus))
                eload_temp=matmul(bmatx_tr,desig)

                ! Version Victor
                do ievab=1,nevab
                    eload(ielem,ievab)=eload(ielem,ievab)+eload_temp(ievab)*dvolu
                end do


!                 Version Owen
!                 do inode=1,nnode
!                     do idofn=1,ndofn
!                         mgash=mgash+1
!                         do istre=1,nstre
!                             eload(ielem,mgash)=eload(ielem,mgash)+bmatx(istre,mgash)*strsg(istre,kgaus)*dvolu
!                         end do
!                     end do
!                 end do


            end do
        end do


    end do !ielem


    return


end


! ****************************************************************************
!       subroutine conver() - Owen. p.213
!       This subroutine checks for convergence of the iteration forces
!       See CONUND p.72
! ****************************************************************************


subroutine conver( )




    use variables
    implicit none

    integer :: itotv, kevab, inode, locno, idofn, nposi, ievab, ielem, ivlib, loca

    double precision :: refor, resid, retot, agash, maxre, ratio, temp, epsilon2
!     double precision, dimension(:),allocatable :: rstfor
    double precision, dimension(nvlib) :: temp_vec
    include 'subroutines/interface.f90'


temp_vec=0.d0
    ncheck=0
    resid=0.d0
    retot=0.d0

    ! Set to zero stfor and tofor
    do itotv=1,ntotv
        stfor(itotv)=0.d0
    end do

    do ielem=1,nelem
        kevab=0
        do inode=1,nnode
            locno=iabs(lnods(ielem,inode))
            do idofn=1,ndofn
                kevab=kevab+1
                nposi=(locno-1)*ndofn+idofn
                stfor(nposi)=stfor(nposi)+eload(ielem,kevab) !Vecteur colonne contenant les forces nodales appliquées a cette iteration
                tofor(nposi)=tofor(nposi)+eload(ielem,kevab) !Vecteur colonne contenant les forces nodales appliquées depuis le debut du processus
            end do
        end do
    end do


    do ivlib=1,nvlib
        loca=gcoord(ivlib)
        rstfor(ivlib)=stfor(loca)
    end do


    ! Calculate new residual vector if no prescribed displacements
    do ivlib=1,nvlib
        vec_residu(ivlib)=( vec_residu(ivlib) - rstfor(ivlib) )
    end do

    norme=0.d0
    maxre=0.d0  !MAXimum REsidual in the residual forces vector

    ! Calculate norm of residual vector
    do ivlib=1,nvlib
        norme = norme + (vec_residu(ivlib)*vec_residu(ivlib))

        if (vec_residu(ivlib) > maxre) then
            maxre = vec_residu(ivlib)
        end if

    end do

    norme=dsqrt(norme)


    if (norme <= epsilon) then !Solution has converged
        ncheck=0
    end if

    if (norme > epsilon) then
        if (iiter == 1) then
            ncheck = 1             !Solution converging
        end if

        if (iiter > 1 .and. norme < pvalue) then
            ncheck = 1              !Solution converging
        end if

        if (iiter > 1 .and. norme >= pvalue) then
            ncheck = 999
        end if

        pvalue=norme
    end if


! Fin Timer
if (ncheck==0 .or. ncheck==999) then
    call cpu_time(finish_inc)

        write(6,30) iiter, ncheck, plast_code, norme, maxre, finish_inc-start_inc
    30 format(5x11hIteration =,i7,3x,18hconvergence code =,i3,3x,17hplasticity code =,i3, &
                & 3x,28hnorm of residual sum ratio =,e14.6,3x,18hmaximum residual =,e14.6,3x,6hTime =,f6.3,9h seconds.)
end if


end


! ****************************************************************************
!       subroutine conver() - Owen. p.213
!       This subroutine checks for convergence of the iteration forces
!       See CONUND p.72
! ****************************************************************************


subroutine conver2( )
 use variables
    implicit none

    integer :: itotv, kevab, inode, locno, idofn, nposi, ievab, ielem, ivlib, loca

    double precision :: refor, resid, retot, agash, maxre, ratio, temp, epsilon2
    double precision, dimension(nvlib) :: diff,pre_vec
    include 'subroutines/interface.f90'



    ncheck=0
    resid=0.d0
    retot=0.d0

    ! Set to zero stfor and tofor
    do itotv=1,ntotv
        stfor(itotv)=0.d0
        !tofor(itotv)=0.d0
    end do

    do ielem=1,nelem
        kevab=0
        do inode=1,nnode
            locno=iabs(lnods(ielem,inode))
            do idofn=1,ndofn
                kevab=kevab+1
                nposi=(locno-1)*ndofn+idofn
                stfor(nposi)=stfor(nposi)+eload(ielem,kevab) !Vecteur colonne contenant les forces nodales appliquées a cette iteration
                tofor(nposi)=tofor(nposi)+eload(ielem,kevab) !Vecteur colonne contenant les forces nodales appliquées depuis le debut du processus
            end do
        end do
    end do


    do ivlib=1,nvlib
        loca=gcoord(ivlib)
        rstfor(ivlib)=stfor(loca)
    end do


    pre_vec=vec_residu

    ! Calculate new residual vector if no prescribed displacements
    do ivlib=1,nvlib
        vec_residu(ivlib)=( vec_residu(ivlib) - rstfor(ivlib) )
    end do


    ! Calculate the variation of the residual vector
    diff=pre_vec-vec_residu

    norme=0.d0
    maxre=0.d0  !MAXimum REsidual in the residual forces vector


    ! Calculate norm of residual vector
    do ivlib=1,nvlib
        norme = norme + (diff(ivlib)*diff(ivlib))

        if (diff(ivlib) > maxre) then
            maxre = diff(ivlib)
        end if

    end do


    norme=dsqrt(norme)


    if (norme <= epsilon) then !Solution has converged
        ncheck=0
    end if

    if (norme > epsilon) then
        if (iiter == 1) then
            ncheck = 1             !Solution converging
        end if

        if (iiter > 1 .and. norme < pvalue) then
            ncheck = 1              !Solution converging
        end if

        if (iiter > 1 .and. norme >= pvalue) then
            ncheck = 999
        end if

        pvalue=norme
    end if



! Fin Timer
call cpu_time(finish_inc)

     write(6,30) iiter, ncheck, plast_code, norme, maxre, finish_inc-start_inc
30 format(5x11hIteration =,i3,3x,18hconvergence code =,i3,3x,17hplasticity code =,i3, &
            & 3x,28hnorm of residual sum ratio =,e14.6,3x,18hmaximum residual =,e14.6,3x,6hTime =,f6.3,9h seconds.)


end

!*******************************************************************************
!   Subroutine linear( ) - Owen p.247
!   This subroutine determines the stresses from given displacements assuming linear elastic behaviour.
!*******************************************************************************

subroutine linear(lprop,kgasp,kgaus)

    use variables
    implicit none

    integer :: idofn, jdofn, inode, istre, jstre, lprop, kgasp,i ,j, kgaus, istr1
    double precision :: bgash
    double precision :: lambda, mu !Coefficient de Lamé

    double precision, dimension(2,2) :: agash

    include 'subroutines/interface.f90'

    poiss=props(lprop,2)
    young=props(lprop,1)


    ! On calcule du_x/dx, du_y/dy, etc... pour avoir les deformations epsilon
    do idofn=1,ndofn
        do jdofn=1,ndofn
            bgash=0.d0

            do inode=1,nnode
                bgash=bgash+cartd(jdofn,inode)*eldis(idofn,inode) !eldis=déplacements de chaque noeud de l'élement calculés dans front( )
            end do

            agash(idofn,jdofn)=bgash
        end do
    end do


    ! Calculate the strains
    stran(1)=agash(1,1)
    stran(2)=agash(2,2)
    stran(3)=agash(1,2)+agash(2,1)
    stran(4)=0.d0

    !Calcul de Epsilon_theta pour le cas exisymmetric
    do inode=1,nnode
        stran(4)=stran(4)+eldis(1,inode)*shap(inode)/gpcod(1,kgasp)
    end do

    if (ntype==2) then
        stran(4)=0.d0
    end if



    ! And the corresponding stresses

    !Version Owen remplacée par matmul
    do istre=1,nstre
       stres(istre)=0.d0
       do jstre=1,nstre
           stres(istre)=stres(istre)+dmatx(istre,jstre)*stran(jstre)
       end do
    end do


    ! Plane stress case: epsilon_z ≠ 0, sigma_z=0
    if (ntype == 1) then
        !stran(4)=(lambda/(lambda+2*mu))*(stran(1)+stran(2))
        stres(4)=0.d0
    end if


    ! Plane strain case: epsilon_z=0, sigma_z ≠ 0
    if (ntype == 2) then
        !stran(4)=0.d0
        stres(4)=poiss*(stres(1)+stres(2))
    end if

    ! Save strain tensor
    do istr1=1,nstr1
        !strsg(istr1,kgaus)=strsg(istr1,kgaus)+desig(istr1)
        defsg(istr1,kgaus)=defsg(istr1,kgaus)+stran(istr1)
    end do

    return

end




!*********************************************************************‡***********
!   Subroutine loadps Owen p.183
!   This subroutine evaluates the consistent nodal forces for each element due to discret point loads, gravity loading and distributed edge loading/unit length of element.
!***********************************************************************************

subroutine loadps( )

    use variables

    !implicit double precision (a-h,o-z)
    implicit none


    integer :: ieval, ievab, ielem2, i, j, idofn, kount, neass, istep, idime, igaus, lnode, iodeg, inode, jgaus
    integer :: jnode, knode, mgash, ncode, nloca, nodeg, ngash, kgaus, lprop, kgasp, ielem, lodpt


    double precision :: dvolu, gycom, gravy, gxcom, radus=0.d0, pxcom, pycom

    double precision    :: etasp, exisp
    integer, dimension(4) :: noprs

    double precision , dimension(4,2) :: press
    double precision , dimension(2) :: point
    double precision , dimension(2) :: pgash
    double precision , dimension(2) :: dgash

    character (len=30)  :: title


    !rload is set to zero
    do ielem=1,nelem
        do ievab=1,nevab
            rload(ielem,ievab)=0.d0
        end do
    end do



    !R12: INTERNAL PRESSURE N°1
    read(5,901) title
    !  901 format(12a6)
    901 format(a17)

    write(6,903) title
    ! 903 format(1h0,12a6)
    903 format(a17)


    ! Read data controling loading types to be inputted
    ! R13
    ! Nonzero values of these respective items (iplod, igrav, iedge) indicate that point loads, gravity loading or distributed edge loading is to be considered
    read(5,919) iplod,igrav,iedge
    write(6,919) iplod,igrav,iedge
    919  format(3i5)



    !Read nodal point loads
    !iplod=0 dans single-element.txt donc cette boucle n'est pas faite
    !iplod=1 dans Example_two_elem.txt donc cette boucle est faite

    !*** Version Victor
    if(iplod /= 0) then
        lodpt=0
        !do while (lodpt < npoin)
            read(5,*) lodpt,(point(idofn),idofn=1,ndofn)
            write(6,931) lodpt,(point(idofn),idofn=1,ndofn)
        931 format(i5,2f10.3)

    !Associate the nodal point loads
            do ielem=1,nelem
                do inode=1,nnode
                    nloca=iabs(lnods(ielem,inode))
                    if(lodpt == nloca) then
                        exit
                    end if
                end do

                do idofn=1,ndofn
                    ngash=(inode-1)*ndofn+idofn
                    rload(ielem,ngash)=point(idofn)
                end do
            end do
        !end do
    end if
    !*** FIN Version Victor

    !*** Version Owen
    !if(iplod /= 0) then

    !    20  read(5,*) lodpt,(point(idofn),idofn=1,ndofn)
    !        write(6,*) lodpt,(point(idofn),idofn=1,ndofn)
    !    931 format(i5,2f10.3)
    !
    !!Associate the nodal point loads
    !        do 30 ielem=1,nelem
    !            do 30 inode=1,nnode
    !                nloca=iabs(lnods(ielem,inode))
    !             30 if(lodpt == nloca) go to 40
    !             40 do 50 idofn=1,ndofn
    !                    ngash=(inode-1)*ndofn+idofn
    !                 50 rload(ielem,ngash)=point(idofn)
    ! if (lodpt < npoin) go to 20
    !end if
    !*** FIN Version Owen

    !write(6,*)
!     do idofn=1,nelem
!        write(6,*) (rload(idofn,ielem),ielem=1,nevab)
!     end do
!     stop

    !igrav=0 dans single-element.txt donc cette boucle n'est pas faite
    if (igrav /= 0) then

        !Gravity loading section
          read(5,*) theta,gravy
      906 format(2f10.3)
          write(6,911)theta,gravy
      911 format(1h0,16h gravity angle =,f10.3,19h gravity constant =,f10.3)

          theta=theta/radian

        !Loop over each element
        do ielem=1,nelem

            !Set up priminilary values
            lprop=matno(ielem)
            thick=props(lprop,3)
            dense=props(lprop,4)

            if(dense == 0.0) exit

            gxcom=dense*gravy*dsin(theta)
            gycom=-dense*gravy*dcos(theta)

            !Compute coordinates
            do inode=1,nnode
                lnode=iabs(lnods(ielem,inode))
                do idime=1,2
                   elcod(idime,inode)=coord(lnode,idime)
                end do
            end do

            !Enter loops for area numerical integration
            kgasp=0
            kgaus=0
            do igaus=1,ngaus
                do jgaus=1,ngaus

                    exisp=posgp(igaus)
                    etasp=posgp(jgaus)
                    !compute the shape functions
                    call sfr2(etasp,exisp)

                    kgasp=kgasp+1
                    kgaus=kgaus+1

                    call jacob2(ielem,kgasp)
                    dvolu=djacb*weigp(igaus)*weigp(jgaus)
                    if(thick /= 0.0) dvolu=dvolu*thick
                    if(ntype == 3) dvolu=dvolu*twopi*gpcod(1,kgasp)

                    !calcuklate loads and associate with element nodal

                    do inode=1,nnode
                        ngash=(inode-1)*ndofn+1
                        mgash=(inode-1)*ndofn+2
                        rload(ielem,ngash)=rload(ielem,ngash)+gxcom*shap(inode)*dvolu
                        rload(ielem,mgash)=rload(ielem,mgash)+gycom*shap(inode)*dvolu
                    end do
                end do
            end do
        end do
    end if

    ! *** Distributed edge loads section
    if(iedge /= 0) then

        !  nedge : number of loaded edges due to external (distributed loads)
        !  nedgel: number of loaded edges due to liquid forces
        !  nedgea: number of loaded edges due to air forces
        !  nedget: number of loaded edges due to temperature forces

        !R14. nedge=2
        read(5,932) nedge
        ! read(5,*) nedge,nedgel,nedgea,nedget
        932 format(i5)

        write(6,912) nedge
        ! write(6,*) nedge,nedgel,nedgea,nedget
        912 format(1h 5x,21hno.of loaded edges = ,i5)

        write(6,915)
        915 format(1h 5x,38hlist of loaded edges and applied loads)

        nodeg=3
        ncode=nnode

        if(nnode == 4) then
            nodeg=2
        end if

        if(nnode == 9) then
            ncode=8
        end if

        ! Loop over each loaded edge

        do iedge=1,nedge !nedge=3

            ! Read data locating the loaded edge and applied load

            ! neass: The element number with which the element edge is associated
            ! noprs: list of nodal points in an anticlockwise sequence, of the nodes forming the element face on which the distributed load acts.
            read(5,902) neass,(noprs(iodeg),iodeg=1,nodeg) !nodeg=3
            902 format(4i5)

            write(6,913) neass,(noprs(iodeg),iodeg=1,nodeg)
            913 format(i10,5x,3i5)

            ! press(i,1): value of normal component distributed load at node noprs(i)
            ! press(i,2): value of tangential component distributed load at node noprs(i)

            read(5,*) ((press(iodeg,idofn),idofn=1,ndofn),iodeg=1,nodeg)  !3 lectures de 2 blocs chacune
            write(6,914) ((press(iodeg,idofn),idofn=1,ndofn),iodeg=1,nodeg)
            914 format(6f15.1)


            etasp=-1.d0 ! Choix abritraire - cf. explication bas de page 184


            ! Calculate the coordinates of the nodes of the element edge
            ! elcod: stores element coordinates
            ! nodeg=3
            do iodeg=1,nodeg
                lnode=noprs(iodeg) !Numero du noeud soumis au chargement
                do idime=1,ndim
                    elcod(idime,iodeg)=coord(lnode,idime) !Coordonnées du noeud en question
                end do
            end do


            ! Enter loop for linear numerical integration
            do igaus=1,ngaus
                exisp=posgp(igaus) !Sampling point coordinates (ksi_p,eta_p) are called exisp and etasp
                !Evaluate the shape functions at the sampling point
                call sfr2(etasp,exisp)

                !calculate components of the equivalent nodal loads

                !      problem with integration on line (boundary)
                !   -shap(i)*(fh1+fh2+fh3+fh4+fh5+fh6)
                !      -shap(i,1)*roda*ha*-1*(nul+nua)
                !      -shap(i,1)*rol*(nul+vv+nua)


                do idofn=1,2
                    pgash(idofn)=0.d0
                    dgash(idofn)=0.d0
                    do iodeg=1,nodeg
                        pgash(idofn)=pgash(idofn)+press(iodeg,idofn)*shap(iodeg)
                        dgash(idofn)=dgash(idofn)+elcod(idofn,iodeg)*deriv(1,iodeg)
                    end do
                end do


                dvolu=weigp(igaus)
                pxcom=dgash(1)*pgash(2)-dgash(2)*pgash(1)  !Equation 6.63 Owen p.184
                pycom=dgash(1)*pgash(1)+dgash(2)*pgash(2)

                ! ntype =3: probleme axisymmetric
                if(ntype == 3) then
                    radus=0.d0
                    do iodeg=1,nodeg
                        radus=radus+shap(iodeg)*elcod(1,iodeg)
                    end do
                    dvolu=dvolu*twopi*radus
                end if

                ! Associate the equivalent nodal edge load with an element

                do inode=1,nnode
                    nloca=iabs(lnods(neass,inode))
                    if(nloca == noprs(1)) then
                        exit
                    end if
                end do


                jnode=inode+nodeg-1
                kount=0

                do knode=inode,jnode
                    kount=kount+1
                    ngash=(knode-1)*ndofn+1  !ndofn=2
                    mgash=(knode-1)*ndofn+2

                    if(knode > ncode) then  !ncode=8
                        ngash=1
                    end if

                    if(knode > ncode) then
                        mgash=2
                    end if

                    rload(neass,ngash)=rload(neass,ngash)+shap(kount)*pxcom*dvolu
                    rload(neass,mgash)=rload(neass,mgash)+shap(kount)*pycom*dvolu

                end do

            end do
        end do
    end if


    write(6,907)
    907 format(1h 5x,36h total nodal forces for each element)

    do ielem=1,nelem
        write(6,905) ielem,(rload(ielem,ievab),ievab=1,nevab)
    end do

    !905 format(1x,i4,5x,8e12.4/(10x,8e12.4))
    905 format(1x,i4,5x,8f16.4/(10x,8f16.4))
    write(6,*) ''
    write(6,*) ''


    return

end



!*******************************************************************
! subroutine nodexy Owen p.178
! This subroutine interpolates the mide side nodes of straight sides of elements
! and the central node of 9 noded elements
!*********************************************************************

! Checks each midside node (a midside node being recognisable from the element topology card).
! If both coordinates of a midside node are found to be zero, its coordinates are linearly interpolated
! between the two adjacent corner nodes.

subroutine nodexy( )

use variables

implicit none

integer :: i, ielem
integer :: nnod1, kount, igash, inode, nodst, midpt, nodmd, nodfn
integer :: lnode, lnod1, lnod3, lnod5, lnod7, drap=0
double precision    :: total

if(nnode == 4) then
    return
end if

! Loop over each element
do ielem=1,nelem

    !Loop over each element edge
    nnod1=9
    if(nnode == 8) then
        nnod1=7
    end if

    do inode=1,nnod1,2

        if(inode == 9) then
            exit
        end if

        !Compute the node number of the first node
        nodst=lnods(ielem,inode)
        igash=inode+2

        if(igash > 8) then
            igash=1
        end if

        !Compute the node number of the last node
        nodfn=lnods(ielem,igash)
        midpt=inode+1

        !Compute the node number of the intermediate node
        nodmd=lnods(ielem,midpt)
        total=abs(coord(nodmd,1))+abs(coord(nodmd,2))

        !If the coordinates of the intermediate node are both zero: interpolate by a straight line
        if (total <= 0.0) then
            do kount=1,2
                coord(nodmd,kount)=(coord(nodst,kount)+coord(nodfn,kount))/2.0
            end do
            write(6,*) "Element: ", ielem, "Noeud:", inode, ": Correction des coordonnées"
            drap=1
        end if

    end do

    if(nnode == 8) then
        inode=7
    end if

    if (inode == 9) then
        lnode=lnods(ielem,inode)
        total=abs(coord(lnode,1))+abs(coord(lnode,2))

        if(total<=0.0) then
            lnod1=lnods(ielem,1)
            lnod3=lnods(ielem,3)
            lnod5=lnods(ielem,5)
            lnod7=lnods(ielem,7)
            kount=1
            do kount=1,2
                coord(lnode,kount)=(coord(lnod1,kount)+coord(lnod3,kount)+coord(lnod5,kount)+coord(lnod7,kount))/4.0
            end do
            write(6,*) "Element: ", ielem, "Noeud:", inode, ": Correction des coordonnées"
            drap=1
        end if
    end if

end do


if (drap /= 0)then
    print *, "Subroutine NODEXY - Status: ERROR"
    stop
end if


return
end



!   *****************************************************************
!   Subroutine zero - Owen p.238
!   This subroutine sets to zero the contents of several arrays employed in the program.
!   These arrays will be employed to accumulate data as the incremental and iterative process continues.
!   ******************************************************************


subroutine zero( )

    use variables
    implicit none

    integer :: ievab, itotv, ivfix, idofn, itotg, istr1, ielem

    tfact=0.d0

    !write(6,*) "nelem",nelem
    !write(6,*) "nevab",nevab
    !write(6,*) "ntotv",ntotv
    !write(6,*) "nvfix",nvfix
    !write(6,*) "ndofn",ndofn
    !write(6,*) "ntotg",ntotg
    !write(6,*) "nstr1",nstr1


    do ielem=1,nelem
        do ievab=1,nevab
            eload(ielem,ievab)=0.d0
            tload(ielem,ievab)=0.d0
        end do
    end do

    tofor=0.d0

    do itotv=1,ntotv
        tdisp(itotv)=0.d0
    end do

    do ivfix=1,nvfix
        do idofn=1,ndofn
            treac(ivfix,idofn)=0.d0
        end do
    end do

    do itotg=1,ntotg
        epstn(itotg)=0.d0
        effst(itotg)=0.d0
        do istr1=1,nstr1
            strsg(istr1,itotg)=0.d0
        end do
    end do


    !Cam Clay
    if (ncrit==5 .or. ncrit==6) then
        ! Initialize yield stress
        do itotg=1,ntotg
            yieldgp(itotg)=props(1,5)
            volspegp(itotg)=0.d0
            epstnp(itotg)=0.d0
            epstnp_pre(itotg)=0.d0
        end do
    end if

    return
end


!****************************************************************************
!   Subroutine sfr2 Owen p.179
!   This subroutine evaluates the shape functions Ni and their derivatives dNi/dksi, dNi,deta at any sampling point (ksi_p,eta_p) within the element for each of the 4-, 8- or 9-noded elements.
!****************************************************************************

!subroutine sfr2(deriv,etasp,exisp,nnode,shap,shap1,derivt)
subroutine sfr2(etasp,exisp)

use variables
implicit none

double precision :: s, t, st, ss, tt, sst, stt, s2, t2, st2, s1, s9, t1, t9
integer :: i,j
double precision :: etasp, exisp

include 'subroutines/interface.f90'
s=exisp  ! ksi_p
t=etasp  ! eta_p

! The evaluated shape functions for each node of an element are stored in the array shape(inode) and their derivatives in deriv(inode, idime).

! Cf. Dhatt & Touzot p. 114
if(nnode == 3) then

    !Shape functions for 3 noded element
    shap(1)=1.d0-s-t
    shap(2)=s
    shap(3)=t


    !Shape functions derivatives
    deriv(1,1)=-1.d0
    deriv(1,2)=1.d0
    deriv(1,3)=0.d0

    deriv(2,1)=-1.d0
    deriv(2,2)=0.d0
    deriv(2,3)=1.d0
end if


if(nnode == 4) then
    st=s*t

    !Shape functions for 4 noded element

    shap(1)=(1.d0-t-s+st)*0.25d0
    shap(2)=(1.d0-t+s-st)*0.25d0
    shap(3)=(1.d0+t+s+st)*0.25d0
    shap(4)=(1.d0+t-s-st)*0.25d0

    !Shape functions derivatives
    deriv(1,1)=(-1.d0+t)*0.25d0
    deriv(1,2)=(+1.d0-t)*0.25d0
    deriv(1,3)=(+1.d0+t)*0.25d0
    deriv(1,4)=(-1.d0-t)*0.25d0
    deriv(2,1)=(-1.d0+s)*0.25d0
    deriv(2,2)=(-1.d0-s)*0.25d0
    deriv(2,3)=(+1.d0+s)*0.25d0
    deriv(2,4)=(+1.d0-s)*0.25d0

end if

! print *,'s=',s,"t=",t
! do i=1,2
!     print *,(deriv(i,j),j=1,4)
! end do

if(nnode == 8) then
    s2=s*2.d0
    t2=t*2.d0
    ss=s*s
    tt=t*t
    st=s*t
    sst=s*s*t
    stt=s*t*t
    st2=s*t*2.d0

    !Shape functions for 8 noded element

!                                     s=+1
!                                     ^
!                                     |
!                                     |
!           7------------6------------5=(1,1)
!           |                         |
!           |                         |         t=eta
!           |                         |         ^
!           8                         4         |
!           |                         |         |
!           |                         |         0-----> s=ksi
!           |                         |
!           1------------2------------3=(1,-1)   ------> t=-1
!


    shap(1)=(-1.d0+st+ss+tt-sst-stt)/4.d0     !=1 en (s,t)=(-1,-1)
    shap(2)=(1.d0-t-ss+sst)/2.d0              !=1 en (s,t)=(0,-1)
    shap(3)=(-1.d0-st+ss+tt-sst+stt)/4.d0     !=1 en (s,t)=(1,-1)
    shap(4)=(1.d0+s-tt-stt)/2.d0
    shap(5)=(-1.d0+st+ss+tt+sst+stt)/4.d0
    shap(6)=(1.d0+t-ss-sst)/2.d0
    shap(7)=(-1.d0-st+ss+tt+sst-stt)/4.d0
    shap(8)=(1.d0-s-tt+stt)/2.d0

    !Shape function derivatives

    deriv(1,1)=(t+s2-st2-tt)/4.d0
    deriv(1,2)=-s+st
    deriv(1,3)=(-t+s2-st2+tt)/4.d0
    deriv(1,4)=(1.d0-tt)/2.d0
    deriv(1,5)=(t+s2+st2+tt)/4.d0
    deriv(1,6)=-s-st
    deriv(1,7)=(-t+s2+st2-tt)/4.d0
    deriv(1,8)=(-1.d0+tt)/2.d0
    deriv(2,1)=(s+t2-ss-st2)/4.d0
    deriv(2,2)=(-1.d0+ss)/2.d0
    deriv(2,3)=(-s+t2-ss+st2)/4.d0
    deriv(2,4)=-t-st
    deriv(2,5)=(s+t2+ss+st2)/4.d0
    deriv(2,6)=(1.d0-ss)/2.d0
    deriv(2,7)=(-s+t2+ss-st2)/4.d0
    deriv(2,8)=-t+st

    ! *** Elkassas version
    !do i=1,nnode
    !    shap1(1,i)=shap(i)
    !end do

    !Change 4 Elkassas
    !do i=1,2
    !    do j=1,nnode
    !        derivt(j,i)=deriv(i,j)
    !    end do
    !end do
    !End of change 4

end if

if(nnode == 9) then
    ss=s*s
    st=s*t
    tt=t*t
    s1=s+1.d0
    t1=t+1.d0
    s2=s*2.d0
    t2=t*2.d0
    s9=s-1.d0
    t9=t-1.d0

    !Shape functions for 9 noded elements

    shap(1)=0.25d0*s9*st*t9
    shap(2)=0.5d0*(1.d0-ss)*t*t9
    shap(3)=0.25d0*s1*st*t9
    shap(4)=0.5d0*s*s1*(1.d0-tt)
    shap(5)=0.25d0*s1*st*t1
    shap(6)=0.5d0*(1.d0-ss)*t*t1
    shap(7)=0.25d0*s9*st*t1
    shap(8)=0.5d0*s*s9*(1.d0-tt)
    shap(9)=(1.d0-ss)*(1.d0-tt)

    !Shape functions

    deriv(1,1)=0.25d0*t*t9*(-1.d0+s2)
    deriv(1,2)=-st*t9
    deriv(1,3)=0.25d0*(1.d0+s2)*t*t9
    deriv(1,4)=0.5d0*(1.d0+s2)*(1.d0-tt)
    deriv(1,5)=0.25d0*(1.d0+s2)*t*t1
    deriv(1,6)=-st*t1
    deriv(1,7)=0.25d0*(-1.d0+s2)*t*t1
    deriv(1,8)=0.5d0*(-1.d0+s2)*(1.d0-tt)
    deriv(1,9)=-s2*(1.0-tt)
    deriv(2,1)=0.25d0*(-1.d0+t2)*s*s9
    deriv(2,2)=0.5d0*(1.d0-ss)*(-1.d0+t2)
    deriv(2,3)=0.25d0*s*s1*(-1.d0+t2)
    deriv(2,4)=-st*s1
    deriv(2,5)=0.25d0*s*s1*(1.d0+t2)
    deriv(2,6)=0.5d0*(1d0-ss)*(1.d0+t2)
    deriv(2,7)=0.25d0*s*s9*(1.d0+t2)
    deriv(2,8)=-st*s9
    deriv(2,9)=-t2*(1.d0-ss)

end if

end



!****************************************************************************
!   Subroutine sfr2_extrap Owen p.179
!   This subroutine evaluates the shape functions Ni and their derivatives dNi/dksi, dNi,deta at any sampling point (ksi_p,eta_p) within the element for each of the 4-, 8- or 9-noded elements.
!****************************************************************************

!subroutine sfr2(deriv,etasp,exisp,nnode,shap,shap1,derivt)
subroutine sfr2_extrap(ksi,eta)

use variables
implicit none

double precision :: s, t, st, ss, tt, sst, stt, s2, t2, st2, s1, s9, t1, t9
integer :: i,j
double precision :: ksi, eta

include 'subroutines/interface.f90'
s=ksi  ! ksi_p
t=eta  ! eta_p

s2=s*2.d0
t2=t*2.d0
ss=s*s
tt=t*t
st=s*t
sst=s*s*t
stt=s*t*t
st2=s*t*2.d0

! print *,"s=",s,"t=",t
! The evaluated shape functions for each node of an element are stored in the array shape(inode) and their derivatives in deriv(inode, idime).

! Cf. Dhatt & Touzot p. 114
if(nnode == 3) then

    !Shape functions for 3 noded element
    shap(1)=1.d0-s-t
    shap(2)=s
    shap(3)=t


    !Shape functions derivatives
    deriv(1,1)=-1.d0
    deriv(1,2)=1.d0
    deriv(1,3)=0.d0

    deriv(2,1)=-1.d0
    deriv(2,2)=0.d0
    deriv(2,3)=1.d0
end if


if(ngaus==2) then
    !Shape functions for 4 noded element

    shap(1)=(1.d0-t-s+st)*0.25d0 ! N1 pour G1
    shap(3)=(1.d0-t+s-st)*0.25d0 ! N2 pour G3
    shap(4)=(1.d0+t+s+st)*0.25d0 ! N3 pour G4
    shap(2)=(1.d0+t-s-st)*0.25d0 ! N4 pour G2

end if

! print *,'s=',s,"t=",t
! do i=1,2
!     print *,(deriv(i,j),j=1,4)
! end do


! If integration only on 8 gp, uncomment this section
! if(ngaus == 3 .and. nnode==8) then
!     s2=s*2.d0
!     t2=t*2.d0
!     ss=s*s
!     tt=t*t
!     st=s*t
!     sst=s*s*t
!     stt=s*t*t
!     st2=s*t*2.d0

!     !Shape functions for 8 noded element

! !                                     s=+1
! !                                     ^
! !                                     |
! !                                     |
! !           7------------6------------5=(1,1)
! !           |                         |
! !           |                         |         t=eta
! !           |                         |         ^
! !           8                         4         |
! !           |                         |         |
! !           |                         |         0-----> s=ksi
! !           |                         |
! !           1------------2------------3=(1,-1)   ------> t=-1
! !


!     shap(1)=(-1.d0+st+ss+tt-sst-stt)/4.d0     !=1 en (s,t)=(-1,-1)
!     shap(4)=(1.d0-t-ss+sst)/2.d0              !=1 en (s,t)=(0,-1)
!     shap(6)=(-1.d0-st+ss+tt-sst+stt)/4.d0     !=1 en (s,t)=(1,-1)
!     shap(7)=(1.d0+s-tt-stt)/2.d0
!     shap(8)=(-1.d0+st+ss+tt+sst+stt)/4.d0
!     shap(5)=(1.d0+t-ss-sst)/2.d0
!     shap(3)=(-1.d0-st+ss+tt+sst-stt)/4.d0
!     shap(2)=(1.d0-s-tt+stt)/2.d0


! end if


if(ngaus == 3) then
    ss=s*s
    st=s*t
    tt=t*t
    s1=s+1.d0
    t1=t+1.d0
    s2=s*2.d0
    t2=t*2.d0
    s9=s-1.d0
    t9=t-1.d0

    !Shape functions for 9 noded elements

    shap(1)=0.25d0*s9*st*t9         !1
    shap(4)=0.5d0*(1.d0-ss)*t*t9    !4
    shap(7)=0.25d0*s1*st*t9         !7
    shap(8)=0.5d0*s*s1*(1.d0-tt)    !8
    shap(9)=0.25d0*s1*st*t1         !9
    shap(6)=0.5d0*(1.d0-ss)*t*t1    !6
    shap(3)=0.25d0*s9*st*t1         !3
    shap(2)=0.5d0*s*s9*(1.d0-tt)    !2
    shap(5)=(1.d0-ss)*(1.d0-tt)     !5

    !Shape functions

end if

end

















