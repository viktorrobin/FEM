subroutine Read_Bound_Cond_Gid()

	use variables
	implicit none

	integer :: ipoin, nb_tab_dyn, ivfix, ifpre, iprop
	integer :: idofn, ifdof, imats, istep, ngash, numat, numel
	integer :: icompt, itotv, itemp

    integer :: ieval, ievab, ielem2, i, j, kount, neass, idime, igaus, lnode, iodeg, inode, jgaus
    integer :: jnode, knode, mgash, ncode, nloca, nodeg, kgaus, lprop, kgasp, ielem, lodpt, iforce


    double precision :: dvolu, gycom, gravy, gxcom, radus=0.d0, pxcom, pycom

    double precision    :: etasp, exisp 
    integer, dimension(4) :: noprs

    double precision , dimension(4,2) :: press
    double precision , dimension(2) :: point
    double precision , dimension(2) :: pgash
    double precision , dimension(2) :: dgash

    character (len=30)  :: title
    character (len=60) :: atemp

	include 'subroutines/interface.f90'

    ! If you want to print console output in a txt file
    !   open(unit=6,file='output_terminal.txt')         !Writing of the variables values

    
    ! Start reading prescribed displacement for this stage


    write(6,917) !7
    ! original: write(6,*) !7
    917 format(//5h node,6x,4hcode,6x,12h fixed values)

    ! Reset to zero
    do itotv=1,ntotv
        iffix(itotv)=0
    end do

    do i=1,8
    	read(5,*)
    end do

    ! Read number of variables blocked (Actually, this already obtained before)
    read(5,*) nvfix

    do i=1,1
    	read(5,*)
    end do

    908 format(1x,i4,5x,i5,5x,f10.4,1x,f10.4)

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



    do i=1,2
    	read(5,*)
    end do

    !Deal with internal pressure (gravity, point load, and distributed loading)

    !rload is set to zero
    do ielem=1,nelem
        do ievab=1,nevab
            rload(ielem,ievab)=0.d0
        end do
    end do

    ! Read data controling loading types to be inputted
    read(5,*) atemp,iplod
    read(5,*) atemp,iedge
    read(5,*) atemp,igrav

  
    !Read nodal point loads
    !*** Version Victor
    if(iplod /= 0) then
    	do i=1,4
    		read(5,*)
    	end do
    	do iforce=1,iplod
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
        end do
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



    ! *** Distributed edge loads section
    if(iedge /= 0) then
    	do i=1,3
    		read(5,*)
    	end do

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

            ! R15,R17
            ! neass: The element number with which the element edge is associated
            ! noprs: list of nodal points in an anticlockwise sequence, of the nodes forming the element face on which the distributed load acts. 
            if (nnode==4) then
            	read(5,*) neass,(noprs(iodeg),iodeg=1,nodeg) !nodeg=3
            end if

            if (nnode==8) then
            	read(5,*) neass,noprs(1),noprs(3),noprs(2)
            end if

            if (nnode==9) then
            	read(5,*) neass,noprs(1),noprs(3),noprs(2)
            end if
            
!             902 format(4i5)

            write(6,913) neass,(noprs(iodeg),iodeg=1,nodeg)
            913 format(i10,5x,3i5)
            
            ! R16,R18
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


    !         write(6,*)posgp
    !         stop

            ! Enter loop for linear numerical integration
            do igaus=1,ngaus
                exisp=posgp(igaus) !Sampling point coordinates (ksi_p,eta_p) are called exisp and etasp
                !Evaluate the shape functions at the sampling point
                !call sfr2(deriv,etasp,exisp,nnode,shap,shap1,derivt)
                call sfr2(etasp,exisp)
            
                !calculate components of the equivalent nodal loads

                !      problem with integration on line (boundry)
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
                !nnode=8

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
        

                !call Val_Variables(pgash,dgash,press,radus,pxcom,pycom)
                end do

            end do
        end do
    end if

  
       !igrav=0 dans single-element.txt donc cette boucle n'est pas faite
    if (igrav /= 0) then
    	do i=1,3
    		read(5,*)
    	end do
        !Gravity loading section
          read(5,*) theta,gravy
      906 format(2f10.3)
          write(6,911)theta,gravy
      911 format(1h0,16h gravity angle =,f10.3,19h gravity constant =,f10.3)
          
          

          !Theta is in degrees, however fortran functions dsin, dcos, etc... takes radians. We need to convert:
          theta=theta*degres
          

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
                    !call sfr2(deriv,etasp,exisp,nnode,shap,shap1,derivt)
                    call sfr2(etasp,exisp)
                    
                    kgasp=kgasp+1
                    kgaus=kgaus+1

                    !call jacob2(xjacm,cartd,deriv,djacb,elcod,gpcod,ielem,kgasp,nnode,shap,cardc,cartdt,kgaus,gpcodg)
                    
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

    write(6,907)
    907 format(1h 5x,36h total nodal forces for each element)

    do ielem=1,nelem
        write(6,905) ielem,(rload(ielem,ievab),ievab=1,nevab)
    end do

    !905 format(1x,i4,5x,8e12.4/(10x,8e12.4))
    905 format(1x,i4,5x,8f16.4/(10x,8f16.4))
    write(6,*) ''
    write(6,*) ''

    !do i=1,nelem
    !    do j=1,nevab
    !    write(6,*) rload(i,j)
    !    end do
    !end do

end

subroutine Read_Bound_Cond_Gid_S2()

    use variables
    implicit none

    integer :: ipoin, nb_tab_dyn, ivfix, ifpre, iprop
    integer :: idofn, ifdof, imats, istep, ngash, numat, numel
    integer :: icompt, itotv, itemp

    integer :: ieval, ievab, ielem2, i, j, kount, neass, idime, igaus, lnode, iodeg, inode, jgaus
    integer :: jnode, knode, mgash, ncode, nloca, nodeg, kgaus, lprop, kgasp, ielem, lodpt, iforce


    double precision :: dvolu, gycom, gravy, gxcom, radus=0.d0, pxcom, pycom

    double precision    :: etasp, exisp 
    integer, dimension(4) :: noprs

    double precision , dimension(4,2) :: press
    double precision , dimension(2) :: point
    double precision , dimension(2) :: pgash
    double precision , dimension(2) :: dgash

    character (len=30)  :: title
    character (len=60) :: atemp

    integer :: nvfix_pre, nvlib_pre

    include 'subroutines/interface.f90'

    ! If you want to print console output in a txt file
    !   open(unit=6,file='output_terminal.txt')         !Writing of the variables values

    !Save the number of fixed variables during stage 1
    nvlib_pre=nvlib
    nvfix_pre=nvfix

!     open(unit=1,form='unformatted',file='temp/estif.tmp')

    print *, ''
    print *, ''
    print *, '-----------------------------------'
    print *, '           STAGE 2'
    print *, '-----------------------------------'

    do i=1,5
        read(5,*)
    end do

    read(5,*) atemp,nincs

    ! Start reading prescribed displacement for this stage
    !Deallocate and Reallocate the vectors and matrix that might be affected by the new conditions and nvfix
    tfact=0.d0


    write(6,917) !7
    ! original: write(6,*) !7
    917 format(//5h node,6x,4hcode,6x,12h fixed values)

    do i=1,3
        read(5,*)
    end do

    ! Read number of variables blocked for the second stage.
    read(5,*) nvfix

    !Deallocate nofix and reallocate it now we have nvfix
    ! In case the number of fixed variables has changed, deallocate and reallocate the concerned vectors.
    if (nvfix /= nvfix_pre) then
        deallocate(nofix,presc0)
        allocate(nofix(nvfix))
        allocate(presc0(nvfix,ndofn))
    end if

     ! Reset to zero
    do itotv=1,ntotv
        iffix(itotv)=0
    end do

    do i=1,1
        read(5,*)
    end do

    908 format(1x,i4,5x,i5,5x,f10.4,1x,f10.4)

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

    ! Deallocage all the arrays affected by the new value of nvfix
    if (nvlib /= nvlib_pre) then
        
        !First deallocate
        deallocate(gcoord, rgstif, rgload, rasdis, vec_residu, vec_residu_pd, rstfor, riffix, rfixed)
        deallocate(ipiv, rgstif_lapack,rasdis_lapack, vec_residu_lapack)
        
        !Then reallocate
        allocate(gcoord(nvlib))    !Tableau 28

!         allocate( rgstif(nvlib,nvlib))
!         allocate( rgload(nvlib))
!         allocate( rasdis(nvlib))
!         allocate( vec_residu(nvlib))
!         allocate( vec_residu_pd(nvlib))
!         allocate(rstfor(nvlib))
!         allocate( riffix(nvlib))
!         allocate( rfixed(nvlib))
        
!         !For lapack
!         allocate( ipiv(nvlib))
!         allocate( rgstif_lapack(nvlib,nvlib))
!         allocate( rasdis_lapack(nvlib))
!          allocate( vec_residu_lapack(nvlib))
     end if

     
    do i=1,2
        read(5,*)
    end do

    !Deal with internal pressure (gravity, point load, and distributed loading)


    !rload is set to zero
    do ielem=1,nelem
        do ievab=1,nevab
            rload(ielem,ievab)=0.d0
        end do
    end do

    ! Read data controling loading types to be inputted
    read(5,*) atemp,iplod
    read(5,*) atemp,iedge
    read(5,*) atemp,igrav

  
    !Read nodal point loads
    !*** Version Victor
    if(iplod /= 0) then
        do i=1,4
            read(5,*)
        end do
        do iforce=1,iplod
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
        end do
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



    ! *** Distributed edge loads section
    if(iedge /= 0) then
        do i=1,3
            read(5,*)
        end do

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

            ! R15,R17
            ! neass: The element number with which the element edge is associated
            ! noprs: list of nodal points in an anticlockwise sequence, of the nodes forming the element face on which the distributed load acts. 
            if (nnode==4) then
                read(5,*) neass,(noprs(iodeg),iodeg=1,nodeg) !nodeg=3
            end if

            if (nnode==8) then
                read(5,*) neass,noprs(1),noprs(3),noprs(2)
            end if

            if (nnode==9) then
                read(5,*) neass,noprs(1),noprs(3),noprs(2)
            end if
            
!             902 format(4i5)

            write(6,913) neass,(noprs(iodeg),iodeg=1,nodeg)
            913 format(i10,5x,3i5)
            
            ! R16,R18
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


    !         write(6,*)posgp
    !         stop

            ! Enter loop for linear numerical integration
            do igaus=1,ngaus
                exisp=posgp(igaus) !Sampling point coordinates (ksi_p,eta_p) are called exisp and etasp
                !Evaluate the shape functions at the sampling point
                !call sfr2(deriv,etasp,exisp,nnode,shap,shap1,derivt)
                call sfr2(etasp,exisp)
            
                !calculate components of the equivalent nodal loads

                !      problem with integration on line (boundry)
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
                !nnode=8

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
        

                !call Val_Variables(pgash,dgash,press,radus,pxcom,pycom)
                end do

            end do
        end do
    end if

  
       !igrav=0 dans single-element.txt donc cette boucle n'est pas faite
    if (igrav /= 0) then
        do i=1,3
            read(5,*)
        end do
        !Gravity loading section
          read(5,*) theta,gravy
      906 format(2f10.3)
          write(6,911)theta,gravy
      911 format(1h0,16h gravity angle =,f10.3,19h gravity constant =,f10.3)
          
          

          !Theta is in degrees, however fortran functions dsin, dcos, etc... takes radians. We need to convert:
          theta=theta*degres
          

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
                    !call sfr2(deriv,etasp,exisp,nnode,shap,shap1,derivt)
                    call sfr2(etasp,exisp)
                    
                    kgasp=kgasp+1
                    kgaus=kgaus+1

                    !call jacob2(xjacm,cartd,deriv,djacb,elcod,gpcod,ielem,kgasp,nnode,shap,cardc,cartdt,kgaus,gpcodg)
                    
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

    write(6,907)
    907 format(1h 5x,36h total nodal forces for each element)

    do ielem=1,nelem
        write(6,905) ielem,(rload(ielem,ievab),ievab=1,nevab)
    end do

    !905 format(1x,i4,5x,8e12.4/(10x,8e12.4))
    905 format(1x,i4,5x,8f16.4/(10x,8f16.4))
    write(6,*) ''
    write(6,*) ''

    !do i=1,nelem
    !    do j=1,nevab
    !    write(6,*) rload(i,j)
    !    end do
    !end do

! if (istage==2) stop


end