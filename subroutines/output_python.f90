!*******************************************************************************
!   Subroutine output_python( )
!   This subroutine outputs displacements, reactions, and stresses
!*******************************************************************************

subroutine output_python( )

use variables
implicit none

integer :: koutp, ipoin, ngash, ngish, kelgs, istre, istr1, ivfix, idofn, igash, igaus, jgaus, kgaus
integer :: ielem, inode,jdim, idime, lnode, knode, kshap, kgaus_rep
double precision :: xgash, xgish, xgesh, xgosh
double precision :: etasp, exisp,alpha_extrap, ksi, eta


double precision, dimension(2) :: centroid_cord
double precision, dimension(nstr1) :: centroid_def

double precision, dimension(3) :: strsp
integer, dimension(nnode) :: ksi_coef, eta_coef
double precision, dimension(nstr1,nnode) :: def_node

double precision, dimension(npoin) :: Vec_Node_Hist
double precision, dimension(npoin,nstr1) :: def_node2
double precision, dimension(npoin) :: plastic_node

!---------------------------------------------------------
! Basic parameters
!---------------------------------------------------------
open(unit=110,file='PostProcessing/Output/parameters.txt')
write(110,*) nelem,nnode,npoin,count_file,ntype,ncrit

!---------------------------------------------------------
!Write connectivity table
!---------------------------------------------------------
! open(unit=111,form='unformatted',file='output/coordinates.txt')
 open(unit=113,file='PostProcessing/Output/conntable.txt')

do ielem=1,nelem
    write(113,*) ielem,(lnods(ielem,inode),inode=1,nnode)
end do

!---------------------------------------------------------
!Write coordinates
!---------------------------------------------------------
! open(unit=111,form='unformatted',file='output/coordinates.txt')
 open(unit=111,file='PostProcessing/Output/coordinates.txt')
write(111,222)
222 format(10x,1hx,18x,1hy)

do inode=1,npoin
    write(111,*) (coord(inode,jdim),jdim=1,ndim)
end do

!---------------------------------------------------------
!Write DISPLACEMENTS
!---------------------------------------------------------
! open(unit=111,form='unformatted',file='output/coordinates.txt')
 open(unit=114,file='PostProcessing/Output/displacements.txt')

do ipoin=1,npoin
    ngish=(ipoin-1)*2+1
    write(114,*) tdisp(ngish),tdisp(ngish+1)
end do

open(unit=115,file='PostProcessing/Output/displacements2.txt')
do ipoin=1,npoin
    ngish=(ipoin-1)*2+1
    write(115,*) (coord(ipoin,jdim),jdim=1,ndim),tdisp(ngish),tdisp(ngish+1)
end do

!---------------------------------------------------------
! write STRAINS
!---------------------------------------------------------
open(unit=116,file='PostProcessing/Output/strains.txt')

if(ntype /= 3) write(116,971)
971 format(7hElement,3x,4hg.p.,3x,3hNo.,9x,1hx,15x,1hy,12x,6hEPS_XX,10x,6hEPS_YY,10x,6hEPS_XY,10x,6hEPS_ZZ,10x,8hmax p.s.,8x,&
    &8hmin p.s.,7x,5hangle)!,3x,6he.p.s.)
if(ntype == 3) write(116,976)
976 format(7hElement,6x,4hg.p.,6x,1hx,6x,1hy,6x,6hEPS_RR,5x,6hEPS_ZZ,5x,6hEPS_RZ,5x,6hEPS_TT,6x,8hmax p.s.,6x, &
	& 8hmin p,s.,3x,5hangle)!,3x),6he.p.s.)

941   format(i7,2x,i3,2x,i4,2x,8e16.6,f10.3)

! Calculate coef to change base of the coordinates system
if (ngaus==2) then
	alpha_extrap=dsqrt(3.d0)
else if (ngaus==3) then
	alpha_extrap=dsqrt(5.d0/3.d0)
end if
! print *,"alpha=",alpha_extrap
! Extrolation of the nodal values from the gauss points

! Coef for new coordinates
if (nnode==4) then
	! Numerotation dans sens antihoraire!!!! Sinon probleme avec les coordonnees
	ksi_coef(1) = -1
	ksi_coef(2) = +1
	ksi_coef(3) = 1
	ksi_coef(4) = -1

	eta_coef(1) = -1
	eta_coef(2) = -1
	eta_coef(3) = 1
	eta_coef(4) = 1
else if (nnode==8) then
	ksi_coef(1) = -1
	ksi_coef(2) = 0
	ksi_coef(3) = 1
	ksi_coef(4) = 1
	ksi_coef(5) = 1
	ksi_coef(6) = 0
	ksi_coef(7) = -1
	ksi_coef(8) = -1

	eta_coef(1) = -1
	eta_coef(2) = -1
	eta_coef(3) = -1
	eta_coef(4) = 0
	eta_coef(5) = 1
	eta_coef(6) = 1
	eta_coef(7) = 1
	eta_coef(8) = 0
else if (nnode==9) then
	ksi_coef(1) = -1
	ksi_coef(2) = 0
	ksi_coef(3) = 1
	ksi_coef(4) = 1
	ksi_coef(5) = 1
	ksi_coef(6) = 0
	ksi_coef(7) = -1
	ksi_coef(8) = -1
	ksi_coef(9) = 0

	eta_coef(1) = -1
	eta_coef(2) = -1
	eta_coef(3) = -1
	eta_coef(4) = 0
	eta_coef(5) = 1
	eta_coef(6) = 1
	eta_coef(7) = 1
	eta_coef(8) = 0
	ksi_coef(9) = 0
end if

! defsg(1,1)=15.0
! defsg(1,2)=10.0
! defsg(1,3)=20.0
! defsg(1,4)=15.0


! defsg(1,1)=10.0
! defsg(1,2)=15.0
! defsg(1,3)=15.0
! defsg(1,4)=20.0


! defsg(1,1)=10.0
! defsg(1,2)=15.0
! defsg(1,3)=20.0
! defsg(1,4)=15.0
! defsg(1,5)=20.0
! defsg(1,6)=25.0
! defsg(1,7)=20.0
! defsg(1,8)=25.0
! defsg(1,9)=30.0

! defsg(1,10)=25.0
! defsg(1,11)=30.0
! defsg(1,12)=35.0
! defsg(1,13)=30.0
! defsg(1,14)=35.0
! defsg(1,15)=40.0
! defsg(1,16)=35.0
! defsg(1,17)=40.0
! defsg(1,18)=45.0

! defsg(1,19)=40.0
! defsg(1,20)=45.0
! defsg(1,21)=50.0
! defsg(1,22)=45.0
! defsg(1,23)=50.0
! defsg(1,24)=55.0
! defsg(1,25)=50.0
! defsg(1,26)=55.0
! defsg(1,27)=60.0

! defsg(1,28)=25.0
! defsg(1,29)=30.0
! defsg(1,30)=35.0
! defsg(1,31)=30.0
! defsg(1,32)=35.0
! defsg(1,33)=40.0
! defsg(1,34)=35.0
! defsg(1,35)=40.0
! defsg(1,36)=45.0


!Extrapolate values from gauss points to nodes
Vec_Node_Hist=0.d0
def_node2=0.d0
do ielem=1,nelem
	! Set to zero matrix where extrapolated values at nodes for each element will be stored
	def_node=0.d0
	kgaus_rep=(ielem-1)*ngaus*ngaus

	do inode=1,nnode

		! Calculate the coordinates of the element nodal points
! 		elcod=0.0d0
		lnode=iabs(lnods(ielem,inode))
! 		print *,'lnode=',lnode

		Vec_Node_Hist(lnode)=Vec_Node_Hist(lnode)+1.d0
! 		do idime=1,ndim
! 			elcod(idime,inode)=coord(lnode,idime)
! 		end do

		! Calculate coordinates of the node for the extrapolation
		ksi=alpha_extrap*ksi_coef(inode)
		eta=alpha_extrap*eta_coef(inode)
! 		print *,"ksi=",ksi,"eta=",eta

		! Calculate shape functions at the node
		call sfr2_extrap(ksi,eta)

! 		print *, (shap(i),i=1,nnode)

		kshap=0
		kgaus=kgaus_rep

		! Extrapolation loop over all the gauss points
		do igaus=1,ngaus
			do jgaus=1,ngaus
				! Increment gauss point number
				kgaus=kgaus+1

				! We only use middle gauss point if 9 noded element
! 				if (nnode==8 .and. igaus==2 .and. jgaus==2) go to 456

				! Increment shape function number
				kshap=kshap+1

				! Extrapolate all the values eps_xx, eps_yy, eps_xy, eps_zz
				do istr1=1,nstr1
					def_node2(lnode,istr1) = def_node2(lnode,istr1) + defsg(istr1,kgaus)*shap(kshap)
				end do
! 				def_node(1,inode) = def_node(1,inode) + defsg(1,kgaus)*shap(kshap)

! 				print *,'value=',defsg(1,kgaus),"shap=",shap(kshap),"prod=",defsg(1,kgaus)*shap(kshap)
				456 continue
	    	end do
		end do
	end do


	if (ngaus==2) then
		!******************************************
		! Calculate value at element centroid
		!******************************************

		centroid_cord=0.d0
		centroid_def=0.d0

		! Calculate coordinates of the centroid for the extrapolation
		ksi=0.d0
		eta=0.d0
! 		print *,"ksi=",ksi,"eta=",eta

		! Calculate coordinates of nodes of the element
		do inode=1,nnode
			lnode=iabs(lnods(ielem,inode))

			do idime=1,ndim
				elcod(idime,inode)=coord(lnode,idime)
			end do
		end do

		!**Calculate coordinates centroid in global coordinates

		!Shape functions
		call sfr2(ksi,eta)

		!Coordinates
		do idime=1,ndim
			centroid_cord(idime)=0.d0
			do inode=1,nnode
				centroid_cord(idime)=centroid_cord(idime)+elcod(idime,inode)*shap(inode)
			end do
		end do

		!Calculate values at centroid

		! Calculate shape functions at the node
		call sfr2_extrap(ksi,eta)

! 		print *, shap(1),shap(2),shap(3),shap(4)

		kshap=0
		kgaus=kgaus_rep

		! Extrapolation loop over all the gauss points

		do igaus=1,ngaus
			do jgaus=1,ngaus
				! Increment gauss point number
				kgaus=kgaus+1

				kshap=kshap+1

				! Extrapolate all the values eps_xx, eps_yy, eps_xy, eps_zz
				do istr1=1,nstr1
					centroid_def(istr1) = centroid_def(istr1) + defsg(istr1,kgaus)*shap(kshap)
				end do
			end do
		end do

		xgash=(centroid_def(1)+centroid_def(2))*0.5d0
		xgish=(centroid_def(1)-centroid_def(2))*0.5d0
		xgesh=centroid_def(3)
		xgosh=dsqrt(xgish*xgish+xgesh*xgesh)
		strsp(1)=xgash+xgosh
		strsp(2)=xgash-xgosh
		if(xgish == 0.d0) xgish=0.1d-20
		strsp(3)=datan(xgesh/xgish)*28.647889757d0  !=90/pi

		!Write values at centroid
		write(116,941) ielem,2,0,(centroid_cord(idime),idime=1,2),(centroid_def(istr1),istr1=1,nstr1),(strsp(istre),istre=1,3)
	end if
end do

! do ipoin=1,npoin
! 	print *,ipoin,Vec_Node_Hist(ipoin)
! end do


! Average values and write into text file
do ipoin=1,npoin
	do istr1=1,nstr1
		def_node2(ipoin,istr1) = def_node2(ipoin,istr1)/Vec_Node_Hist(ipoin)
	end do

	xgash=(def_node2(ipoin,1)+def_node2(ipoin,2))*0.5d0
	xgish=(def_node2(ipoin,1)-def_node2(ipoin,2))*0.5d0
	xgesh=def_node2(ipoin,3)
	xgosh=dsqrt(xgish*xgish+xgesh*xgesh)
	strsp(1)=xgash+xgosh
	strsp(2)=xgash-xgosh
	if(xgish == 0.d0) xgish=0.1d-20
	strsp(3)=datan(xgesh/xgish)*28.647889757d0  !=90/pi

	write(116,941) 0,0,ipoin,(coord(ipoin,idime),idime=1,2),(def_node2(ipoin,istr1),istr1=1,nstr1),(strsp(istre), &
	    	& istre=1,3)

end do



! Write coordinates of gp and the values
kgaus=0
do 61 ielem=1,nelem
	kelgs=0
! 	write (6,931) ielem
! 931 format(6x,13helement no. =,i5)
	! Evaluate the coordinates of the element nodal points
	do inode=1,nnode
		lnode=iabs(lnods(ielem,inode))

		do idime=1,ndim
			elcod(idime,inode)=coord(lnode,idime)
		end do
	end do

	do 61 igaus=1,ngaus
		exisp=posgp(igaus)
		do 61 jgaus=1,ngaus
			etasp=posgp(jgaus)
			kgaus=kgaus+1
			kelgs=kelgs+1

			!**Calculate coordinates Gauss Points in global coordinates

			!Shape functions
			call sfr2(etasp,exisp)

			!Coordinates
			do idime=1,ndim
				gpcod(idime,kelgs)=0.d0
				do inode=1,nnode
					gpcod(idime,kelgs)=gpcod(idime,kelgs)+elcod(idime,inode)*shap(inode)
				end do
			end do


			xgash=(defsg(1,kgaus)+defsg(2,kgaus))*0.5d0
			xgish=(defsg(1,kgaus)-defsg(2,kgaus))*0.5d0
			xgesh=defsg(3,kgaus)
			xgosh=dsqrt(xgish*xgish+xgesh*xgesh)
			strsp(1)=xgash+xgosh
			strsp(2)=xgash-xgosh
			if(xgish == 0.d0) xgish=0.1d-20
	        strsp(3)=datan(xgesh/xgish)*28.647889757d0  !=90/pi

61 write(116,941) ielem,1,kelgs,(gpcod(idime,kelgs),idime=1,2),(defsg(istr1,kgaus),istr1=1,4),(strsp(istre),istre=1,3)!,epstn(kgaus)


!---------------------------------------------------------
! write STRESSES
!---------------------------------------------------------
open(unit=117,file='PostProcessing/Output/stresses.txt')

if(ntype /= 3) write(117,871)
871 format(7hElement,3x,4hg.p.,3x,3hNo.,9x,1hx,15x,1hy,12x,8hSIGMA_XX,8x,8hSIGMA_YY,8x,8hSIGMA_XY,8x,8hSIGMA_ZZ,8x,&
	&8hmax p.s.,8x,8hmin p.s.,7x,5hangle)!,3x,6he.p.s.)
if(ntype == 3) write(117,876)
876 format(7hElement,3x,4hg.p.,3x,3hNo.,9x,1hx,15x,1hy,12x,8hSIGMA_RR,8x,8hSIGMA_ZZ,8x,8hSIGMA_RZ,8x,8hSIGMA_TT,8x,&
	&8hmax p.s.,8x,8hmin p.s.,7x,5hangle)!,3x,6he.p.s.)

! 841   format(i7,2x,i3,2x,i4,2x,8e16.6,f10.3)
841   format(i7,2x,i3,2x,i4,2x,8e16.6,f10.3,2x,8e16.6)
! Calculate coef to change base of the coordinates system
if (ngaus==2) then
	alpha_extrap=dsqrt(3.d0)
else if (ngaus==3) then
	alpha_extrap=dsqrt(5.d0/3.d0)
end if
! print *,"alpha=",alpha_extrap
! Extrolation of the nodal values from the gauss points


! defsg(1,1)=15.0
! defsg(1,2)=10.0
! defsg(1,3)=20.0
! defsg(1,4)=15.0


! defsg(1,1)=10.0
! defsg(1,2)=15.0
! defsg(1,3)=15.0
! defsg(1,4)=20.0


! defsg(1,1)=10.0
! defsg(1,2)=15.0
! defsg(1,3)=20.0
! defsg(1,4)=15.0
! defsg(1,5)=20.0
! defsg(1,6)=25.0
! defsg(1,7)=20.0
! defsg(1,8)=25.0
! defsg(1,9)=30.0

! defsg(1,10)=25.0
! defsg(1,11)=30.0
! defsg(1,12)=35.0
! defsg(1,13)=30.0
! defsg(1,14)=35.0
! defsg(1,15)=40.0
! defsg(1,16)=35.0
! defsg(1,17)=40.0
! defsg(1,18)=45.0

! defsg(1,19)=40.0
! defsg(1,20)=45.0
! defsg(1,21)=50.0
! defsg(1,22)=45.0
! defsg(1,23)=50.0
! defsg(1,24)=55.0
! defsg(1,25)=50.0
! defsg(1,26)=55.0
! defsg(1,27)=60.0

! defsg(1,28)=25.0
! defsg(1,29)=30.0
! defsg(1,30)=35.0
! defsg(1,31)=30.0
! defsg(1,32)=35.0
! defsg(1,33)=40.0
! defsg(1,34)=35.0
! defsg(1,35)=40.0
! defsg(1,36)=45.0


!Extrapolate values from gauss points to nodes
Vec_Node_Hist=0.d0
def_node2=0.d0
plastic_node=0.d0

do ielem=1,nelem
	! Set to zero matrix where extrapolated values at nodes for each element will be stored
	def_node=0.d0
	kgaus_rep=(ielem-1)*ngaus*ngaus

	do inode=1,nnode

		! Calculate the coordinates of the element nodal points
! 		elcod=0.0d0
		lnode=iabs(lnods(ielem,inode))
! 		print *,'lnode=',lnode

		Vec_Node_Hist(lnode)=Vec_Node_Hist(lnode)+1.d0
! 		do idime=1,ndim
! 			elcod(idime,inode)=coord(lnode,idime)
! 		end do

		! Calculate coordinates of the node for the extrapolation
		ksi=alpha_extrap*ksi_coef(inode)
		eta=alpha_extrap*eta_coef(inode)
! 		print *,"ksi=",ksi,"eta=",eta

		! Calculate shape functions at the node
		call sfr2_extrap(ksi,eta)

! 		print *, (shap(i),i=1,nnode)

		kshap=0
		kgaus=kgaus_rep

		! Extrapolation loop over all the gauss points
		do igaus=1,ngaus
			do jgaus=1,ngaus
				! Increment gauss point number
				kgaus=kgaus+1

				! We only use middle gauss point if 9 noded element
! 				if (nnode==8 .and. igaus==2 .and. jgaus==2) go to 457

				! Increment shape function number
				kshap=kshap+1

				! Extrapolate all the values eps_xx, eps_yy, eps_xy, eps_zz
				do istr1=1,nstr1
					def_node2(lnode,istr1) = def_node2(lnode,istr1) + strsg(istr1,kgaus)*shap(kshap)
				end do
! 				def_node(1,inode) = def_node(1,inode) + strsg(1,kgaus)*shap(kshap)

				!Extrapolate plastic deformations
				plastic_node(lnode) = plastic_node(lnode) + epstn(kgaus)*shap(kshap)

! 				print *,'value=',strsg(1,kgaus),"shap=",shap(kshap),"prod=",strsg(1,kgaus)*shap(kshap)
				457 continue
	    	end do
		end do
	end do


	if (ngaus==2) then
		!******************************************
		! Calculate value at element centroid
		!******************************************

		centroid_cord=0.d0
		centroid_def=0.d0

		! Calculate coordinates of the centroid for the extrapolation
		ksi=0.d0
		eta=0.d0
! 		print *,"ksi=",ksi,"eta=",eta

		! Calculate coordinates of nodes of the element
		do inode=1,nnode
			lnode=iabs(lnods(ielem,inode))

			do idime=1,ndim
				elcod(idime,inode)=coord(lnode,idime)
			end do
		end do

		!**Calculate coordinates centroid in global coordinates

		!Shape functions
		call sfr2(ksi,eta)

		!Coordinates
		do idime=1,ndim
			centroid_cord(idime)=0.d0
			do inode=1,nnode
				centroid_cord(idime)=centroid_cord(idime)+elcod(idime,inode)*shap(inode)
			end do
		end do

		!Calculate values at centroid

		! Calculate shape functions at the node
		call sfr2_extrap(ksi,eta)

! 		print *, shap(1),shap(2),shap(3),shap(4)

		kshap=0
		kgaus=kgaus_rep

		! Extrapolation loop over all the gauss points

		do igaus=1,ngaus
			do jgaus=1,ngaus
				! Increment gauss point number
				kgaus=kgaus+1

				kshap=kshap+1

				! Extrapolate all the values eps_xx, eps_yy, eps_xy, eps_zz
				do istr1=1,nstr1
					centroid_def(istr1) = centroid_def(istr1) + strsg(istr1,kgaus)*shap(kshap)
				end do
			end do
		end do

		xgash=(centroid_def(1)+centroid_def(2))*0.5d0
		xgish=(centroid_def(1)-centroid_def(2))*0.5d0
		xgesh=centroid_def(3)
		xgosh=dsqrt(xgish*xgish+xgesh*xgesh)
		strsp(1)=xgash+xgosh
		strsp(2)=xgash-xgosh
		if(xgish == 0.d0) xgish=0.1d-20
		strsp(3)=datan(xgesh/xgish)*28.647889757d0  !=90/pi

		!Write values at centroid
		write(117,841) ielem,2,0,(centroid_cord(idime),idime=1,2),(centroid_def(istr1),istr1=1,nstr1),(strsp(istre),istre=1,3)
	end if
end do


! Average values and write into text file
do ipoin=1,npoin
	do istr1=1,nstr1
		def_node2(ipoin,istr1) = def_node2(ipoin,istr1)/Vec_Node_Hist(ipoin)
	end do

	!Plastic deformations
	plastic_node(ipoin) = plastic_node(ipoin)/Vec_Node_Hist(ipoin)

	xgash=(def_node2(ipoin,1)+def_node2(ipoin,2))*0.5d0
	xgish=(def_node2(ipoin,1)-def_node2(ipoin,2))*0.5d0
	xgesh=def_node2(ipoin,3)
	xgosh=dsqrt(xgish*xgish+xgesh*xgesh)
	strsp(1)=xgash+xgosh
	strsp(2)=xgash-xgosh
	if(xgish == 0.d0) xgish=0.1d-20
	strsp(3)=datan(xgesh/xgish)*28.647889757d0  !=90/pi

	write(117,841) 0,0,ipoin,(coord(ipoin,idime),idime=1,2),(def_node2(ipoin,istr1),istr1=1,nstr1),(strsp(istre), &
	    	& istre=1,3),plastic_node(ipoin)

end do



! Write coordinates of gp and the values
kgaus=0
do 62 ielem=1,nelem
	kelgs=0
! 	write (6,931) ielem
! 931 format(6x,13helement no. =,i5)
	! Evaluate the coordinates of the element nodal points
	do inode=1,nnode
		lnode=iabs(lnods(ielem,inode))

		do idime=1,ndim
			elcod(idime,inode)=coord(lnode,idime)
		end do
	end do

	do 62 igaus=1,ngaus
		exisp=posgp(igaus)
		do 62 jgaus=1,ngaus
			etasp=posgp(jgaus)
			kgaus=kgaus+1
			kelgs=kelgs+1

			!**Calculate coordinates Gauss Points in global coordinates

			!Shape functions
			call sfr2(etasp,exisp)

			!Coordinates
			do idime=1,ndim
				gpcod(idime,kelgs)=0.d0
				do inode=1,nnode
					gpcod(idime,kelgs)=gpcod(idime,kelgs)+elcod(idime,inode)*shap(inode)
				end do
			end do


			xgash=(strsg(1,kgaus)+strsg(2,kgaus))*0.5d0
			xgish=(strsg(1,kgaus)-strsg(2,kgaus))*0.5d0
			xgesh=strsg(3,kgaus)
			xgosh=dsqrt(xgish*xgish+xgesh*xgesh)
			strsp(1)=xgash+xgosh
			strsp(2)=xgash-xgosh
			if(xgish == 0.d0) xgish=0.1d-20
	        strsp(3)=datan(xgesh/xgish)*28.647889757d0  !=90/pi

62 write(117,841) ielem,1,kelgs,(gpcod(idime,kelgs),idime=1,2),(strsg(istr1,kgaus),istr1=1,4),(strsp(istre),istre=1,3)!,epstn(kgaus)


!---------------------------------------------------------
! write STRESSES
!---------------------------------------------------------
! open(unit=115,file='PostProcessing/Output/stresses.txt')

! if(ntype /= 3) write(115,970)
!   970 format(2x,4hg.p.,6x,9hxx-stress,5x,9hyy-stress,5x,9hxy-stress,5x,9hzz-stress,6x,8hmax p.s.,6x,8hmin p.s.,3x, &
!                     5hangle,3x,6he.p.s.)
!       if(ntype == 3) write(115,975)
!   975 format(2x,4hg.p.,6x,9hrr-stress,5x,9hzz-stress,5x,9hrz-stress,5x,9htt-stress,6x,8hmax p.s.,6x,8hmin p,s.,3x, &
!                     5hangle,3x,6he.p.s.)


!  kgaus=0
!       do 60 ielem=1,nelem
!         kelgs=0
!         write (115,930) ielem
!   930   format(6x,13helement no. =,i5)

!         do 60 igaus=1,ngaus
!             do 60 jgaus=1,ngaus
!                 kgaus=kgaus+1
!                 kelgs=kelgs+1
!                 xgash=(strsg(1,kgaus)+strsg(2,kgaus))*0.5d0
!                 xgish=(strsg(1,kgaus)-strsg(2,kgaus))*0.5d0
!                 xgesh=strsg(3,kgaus)
!                 xgosh=dsqrt(xgish*xgish+xgesh*xgesh)
!                 strsp(1)=xgash+xgosh
!                 strsp(2)=xgash-xgosh
!                 if(xgish == 0.d0) xgish=0.1d-20
!                 strsp(3)=datan(xgesh/xgish)*28.647889757d0  !=90/pi

!  60   write(115,940) kelgs,(strsg(istr1,kgaus),istr1=1,4),(strsp(istre),istre=1,3),epstn(kgaus)
! 940   format(i5,2x,6e14.6,f8.3,e14.6)


! ! jacob2(ielem,kgasp)
! ! integer ::  idime, i, inode, j, jdime, kgasp, ielem
! ! calculate coordinates of sampling point
! ! Eq. 6.43 p. 171
! do idime=1,ndim
!     gpcod(idime,kgasp)=0.d0
!     do inode=1,nnode
!         gpcod(idime,kgasp)=gpcod(idime,kgasp)+elcod(idime,inode)*shap(inode)
!     end do
! end do


end subroutine


