subroutine output_vtk( )
	use variables
	implicit none
	character(len=50) :: title_file, numfile_char,part4,part1, npoin_char, nelem_char,x_char, y_char, z_char, nelem_tot_char
	character(len=50) :: n1_c, n2_c, n3_c, n4_c, n5_c, n6_c, n7_c, n8_c, n9_c, nnode_char
	integer :: i, j, ipoin, ngash,ngish,igish,igash, nelem_tot
	double precision :: z=0.d0

	!http://people.sc.fsu.edu/~jburkardt/data/vtk/vtk.html
	!http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html

	!Creating name file
	part1='output/'
	title_file='output_vtk_'
	write(numfile_char,'(i10)') iincs
    part4='.vtk'

    ! Concatenation - trim supprime tous les espaces de fin
    title_file= trim(part1) // trim(title_file) // trim(numfile_char) // trim(part4) 

  	!Creating output file
	open(13,file=title_file)
	
	!Header
	write(13,900) '# vtk DataFile Version 1.0'
	write(13,900) 'Unstructured Grid Example'
	write(13,900) 'ASCII'
	900 format(a)

	write(13,901) 'DATASET UNSTRUCTURED_GRID'
	901 format(/a)

	!Pointdata
		!header
	write(npoin_char,'(i10)') npoin
	npoin_char=adjustl(npoin_char)
	write(13,910) trim(npoin_char)
	910 format(7hPOINTS ,a,6h float)

		!Coordonnees
	z_char='0'
	do i=1,npoin
		write(x_char,'(f10.6)') coord(i,1)
		write(y_char,'(f10.6)') coord(i,2)
		write(z_char,'(f10.6)') z
		x_char=adjustl(x_char)
		y_char=adjustl(y_char)
		z_char=adjustl(z_char)
		write(13,920) trim(x_char), trim(y_char), trim(z_char)
	end do
	920 format(a,1x,a,1x,a)
	
	
	!Celldata
		!header
	nelem_tot=nelem*(nnode+1)
	write(nelem_char,'(i10)') nelem
	write(nelem_tot_char,'(i10)') nelem_tot
	nelem_char=adjustl(nelem_char)
	nelem_tot_char=adjustl(nelem_tot_char)
	write(13,930) trim(nelem_char), trim(nelem_tot_char)
	930 format(/6hCELLS ,a,1x,a)

		!Connectivity table
	write(nnode_char,'(i10)') nnode
	nnode_char=adjustl(nnode_char)

	do i=1,nelem
		
! 		write(n1_c,'(i10)') (iabs(lnods(i,1))-1)
! 		write(n2_c,'(i10)') (iabs(lnods(i,2))-1)
! 		write(n3_c,'(i10)') (iabs(lnods(i,3))-1)
! 		write(n4_c,'(i10)') (iabs(lnods(i,4))-1)
! 		write(n5_c,'(i10)') (iabs(lnods(i,5))-1)
! 		write(n6_c,'(i10)') (iabs(lnods(i,6))-1)
! 		write(n7_c,'(i10)') (iabs(lnods(i,7))-1)
! 		write(n8_c,'(i10)') (iabs(lnods(i,8))-1)
! 		if (nnode == 9) then
! 			write(n9_c,'(i10)') (iabs(lnods(i,9))-1)
! 		end if

! 		n1_c=adjustl(n1_c)
! 		n2_c=adjustl(n2_c)
! 		n3_c=adjustl(n3_c)
! 		n4_c=adjustl(n4_c)
! 		n5_c=adjustl(n5_c)
! 		n6_c=adjustl(n6_c)
! 		n7_c=adjustl(n7_c)
! 		n8_c=adjustl(n8_c)
! 		if (nnode == 9) then
! 			n9_c=adjustl(n9_c)
! 		end if
		select case (nnode)
			case (3)
				write(n1_c,'(i10)') (iabs(lnods(i,1))-1)
				write(n2_c,'(i10)') (iabs(lnods(i,2))-1)
				write(n3_c,'(i10)') (iabs(lnods(i,3))-1)
				n1_c=adjustl(n1_c)
				n2_c=adjustl(n2_c)
				n3_c=adjustl(n3_c)
				write(13,940) trim(nnode_char), trim(n1_c), trim(n2_c), trim(n3_c)
				940 format(a,1x,a,1x,a,1x,a)

			case(4)
				write(n1_c,'(i10)') (iabs(lnods(i,1))-1)
				write(n2_c,'(i10)') (iabs(lnods(i,2))-1)
				write(n3_c,'(i10)') (iabs(lnods(i,3))-1)
				write(n4_c,'(i10)') (iabs(lnods(i,4))-1)
				n1_c=adjustl(n1_c)
				n2_c=adjustl(n2_c)
				n3_c=adjustl(n3_c)
				n4_c=adjustl(n4_c)
				write(13,944) trim(nnode_char), trim(n1_c), trim(n2_c), trim(n4_c), trim(n3_c)
				944 format(a,1x,a,1x,a,1x,a,1x,a,1x)

			case(8)
				write(n1_c,'(i10)') (iabs(lnods(i,1))-1)
				write(n2_c,'(i10)') (iabs(lnods(i,2))-1)
				write(n3_c,'(i10)') (iabs(lnods(i,3))-1)
				write(n4_c,'(i10)') (iabs(lnods(i,4))-1)
				write(n5_c,'(i10)') (iabs(lnods(i,5))-1)
				write(n6_c,'(i10)') (iabs(lnods(i,6))-1)
				write(n7_c,'(i10)') (iabs(lnods(i,7))-1)
				write(n8_c,'(i10)') (iabs(lnods(i,8))-1)
				n1_c=adjustl(n1_c)
				n2_c=adjustl(n2_c)
				n3_c=adjustl(n3_c)
				n4_c=adjustl(n4_c)
				n5_c=adjustl(n5_c)
				n6_c=adjustl(n6_c)
				n7_c=adjustl(n7_c)
				n8_c=adjustl(n8_c)
				write(13,941) trim(nnode_char), trim(n1_c), trim(n3_c), trim(n5_c), trim(n7_c), &
					& trim(n2_c), trim(n4_c), trim(n6_c), trim(n8_c)
				941 format(a,1x,a,1x,a,1x,a,1x,a,1x,a,1x,a,1x,a,1x,a)

		end select
	end do
		

	!nnode=nombre de noeuds par elements
	!lnods=table de connectivité (ligne 1=element 1, ligne 2=element 2, etc...)
	!nelem=nombre d'elements total dans la structure
	!npoin=nombre total de points dans la structure
	!coord=coordonnées des points (attention ici c'est en 2D)
	!ndim=nombre de dimensions

    
	!Cell_Type
		!header
	write(13,950) trim(nelem_char)
	950 format(/11hCELL_TYPES ,a)

	select case (nnode)
		case(3)
			do i=1,nelem
				write(13,'(a)') '5'
			end do

		case(4)
			do i=1,nelem
				write(13,'(a)') '8'
			end do

		case(8)
			do i=1,nelem
				write(13,'(a)') '23'
			end do
	end select

	! Point data
		!header
	write(13,970) trim(npoin_char)
	970 format(/11hPOINT_DATA ,a)
	write(13,'(a)') 'SCALARS Fx float'
	write(13,'(a)') 'LOOKUP_TABLE default'
	
	do ipoin=1,npoin
	    ngash=ipoin*2
	    ngish=ngash-2+1

	    write(x_char,'(f12.6)') tofor(ngish)
		write(y_char,'(f10.6)') tofor(ngash)*100
		write(z_char,'(f10.6)') z

		x_char=adjustl(x_char)
		y_char=adjustl(y_char)
		z_char=adjustl(z_char)

		write(13,920) trim(y_char)!, trim(y_char), trim(z_char)
		!write(13,*) (tdisp(igash),igash=ngish,ngash),0.0
	end do

	

	close(13)
end
