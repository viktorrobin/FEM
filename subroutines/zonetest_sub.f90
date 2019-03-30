subroutine affiche(t)

	implicit none
	integer, dimension(:,:), allocatable  :: t
	integer :: i, j, n

	interface
		subroutine affiche2(t)
			integer, dimension(:,:), allocatable  :: t
		end subroutine
	end interface


	do i=1,size(t,1)
		print '(1x,160i5)', t(i,:)
	end do


end subroutine

subroutine affiche2(t)

	implicit none
	integer, dimension(:,:), allocatable  :: t
	integer :: i, j, n

	do i=1,size(t,1)
		print '(1x,160i5)', t(i,:)
	end do


end subroutine

subroutine zero(t,nblignes,nbcol)

	implicit none
	integer, dimension(nblignes,nbcol)  :: t
	integer :: i, j, n,nbcol,nblignes


		do i=1,nblignes
			do j=1,nbcol
				t(i,j)=0.d0
			end do
		end do

	return
end subroutine

