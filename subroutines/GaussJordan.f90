! Pivot de Gauss

subroutine GaussJordan()
use variables

implicit none 
double precision, dimension(nvlib) :: fresv
integer :: i, j, kount, ieqns, icols, ieqn1, irows, jevab

double precision :: factr, pivot

integer :: nback,nbac1, ipoin, idofn, neqn1

double precision, dimension(nvlib) :: xdisp, react
double precision, dimension(nvlib,2) :: gtdisp, gtreac
double precision, dimension(nvlib,nvlib) :: grgstif
double precision :: resid, k1, k2, k3, k4

!Gaussian Reduction
kount=0

do i=1,nvlib
    do j=1,nvlib
    grgstif(i,j)=rgstif(i,j)
end do
end do

call disp_rgstif()

print *,'gload=',rgload
print *,'rfixed=',rfixed
print *,'riffix=',riffix
! stop

react=0.d0
xdisp=0.d0

do i=1,size(rgload,1)
    rgload(i)=vec_residu(i)
end do


do 70 ieqns=1,nvlib
    if (riffix(ieqns) == 1) go to 40

        ! Get the pivot
        pivot=rgstif(ieqns,ieqns)
        print *,'pivot=',pivot

        if (abs(pivot) < 1.0d-9) then
            write(6,*) 'Incorrect pivot=',pivot,'line=',ieqns
            stop
        end if

        if (ieqns==nvlib) go to 70

        ieqn1=ieqns+1

        do 30 irows=ieqn1,nvlib
        	kount=kount+1
!             if (rgstif(irows,ieqns) /= 0.d0) then
                factr=rgstif(irows,ieqns)/pivot
                fresv(kount)=factr

                if (factr==0.d0) go to 30

                ! Modify matrix [K]=rgstif
                do icols=ieqns,nvlib
                    rgstif(irows,icols)=rgstif(irows,icols) - factr*rgstif(ieqns,icols)
                end do

                ! Modify rgload
                rgload(irows)=rgload(irows) - factr*rgload(ieqns)
        30 continue
!             end if
        
                go to 70

!     if (riffix(ieqns) == 1) then

        !Displacement is prescribed; we substitute the prescribed value into the remaining equations
        40 continue
        do irows=ieqns,nvlib        
!             if (rgstif(irows,ieqns) /= 0.d0) then
                rgload(irows)=rgload(irows)-rgstif(irows,ieqns)*rfixed(ieqns)
!                 rgstif(irows,ieqns)=0.d0
!             end if
        end do
!     end if
70 continue


! print *,'rgload=',rgload

! do ieqns=1,nvlib
!     print *,(rgstif(ieqns,jevab),jevab=1,nvlib)
! end do


!Back-substitution routine

do 100 ieqns=1,nvlib
	react(ieqns)=0.d0
100 continue

neqn1=nvlib+1

do 400 ieqns=1,nvlib
	nback=neqn1-ieqns

	pivot=rgstif(nback,nback)
	resid=rgload(nback)

	if (nback == nvlib) go to 300

	nbac1=nback+1

	do 200 icols=nbac1,nvlib
		resid=resid - rgstif(nback,icols)*xdisp(icols)
	200 continue

300	if (riffix(nback) == 0) xdisp(nback)=resid/pivot
	if (riffix(nback) == 1) xdisp(nback)=rfixed(nback)
	if (riffix(nback) == 1) react(nback) = -resid



400 continue

kount=0
do ipoin=1,npoin
	do idofn=1,ndofn
		kount=kount+1
		gtdisp(ipoin,idofn)=xdisp(kount)
		gtreac(ipoin,idofn)=react(kount)
	end do
end do

print *, 'react=',react
print *, 'disp=',xdisp

do i=1,nvlib
    do j=1,nvlib
    rgstif(i,j)=grgstif(i,j)
end do
end do

! asdis=xdisp
! vec_residu_big=rgload
rasdis=xdisp
print *,vec_residu
! stop
end subroutine

