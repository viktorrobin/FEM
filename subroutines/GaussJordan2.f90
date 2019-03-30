! Pivot de Gauss

subroutine GaussJordan2()
use variables

implicit none 
double precision, dimension(ntotv) :: fresv
integer :: i, j, kount, ieqns, icols, ieqn1, irows, jevab

double precision :: factr, pivot

integer :: nback,nbac1, ipoin, idofn, neqn1

double precision, dimension(ntotv) :: xdisp, react
double precision, dimension(ntotv,2) :: gtdisp, gtreac

double precision :: resid, k1, k2, k3, k4

!Gaussian Reduction
kount=0

! call disp_gstif()

print *,'gload=',gload
print *,'fixed=',fixed
print *,'iffix=',iffix
! stop

react=0.d0
xdisp=0.d0

! do i=1,size(gload,1)
!     gload(i)=vec_residu(i)
! end do

! if (iiter>1) then
!     gload=vec_residu
! end if

do 70 ieqns=1,ntotv
    if (iffix(ieqns) == 1) go to 40

        ! Get the pivot
        pivot=gstif(ieqns,ieqns)
        print *,'pivot=',pivot

        if (abs(pivot) < 1.0d-9) then
            write(6,*) 'Incorrect pivot=',pivot,'line=',ieqns
            stop
        end if

        if (ieqns==ntotv) go to 70

        ieqn1=ieqns+1

        do 30 irows=ieqn1,ntotv
        	kount=kount+1
!             if (gstif(irows,ieqns) /= 0.d0) then
                factr=gstif(irows,ieqns)/pivot
                fresv(kount)=factr

                if (factr==0.d0) go to 30

                ! Modify matrix [K]=gstif
                do icols=ieqns,ntotv
                    gstif(irows,icols)=gstif(irows,icols) - factr*gstif(ieqns,icols)
                end do

                ! Modify gload
                gload(irows)=gload(irows) - factr*gload(ieqns)
        30 continue
!             end if
        
                go to 70

!     if (iffix(ieqns) == 1) then

        !Displacement is prescribed; we substitute the prescribed value into the remaining equations
        40 continue
        do irows=ieqns,ntotv        
!             if (gstif(irows,ieqns) /= 0.d0) then
                gload(irows)=gload(irows)-gstif(irows,ieqns)*fixed(ieqns)
!                 gstif(irows,ieqns)=0.d0
!             end if
        end do
!     end if
70 continue


! print *,'gload=',gload

! do ieqns=1,ntotv
!     print *,(gstif(ieqns,jevab),jevab=1,ntotv)
! end do


!Back-substitution routine

do 100 ieqns=1,ntotv
	react(ieqns)=0.d0
100 continue

neqn1=ntotv+1

do 400 ieqns=1,ntotv
	nback=neqn1-ieqns

	pivot=gstif(nback,nback)
	resid=gload(nback)

	if (nback == ntotv) go to 300

	nbac1=nback+1

	do 200 icols=nbac1,ntotv
		resid=resid - gstif(nback,icols)*xdisp(icols)
	200 continue

300	if (iffix(nback) == 0) xdisp(nback)=resid/pivot
	if (iffix(nback) == 1) xdisp(nback)=fixed(nback)
	if (iffix(nback) == 1) react(nback) = -resid



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



! asdis=xdisp
gload=react
asdis=xdisp
! print *,vec_residu
! stop
end subroutine