
!******************************************************************
!   Subroutine jacob2 - Owen p.182
!   Calculates, for any sampling position (ksi_p,eta_p) (usually the gauss point) the following quantity:
!       * The cartesian coordinates of the gauss point which are stored in the array gpcod()
!       * The Jacobian matrix which is stored in xjcam()
!       * The determinant of the Jacobian matrix, djacb
!       * The inverse of the Jacobian matrix which is stored in xjaci()
!       * The cartesian derivatives dNi/dx, dNi/dy (or dNi/dr, dNi/dz) of the element shape functions
! *******************************************************

!subroutine jacob2(xjacm,cartd,deriv,djacb,elcod,gpcod,ielem,kgasp,nnode,shap,cardc,cartdt,kgaus,gpcodg)
!subroutine jacob2(ielem, kgasp)
subroutine jacob2(ielem,kgasp)

use variables

implicit none

integer ::  idime, i, inode, j, jdime, kgasp, ielem



!if (nnode>3) then
! calculate coordinates of sampling point
! Eq. 6.43 p. 171
do idime=1,ndim
    gpcod(idime,kgasp)=0.d0
    do inode=1,nnode
        gpcod(idime,kgasp)=gpcod(idime,kgasp)+elcod(idime,inode)*shap(inode)
    end do
end do

! print *,gpcod(1,kgasp),gpcod(2,kgasp)
! pause

! write(6,*) elcod
! write(6,*) shap
! stop
!end if

!if (nnode==3) then
!do idime=1,ndim
!    gpcod(idime,kgasp)=0.0
!    do inode=1,nnode
!        gpcod(idime,kgasp)=elcod(idime,inode)
!    end do
!end do
!
!end if

!do idime=1,2
!    write(6,*) (gpcod(idime,inode),inode=1,1)
!end do
!stop


!       kgaus=1 !ajouter par victor

!        write(*,*) 'ok jacob2',kgaus, kgasp
!        pause

!do idime=1,2
!    gpcodg(idime,kgaus)=gpcod(idime,kgasp)
!end do



!***    create jacobian matrix , xjacm
!
do idime=1,2
    do jdime=1,2
        xjacm(idime,jdime)=0.d0
        do inode=1,nnode
            xjacm(idime,jdime)=xjacm(idime,jdime)+deriv(idime,inode)*elcod(jdime,inode)
        end do
    end do
end do


!***    calculates determinant and inverse of jacobian matrix
djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)

! write(6,*) djacb
! stop

if (djacb<0.d0) then
    write(6,600) ielem
    stop
else if (djacb==0.d0) then
    write(6,600) ielem
    stop
end if

600 format(//,36h programm halted in subroutin jacob2,/,10x,22h zero or negative area,/,10x,16h element number ,i5)    

xjaci(1,1)=xjacm(2,2)/djacb
xjaci(2,2)=xjacm(1,1)/djacb
xjaci(1,2)=-xjacm(1,2)/djacb
xjaci(2,1)=-xjacm(2,1)/djacb

!do idime=1,2
!    write(6,*) (xjaci(idime,jdime),jdime=1,2)
!end do
!stop

!***    calculate cartesian derivatives
!
do idime=1,2
    do inode=1,nnode
        cartd(idime,inode)=0.d0
        do jdime=1,2
            cartd(idime,inode)=cartd(idime,inode)+xjaci(idime,jdime)*deriv(jdime,inode)
        end do
    end do
end do

! do idime=1,ndim
!    write(6,*) (cartd(idime,jdime),jdime=1,nnode)
! end do
!stop


! change 5
!do i=1,nnode
!    do j=1,2
!        cartdt(i,j)=cartd(j,i)
!    end do
!end do
!  end of change 5




!! store in a vector form
!do j=1,nnode
!    cardc(2*j-1,1)=cartd(1,j)
!    cardc(2*j,1)=cartd(2,j)
!end do

return

end




!   *****************************************************************
!   This subroutine sets up the gauss-legendre integration constants Owen p.179
!   Its function is to set up the sampling point positions and weighting factors for numerical integration.
!   These processes are rectricted to either two or three point integration rules. cf. def. ngaus dans Module.f90
!   ******************************************************************

!subroutine gaussq(ngaus,posgp,weigp)
subroutine gaussq( )

    use variables
        
    integer :: igash, kgaus
    !ngaus=2 avec exemple de Owen
    ! cf. http://fr.wikipedia.org/wiki/M%C3%A9thodes_de_quadrature_de_Gauss
    ! weigp: poids de l'interpolation
    ! posgp: Coordonnées des points du support pour l'interpolation 



    !if (ngaus > 2) then
    !        posgp(1)=-0.77459666924148337703585307995648d0  ! = -(3/5)^(1/2)
    !        posgp(2)=0.d0
    !        weigp(1)=0.74535599249992989880305788957709d0  ! = 5/9
    !        weigp(2)=0.8888888888888889d0  ! = 8/9 
    !    else ! ie if ngaus<=2
    !        posgp(1)=-0.9428090415820633658677924828065d0  !=-(1/3)^(1/2)
    !        weigp(1)=1.d0
    !end if
    
        if (ngaus == 3) then
                posgp(1)=-0.77459666924148337703585307995648d0  ! = -(3/5)^(1/2)
                posgp(2)=0.d0
                weigp(1)=0.55555555555555555555555555555556d0  ! = 5/9
                weigp(2)=0.8888888888888888888888888888889d0  ! = 8/9

        elseif (ngaus<=2) then! ie if ngaus<=2
                posgp(1)=-0.57735026918962576450914878050196d0  !=-(1/3)^(1/2)
                weigp(1)=1.d0
        end if

          
     
    kgaus=ngaus/2

    do igash=1,kgaus
        jgash=ngaus+1-igash
        posgp(jgash)=-posgp(igash)
        weigp(jgash)=weigp(igash)
    end do



!     do igash=1,ngaus
!     write(6,*) posgp(igash)
!     end do
!     stop





return
end







! -------------------------- MODULE cg.f90 ----------------------------

!************************************************************************
!*                                                                      *
!* Conjugate Gradient Method (CG Method)                                *
!* -------------------------------------                                *
!*                                                                      *
!* Programming language: ANSI C                                         *
!* Compiler:             Turbo C 2.0                                    *
!* Computer:             IBM PS/2 70 with 80387                         *
!* Sources:              [BUNS85], [SCHW], [MAES84]                     *
!* Author:               Juergen Dietel, Computer Center, RWTH Aachen   *
!* Date:                 7.31.1992                                      *
!*                                                                      *
!*             F90 version by J-P Moreau (without dynamic allocations). *
!************************************************************************
Subroutine cg_method (A,X,B,n)!  &     ! Conjugate Gradient Method
                       !n, &     ! Size of the linear system
                       !a, &     ! System matrix
                       !y, &     ! right hand side
                       !x, &     ! solution vector
                       !fehler & ! error code
                     !)
integer n
  double precision, dimension(n,n) :: A
  double precision, dimension(n) :: B, X
! original name: cg_verfahren()
integer, parameter :: SIZE=24
double precision, parameter :: ZERO=0.d0, MACH_EPS=1.0d-15
!real*8  a(0:SIZE,0:SIZE),x(0:SIZE),y(0:SIZE), delta(0:SIZE)
integer fehler

!************************************************************************
!* CG solves the linear system                                          *
!*                         A * X = Y                                    *
!* for a symmetric, positive definite matrix A via the conjugate        *
!* gradient method.                                                     *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* n  Size of the linear system                                         *
!* a  [0..n-1,0..n-1] system matrix A. Only the upper triangle of A is  *
!*    used.                                                             *
!* y  [0..n-1] vector of the right hand side                            *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* x  [0..n-1] vector giving the solution                               *
!*                                                                      *
!* Return value:                                                        *
!* =============                                                        *
!* = 0: all is ok                                                       *
!* = 1: n < 2 or other disallowed input parameters                      *
!* = 2: memory exceeded                                                 *
!*                                                                      *
!************************************************************************
double precision d(n), &   ! (0..n-1) auxiliary vectors d and g
         g(n), &
         AmalD(n)  ! (0..n-1) auxiliary vector A * d

double precision  alpha,    &   ! coefficient
         beta,     &   ! coefficient
         dividend, &   ! numerator and denominator of a fraction
         divisor,  &   ! respectively, used to compute alpha, beta
         hilf,     &   ! auxiliary variables
         hilf2,    &
         abstand,  &   ! distance of two successive approximations
                       ! for the solution vector x (taken in the
                       ! euclidean norm)
         xnorm         ! euklidean norm of x
  integer k, i, j      ! loop variables

  if (n < 2) then      ! invalid parameter?
    fehler=1
    print *,'n<2'
    return
  end if



    !------------------------------------------------------------------
  ! start with x at the origin                                       
  !------------------------------------------------------------------
  do i = n , 1, -1
    X(i) = ZERO
  end do

  !------------------------------------------------------------------
  ! initialize  d and g :                                            
  ! d = -g = -(a*x - y) = y (since x = 0)                            
  !------------------------------------------------------------------
  do i = n, 1, -1
    hilf = B(i)
    d(i) = hilf
    g(i) = -hilf
  end do


  !------------------------------------------------------------------
  ! perform at most n steps of the CG Method                         
  !------------------------------------------------------------------
  do k = n, 0, -1
    !----------------------------------------------------------------
    ! compute new alpha:                                             
    ! alpha = -(d(transp) * g) / (d(transp) * (a * d))               
    !----------------------------------------------------------------

    dividend = ZERO
    divisor  = ZERO

    do i = n, 1, -1
      dividend = dividend + d(i) * g(i)
      hilf = ZERO
      do j = 1, i-1
        hilf = hilf + A(j,i) * d(j)
      end do
      do j = i, n
        hilf = hilf + A(i,j) * d(j)
      end do
      AmalD(i) = hilf
      divisor = divisor + d(i) * hilf
    end do

    if (divisor.eq.ZERO) then
      fehler=0
!       print *,'divisor=0'
      return
    end if

    alpha = -dividend / divisor

    !----------------------------------------------------------------
    ! compute the norm of x und  alpha * d  and find a new x:        
    ! x  =  x + alpha * d, then check whether x is close enough,     
    ! in order to stop the process before n complete steps           
    !----------------------------------------------------------------
    xnorm   = ZERO
    abstand = ZERO

    do i = n , 1, -1
      hilf =  X(i)
      xnorm   = xnorm + hilf*hilf
      hilf2   =  alpha * d(i)
      abstand = abstand + hilf2*hilf2
      X(i)    =  hilf + hilf2
    end do

    if (abstand < MACH_EPS * xnorm) then
      fehler=0
      return
    end if


    !----------------------------------------------------------------
    ! compute new g:   g  =  g + alpha * (a * d)                     
    !----------------------------------------------------------------
    do i = n, 1, -1
      g(i) = g(i) + alpha * AmalD(i)
    end do

    !----------------------------------------------------------------
    ! compute new beta :                                             
    ! beta = (g(transp) * (a * d)) / (d(transp) * (a * d))           
    !----------------------------------------------------------------
    dividend = ZERO
    do i = n, 1, -1
      dividend = dividend + g(i) * AmalD(i)
    end do

    beta = dividend / divisor

    !----------------------------------------------------------------
    ! compute new d :   d  =  - g + beta * d                         
    !----------------------------------------------------------------
    do i = n, 1, -1
      d(i) = -g(i) + beta * d(i)
    end do

  end do  !k loop

  fehler=0
  return
end




