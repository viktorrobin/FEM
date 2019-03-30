program cong_grad

  double precision, dimension(:,:), allocatable :: A
  double precision, dimension(:), allocatable :: B, X

  integer :: i, j, n

  interface
    subroutine cg_method(A,B,X,n)
        double precision, dimension(:,:), allocatable :: A
        double precision, dimension(:), allocatable :: B, X    
    end subroutine
  end interface

    n=3
    allocate(A(n,n))
    allocate(B(n))
    allocate(X(n))

do i=1,n
    X(i)=0.d0
end do

!do i=1,3
!write(6,*) xx(i)
!end do

! A(0,0)=110.d0
! A(0,1)=-30.d0
! A(0,2)=-30.d0
! A(1,1)=110.d0
! A(1,2)=50.d0
! A(2,2)=110.d0

! A(1,0)=a(0,1)
! A(2,0)=a(0,2)
! A(2,1)=a(1,2)

! B(0)=0.d0
! B(1)=0.d0
! B(2)=0.7d0


A(1,1)=110.d0
A(1,2)=-30.d0
A(1,3)=-30.d0
A(2,2)=110.d0
A(2,3)=50.d0
A(3,3)=110.d0

A(2,1)=a(0,1)
A(3,1)=a(0,2)
A(3,2)=a(1,2)

B(1)=0.d0
B(2)=0.d0
B(3)=0.7d0


do i=1,3
    write(6,*) (A(i,j),j=1,3)
end do

do i=1,3
    write(6,*) (B(i))
end do

call cg_method(A,B,X,n)

print *,X
end program


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
Subroutine cg_method (A,B,X,nvlib)!  &     ! Conjugate Gradient Method
                       !n, &     ! Size of the linear system
                       !a, &     ! System matrix
                       !y, &     ! right hand side
                       !x, &     ! solution vector
                       !fehler & ! error code
                     !)
  double precision, dimension(:,:), allocatable :: A
  double precision, dimension(:), allocatable :: B, X
! original name: cg_verfahren()
integer, parameter :: SIZE=24
real*8, parameter :: ZERO=0.d0, MACH_EPS=2.d-10
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
double precision d(nvlib), &   ! (0..n-1) auxiliary vectors d and g
         g(nvlib), &
         AmalD(nvlib)  ! (0..n-1) auxiliary vector A * d

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
  integer k, i, j,n      ! loop variables

  if (n < 2) then      ! invalid parameter?
    fehler=1
    return
  end if
n=nvlib


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
      print *,'divisor=0'
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
