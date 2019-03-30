!******************************************************
!
! This file contains subroutines used for plasticity
!
!******************************************************


! **********************************************************************
!   Subroutine invar - Owen p. 239
!   This subroutine evaluates the stress invariants and the current value of
!       the yield function
! ***********************************************************************

subroutine invar(lprop,stemp)

use variables
implicit none

double precision, dimension(4) :: stemp
double precision :: varj1
integer :: lprop

! write(6,*) 'stemp dans invar',stemp
! stop


! First invariant I1=Tr(sigma)
I1=stemp(1)+stemp(2)+stemp(4)
smean=(stemp(1)+stemp(2)+stemp(4))/3.d0

!Compute deviatoric stresses according to (7.7): s=sigma -(1/3)*I1
devia(1)=stemp(1)-smean         ! sigma_xx
devia(2)=stemp(2)-smean         ! sigma_yy
devia(3)=stemp(3)               ! sigma_xy
devia(4)=stemp(4)-smean         ! sigma_ZZ


!On vérifie qu'il n'y a pas d'incohérences dans le tenseur déviatorique
!We check there are no inconstitencies with the deviatoric tensor, i.e. if J1=0 (First deviatoric stress invariant J1')
	varj1=devia(1)+devia(2)+devia(4)

! 	write(6,*) "smean, devia(1),devia(2),devia(3),devia(4)="
! 	write(6,*) smean,devia(1),devia(2),devia(3),devia(4)
	!stop

! 	if (abs(varj1) > 1.d-8) then
    if (abs(varj1) > epsilon) then
    	write(6,*) 'Error in the deviatoric stress tensor: J1 ≠ 0 =', varj1
    	stop
	end if
! End test

!Calculate second deviatoric stress invariant J2'
varj2=devia(3)*devia(3)+0.5d0*(devia(1)*devia(1) + devia(2)*devia(2) + devia(4)*devia(4))
! write(6,*) "varj2=",varj2

!Calculate third deviatoric stress invariant J3'
varj3=devia(4)*(devia(4)*devia(4)-varj2)

!Compute √(J2')
steff=dsqrt(varj2)

if (steff == 0.d0) go to 10
sint3=-3.d0*root3*varj3/(2.d0*varj2*steff)
if(sint3.GT.1.d0) sint3=1.d0
go to 20
10 sint3=0.d0
20 continue
if(sint3.LT.-1.d0) then
    sint3=-1.d0
end if

if(sint3.GT.1.d0) then
    sint3=1.d0
end if
theta=dasin(sint3)/3.d0


! Evaluate yield function according to the criterion
select case(ncrit)
    case(1)     !ncrit=1: Tresca yield criterion
        yield=2.d0*dcos(theta)*steff
        return
    case(2)     !ncrit=2: Von Mises
        yield=root3*steff
!         print *,"yield=",yield
        return
    case(3)     !ncrit=3: Mohr-Coulomb
        phira=props(lprop,7)*degres
        snphi=dsin(phira)
        yield=smean*snphi+steff*(dcos(theta)-dsin(theta)*snphi/root3)
        return
    case(4)     !ncrit=4: Drucker-Prager
        phira=props(lprop,7)*degres
        snphi=dsin(phira)
        yield=6.d0*smean*snphi/(root3*(3.d0-snphi))+steff
        return
    case(5)     !ncrit=5: Cam-Clay
        csl=props(lprop,11)
!         yield=(3.d0*varj2/(smean*csl*csl))+smean
        yield=-(I1/3.d0)-(9.d0*varj2)/(I1*csl*csl)
!         print *,'yield ss=',yield
        return
    !Change  44
   case(6)     !ncrit=6: Lime treated soils
!        q=1.732050808*steff
        csl=props(lprop,11)
        pb=props(lprop,16)
       yield=-(I1/3.d0)-(9.d0*varj2)/((I1+3.d0*pb)*csl*csl)
end select


end

!********************************************************************
!   Subroutine yieldf - Owen p.241
!   This subroutine evaluates the flow vector a defined in (7.74)
!********************************************************************

subroutine yieldf(lprop,kgaus)

use variables

implicit none

integer :: istr1, lprop, kgaus


if(steff == 0.d0) then
    return
end if


frict=props(lprop,7)
! snphi=dsin(frict*degres) ! ?
tanth=dtan(theta)
tant3=dtan(3.d0*theta)
sinth=dsin(theta)
costh=dcos(theta)
cost3=dcos(3.d0*theta)
!root3=1.73205080757

! Calculate vector a1

veca1(1)=1.d0
veca1(2)=1.d0
veca1(3)=0.d0
veca1(4)=1.d0

! Calculate vector a2

do istr1=1,nstr1
    veca2(istr1)=devia(istr1)/(2.d0*steff)
end do

veca2(3)=devia(3)/steff

! Calculate vector a3

veca3(1)=devia(2)*devia(4)+varj2/3.d0
veca3(2)=devia(1)*devia(4)+varj2/3.d0
veca3(3)=-2.d0*devia(3)*devia(4)
veca3(4)=devia(1)*devia(2)-devia(3)*devia(3)+varj2/3.d0

select case(ncrit)
    case(1)     ! ncrit=1: Tresca
        cons1=0.d0
        abthe=dabs(theta*radian)
        if(abthe >= 29.d0) then
            cons2=root3
            cons3=0.d0
        end if

        if(abthe < 29.d0) then
            cons2=2.d0*(costh+sinth*tant3)
            cons3=root3*sinth/(varj2*cost3)
        end if

    case(2)     ! ncrit=2: Von Mises
        cons1=0.d0
        cons2=root3
        cons3=0.d0

    case(3)     ! ncrit=3: Mohr-Coulomb
        cons1=dsin(frict*degres)/3.d0
        abthe=dabs(theta*radian)


        if(abthe.LT.29.d0) go to 30
        cons3=0.d0
        plumi=1.d0
        if (theta.GT.0.d0) plumi=-1.d0
        cons2=0.5*(root3+plumi*cons1*root3)
        go to 40
        30 cons2=costh*((1.d0+tanth*tant3)+cons1*(tant3-tanth)*root3)
        cons3=(root3*sinth+3.d0*cons1*costh)/(2.d0*varj2*cost3)
        go to 40

    case(4)     ! ncrit=4: Drucker-Prager
        snphi=dsin(frict*degres)
        cons1=2.d0*snphi/(root3*(3.d0-snphi))
        cons2=1.d0
        cons3=0.d0

    case(5)     ! ncrit=5: Cam Clay
!
        cons1=(1.d0/9.d0)*(2.d0*I1+3.d0*yieldgp(kgaus))
!         print *,"iincs=",iincs,'kgaus=',kgaus,'yieldgp=',yieldgp(kgaus)
        cons2=(6.d0*steff)/(csl*csl)
        cons3=0.d0
    case(6)  ! ncrit=6: MASS
        cons1=(1.d0/9.d0)*(2.d0*I1+3.d0*(yieldgp(kgaus)+pb))
!         print *,"iincs=",iincs,'kgaus=',kgaus,'yieldgp=',yieldgp(kgaus)
        cons2=(6.d0*steff)/(csl*csl)
        cons3=0.d0

end select

40 continue
do istr1=1,nstr1
    avect(istr1)=cons1*veca1(istr1)+cons2*veca2(istr1)+cons3*veca3(istr1)
end do

return
end


! *********************************************************************
!   Subroutine flowpl - Owen p.243
!   This subroutine evaluates the plastic d vector.
! ***********************************************************************


subroutine flowpl(lprop,kgaus)

use variables

implicit none

integer :: istr1, lprop, kgaus

double precision :: denom, fmul1, fmul2, fmul3, hards, part1, part2, part3, pyIIo,betas


young=props(lprop,1)
poiss=props(lprop,2)
hards=props(lprop,6)


fmul1=young/(1.d0+poiss)

if (ntype == 2 .or. ntype == 3) then
    !Plane strain: ntype = 2 / Axial symmetry: ntype = 3
    fmul2=young*poiss*(avect(1)+avect(2)+avect(4))/((1.d0+poiss)*(1.d0-2.d0*poiss))
    dvect(1)=fmul1*avect(1)+fmul2
    dvect(2)=fmul1*avect(2)+fmul2
    dvect(3)=0.5d0*avect(3)*young/(1.d0+poiss)
    dvect(4)=fmul1*avect(4)+fmul2
else
    !Plane stress: ntype = 1
    fmul3=young*poiss*(avect(1)+avect(2))/(1.d0-poiss*poiss)
    dvect(1)=fmul1*avect(1)+fmul3
    dvect(2)=fmul1*avect(2)+fmul3
    dvect(3)=0.5d0*avect(3)*young/(1.d0+poiss)
    dvect(4)=fmul1*avect(4)+fmul3
end if


! Compute 1/(H'+dvect*avect) for later evaluation of the elasto-plastic matrix Dep
! /!\ matrix D is constant in this version, ie not re-evaluated even if yield occured

!For Cam Clay, hards correspond au parametre A et n'est pas constant -> A CALCULER

if (ncrit/=5) then
    70 denom=hards
end if

if (ncrit==5) then
    denom = ( yieldgp(kgaus)*volspegp(kgaus)*I1*( 3.d0*yieldgp(kgaus)+2.d0*I1) )/( 9.d0*(CCkappa-CClambda) )
end if

!Model MASS
if (ncrit==6) then
    part1 = volspegp(kgaus)*(I1+3.d0*pb)*(2.d0*I1 + 3.d0*(yieldgp(kgaus)+pb))
    part2 = exp(beta*yieldgp(kgaus))*beta*(exp(pyI*beta) + exp(pyII*beta))*(Dec-Dei)
    part3 = (exp(yieldgp(kgaus)*beta)+exp(pyII*beta))*(exp(yieldgp(kgaus)*beta)+exp(pyII*beta))
    denom = part1/( 9.d0*( -(part2/part3) + (CClambda-CCkappa)/yieldgp(kgaus)) )
end if


do istr1=1,nstr1
    denom=denom+avect(istr1)*dvect(istr1)
end do

abeta=1.d0/denom


return
end

