!*******************************************************************************
!   Subroutine output_python( )
!   This subroutine outputs displacements, reactions, and stresses
!*******************************************************************************

subroutine output_deviator( )

use variables
implicit none

integer :: koutp, ipoin, ngash, ngish, kelgs, istre, istr1, ivfix, idofn, igash, igaus, jgaus, kgaus
integer :: ielem, inode,jdim, idime, lnode, knode, kshap, kgaus_rep
double precision :: xgash, xgish, xgesh, xgosh
double precision :: etasp, exisp,alpha_extrap, ksi, eta
double precision :: s_xx, s_yy, s_xy, s_zz, varj2_triax, eps_a
! (strsg(istr1,kgaus),istr1=1,4)

!We work on first gauss point to begin with
s_xx=strsg(1,1)
s_yy=strsg(2,1)
s_xy=strsg(3,1)
s_zz=strsg(4,1)

!Calculate second deviatoric stress invariant J2'
varj2_triax=s_xy*s_xy+0.5d0*(s_xx*s_xx + s_yy*s_yy + s_zz*s_zz)

eps_a=defsg(1,1)

open(unit=152,file='PostProcessing/Output/Deviatoric/deviatoric.txt')

if (iincs==1) then
	write(152,*) "Increment		eps_a 		varj2 		sqrt(varj2)"
end if
write(152,*) iincs,eps_a*100.d0,varj2_triax,dsqrt(varj2_triax)


end subroutine


