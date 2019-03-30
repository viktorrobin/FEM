

subroutine outputvalidationCC( )

    use variables
implicit none

! integer :: koutp, ipoin, ngash, ngish, kelgs, istre, istr1, ivfix, idofn, igash, igaus, jgaus, kgaus
! integer :: ielem, inode,jdim, idime, lnode, knode, kshap, kgaus_rep,i
double precision :: p_eff_mean,temp,q_devia, depspp

! double precision :: etasp, exisp,alpha_extrap, ksi, eta

integer :: num_gaus

double precision, dimension(4) :: deviateur
! double precision, dimension(npoin) :: yield_node,spe_vol_node,epstnp_node,deviatoric_stress_node,effective_mean_stress_node

! character(len=200) :: filename,part4,part1, numfile_char,npoin_char, nelem_char,x_char, y_char, z_char, nelem_tot_char

!-----------------------------------------------------------------------------------------
!Write cam clay parameters
!-----------------------------------------------------------------------------------------
if (istage==1 .and. iincs==nincs) then
		open(unit=133,file='PostProcessing/Output/Validation_CamClay/para_camclay.txt')
        write(133,*) 'Nkappa     Nlambda      lambda     kappa     M     py      pi     '
        write(133,*) Nkappa,Nlambda,CClambda,CCkappa,csl,props(1,5),strsg(1,1)
        close(133)
end if

!-----------------------------------------------------------------------------------------
!Write axial displacement, specific volume, effective mean stress, and deviatoric stress
!-----------------------------------------------------------------------------------------
if (istage==1 .and. iincs==1) then
        open(unit=132,file='PostProcessing/Output/Validation_CamClay/validation_camclay.txt')
        write(132,*) 'p     q     Eps_a     v'
end if

! we choose gauss point number
    num_gaus=1

! strsg(istr1,kgaus)

depspp=epstnp(num_gaus)-epstnp_pre(num_gaus)

epstnp_pre=epstnp

p_eff_mean=(strsg(1,num_gaus)+strsg(2,num_gaus)+strsg(4,num_gaus))/3.d0 

! print *, (strsg(i,num_gaus),i=1,4)
!Compute deviatoric stresses according to (7.7): s=sigma -(1/3)*I1
deviateur(1)=strsg(1,num_gaus)-p_eff_mean         ! sigma_xx / rr
deviateur(2)=strsg(2,num_gaus)-p_eff_mean         ! sigma_yy / zz
deviateur(3)=strsg(3,num_gaus)               		! sigma_xy / rz
deviateur(4)=strsg(4,num_gaus)-p_eff_mean         ! sigma_ZZ / tt


!Calculate second deviatoric stress invariant J2'
temp=deviateur(3)*deviateur(3)+0.5d0*(deviateur(1)*deviateur(1) + deviateur(2)*deviateur(2) + deviateur(4)*deviateur(4))
q_devia=dsqrt(3.d0*temp)
! print *, p_eff_mean,dsqrt(3.d0*varj22),volspegp(num_gaus)

! if (ntype==3) then
! 	print *,p_eff_mean,q_devia,defsg(2,num_gaus), volspegp(num_gaus)
! end if

	write(132,*) p_eff_mean,q_devia,defsg(2,num_gaus), volspegp(num_gaus),depspp
! end if

! write(1304,*) istage,istage,iincs,iincs
return

end subroutine