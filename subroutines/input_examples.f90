subroutine initial_cond(num_example)
	use variables
	implicit none

	integer :: num_example

	if (num_example==0) then

		!*** EXEMPLE OWEN
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=12               !Maximum of ELEMents
		npoin=51                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=12                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE OWEN

		open(5,file='Examples/Example_2D_p542.txt')		


	else if (num_example==1) then  

		!*** EXEMPLE EXEMPLE 2 elements carres
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=2               !Maximum of ELEMents
		npoin=6                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=2                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 2 elements carres

		open(5,file='Examples/Example_two_elem.4d.txt')
	

	else if (num_example==2) then

		!*** EXEMPLE EXEMPLE 3  1 element Q4
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/Example_one_elem_4d.txt')
		    
		
	else if (num_example==3) then
		! *** EXEMPLE EMAD ENSM 2T3
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=3                 !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=1                 !?? Maximum number of material sets??
		mpoin=4               !Maximum permissible number of nodal points in the program
		nnode=3                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		                           ! =4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=2                !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE EMAD ENSM
		open(5,file='Examples/Example_two_elem.txt')
	
else if (num_example==4) then

		!*** EXEMPLE EXEMPLE 3  1 element carres Q8
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=8                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/1elemQ8.txt')

else if (num_example==5) then

		!*** EXEMPLE EXEMPLE 3  4 elements Q8 POINT LOAD
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=21                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/4elemQ8.txt')

else if (num_example==6) then

		!*** EXEMPLE EXEMPLE 3  4 elements Q8 DISTRIBUTED LOAD
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=21                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/4elemQ8_DistLoading.txt')

else if (num_example==7) then

		!*** EXEMPLE EXEMPLE 3  4 elements Q8 Deplacements imposes
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=21                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/4elemQ8_DispImp.txt')

else if (num_example==8) then

		!*** EXEMPLE EXEMPLE 3  1 elements Q4 Deplacements imposes
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/1elemQ4_DispImp.txt')

else if (num_example==9) then

		!*** EXEMPLE EXEMPLE 3  4 elements Q8 Distributed loading en axisymetrie
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=21                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=4                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure

			open(5,file='Examples/4elemQ8_DistLoading_Axisym.txt')

else if (num_example==10) then

		!*** EXEMPLE EXEMPLE 3  1 elements Q8 Point Load en axisymetrie
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=8                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=8                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/1elemQ8_Axis.txt')

	else if (num_example==11) then

		!*** EXEMPLE EXEMPLE 3  1 element Q4 Axisymmetric case Point Load.
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/1Q4_Axi.txt')

else if (num_example==12) then
		! *** EXEMPLE EMAD ENSM 2T3
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=3                 !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=1                 !?? Maximum number of material sets??
		mpoin=4               !Maximum permissible number of nodal points in the program
		nnode=3                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		                           ! =4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode


		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=2                !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE EMAD ENSM
		open(5,file='Examples/Example_two_elem_dispImp.txt')


else if (num_example==13) then

		!*** EXEMPLE EXEMPLE 3  1 element Q4 Plane strain Plasticite parfaite
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/1Q4_Plast.txt')

	else if (num_example==14) then

		!*** EXEMPLE EXEMPLE 3  1 element Q4 Plane strain Plasticite parfaite
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/1Q4_Plast_DispImp.txt')

	else if (num_example==15) then

		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=3                 !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=1                 !?? Maximum number of material sets??
		mpoin=4               !Maximum permissible number of nodal points in the program
		nnode=3                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		                           ! =4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=2                !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure

			open(5,file='Examples/2T3_DispImp.txt')

else if (num_example==16) then

		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=4                 !Maximum of ELEMents
		npoin=9                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=1                 !?? Maximum number of material sets??
		mpoin=9               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		                           ! =4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=4                !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure

			open(5,file='Examples/4Q4_Plast.txt')
else if (num_example==17) then

		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=4                 !Maximum of ELEMents
		npoin=9                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=1                 !?? Maximum number of material sets??
		mpoin=9               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=8                 !Number of material parameters required to define characteristics of a material completely
		                           ! =4 for elasto-plastic problem; =2 for other applications
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=4                !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure

			open(5,file='Examples/4Q4_Plast_PrescDisp.txt')

else if (num_example==18) then

		!*** EXEMPLE EXEMPLE 3  1 element Q4 Axi-symmetry Plasticite Modified Cam Clay
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=11                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications; =11 for Cam Clay
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/1Q4_CamClay.txt')

else if (num_example==19) then

		!*** EXEMPLE EXEMPLE 3  1 element Q4 Axi-symmetry Plasticite Modified Cam Clay
		mbufa=10                !Variables pour buffer; pas d'interet au XXIe siecle

		melem=1               !Maximum of ELEMents
		npoin=4                 !Total number of nodal points in the structure/problem
		ndofn=2                 !Number of Degrees Of Freedom per Node. =2: X et Y
		mfron=ndofn*npoin       !Maximum front width
		mmats=5                 !?? Maximum number of material sets??
		mpoin=npoin               !Maximum permissible number of nodal points in the program
		nnode=4                 !Number of NODes par Element

		mstif=(mfron*mfron-mfron)/2.0 +mfron

		mgaus=melem*nnode           !Utilisée dans problemes transitoire cf. Owen p.415
		mtotg=melem*nnode           !Egal a mgaus; parametre original dans Owen > subroutine DIMEN. Utilite: dimension de certaines matrices

		mtotv=ndofn*npoin       !Maximum total number of degrees of freedom =16
		mvfix=npoin				!Maximum number of prescribed boundary nodes


		mevab=nnode*ndofn           !Maxium of Element VAriaBles per element

		nprop=16                 !Number of material parameters required to define characteristics of a material completely
		!=4 for elasto-plastic problem; =2 for other applications; =11 for Cam Clay
		nevab=ndofn*nnode

		mdouble=2*nnode


		ndim=2                  !nombre de dimensions du problème: ici 2 => ndim=2
		nelem=1                 !Total number of elements in the structure. =1 pour le moment

		nmats=1                 !Number of different MATerialS in the structure
		! *** FIN EXEMPLE 3 1 element carres

			open(5,file='Examples/LimeTreatedSoils.txt')



end if



	return

end subroutine