
!******************************************************************
!   PROGRAM FOR THE ELASTO PLASTIC ANALYSIS OF PLANE STRESSES,
!   PLAN STRAN,AND AXISYMETRIC SOLIDS (FOR SATURATED SOILS)
!******************************************************************

!   Based on Finite Elements in Plasticity: "Theory and Practice by Owen & Hinton"
!   Fortran 90 version by Victor ROBIN - 22.09.12, Nancy

!     Number of degrees of freedom per node: 2
!     1-X-DISPLACENETS
!     2-Y-DISPLACEMENTS

!   Version 2 for ELASTIC & PLASTIC behaviour: Von Mises, Mohr-Coulomb, Tresca, Drucker-Prager
!   Started: 10.08.2013, Exeter
!   Finished: 18.08.2013
!   By: Victor ROBIN
!******************************************************************

include 'workspace/Module.f90' !Compile global variables

program elast

    use variables     ! Use gobal variables

    implicit none

    ! Declare variables
    integer :: i, j, loca, itotv, rep_a, rep_c
    double precision :: prev_res

    include 'subroutines/interface.f90'

    real, dimension(2) :: tarray
    real :: e
    INTEGER :: clock_start,clock_end,clock_rate, ielapsed,ielapsed1,ielapsed2,ielapsed3,ielapsed4
    double precision:: elapsed_time,cstart,cend,crate, ave1,ave2,ave3,ave4
    integer :: info,n

    ! To deactivate the ouput in terminal (time optimization), uncomment this line:
!     open(unit=6,file='temp/output_terminal.txt')


    CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
    CALL SYSTEM_CLOCK(COUNT=clock_start) ! Start timing

    e=ETIME(tarray)

    ! Timer
    call cpu_time(start)

    !***************************************
    ! Version without GID:

    ! Preset variables associated with dynamic dimensioning
!      call dimen( )

!     ! Call the subroutine which reads most of the problem data
!      call input( )

!     ! Call the subroutine which calculates the consistent load vector for each element after reading the relevant input data
!      call loadps( )
    !***************************************

    call input_data_GiD( )
    call Read_Bound_Cond_Gid( )

    ! Initialise certain arrays
     call zero( )

    ! Output initial state
    iincs=0
    count_file=0

    call output_increm( )

    do istage=1,nstage

        ! In case there are more than one stage, read the boundaries conditions for the second stage
        if (istage == 2) then
            call Read_Bound_Cond_Gid_S2()
        end if

        ! Loop over each increment
        do iincs=1,nincs

            count_file=count_file+1
            ! Start Timer
            call cpu_time(start_inc)

            ! Read data for current increment
            call increm( )

            ! Evaluation of the element stiffness matrices Ke
            call stiffp( )

            ! Assembly global stiffness matrix [K] and global load vector
            call assemb( )

    !         If required, correct global stiffness matrix [K] and global load vector to account for prescribed displacements
    !         if (drap_pres==1) then
    !             call greduc( )
    !         end if

            ! Reduction of the global stiffness matrix and global load vector using boundary conditions
            call reduc( )


    !         Initialization of the residual vector
            do i=1,size(rgload,1)
                vec_residu(i)=rgload(i)
            end do

            iiter=0
            norme=epsilon+1
            prev_res=10000

            gload_backup=gload

    !     Display Reduced Stiffness Matrix
    !     call disp_rgstif( )
    !     stop

            ! Start iteration process to solve [K].{Un}={F}
            do while (norme > epsilon)


                ! Counter for number of iterations
                iiter=iiter+1

                ! If prescribed displacements, modify load vectors for first iteration
                if (drap_pres ==1 .and. iiter==1) then

    !                ! Add prescribe displacements to load vectors
                    call prescdisp()

                    ! Solve linear system A.X=B using conjugate gradient: [K].{Un}={F}
                        n=size(rgstif,1)
                        rgstif_lapack=rgstif
                        rasdis_lapack=vec_residu_pd

                        CALL DGESV(n,1,rgstif_lapack,n,ipiv,rasdis_lapack,n,info)
    !                     CALL dposv('Upper',n,1,rgstif_lapack,n,rasdis_lapack,n,info)


                        rasdis=rasdis_lapack

                else
                    ! No prescribed displacements:
                    ! Solve linear system A.X=B using conjugate gradient: [K].{Un}={F}
                    n=size(rgstif,1)

                        rgstif_lapack=rgstif
                        rasdis_lapack=vec_residu

                        CALL DGESV(n,1,rgstif_lapack,n,ipiv,rasdis_lapack,n,info)

                        rasdis=rasdis_lapack

                end if

                  asdis=0.d0

    !             Inject displacement values in asdis(ntotv) (global vector)
                do i=1,nvlib
                   loca=gcoord(i)
                   asdis(loca)=rasdis(i)
                end do


    !             Only for the first iteration, add prescribed displacements
                if (drap_pres==1 .and. iiter==1) then
                    do i=1,ntotv
                        if (iffix(i)==1 .and. fixed(i)/=0.d0) then
                            asdis(i)=fixed(i)
                        end if
                    end do
                end if

                ! Add displacement to previous total displacements values
                do itotv=1,ntotv
                        tdisp(itotv)=tdisp(itotv)+asdis(itotv)
                end do

                ! Calculate strain tensor, stress tensor, and global load forces
                call residu( )

                ! Check for convergence using load forces vector previously calculated
                call conver( )


            end do

            !Calculate average
            ave1=0.d0
            ave2=0.d0
            ave3=0.d0

            ! output results Owen & Hilton Style
            call output( )

            ! output results Owen & Hilton Style
            call output_increm( )

            call output_deviator( )

        end do !nincs

    end do !nstage

    ! Output results of this increment in vtk files
    print *,' '
    print *,'-----------------------------------'
    print *,'Exporting results...'
    print *,'-----------------------------------'
    print *,' '
    call output_python( )


    close(5)
    close(1)
    close(55)
    close(56)
    close(132)

    print *, " "  !Mise en page de la sortie console
    print *, "FINITE ELEMENT ANALYSIS COMPLETED"

    e=ETIME(tarray)
                  print *, 'elapsed:',e
                  print *, 'user:',tarray(1)
                  print *, 'sys:',tarray(2)


    print *, ''
    print *,'Total time execution:'

    CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing
       ! Calculate the elapsed time in seconds:
       cstart=clock_start
       cend=clock_end
       crate=clock_rate
    !    elapsed_time=((clock_end-clock_start)/clock_rate)
        elapsed_time=(cend-cstart)/crate

    if (elapsed_time < 60) then
        !Seconds

        write(6,222) elapsed_time
        222 format(f5.2,1x,3hsec)

    else if (elapsed_time > 60 .and. elapsed_time < 3600) then
        !Min - Sec
        ielapsed=elapsed_time/60
        ielapsed2=(elapsed_time/60 - ielapsed)*60

        write(6,223) ielapsed, ielapsed2
        223 format(i2,1x,3hmin,1x,i2,1x,3hsec)

    else if (elapsed_time >= 3600 .and. elapsed_time < 86400) then
        !Hours-Min-Seconds
        ielapsed3=elapsed_time/3600
        ielapsed=(elapsed_time/3600 - ielapsed3)*60
        ielapsed2=(elapsed_time - ielapsed3*3600-ielapsed*60)

        write(6,224) ielapsed3,ielapsed, ielapsed2
        224 format(i2,1x,5hhours,1x,i2,1x,3hmin,1x,i2,1x,3hsec)

    else if (elapsed_time >= 86400) then
        !Days-Hours-Min-Seconds
        ielapsed1=elapsed_time/86400
        ielapsed2=elapsed_time/3600 - ielapsed1*24
        ielapsed3=elapsed_time/60 - ielapsed1*1440-ielapsed2*60
        ielapsed4=elapsed_time - ielapsed1*86400-ielapsed2*60*60-ielapsed3*60
        write(6,225) ielapsed1,ielapsed2, ielapsed3,ielapsed4
        225 format(i2,1x,4hDays,1x,i2,1x,5hhours,1x,i2,1x,3hmin,1x,i2,1x,3hsec)


    end if

end program elast

! Load subroutines files
include 'subroutines/geometry.f90'              ! Import mesh, geometry (nnode, convectivity table), and boundary conditions
include 'subroutines/input_examples.f90'        ! Different simulations
include 'subroutines/FEM.f90'                   ! Subroutines inherent to FEM analysis: shape functions, element stiffness matrix, assembling, etc...
include 'subroutines/numerical_analysis.f90'    ! Classic numerical methods: conjugate gradient, jacobian matrix, etc...
include 'subroutines/control.f90'               ! Procedures to check geometry is consistent
include 'subroutines/output.f90'                ! Subroutines to print output on the terminal
include 'subroutines/output_vtk.f90'            ! Subroutines to generate .vtk file for Paraview
include 'subroutines/plasticity.f90'
include 'subroutines/GaussJordan2.f90'
include 'subroutines/output_python.f90'
include 'subroutines/output_deviator.f90'
include 'subroutines/input_bound_cond_GiD.f90'
include 'subroutines/input_data_GiD.f90'
include 'subroutines/validation.f90'


