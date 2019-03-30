
!*************************************************************************
! check1: this subroutine checks the main control data
! Owen p.200
!**************************************************************************

!subroutine check1(ndofn,nelem,ngaus,nmats,nnode,npoin,nstre,ntype,nvfix,ncrit,nalgo,nincs)
subroutine check1( )
! On verifie que tous les parametres d'entrée type nombre de noeuds, etc... sont corrects

    use variables 

    implicit none
    
    integer :: ieror, keror


    do ieror=1,12
      neror(ieror)=0
      !read (*,'()')
    end do

   ! create the diagostic messages

!      if ( npoin <= 0)                    neror(1)=1
!      if ( nelem*nnode < npoin)           neror(2)=1
!      if ( nvfix < 2 .or. nvfix>npoin)    neror(3)=1
!      if ( nincs < 1)                     neror(4)=1
!      if ( ntype < 1 .or. ntype > 3)      neror(5)=1
!      if ( nnode < 4 .or. nnode > 9)      neror(6)=1
!      if ( ndofn < 2 .or. ndofn > 6)      neror(7)=1
!      if ( nmats < 1 .or. nmats > nelem)  neror(8)=1
!      if ( ncrit < 1 .or. ncrit > 6)      neror(9)=1
!      if ( ngaus < 2 .or. ngaus > 3)      neror(10)=1
!      if ( nalgo < 1 .or. nalgo > 4)      neror(11)=1
!      if ( nstre < 3 .or. nstre > 5)      neror(12)=1

    if ( npoin <= 0)                    neror(1)=1
      if ( nelem*nnode < npoin)           neror(2)=1
      if ( nvfix < 2 .or. nvfix>npoin)    neror(3)=1
      if ( nincs < 1)                     neror(4)=1
      if ( ntype < 1 .or. ntype > 3)      neror(5)=1
      if ( nnode < 3 .or. nnode > 9)      neror(6)=1
      if ( ndofn < 2 .or. ndofn > 6)      neror(7)=1
      if ( nmats < 1 .or. nmats > nelem)  neror(8)=1
      if ( ncrit < 1 .or. ncrit > 6)      neror(9)=1
      if ( ngaus < 2 .or. ngaus > 3)      neror(10)=1
      if ( nalgo < 1 .or. nalgo > 10)     neror(11)=1
      if ( nstre < 3 .or. nstre > 5)      neror(12)=1
      
    ! either return or else print the errors

    keror=0 !Indique si il y a un probleme dans les variables d'entrée si keror=1
    do ieror=1,12
      if(neror(ieror) /= 0) then
        keror=1
        write(6,900) ieror
      end if 
    end do

900 format(//31h *** diagnosis by check1,error,i3)


    if (keror /= 0) then
        call echo() ! otherwise ech all the remaining data without further comment
      else
        return
    end if

end


!****************************************************************************
!   This subroutine checks the remainder of the input data Owen p.202
!   Check geometric data, boundary conditions, and material properties
!****************************************************************************

!subroutine check2(coord,iffix,lnods,matno,melem,mfron,mpoin,mtotv,mvfix,ndfro,ndofn,nelem,nmats,nnode,nofix,npoin,nvfix)
 
 subroutine check2( )

use variables

implicit none


integer :: kpoin, drap=0,ipoin,jpoin,idime, ielem
integer :: i, idofn, ieror, inode, ivfix, jvfix, keror, kfron, klast
integer :: kount, kstar, kvfix, kzero, nfron, nloca, nlast
double precision    :: sigma2


! check against two identical nonzero nodal coordinates
do ieror=13,24 !A partir de 13 car c'est la suite de la subroutine check1
    neror(ieror)=0
end do
      
do ielem=1,nelem
    ndfro(ielem)=0
end do


!Verifie qu'il n'y pas des noeuds avec les même coordonnées ie tous les noeuds sont uniques

!Version originale, meilleure sur le principe mais moins dans l'écriture
! do 40 ipoin=2,npoin
!     kpoin=ipoin-1
!     do 30 jpoin=1,kpoin
!         do 20 idime=1,2
!             if(coord(ipoin,idime) /= coord(jpoin,idime)) then
!                 goto 30
!             end if
!         20 continue
!         neror(13)=neror(13)+1
!     30 continue
! 40 continue


do ipoin=2,npoin
    kpoin=ipoin-1
    do jpoin=1,kpoin
        drap=0
        do idime=1,ndim
            if ( coord(ipoin,idime) == coord(jpoin,idime) ) then
                drap=drap+1
            end if
        end do
        if (drap == ndim) then
            neror(13)=neror(13)+1  
        end if   
    end do
end do




!Check the list of element property numbers
do ielem=1,nelem
    if(matno(ielem) <= 0.or.matno(ielem) > nmats) then
        neror(14)=neror(14)+1
    end if
end do



!Check for impossible node numbers
do ielem=1,nelem
    do inode=1,nnode
        if(lnods(ielem,inode) == 0) then
            neror(15)=neror(15)+1
        end if
        if(lnods(ielem,inode) < 0 .or. lnods(ielem,inode) > npoin) then
            neror(16)=neror(16)+1
        end if
    end do
end do


!Check for any repetition of a node number within an element
do ipoin=1,npoin
    kstar=0
    do ielem=1,nelem
        kzero=0
        do inode=1,nnode
            if(lnods(ielem,inode) /= ipoin) then
                cycle
            end if
            
            kzero=kzero+1
            
            if(kzero > 1) neror(17)=neror(17)+1
    ! Seek first, last and intermediate appearance of node ipoin
            if(kstar == 0) then
                kstar=ielem
                ! Calculate increase or decrease in frontwidth at each element stage...
                ndfro(ielem)=ndfro(ielem)+ndofn
            end if
            ! ... and change the signe of the last appearance of each node
            klast=ielem
            nlast=inode
        end do
    end do


    if(klast < nelem) then
        ndfro(klast+1)=ndfro(klast+1)-ndofn
        lnods(klast,nlast)=-ipoin
        cycle
    end if
    
    if(kstar == 0) then
        
        !Check that coordinates for an unused node have not been specified
        
        write(6,900) ipoin
    900 format(/15h check why node,i4,14h never appears)
        neror(18)=neror(18)+1
        sigma2=0.0
        
        do idime=1,2
         sigma2=sigma2+abs(coord(ipoin,idime))
        end do
      
        if(sigma2 /= 0.0) then
            neror(19)=neror(19)+1
        end if
        
        !Check that an unused node number is not restrained node
        
        do ivfix=1,nvfix
            if(nofix(ivfix) == ipoin) then
                neror(20)=neror(20)+1
            end if
        end do
    
    end if
end do



!Calculate the largest front width

nfron=0
kfron=0
do ielem=1,nelem
    nfron=nfron+ndfro(ielem)
    if(nfron > kfron) then
        kfron=nfron
    end if
end do


      
write(6,905) kfron
905 format(//33h maximum frontwidth encountered =,i5)
      
if(kfron > mfron) then
    neror(21)=1
end if



!Continue checking the data for the fixed values

do ivfix=1,nvfix   
      if(nofix(ivfix) <= 0 .or. nofix(ivfix) > npoin) then
          neror(22)=neror(22)+1
      end if
      kount=0
      nloca=(nofix(ivfix)-1)*ndofn

      do idofn=1,ndofn
         nloca=nloca+1
         if(iffix(nloca) > 0) then
             kount=1
         end if
      end do
    

      if(kount == 0) then
          neror(23)=neror(23)+1
        !write(6,*) ivfix
      end if
      kvfix=ivfix-1
      
      do jvfix=1,kvfix
         if(ivfix /= 1 .and. nofix(ivfix) == nofix(jvfix)) then
             neror(24)=neror(24)+1
         end if
      end do
      
end do



keror=0
do ieror=13,24
    if(neror(ieror) == 0) then
        cycle
    end if
    keror=1
    write(6,910) ieror,neror(ieror)
910 format(//31h ***diagnosis by check2, error,i3,6x,18h associated number,i5)
end do


if(keror == 0) then 
!Return all nodal connections numbers to positive values
    do ielem=1,nelem
        do inode=1,nnode
           lnods(ielem,inode)=iabs(lnods(ielem,inode))
        end do
    end do
    
    return
end if

if(keror /= 0) then 
    call echo
end if      
 

do i=1,1
    print *, " "  !Mise en page de la sortie console
    print *, "Subroutine CHECK2 - Status: ERROR"
    print *, " "  !Mise en page de la sortie console
end do


end
!END CHECK2





! ************************************************************************
! Subroutine echo Owen p.201
! if data errors have been detected by subroutines check1 or
! check2,this subroutine reads and writes the remaining data cards
! **************************************************************************

subroutine echo()

    dimension ntitl(80)
    print *, "INFORMATION PROVIDED BY SUBROUTINE ECHO"
    
    write(6,900)
900 format(//50h now follows a listing of post-disaster data cards/)
10  read(5,905) ntitl
905 format(80a1)
    write(6,910) ntitl
910 format(20x,80a1)
    goto 10

end