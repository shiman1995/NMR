module comanmr
implicit none
integer,parameter::MAXAT=25000,MAXRES=6000
integer lastprotatom,firstprotres,lastprotres,resno(MAXAT),charge(MAXRES)
integer natom,nres,modum,molatom
character(len=4)::atomname(MAXAT),dlabel(MAXAT)
character(len=3)::residue(MAXAT),residuename(MAXAT)
character(len=2)::element(MAXAT)
real(kind=8)::coord(3,MAXAT)
real(kind=8)::qmcharge(MAXAT),rad(MAXAT)
logical::connect(MAXRES,MAXRES)
logical::atomsign(MAXAT)
character*80 filek
character(len=4) basename
logical::ter(0:MAXRES)
integer selectC3(MAXRES),selectO3(MAXRES)
end module comanmr

Program main
use comanmr
lastprotres=500
molatom=33
basename = 'CITD'
call fileRead
call jconnect
!call fragselect
call fragCreate
end

subroutine fileRead
use comanmr
implicit none
integer i,j,k,l,m
integer nptemp
integer ist
integer ttnumber
character*80 line
character*6 sn(MAXAT)
double precision chargef(MAXRES)

open(10,file='newcell.pdb')
do i=1,MAXRES
    ter(i) = .false.
    chargef(i) = 0.d0
enddo
ter(0) = .false.
lastprotatom = 0
nptemp = 0
i = 0
open(40,file='fitting_charge.txt')
do l = 1,molatom
    read(40,'(64x,f12.6)') qmcharge(l)
enddo
l = 0
m = 0
do k=1,MAXAT
read(10,'(a80)',end = 101) line
    if(line(1:4) .eq. 'ATOM') then
        i = i + 1
        l = l + 1
        if(mod(i,molatom)==1) then
          m = m + 1
          l = 1
        endif
          read(line,100)sn(i),ttnumber,atomname(i),residue(i),resno(i),&
            (coord(j,i),j=1,3),element(i)
100 format(a4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3,23x,a2)
          resno(i) = m
          nptemp=m
          qmcharge(i) = qmcharge(l)
          if(nptemp.gt.MAXRES) then
            write(0,*) 'too many residues: ', nptemp, MAXRES
            stop
          endif
    endif
        residuename(nptemp)=residue(i)
enddo
101 natom = i
nres = max(1,nptemp)
firstprotres=resno(1)
close(10)

do i=1,natom
    if(element(i).eq.'  ') then
        element(i)(1:1) = atomname(i)(2:2)
    endif
    write(*,*) i,atomname(i),resno(i),qmcharge(i)
enddo

do i=firstprotres, nres
    do j=firstprotres, nres
        connect(i,j)=.false.
    enddo
enddo

do i=firstprotres, nres
    connect(i,i)=.true.
enddo
end subroutine fileRead

subroutine jconnect
use comanmr
implicit none
integer i,j,k
double precision dis,nbcut,nhcut
nbcut = 3.5d0
nhcut = 2.5d0
do i=1,natom
    do j=i+1,natom
        if(resno(i).eq.resno(j)) cycle
        dis=dsqrt((coord(1,i)-coord(1,j))**2+ &
            (coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
        if(dis.le.nbcut .and. &
            (element(i)(1:1).ne.'H'.or.element(j)(1:1).ne.'H')) then
            connect(resno(i),resno(j)) = .true.
            connect(resno(j),resno(i)) = .true.
        endif
        if(dis.le.nhcut.and.element(i)(1:1).eq.'H'.and.  &
            element(j)(1:1).eq.'H') then
            connect(resno(i),resno(j)) = .true.
            connect(resno(j),resno(i)) = .true.
        endif
    enddo
enddo

end subroutine jconnect

!subroutine fragselect
!use comanmr
!implicit none
!integer i
!
!do i=1,natom
!    if(resno(i).le.lastprotres) then
!        if(atomname(i).eq." C3'") selectC3(resno(i))=i
!        if(atomname(i).eq." O3'") selectO3(resno(i))=i
!    endif
!enddo
!end subroutine fragselect

subroutine fragCreate
use comanmr
implicit none
integer i,j,k,iqm,cfrag,kk,iqmca,jj,kkatom,n1,n2,iitemp
integer kstart,kfinal,ktemp,iresm1,kbstart,kbfinal
double precision x,y,z,a,b,c,d
do k=firstprotres,lastprotres
    iqm=0
    modum=0
    filek = basename(1:4)//char(48+k/100) &
        //char(48+(k-k/100*100)/10)//char(48+(k-k/10*10))
    open(30,file=filek(1:7)//'.com')
    open(31,file=filek(1:7)//'.pqr')
    write(30,'(a)') '%mem=16GB'
    write(30,'(a)') '%nprocshared=8'
    write(30,'(a)') &
        '#B3LYP/6-31g** charge nmr(printeigenvectors) nosymm &
        integral(grid=ultrafine)'
    write(30,*)
    write(30,'(a,i4)') ' AF-NMR fragment for residue ',k
    write(30,*)
    
    do i=1,natom
        atomsign(i)=.false.
    enddo

    cfrag=0
    do ktemp=firstprotres,nres
        if(connect(k,ktemp))then
            cfrag = cfrag + charge(ktemp)
        endif
     enddo
    write(30,'(1x,i3,2x,i2)') cfrag,1
!
! core region
!
    do kk=1,natom
        if(resno(kk).eq.k) then
        iqm = iqm + 1
        atomsign(kk)=.true.
        call addatom(kk,iqm)
        endif
    enddo
!
! buffer region
!
    do ktemp=firstprotres,nres
        if((ktemp.ne.k).and.(connect(k,ktemp)))then
            do jj=1,natom
                if(resno(jj).eq.ktemp) then
                iqm = iqm + 1
                atomsign(jj)=.true.
                call addatom(jj,iqm)
                endif
            enddo
        endif
    enddo
!
! add H
!
!           if(ktemp.ne.firstprotres.and.ktemp.le.lastprotres &
!               .and. .not.ter(ktemp-1)) then
!               if(.not.atomsign(kbstart-1))then
!                   n1=kbstart
!                   n2=selectC3(ktemp-1)
!                   call xyzchangeO(coord(1,n2),coord(2,n2),coord(3,n2), &
!                       coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
!                   iqm = iqm+1
!                   call addH(iqm,x,y,z)
!               endif
!           endif
!           if(ktemp.lt.lastprotres.and..not.ter(ktemp)) then
!               if(.not.atomsign(kbfinal+3))then
!                   n1=selectC3(ktemp)
!                   n2=selectO3(ktemp)
!                   call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
!                       coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
!                   iqm = iqm+1
!                   call addH(iqm,x,y,z)
!               endif
!           endif
!       endif
!   enddo
    close(31)
    write(30,*)
!
! mm region
!
    do kk=1,natom
        if(.not.atomsign(kk))then
            write(30,'(3f10.4,2x,f12.8)') (coord(j,kk),j=1,3),qmcharge(kk)
        endif
    enddo
!   open(23,file='qvalue.dat')
!   read(23,*)
!       do iitemp=999999
!           read(23,*,end=60)a,b,c,d
!           write(30,'3f15.6,2x,f12.8)')a,b,c,d
!       enddo
!   60 continue
!   close(23)
    write(30,*)
    close(30)
enddo
end subroutine fragCreate

subroutine xyzchange(xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew)
implicit none
double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew,grad

grad=dsqrt(1.09d0**2/((xold-xzero)**2+(yold-yzero)**2 &
    +(zold-zzero)**2))
xnew=xzero+grad*(xold-xzero)
ynew=yzero+grad*(yold-yzero)
znew=zzero+grad*(zold-zzero)
return   
end subroutine xyzchange

subroutine xyzchangeO(xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew)
implicit none
double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew,grad

grad=dsqrt(0.96d0**2/((xold-xzero)**2+(yold-yzero)**2 &
    +(zold-zzero)**2))
xnew=xzero+grad*(xold-xzero)
ynew=yzero+grad*(yold-yzero)
znew=zzero+grad*(zold-zzero)
return   
end subroutine xyzchangeO

subroutine addH(iqm,x,y,z)
    use comanmr
    implicit none
    integer iqm
    real(kind=8) x,y,z
    dlabel(iqm)(1:1) = 'H'
    if(iqm.lt.10) then
        write(dlabel(iqm)(2:2),'(i1)') iqm
        dlabel(iqm)(3:4) = '  '
    elseif(iqm.lt.100) then
        write(dlabel(iqm)(2:3),'(i2)') iqm
        dlabel(iqm)(4:4) = ' '
    elseif(iqm.lt.1000) then
        write(dlabel(iqm)(2:4),'(i3)') iqm
    else
        write(0,*) 'Error: more than 10000 atoms'
        stop
    endif
    write(30,'(a,2x,3f12.5)') dlabel(iqm)(1:1),x,y,z
    modum = modum+1
    if(modum.lt.10) then
        write(31,'(a,i5,2x,a,i1,a,3f8.3,f8.4,f8.3)') 'ATOM  ', &
            iqm,'H',modum,'  MOD  9999    ',x,y,z,0.0,1.2
    else
        write(31,'(a,i5,2x,a,i2,a,3f8.3,f8.4,f8.3)') 'ATOM  ', &
            iqm,'H',modum,'  MOD  9999    ',x,y,z,0.0,1.2
    endif
    return 
end subroutine addH

!subroutine get_atom_range(jstart,jfinal,jtemp,jselectO3,jter,lastpro)
!implicit none
!integer jstart,jfinal,jtemp,jselectO3(6000),lastpro
!logical jter(0:6000)
!    if(jtemp.le.lastpro) then
!        jstart=jselectC3(jtemp-1)
!        jfinal=jselectO3(jtemp+1)-1
!        if(jter(jtemp))then
!            jfinal=jselectO3(jtemp)+1
!        endif
!        if(jtemp==1)then
!            jstart=1
!        endif
!    else
!        jstart=jselectO3(jtemp)
!        jfinal=jselectO3(jtemp+1)-1
!    endif
!
!end subroutine get_atom_range

subroutine addatom(jj,jiqm)
use comanmr
implicit none
integer jj,jiqm,j
    if(jiqm.lt.10)then
        dlabel(jiqm)(1:1) = element(jj)(1:1)
        write(dlabel(jiqm)(2:2),'(i1)') jiqm
        dlabel(jiqm)(3:4) = '  '
    elseif(jiqm.lt.100) then
        dlabel(jiqm)(1:1) = element(jj)(1:1)
        write(dlabel(jiqm)(2:3),'(i2)') jiqm
        dlabel(jiqm)(4:4) = ' '
    elseif(jiqm.lt.1000) then
        dlabel(jiqm)(1:1) = element(jj)(1:1)
        write(dlabel(jiqm)(2:4),'(i3)') jiqm
    else
        write(0,*) 'Error: more than 10000 atoms'
        stop
    endif
    write(30,'(a,2x,3f12.5)') dlabel(jiqm)(1:1),(coord(j,jj),j=1,3)
    write(31,'(a,i5,1x,a4,1x,a3,i6,4x,3f8.3,f8.4,f8.3)') 'ATOM  ', &
        jj,atomname(jj),residue(jj),resno(jj),(coord(j,jj),j=1,3), &
        qmcharge(jj),rad(jj)
end subroutine addatom

