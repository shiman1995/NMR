!
! @Author: hanyanqiang
! @LastEditors: hanyanqiang
! @Date: 2019-03-13 20:55:13
! @LastEditTime: 2019-03-14 11:19:44
!
! Produce supercell(pdb) by pdb
! By HYQ
program main
  implicit none
  integer::ist,i,j,k,m,n,natom,nres,na,nxl,nyl,nzl
  integer::un,ia,ua,unres
  character(len=4)::molename
  character(len=80)::fline
  real::clata,clatb,clatc,angle1,angle2,angle3
  real::sin1,sin2,sin3,cos1,cos2,cos3,fy,fz,ff
  real,parameter::PI=3.1415926
  integer,allocatable::resnum(:)
  character(len=4),allocatable::atomname(:)
  real,allocatable::coord(:,:),cellcoord(:,:,:,:,:)

  open(11,file='parameter.in')
    read(11,*)
    read(11,*)nxl, nyl, nzl
    read(11,*)
    read(11,*)natom
    read(11,*)
    read(11,*)clata, clatb, clatc
    read(11,*)
    read(11,*)angle1,angle2,angle3
  close(11)

  allocate(coord(natom,3))
  allocate(resnum(natom))
  allocate(atomname(natom))
  allocate(cellcoord(-nxl:nxl,-nyl:nyl,-nzl:nzl,natom,3))

  ! read pdb
  open(5,file='unitcell.pdb')
    do while(.true.)
    read(5,'(A80)',iostat=ist)fline
    if(ist/=0)exit
    if(fline(1:4).eq.'ATOM')then
      read(fline,'(6x,I5)')i
      read(fline,'(12x,a4)')atomname(i)
      read(fline,'(16x,a4)')molename
      read(fline,'(22x,I4)')resnum(i)
      read(fline,'(30x,f8.3)')coord(i,1)
      read(fline,'(38x,f8.3)')coord(i,2)
      read(fline,'(46x,f8.3)')coord(i,3)
    endif
    enddo
  close(5)
 
  write(*,*)molename,atomname(6),atomname(2),coord(1,3),coord(3,2)
  !molecular per cell
  nres=resnum(natom)
  !atoms per molecular
  na=natom/nres
  write(*,*)natom,nres,na

  sin1=sin(angle1*PI/180.0)
  sin2=sin(angle2*PI/180.0)
  sin3=sin(angle3*PI/180.0)
  cos1=cos(angle1*PI/180.0)
  cos2=cos(angle2*PI/180.0)
  cos3=cos(angle3*PI/180.0)

  fy = (cos1-cos2*cos3)/sin3
  fz = sqrt((sin2)**2-fy**2)
  ff = sin3*fz

  ! produce supercell
  do i=-nxl,nxl
    do j=-nyl,nyl
      do k=-nzl,nzl
        do m=1,natom
           cellcoord(i,j,k,m,1)=coord(m,1)+real(i)*clata&
                                +real(j)*clatb*cos3&
                                +real(k)*clatc*cos2
           cellcoord(i,j,k,m,2)=coord(m,2)+real(j)*clatb*sin3&
                                +real(k)*clatc*fy
           cellcoord(i,j,k,m,3)=coord(m,3)+real(k)*clatc*fz
         enddo
       enddo
     enddo
   enddo
  write(*,*) nxl,nyl,nzl

  !num of unitcell
  un=0
  !num of atoms in one cell
  ia=0
  !num of atoms in supper cell
  ua=0
  !num of mole in supper cell
  unres=0

  ! print supercell with pdb 
  open(5,file='newcell.pdb')
  do i=-nxl,nxl
    do j=-nyl,nyl
      do k=-nzl,nzl
        un=un+1
          do n=1,natom
            ua=(un-1)*natom+n
            if(mod(ua,na)==1)then
              unres=unres+1
            endif
            write(5,100)'ATOM',ua,atomname(n),molename,unres,cellcoord(i,j,k,n,1),&
                        cellcoord(i,j,k,n,2),cellcoord(i,j,k,n,3)
 !           if(mod(ua,na)==0)then 
 !             write(5,'(A3)')'TER'
 !           endif
          enddo
      end do
    enddo
  enddo
  close(5)
100   format(A4,I7,1X,A4,A4,I6,4X,3F8.3)

      ! print supercell with inpcrd
      open(5,file='newcell.inpcrd')
      write(5,*)molename
      write(5,*)natom
      write(5,'(6f12.7)')((cellcoord(0,0,0,m,n),n=1,3),m=1,natom)
      write(5,'(3f12.7)')clata,clatb,clatc
      write(5,'(3f12.6)')angle1,angle2,angle3
      close(5)

  deallocate(coord)
  deallocate(resnum)
  deallocate(atomname)
  deallocate(cellcoord)
end

