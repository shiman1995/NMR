      module comanmr
      implicit none
      integer,parameter::MAXAT=250000,MAXRES=6000
      integer lastprotatom,firstprotres,lastprotres,resno(MAXAT)
      integer natom,nres,modum,molatom,n_a,i,j,n,k,g1(40),g2(40),g3(40)
	integer ng1,ng2,ng3,m
      character(len=4)::atomname(MAXAT),dlabel(MAXAT)
      character(len=3)::residue(MAXAT),residuename(MAXAT)
      character(len=2)::element(MAXAT)
      real(kind=8)::coord(3,MAXRES,100)
      real(kind=8)::qmcharge(MAXAT),rad(MAXAT),charge(MAXRES)
      logical::connect1(3,MAXRES,MAXRES),connect2(MAXRES,MAXRES)
      logical::atomsign(MAXAT),connect3(MAXRES,MAXRES)
      logical::connect(MAXRES,MAXRES)
      character*80 filek
      character(len=4) basename
      logical::ter(0:MAXRES)
      integer selectC3(MAXRES),selectO3(MAXRES)
      end module comanmr


      Program main
      use comanmr
	real*8 nbcut , nhcut
      lastprotres=500
      molatom=33
      basename = 'CITD'
      nbcut = 3.5
      nhcut = 2.5
	  
      call fileRead

      open(11,file='lgroup')

      read(11,*) ng1
      read(11,*) ng2
      read(11,*) ng3 
  
      do i=1,ng1
         read(11,*) g1(i)
      end do
      do i=1,ng2
         read(11,*) g2(i)
      end do     
      do i=1,ng3
         read(11,*) g3(i)
      end do
    	     
      close(11)  
	
	do i =  1 , lastprotres

	   do j = 1 , lastprotres

            if(i/=j)then		  
		    
			do n = 1 , molatom

                 do k = 1 , ng1

        dis=dsqrt((coord(1,j,n)-coord(1,i,g1(k)))**2+(coord(2,j,n) 
     &-coord(2,i,g1(k)))**2+(coord(3,j,n)-coord(3,i,g1(k)))**2)

        if(dis.le.nhcut.and.element(n)(1:1).eq.'H'.and.  
     &       element(g1(k))(1:1).eq.'H') then
            connect1(1,i,j) = .true.
	      goto 101
	  elseif(dis.le.nbcut.and.(element(n)(1:1)/='H'.or.  
     &  element(g1(k))(1:1)/='H'))then
            connect1(1,i,j) = .true.
	      goto 101
        endif	               
		      
			 enddo	  
		              
101	         do k = 1 , ng2  
        dis=dsqrt((coord(1,j,n)-coord(1,i,g2(k)))**2+(coord(2,j,n) 
     &-coord(2,i,g2(k)))**2+(coord(3,j,n)-coord(3,i,g2(k)))**2)

        if(dis.le.nhcut.and.element(n)(1:1).eq.'H'.and.  
     &element(g2(k))(1:1).eq.'H') then       
            connect1(2,i,j) = .true.
	      goto 102
	  elseif(dis.le.nbcut.and.(element(n)(1:1)/='H'.or.  
     &  element(g2(k))(1:1)/='H'))then
            connect1(2,i,j) = .true.
	      goto 102
        endif	               
		      
			 enddo	
			 
102	         do k = 1 , ng3

        dis=dsqrt((coord(1,j,n)-coord(1,i,g3(k)))**2+(coord(2,j,n) 
     &-coord(2,i,g3(k)))**2+(coord(3,j,n)-coord(3,i,g3(k)))**2)

        if(dis.le.nhcut.and.element(n)(1:1).eq.'H'.and.  
     &   element(g3(k))(1:1).eq.'H') then
            connect1(3,i,j) = .true.
	      goto 103
	  elseif(dis.le.nbcut.and.(element(n)(1:1)/='H'.or.  
     &  element(g3(k))(1:1)/='H'))then
            connect1(3,i,j) = .true.
	      goto 103
        endif	               
		      
			 enddo	 
103      continue
         enddo

           endif
	
	   enddo

	enddo
      call fragCreate
      end	
			 			  	
      subroutine fileRead
      use comanmr
      implicit none
      integer l
      integer nptemp
      integer ist
      integer ttnumber
      character*80 line
      character*6 sn(MAXAT)
	real*8 rad2,rad1
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
      m = 1
      n_a = 0

      do k=1,MAXAT
      read(10,'(a80)',end = 101) line

          if(line(1:4) .eq. 'ATOM') then

              i = i + 1
              l = l + 1
      		n_a = n_a + 1

              if(i>1.and.mod(i,molatom)==1) then
                m = m + 1
                l = 1
      		  n_a = 1
              endif
            
                read(line,100)sn(i),ttnumber,atomname(i),residue(i),
     & resno(i),(coord(j,m,n_a),j=1,3),element(n_a)
  100 format(a4,1x,i6,1x,a4,1x,a3,2x,i4,4x,3f8.3,23x,a2)
                resno(i) = m
                nptemp = m
                qmcharge(i) = qmcharge(l)
                if(nptemp.gt.MAXRES) then
                  write(0,*) 'too many residues: ', nptemp, MAXRES
                  stop
                endif
              residuename(nptemp)=residue(i)

          endif  
 
      enddo
     
  101 natom = i
      nres = max(1,nptemp)
      firstprotres=1
	if(lastprotres/=nres)write(*,*)'wrong number of residues: ',
     &lastprotres , 'true number of residues: ',nres
	if(n_a/=molatom)write(*,*)'wrong number of molatom: ',
     &molatom , 'true number of residues: ',n_a
      close(10)
  
      do i=1,natom
          if(element(i).eq.'  ') then
              element(i)(1:1) = atomname(i)(2:2)
          endif
      enddo
      do m = 1 , 3
      do i=1, lastprotres
          do j=1, lastprotres
              connect1(m,i,j)=.false.
          enddo
      enddo
      enddo

      do m = 1 , 3
      do i=1, lastprotres	
              connect1(m,i,i)=.true.
      enddo
      enddo

      do i=1, lastprotres
          connect(i,i)=.true.
      enddo

      end subroutine fileRead	       	      


      subroutine fragCreate
      use comanmr
      implicit none
      integer iqm,cfrag,kk,iqmca,jj,kkatom,n1,n2,iitemp
      integer kstart,kfinal,ktemp,iresm1,kbstart,kbfinal
	integer nlowatom
      double precision x,y,z,a,b,c,d
      
      do m=1,3
	 
      do k=1,lastprotres
          iqm=0
          modum=0
           filek = basename(1:4)//char(48+k/1000)//
     &		char(48+(k-k/1000*1000)/100)    
     &             //char(48+(k-k/100*100)/10)//char(48+(k-k/10*10))
          open(30,file=filek(1:8)//'_'//char(48+(m))//'.com')
          open(31,file=filek(1:8)//'_'//char(48+(m))//'.pqr')
          write(30,'(a)') '%mem=10GB'
          write(30,'(a)') '%nprocshared=4'
          write(30,'(a)') 
     &         '#B3LYP/6-31g** charge nmr(printeigenvectors) nosymm 
     &         integral(grid=ultrafine)'
          write(30,*)
          write(30,'(a,i4)') ' AF-NMR fragment for residue ',k
          write(30,*)
          
          do i=1,natom
              atomsign(i)=.false.
          enddo
      
          cfrag=0

          write(30,'(1x,i3,2x,i2)') cfrag,1
      !
      ! core region
      !
          do kk=1,natom
              if(resno(kk).eq.k) then
              iqm = iqm + 1
              atomsign(kk)=.true.
              call addatom(iqm,iqm,k)
              endif
          enddo
      !
      ! buffer region
      !
          do ktemp=1,lastprotres
              if(ktemp.ne.k)then
              if(connect1(m,k,ktemp))then
                  do jj=1,molatom
                      iqm = iqm + 1
                      call addatom(jj,iqm,ktemp)
                  enddo
              endif
	        endif
          enddo

          close(31)
          write(30,*)

          do ktemp=1,lastprotres
              if(.not.(connect1(m,k,ktemp)))then	
                  do jj=1,molatom		       
       write(30,'(3f10.4,2x,f12.8)') (coord(j,ktemp,jj),j=1,3)
     &,qmcharge(jj)
                  enddo
	        endif
	    enddo
      !
      ! mm region
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
	enddo

      end subroutine fragCreate

      subroutine xyzchange(xold,yold,zold,xzero,yzero,zzero,xnew,
     &ynew,znew)

        implicit none
        double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,
     &znew,grad

         grad=dsqrt(1.09d0**2/((xold-xzero)**2+(yold-yzero)**2   
     &  +(zold-zzero)**2))
         xnew=xzero+grad*(xold-xzero)
         ynew=yzero+grad*(yold-yzero)
         znew=zzero+grad*(zold-zzero)
         return
      end subroutine xyzchange

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
          write(31,'(a,i5,2x,a,i1,a,3f8.3,f8.4,f8.3)') 'ATOM  ', 
     & iqm,'H',modum,'  MOD  9999    ',x,y,z,0.0,1.2
      else
          write(31,'(a,i5,2x,a,i2,a,3f8.3,f8.4,f8.3)') 'ATOM  ', 
     & iqm,'H',modum,'  MOD  9999    ',x,y,z,0.0,1.2
      endif
      return 
      end subroutine addH

      subroutine addatom(jj,jiqm,k_n)
      use comanmr
      implicit none
      integer jj,jiqm,k_n
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
      write(30,'(a,2x,3f12.5)')element(jj),(coord(j,k_n,jj),j=1,3)
      write(31,'(a,i5,1x,a4,1x,a3,i6,4x,3f8.3,f8.4,f8.3)') 'ATOM  ', 
     & jj,atomname(jj),residue(jj),resno(jj),(coord(j,k_n,jj),j=1,3), 
     &  qmcharge(jj),rad(jj)

      end subroutine addatom
