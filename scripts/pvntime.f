



      program pvn
      implicit none
      integer argmax,nmax,tmax
      parameter (argmax=10,nmax=15000,tmax=21000)
      real pi,temp,gasconst,mu,kb
      integer narg,nskip,nuse,ndone,nframe
      character*200 arg(argmax),xtcfile,pdbfile
      integer uz,ret,xtcm,istep,natoms,i
      real gbox(3,3),time,prec,coord(nmax*3)
C========================================================================
C Declarations specific to this code: -
      character*200 outfile, out
      integer j,k,indice, indx, bin
      real distance, del, volume, count, rad, water
      real xdist,ydist,zdist,a(500),pv(500)
      real coordinates(nmax,tmax)
      real list(500), mindist, nw(100000)
      list = 0
      rad = 0.22
      temp = 298
      kb = 0.00831

C=======================================================================




      write(*,*)' '
      write(*,*)'PvN calculation started'
      write(*,*)' '

C     adding defaults---

      outfile='pvn.dat'
	out='NvTt.dat'
      nuse=100000
      nskip=0
C     natoms =0



C     reading input---
      narg = iargc()
      do i=1,narg
         call getarg(i,arg(i))
      enddo
C     reading input---

C MAKE NECESSARY CHANGES!!!!=============================================
      if(narg.eq.0)then
         write(*,*)'Use flags: -xtc -pdb -nskip -nuse'
         stop
      endif

      do i=1,narg
         if(arg(i).eq.'-xtc') xtcfile=arg(i+1)
         if(arg(i).eq.'-pdb') pdbfile=arg(i+1)
         if(arg(i).eq.'-o') outfile=arg(i+1)
         if(arg(i).eq.'-nskip') read(arg(i+1),*) nskip
         if(arg(i).eq.'-nuse') read(arg(i+1),*) nuse
C========================================================================
      enddo

      nframe=0
      ndone=0

C OPEN WHATEVER FILES YOU NEED TO HERE:=================================

      open(unit=20,file=outfile)
      open(unit=30,file=out)
      

c=======================================================================


C     opening xtc file---
      call xdrfopen(uz, xtcfile, "r",ret)
      if(ret.eq.1)then
         do while (ret.eq.1.and.nframe.le.nskip+nuse)
            call readxtc(uz,natoms,istep,time,gbox,coord,prec,ret)
		
C     opening xtc file---
            nframe = nframe+1
            if(nframe.gt.nskip.and.nframe.le.(nskip+nuse))then
               ndone=ndone+1

	
C========================================================================

C Part specific to this code starts here

C------------------------------------------------------------------------
                   
                   indice=0
                   water = 0
                   kloop: do k=1, natoms*3, 9
                      xdist = abs(1 - coord(k))
                      ydist = abs(1 - coord(k+1))
                      zdist = abs(1 - coord(k+2))
		if (xdist .ge. 2.5) then 
                         xdist = 5-xdist
                      end if
                      if (ydist .ge. 2.5) then 
                         ydist = 5-ydist 
                      end if
                      if (zdist .ge. 2.5) then 
                         zdist = 5-zdist
                      end if

                      
                      distance = sqrt((xdist**2)+(ydist**2)+(zdist**2))
                      if(distance .le. rad) then
                      indice = indice +1
		   water = water +1
                      end if
                      end do kloop
      list(indice)= list(indice)+1
      nw(ndone)= water
	write(30,*)nw(ndone)		
                

                   
                         
C------------------------------------------------------------------------

C Part specific to this code ends here

C==========================================================================

            endif
         enddo
      else
         write(*,*)'Something wrong'
         stop
      endif



C========================================================================

C Part specific to this code starts here

C------------------------------------------------------------------------
	

count = 0
      do j=1, 500
         
        a(j) = j
        pv(j) = list(j)/real(ndone)
        count = count + pv(j) 
         
      end do
	if (count .ne. 0 .and. pv(1) .ne. 0) then
	do while(j .le. 50)
	count = count + list(j)
	end do
	pv(0) = 1 - count
	end if
	
	do k=0, 100
	mu = -kb*temp*log(pv(0))
         write(20,*)a(k),pv(k),mu   
	end do
C------------------------------------------------------------------------

C Part specific to this code ends here

C==========================================================================
      write(*,*)' '
      write(*,*)'Free Energy' 
      write(*,*)' '
      
	

      write(*,*)real(ndone)
      write(*,*)'pvn calculation ended'
      write(*,*)' '



      stop
      end
