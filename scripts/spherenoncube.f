	subroutine spherenoncube(box_dims,v,dr)
	! program spherenoncube
************************************************************************
* volume of a sphere of radius r inside a non-cube                     *
*                                                                      *
* Lu Yang, April 29, 2000. RPI.                                        *
************************************************************************
        
	implicit none
	real*8 box_dims(1:3)
	real*8 r,dr,r1,r2,r3,pi,rr,box,boy,boz,b(4)
	real*8 L1,L2,L3,L4,L5,L6,L7,ss,ss1,ss2,ss3
        real*8 limit0,limit1,b_func,func,dv
	integer i,nr,nrmax
	parameter (nrmax=5000)
	real*8 v(0:nrmax)
        common /sphere/r,b_func
        external func
	

	box=box_dims(1)
	boy=box_dims(2)
	boz=box_dims(3)
	! box=1
	! boy=2
	! boz=3
	! dr=0.0005
	nr=int(sqrt((box/2)**2 + (boy/2)**2 + (boz/2)**2)/dr)+1
	v(0)=0.
	do i=1,nr
	   r=i*dr
	   
c Compare box,boy,boz and sort them from minize to maximum
	   b(1)=box
	   b(2)=boy
	   b(3)=boz
	   rr = r**2
	   if(b(1).LT.b(2).and.b(2).LT.b(3))then
	      b(1) = b(1)
	      b(2) = b(2)
	      b(3) = b(3)
	   else if(b(1).LT.b(3).and.b(3).LT.b(2))then
	      b(1) = b(1)
	      b(4) = b(2)
	      b(2) = b(3)
	      b(3) = b(4)
	   else if(b(2).LT.b(1).and.b(1).LT.b(3))then
	      b(4) = b(1)
	      b(1) = b(2)
	      b(2) = b(4)
	      b(3) = b(3)
	   else if(b(2).LT.b(3).and.b(3).LT.b(1))then
	      b(4) = b(1)
	      b(1) = b(2)
	      b(2) = b(3)
	      b(3) = b(4)
	   else if(b(3).LT.b(2).and.b(2).LT.b(1))then
	      b(4) = b(1)
	      b(1) = b(3)
	      b(2) = b(2)
	      b(3) = b(4)
	   else
	      b(4) = b(1)
	      b(1) = b(3)
	      b(3) = b(2)
	      b(2) = b(4)
	   end if

c Calculate the half length of the box
	   L1=b(1)/2.0d0
	   L2=b(2)/2.0d0
	   L3=b(3)/2.0d0
	   L4=sqrt(L1**2+L2**2)
	   L5=sqrt(L1**2+L3**2)
	   L6=sqrt(L2**2+L3**2)
	   L7=sqrt(L1**2+L2**2+L3**2)
	   pi=4.d0*atan(1.0)
	   

c Compare L3 and L4

	   if(L3.lt.L4) then
	      if(r.le.0.)then
		 v(i)=0.d0
	      elseif(r.lt.L1)then
		 v(i)=4.d0*pi*r**3/3.d0
	      elseif(r.lt.L2)then
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L1*pi*rr+L1**3*pi/3.d0)
	      elseif(r.lt.L3)then
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L1*pi*rr+L1**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-L2*pi*rr+L2**3*pi/3.d0)
	      elseif(r.lt.L4)then
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L1*pi*rr+(L1**3)*pi/3.d0)
	2	      -2.d0*(2.d0*pi*(r**3)/3.d0-L2*pi*rr+(L2**3)*pi/3.d0)
	3	      -2.d0*(2.d0*pi*r**3/3.d0-L3*pi*rr+(L3**3)*pi/3.d0)
	      elseif(r.lt.L5)then
		 limit0=L2
		 limit1=sqrt(rr-L1**2)
		 b_func=2.d0*L1
		 call qromb(func,limit0,limit1,ss)
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L3*pi*rr+L3**3*pi/3.d0)-
	2	      2.d0*(2.d0*pi*r**3/3.d0-L1*pi*rr+L1**3*pi/3.d0)-
	3	      2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L1**2)*pi*rr+
	4	      (sqrt(rr-L1**2))**3*pi/3.d0)-
	5	      ss
		 
	      elseif(r.lt.L6)then
		 limit0=L1
		 limit1=sqrt(rr-L3**2)
		 b_func=2.d0*L3
		 call qromb(func,limit0,limit1,ss)
		 ss1=ss
		 limit0=L2
		 limit1=sqrt(rr-L1**2)
		 b_func=2.d0*L1
		 call qromb(func,limit0,limit1,ss)
		 ss2=ss
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L3*pi*rr+L3**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L3**2)*pi*rr
	3	      +(sqrt(rr-L3**2))**3*pi/3.d0)
	4	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L1**2)*pi*rr
	5	      +(sqrt(rr-L1**2))**3*pi/3.d0)
	6	      -ss1-ss2
	      elseif(r.lt.L7)then
		 limit0=L1
		 limit1=sqrt(rr-L3**2)
		 b_func=2.d0*L3
		 call qromb(func,limit0,limit1,ss)
		 ss1=ss
		 limit0=L2
		 limit1=sqrt(rr-L1**2)
		 b_func=2.d0*L1
		 call qromb(func,limit0,limit1,ss)
		 ss2=ss
		 limit0=L2
		 limit1=sqrt(rr-L3**2)
		 b_func=2.d0*L3
		 call qromb(func,limit0,limit1,ss)
		 ss3=ss
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L3*pi*rr+L3**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L3**2)*pi*rr
	3	      +(sqrt(rr-L3**2))**3*pi/3.d0)
	4	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L1**2)*pi*rr
	5	      +(sqrt(rr-L1**2))**3*pi/3.d0)
	6	      +2.d0*pi*(rr*sqrt(rr-L3**2)-(sqrt(rr-L3**2))**3/3.d0
	1	      -rr*L2+L2**3/3.d0)
	2	      -ss1-ss2-ss3
	      else
		 v(i)=8.d0*L1*L2*L3
	      endif
	   else
	      if(r.le.0.)then
		 v(i)=0.d0
	      elseif(r.lt.L1)then
		 v(i)=4.d0*pi*r**3/3.d0
	      elseif(r.lt.L2)then
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L1*pi*rr+L1**3*pi/3.d0)
	      elseif(r.lt.L4)then
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L1*pi*rr+L1**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-L2*pi*rr+L2**3*pi/3.d0)
	      elseif(r.lt.L3)then
		 limit0=L2
		 limit1=sqrt(rr-L1**2)
		 b_func=2.d0*L1
		 call qromb(func,limit0,limit1,ss)
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L1*pi*rr+L1**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L1**2)*pi*rr
	3	      +(sqrt(rr-L1**2))**3*pi/3.d0)-ss
	      elseif(r.lt.L5)then
		 limit0=L2
		 limit1=sqrt(rr-L1**2)
		 b_func=2.d0*L1
		 call qromb(func,limit0,limit1,ss)
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L3*pi*rr+L3**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-L1*pi*rr+L1**3*pi/3.d0)
	3	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L1**2)*pi*rr
	4	      +(sqrt(rr-L1**2))**3*pi/3.d0)
	5	      -ss
	      elseif(r.lt.L6)then
		 limit0=L1
		 limit1=sqrt(rr-L3**2)
		 b_func=2.d0*L3
		 call qromb(func,limit0,limit1,ss)
		 ss1=ss
		 limit0=L2
		 limit1=sqrt(rr-L1**2)
		 b_func=2.d0*L1
		 call qromb(func,limit0,limit1,ss)
		 ss2=ss
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L3*pi*rr+L3**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L3**2)*pi*rr
	3	      +(sqrt(rr-L3**2))**3*pi/3.d0)
	4	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L1**2)*pi*rr
	5	      +(sqrt(rr-L1**2))**3*pi/3.d0)
	6	      -ss1-ss2
	      elseif(r.lt.L7)then
		 limit0=L1
		 limit1=sqrt(rr-L3**2)
		 b_func=2.d0*L3
		 call qromb(func,limit0,limit1,ss)
		 ss1=ss
		 limit0=L2
		 limit1=sqrt(rr-L1**2)
		 b_func=2.d0*L1
		 call qromb(func,limit0,limit1,ss)
		 ss2=ss
		 limit0=L2
		 limit1=sqrt(rr-L3**2)
		 b_func=2.d0*L3
		 call qromb(func,limit0,limit1,ss)
		 ss3=ss
		 v(i)=4.d0*pi*r**3/3.d0-2.d0*(2.d0*pi*r**3/3.d0-
	1	      L3*pi*rr+L3**3*pi/3.d0)
	2	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L3**2)*pi*rr
	3	      +(sqrt(rr-L3**2))**3*pi/3.d0)
	4	      -2.d0*(2.d0*pi*r**3/3.d0-sqrt(rr-L1**2)*pi*rr
	5	      +(sqrt(rr-L1**2))**3*pi/3.d0)
	6	      +2.d0*pi*(rr*sqrt(rr-L3**2)-(sqrt(rr-L3**2))**3/3.d0
	7	      -rr*L2+L2**3/3.d0)
	8	      -ss1-ss2-ss3
	      else
		 v(i)=8.d0*L1*L2*L3
	      endif
	   end if
	enddo
	! do i=1,nr-1
	!    dv=4*pi*((i*dr)**3-((i-1)*dr)**3)/3
	!    write(*,*)i*dr,(v(i)-v(i-1))/dv
	! enddo

	end
	
	function func(y)
	real*8 func,r,y,b_func
	common /sphere/r,b_func
	
c       write(*,*)r,b,y,r**2-y**2-(b/2.)**2
	if((r**2-y**2-(b_func/2.)**2).lt.0.0)then
	   if(abs(r**2-y**2-(b_func/2.)**2).lt.0.0001)then
	      func = 0.0
	   endif
	   return
	endif
c       stop'problem evaluating func'
	pi = 4.d0*datan(1.0d0)
	func = 2.d0*(b_func*sqrt(r**2-y**2-(b_func/2.)**2)
	1    +2.d0*(r**2-y**2)*asin(b_func/2./sqrt(r**2-y**2)))
	return
	end
	
	SUBROUTINE QROMB(FUNC,A,B,SS)
	implicit real*8 (a-h,o-z)
	integer n
	parameter (n=20)
	it = n
	TNM=real(it)
	DEL=(B-A)/TNM
	X=A+0.5*DEL
	SUM=0.0d0
	DO 11 J=1,IT
	   SUM=SUM+FUNC(X)
	   X=X+DEL
 11	CONTINUE
	SS=DEL*SUM

	RETURN
	END






