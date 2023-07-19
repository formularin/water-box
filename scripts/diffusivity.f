C Mandatory arguments: -xtc
C =============================================================================
C General declarations
C uz: ID value for internal use in libxdr
C ret: 1 indicates success in opening XTC file, 0 indicates error

      program pvn
      implicit none
      integer argmax, nmax, tmax
      parameter (argmax=10, nmax=15000, tmax=21000)
      real pi, temp, gasconst, mu, kb
      integer narg, nskip, nuse, ndone, nframe
      character*200 arg(argmax), xtcfile, pdbfile
      integer uz, ret, xtcm, istep, natoms, i
      real gbox(3, 3), time, prec, coord(nmax*3)
      character*200 outfile

C =============================================================================
C Declarations specific to this code

      integer num_samples, sample_size
      parameter (num_samples=125, sample_size=1250)

C     f0 is the starting frame of the current sample
      integer f0, f0_interval
      real vacf(num_samples), MSD(num_samples)
      integer num_part
C     Can't have a dynamic array size, so unfortunately the program has to be edited with the number of water molecules and number of frames.
      parameter (num_part=392, nuse=2500)
      real all_coords(nuse, num_part*3), velocities(nuse, num_part*3)

      f0_interval = nuse / num_samples

C =============================================================================
C Reading command-line arguments

      outfile="out.dat"

      narg = iargc()
      do i=1,narg
         call getarg(i,arg(i))
      enddo

      if(narg.eq.0)then
         write(*,*)'Use flags: -xtc -pdb -nskip -nuse'
         stop
      endif

      do i=1,narg
         if(arg(i).eq.'-xtc') xtcfile=arg(i+1)
         if(arg(i).eq.'-pdb') pdbfile=arg(i+1)
         if(arg(i).eq.'-o') outfile=arg(i+1)
         if(arg(i).eq.'-nskip') read(arg(i+1),*) nskip
C         if(arg(i).eq.'-nuse') read(arg(i+1),*) nuse
      enddo

      f0 = nskip + 1

C =============================================================================
C Opening files

      open(unit=20, file=outfile)

      call xdrfopen(uz, xtcfile, "r", ret)
      if (ret .eq. 1) then

C =============================================================================

         do while (ret .eq. 1 .and. nframe .le. nskip+nuse)
            call readxtc(uz, natoms, istep, time, gbox, coord, prec, ret)
            nframe = nframe + 1
            if (nframe .gt. nskip .and. nframe .le. (nskip + nuse)) then
               ndone = ndone + 1

C =============================================================================
C Diffusivity/VACF calculation

               if (ndone - f0 .eq. sample_size) then
                  f0 = f0 + f0_interval 
               endif

               

C =============================================================================

            endif
        enddo
      else
        write(*,*) "Something's wrong"
      endif
      stop
      end program pvn