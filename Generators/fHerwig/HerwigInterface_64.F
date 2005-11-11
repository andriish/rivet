clhwout Herwig version 6.1, June-July-2000
c  6.3 versions Sept 2001
      subroutine herwiginterface
c----------------------------------------------------------------------
c   program initialization
c----------------------------------------------------------------------
c
c  HERWIG common blocks (`external` loads in  block data)
c-----------------------------------------------------------------------
#include  "herwig6400.inc"
      external hwudat
c-----------------------------------------------------------------------
c  local variables
c-----------------------------------------------------------------------
      integer       status, luherc, ierr, i, mm
cjl
      integer j2,j3,ii
      INTEGER IPDG,IWIG
      CHARACTER*8 NAME
c------------------------------------------------------------------------------ 
c  set defaults
c------------------------------------------------------------------------------

      call hwigin                       ! set HERWIG defaults
c
      maxpr         = 0                 ! Turn off screen data output
      lwsud         = 0
      lwevt         = 0
      iprint        = 1
c      PART1         = 'P       '
c      PART2         = 'P       '
c      IPROC         = 1500
      maxev         = 10000
       maxer         = 100
c      pbeam1        = 7000.0
c      pbeam2        = 7000.0
c      ptmin         = 100.0      
      lhwout        = 6
cjl sudord=2 for v 6.1
      sudord        = 2
cjl set default tlout to large negative, so timing-out feature is inoperative;
cjl   on systems where herwig's call timel works, user can reset tlout to 10
cjl   (or whatever) in talkto if desired (see Cern KERNLIB short writup Z007).
      tlout         = -10000000.0
cjl turn on mixing in Herwig - not using extenal b decay package
      mixing = .true.
cjl set stable all particles with b or c quarks - let QQ do them
c commented out for atlas
c      do ii=1,nres
cc        jj=abs(iflav(ii))
c        j1=mod(jj,10)
c        jj=jj/10
c        j2=mod(jj,10)
c        jj=jj/10
c        j3=mod(jj,10)
c        if(j1.eq.4.or.j2.eq.4.or.j3.eq.4) then
cjl  copied from hwusta.f
c          CALL HWUIDT(3,IPDG,IWIG,rname(ii))
c          RSTAB(IWIG)=.TRUE.
c        else if(j1.eq.5.or.j2.eq.5.or.j3.eq.5) then
c          CALL HWUIDT(3,IPDG,IWIG,rname(ii))
c          RSTAB(IWIG)=.TRUE.
c        end if
c      end do
c
c      WRITE (LHWOUT,10)
   10 FORMAT(/10X,
     &    'HERWIG: ALL PARTICLE TYPE WITH A B QUARK SET STABLE, USE QQ')
c      WRITE (LHWOUT,12)
 12   FORMAT(/10X,
     &    'HERWIG: ALL PARTICLE TYPE WITH A C QUARK SET STABLE, USE QQ')
cjl set pipsmr false, so primary vertex at coors (0.,0.,0.)
      pipsmr=.false. 
c
      return
      end

C
C     The following two dummy routines are needed by Herwig...

      Subroutine Hvcbvi
#include  "herwig6400.inc"
      Write(lhwout,10)
 10   Format(/10x,'HVCBVI called but not linked')
      end


      Subroutine Hwaend
      end

c
      Subroutine cdfreadsusy(filnam,iunit)
      character*(*)  filnam
      integer status
c      integer sigtlu, sifrlu
      integer iunit
c-----------------------------------------------------------------------
c  open the SUSY particle properties file
c-----------------------------------------------------------------------

c      status        = sigtlu(iunit)
      write(*,*) ' '
      write (*,6660)
 6660 format (' SUSY filename is ')
      write (*,*) filnam
      write (*,6664) iunit
 6664 format ('      and unit number is ',i6)
C
      open(unit=iunit,file=filnam,status='old',err=101
     &     ,FORM='FORMATTED',ACCESS='SEQUENTIAL')
      goto 102
c-----------------------------------------------------------------------
c  error: can't open SUSY file
c-----------------------------------------------------------------------
101   continue
      write(*,*) ' '
      write(*,*) ' ********************************************* '
      write(*,*) ' '
      write(*,*) 
     $     ' Error in opening SUSY file for READING'
      write(*,*)
     $     ' Program will abort.'
      write(*,*) ' '
      write(*,*) ' ********************************************* '
      write(*,*) ' '
102   continue
      end
c
c-----------------------------------------------------------------------
c  dummy hwanal, called from HerwigModule.cxx, can be useful
      subroutine hwanal(ihwcod)
c  input ihwcod = 1, set to zero if want this event skipped
#include  "herwig6400.inc"
      integer ihwcod
c
c
c skip if ierror ne 0
      if(ierror.ne.0) return
c
cc  do whatever here
cc  if want to skip this event, say  ihwcod = 0
c
      end
