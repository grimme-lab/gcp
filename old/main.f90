program gcp_main
use strings
use gcp
implicit none
integer n,maxat
integer nn,nb,i,zz,j
parameter (maxat =500000)
! coordinates in au
real*8 xyz(3,maxat)
! gradient
real*8 gcp_g(3,maxat)
! cardinal numbers of elements
integer   iz(maxat)
! frozen coordinate
integer ifrez(maxat)
real*8  xx(10),gcp_e,r
real*8 t0,t1
logical ex,echo,dograd,dohess,helpme,warn,lib,parm,test,verb,tmen,parfile
character*20 method,ctmp
character*80 arg(10)
character*80 atmp,ftmp,dtmp
!pbc stuff
logical vasp, stress,pbc
character*80 args(90), comment, line
real*8 lat(3,3), gcp_glat(3,3)
integer nargs, ios

!** subroutine ?? **
lib=.false.

!************
!* Defaults *
!************
warn=.false.
parm=.false.
dohess=.false.
helpme=.false.
echo=.true.
dograd=.false.
verb=.false.
tmen=.false.
pbc=.false.
vasp=.false.
stress=.false.
parfile=.false.
arg=''
method=''


!****************************************
!* Get and process command line options *
!****************************************
if(.not.lib) then
  do i=1,10
     call getarg(i,arg(i))
  enddo
endif

do i=1,10
ftmp=arg(i)
if(index(ftmp,'-h ').ne.0) helpme=.true.
if(index(ftmp,'-noprint ').ne.0) echo=.false.
if(index(ftmp,'-level ').ne.0) method=arg(i+1)
if(index(ftmp,'-func ').ne.0) method=arg(i+1)
if(index(ftmp,'-l ').ne.0) method=arg(i+1)
if(index(ftmp,'-grad ').ne.0) dograd=.true.
if(index(ftmp,'-hess ').ne.0) dohess=.true.
if(index(ftmp,'-v ').ne.0) verb=.true.
if(index(ftmp,'-pbc ').ne.0) pbc=.true.
if(index(ftmp,'-vasp ').ne.0) vasp=.true.
if(index(ftmp,'-stress ').ne.0) stress=.true.
if(index(ftmp,'-tme'   ).ne.0)tmen=.true.
if(index(ftmp,'-parfile ').ne.0) parfile=.true.
enddo
method=lowercase(method)

!************
!* Header *
!************
if(echo)then
write(*,*)' ___________________________________________ '
write(*,*)'                                             '
write(*,*)'|              **  g C P  **                |'
write(*,*)'|  a geometrical counterpoise correction    |'
write(*,*)'|     H.Kruse J.G.Brandenburg S.Grimme      |'
write(*,*)'|          Version 2.02 Nov 2016            |'
write(*,*)'|                                           |'
write(*,*)' ___________________________________________ '
write(*,*)
write(*,*)'Please cite work done with this code as:'
write(*,*)'H. Kruse, S. Grimme J. Chem. Phys. 136, 154101 (2012)'
write(*,*)'DOI: 10.1063/1.3700154 '
write(*,*)'For the periodic version, please also cite:'
write(*,*)'J. G. Brandenburg, M. Alessio, B. Civalleri, M. F. Peintinger'
write(*,*)'T. Bredow, S.Grimme J. Phys. Chem. A 117, 9282-9292 (2013).'
write(*,*)'DOI: 10.1021/jp406658y '
write(*,*)
write(*,*)
endif

if (helpme) call help(.true.)

!***************************************
!* read coordinates                    *
!* Turbomole format or xmol format     *
!***************************************
if(.not.lib) then
   if(pbc) then
      atmp=arg(1)
      inquire(file=atmp,exist=ex)
      if(.not.ex) then
         write(*,'(a)')   ''
         write(*,'(a)')  '** ERROR:  no <coordinates> found! **'
         write(*,'(a)')   ''
         write(*,'(a)')   'gcp <coordintes> -level <method>'
         write(*,'(a)')   'see gcp -h for further info '
         write(*,'(a)')   ''
         stop 'input error'
      endif
      !Determine file format
      open(unit=3,file=atmp)
      read(3,'(a)') line
      call parse(line,' ',args,nargs)
     ! call value(args(1),comment,ios)
      if(args(1).eq.'$cell'.and.args(2).eq.'vectors') then
         vasp=.false.
      end if
      close(3)
      if(vasp) then
         call rdatomnumbervasp(atmp,n,echo)
         call rdcoordvasp(xyz,lat,iz,n,atmp,echo)
!         call rdcoordvasp(xyz,lat,iz,n,atmp,echo)
         else
         if(atmp.eq.'POSCAR') then
            call rdatomnumbervasp(atmp,n,echo)
            call rdcoordvasp(xyz,lat,iz,n,atmp,echo)
         end if
      end if
   else
      atmp=arg(1)
      inquire(file=atmp,exist=ex)
      if(.not.ex) then
         write(*,'(a)')   ''
         write(*,'(a)')  '** ERROR:  no <coordinates> found! **'
         write(*,'(a)')   ''
         write(*,'(a)')   'gcp <coordintes> -level <method>'
         write(*,'(a)')   'see gcp -h for further info '
         write(*,'(a)')   ''
         stop 'input error'
      endif
      call tmolrd(maxat,xyz,iz,ifrez,n,atmp,echo)
   end if
endif

if(n.lt.1)     stop 'no atoms'
if(n.gt.maxat) stop 'too many atoms'


call gcp_call(n,xyz,lat,iz,gcp_e,gcp_g,gcp_glat,dograd,dohess,pbc,method,echo,parfile)


! print gradient to std outout
  if(verb) then
   write(*,'(''gradient: Ggcp'')')
   do i=1,n
    write(*,'(3E22.13)')gcp_g(1,i),gcp_g(2,i),gcp_g(3,i)
   enddo
  endif

    if(echo)then
      write(*,*) '                '
      write(*,*) '** gCP correction ** '
      write(*,'(2x,a7,F18.10,'' / (a.u.) || '',x,F11.4,'' / (kcal/mol)'')') 'Egcp:  ', gcp_e,gcp_e*627.5099d0
      write(*,*) '                '
       if(dograd)then
       write(*,*)'|G|=',sum(abs(gcp_g(1:3,1:n)))
       endif
    endif

! write gradient
    if(dograd.and..not.verb) then
       call wregrad_tm(maxat,n,xyz,iz,ifrez,gcp_e,gcp_g,echo)
    endif

    call cpu_time(t1)
    if(echo) call prtim(6,t1-t0,'t','gCP ')
    if(warn.and.echo) write(*,*) 'Carefully read the notifications/warnings given at loadup'
    if(echo) call done('gCP',6)
 end program gcp_main
