! This file is part of mctc-gcp.
! SPDX-Identifier: GPL-3.0-or-later
!
! mctc-gcp is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! mctc-gcp is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with mctc-gcp.  If not, see <https://www.gnu.org/licenses/>.

module gcp
   use gcp_version, only : get_gcp_version
contains
!subroutine for interfacing with orca/turbomole/crystal
subroutine gcp_call(n,xyz,lat,iz,gcp_e,gcp_g,gcp_glat,dograd,dohess,pbc,method,echo,parfile)
implicit none
integer n   !number of atoms
real*8 xyz(3,n) !xyzcoordinates
real*8 lat(3,3) !lattice matrix
integer iz(n) !element number
real*8 gcp_g(3,n) !xyz gradient
real*8 gcp_glat(3,3)
logical dograd,dohess,echo !calculate gradient, calculate hessian, print informations
logical pbc !periodic system
character*20 method !method string
integer max_elem,max_para
integer nn,nb,i,zz,j
! max elements possible
parameter (max_elem=94)
! actually parameterized
parameter (max_para=36)
real*8 autoang
parameter (autoang =0.5291772083d0)
! restore cardinal numbers of elements
integer   iz2(n)
! frozen coordinate
integer ifrez(n)
! #electrons of elements
integer   nel(max_elem)
! # basis functions
integer nbas(max_elem)
real*8 emiss(max_elem)
real*8 tmp  (max_elem)
! parameters
real*8 p(4)
real*8 xva(n),xvb(n)
real*8  xx(10),gcp_e,r
real*8 t0,t1
logical ex,warn,lib,parm,test,verb
character*20 ctmp
character*80 arg(10)
character*80 atmp,ftmp,dtmp
character*80 basname
!v2.0 stuff
logical tmen,base,damp,stress,parfile,srb,newgcp
character*80 args(90), comment, line
real*8 er,el,ebas,gbas(3,n)
integer nargs, ios
real*8 dmp_scal,dmp_exp,rscal,qscal


!default no srvb
srb=.false.
!damping params
dmp_scal=4.0d0  ! 4
dmp_exp=6.0d0   ! 6
!default
base=.false.
damp=.true.

!************************
!* modify input string  *
!* 1. use lower case    *
!* 2. delete hyphends   *
!************************
do while (index(method,'-').ne.0)
  i=scan(method,'-')
  ctmp=trim(method(:(i-1)))//trim(method((i+1):))
  method=trim(ctmp)
enddo

select case(method)
case('hf3c')
!JGB special method combination 1: HF-3c
if(echo) write(*,*) ' Perform special gCP-SRB correction for HF-3c '
base=.true.
case('pbeh3c','pbeh3cmod')
!JGB special method combination 2: PBEh-3c
if(echo) write(*,*) ' Perform special gCP(damped) correction for PBEh-3c'
damp=.true.
case('hse3c')
!JGB special method combination 3: HSE-3c
if(echo) write(*,*) ' Perform special gCP(damped) correction for HSE-3c'
damp=.true.
!case('b3pbe3c')
!JGB special method combination 3: B3PBE-3c
!if(echo) write(*,*) ' Perform special gCP(damped) correction for B3PBE-3c'
!damp=.true.
case('b973c')
!JGB special method combination 4: B97-3c
!JGB no gCP needed, switch to SRB correction
if(echo) write(*,*) ' Perform special SRB correction for B97-3c '
srb=.true.
rscal=10.00d0
qscal=0.08d0
case default
end select


!skip for B97-3c
if(.not.srb) then
!get parameters
!****************************************************************************
!* read parameter file                                                      *
!* format:                                                                  *
!* basis sigma eta alpha beta                                               *
!* To read in nbas and emiss, generate param file with                      *
!* with '-parfile', edit copy it to .gcppar:                              *
!* A '#' and the beginning will tell the program to read a different format *
!****************************************************************************
if(method.eq.''.or.method.eq.'file')then
!* Fortran 2003 *
  call get_environment_variable('HOSTNAME',ftmp)
  call get_environment_variable('HOME',atmp)

!* uncomment the following line                  *
!* if your compiler wont support above features  *
!
!  call system('hostname > .tmpx')
!  open(unit=43,file='.tmpx')
!  read(43,'(a)')ftmp
!  close(43,status='delete')
!  call system('echo $HOME > .tmpx')
!  open(unit=43,file='.tmpx')
!  read(43,'(a)')atmp
!  close(43,status='delete')
  write(dtmp,'(a,''/.gcppar.'',a)') trim(atmp), trim(ftmp)
  inquire(file=trim(dtmp),exist=ex)
  if(ex)then
     if(echo)write(*,*) 'reading ',trim(dtmp)
     open(unit=42,file=dtmp)
     read(42,'(a)',end=9)ftmp
     call charXsplit(ftmp,basname,1)
     call lower_case(basname)
     call readl(ftmp,xx,nn)
     p(1:4)=xx(1:4)
      if(adjustl(trim(basname)).eq.'#') then
       if(echo) then
       write(*,*) ' found extended param file ! '
       endif
       parm=.true.
       read(42,'(x,F8.4)',end=9) p(1:4)
       do i=1,36
          read(42,'(x,I3,x,F9.6)',end=9) nbas(i), emiss(i)
       enddo
!additional read for new param file
      elseif(adjustl(trim(basname)).eq.'##') then
      if(echo) then
       write(*,*) ' found extended param file v2! '
      endif
      call readl(ftmp,xx,nn)
      parm=.true.
      damp=.true.
      dmp_scal=xx(5)
      dmp_exp=xx(6)
      do i=1,36
          read(42,'(x,I3,x,F9.6)',end=9) nbas(i), emiss(i)
       enddo
     endif
9   close(42)
     if(.not.parm) then
      if(echo) write(*,*) 'loading ',trim(basname),' params'
      call setparam(emiss,nbas,p,basname)
     endif
  else
   write(*,'(a)') ''
   write(*,'(2a)')  'found not param file in ',trim(dtmp)
   error stop '* ERROR *'
  endif
method=trim(basname)
!*******************
!* built-in method *
!*******************
else
    call setparam(emiss,nbas,p,method)
endif

!open(unit=42,file='~/.gcppar')
!read(42,*) p(1:4),dmp_scal,dmp_exp
!do i=1,36
!read(42,*) emiss(i),tmp(i)
!enddo
!close(42)


!**********************************************************
!* check for missing parameters for current elements      *
!* elements >36 will be treated as their lower homologues *
!* eg: Ag -> Cu, In -> Ga                                 *
!* prepare virtual bf array xva,xvb                       *
!**********************************************************
! backup original
iz2=iz
! re-arrange elements
do i=1,n
zz=iz(i)
select case (zz)
 case(37:54)
  iz(i)=iz(i)-18
  write(*,'(4a)') '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
 case(55:57)
  iz(i)=iz(i)-18*2
  write(*,'(4a)') '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
 case(58:71,90:94)
  iz(i)=21
  write(*,'(4a)') '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
 case(72:89)
  iz(i)=iz(i)-2*18-14
  write(*,'(4a)') '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
end select
enddo

! set electrons
do i=1,maxval(iz)
  nel(i)=i
enddo


!* handle ECP basis set
if(index(method,'vmb').ne.0) then
write(*,'(a)') 'VMB basis has ECPs. Adjusting electrons'
do i=5,10
nel(i)=nel(i)-2
enddo
do i=11,18
nel(i)=nel(i)-10
enddo
endif

! count bf
nb=0
do nn=1,n
   nb=nb+nbas(iz(nn))
enddo
if(nb.lt.1) error stop 'Nbf setup gone wrong'

! set virtuals and look for missing parameters
xva=0
do i=1,n
  if(method.eq.'def2mtzvpp')then  !SG
   xva(i)=1.0d0                   !SG
   if(iz(i).eq.6) xva(i)=3.0d0    !SG
   if(iz(i).eq.7) xva(i)=0.5d0    !SG
   if(iz(i).eq.8) xva(i)=0.5d0    !SG
  else
   xva(i)=(nbas(iz(i))-0.5d0*nel(iz(i)))
  endif
  if(emiss(iz(i)).gt.01e-6.and.xva(i).lt.0.0d0) then
  write(*,'(a)') 'element/emiss/nvirt/nel/nbas'
  write(*,'(a4,f8.5,f8.5,i3,i3)') esym(iz(i)),emiss(iz(i)),xva(i),nel(iz(i)),nbas(iz(i))
   error stop 'negative number of virtual orbitals. Something is wrong in the parameters !'
  endif
  if( emiss(iz(i)).lt.0.1e-6) then
!    print*, '** WARNING ** -> element ',esym(iz(i)),' has no parameters       (no contribution)!'
!    xva(i)=0
  cycle
  warn=.true.
  endif
  if(xva(i).lt.0.5d0) then
  write(*,'(3a)') '** WARNING ** -> element ',esym(iz(i)),' has no virtual orbitals (no contribution)!'
  warn=.true.
  endif
enddo
xvb=xva

endif !B97-3c

!*********************************
!* Parameter handling/setup done *
!*********************************
  if(echo)then
     write(*,*)''
     write(*,'(2x,''level '',3x,a12,2x,a)') adjustr(trim(method)),adjustr(trim(getlevel(method)))
     write(*,'(2x,''Nbf   '',3x,I12)')  nb
     write(*,'(2x,''Atoms '',3x,I12)')  n
     write(*,*)''
     if(.not.srb) then
        write(*,'(2x,a)') 'Parameters: '
        write(*,'(2x,a,2x,f10.4)') 'sigma   ',p(1)
        write(*,'(2x,a,2x,f10.4)') 'eta     ',p(2)
        write(*,'(2x,a,2x,f10.4)') 'alpha   ',p(3)
        write(*,'(2x,a,2x,f10.4)') 'beta    ',p(4)
        if(damp) then
           write(*,'(2x,a,2x,f10.4)') 'dmp_scal',dmp_scal
           write(*,'(2x,a,2x,f10.4)') 'dmp_exp ',dmp_exp
        endif
     else
        write(*,'(/2x,a22)')'Parameters for SRB: '
        write(*,'(2x,a,2x,f10.4)') 'rscal',rscal
        write(*,'(2x,a,2x,f10.4)') 'qscal',qscal
     endif
     write(*,*)''
   endif



! **************************
! * print parameter table  *
! **************************

if(.not.srb) then
if(echo) then
  write(*,*) ' '
  write(*,*) 'element parameters ',trim(method)
  write(*,'(2x,3(a4,2x,a6,3x,a4,4x))') 'elem','emiss','nbas','elem','emiss','nbas','elem','emiss','nbas'
  do i=1,36,3
      write(*,'(2x,3(A4,2x,F8.5,2x,I3,4x))') esym(i), emiss(i),nbas(i), &
                                           esym(i+1), emiss(i+1),nbas(i+1),&
                                           esym(i+2), emiss(i+2),nbas(i+2)
  enddo
  write(*,*) ' '
endif
endif


! **************************
! * write parameter file   *
! **************************
if(parfile) then
  write(*,'(a)')  'printing extended parameter file ...'
  write(*,'(a)')  'modfiy & copy gcp.param to $HOME/.gcppar.$HOSTNAME'
  open(unit=22,file='gcp.param')
    write(22,*) ' # ',trim(method),' (comment line)'
    write(22,'(x,F8.5)') p(1:4)
    do i=1,36
      write(22,'(x,I3,x,F8.5)') nbas(i), emiss(i)
    enddo
  close(22)
endif

!**************************
!* calc energy & gradient *
!**************************
gcp_g(1:3,1:n)=0.0d0
gcp_e=0.0d0
if(srb) then
   call srb_egrad2(xyz,iz,lat,n,gcp_e,gcp_g,dograd,rscal,qscal,echo,pbc)
else
   call gcp_egrad(n,max_elem,emiss,xyz,iz,p,gcp_e,gcp_g,dograd,echo,xva,xvb,pbc,lat,damp,base,dmp_scal,dmp_exp)
endif
!*************************
!* calc num. cell stress *
!*************************
if(dograd.and.pbc) then
call sgcp(n,max_elem,emiss,xyz,iz,p,gcp_e,gcp_g,dograd,echo,xva,xvb,lat,gcp_glat,base,damp,dmp_scal,dmp_exp)

! calc hessian
if(dohess) then
!     call  hessian(n,max_elem,emiss,xyz,iz,p,gcp_e,xva,xvb,dohess,damp,base,dmp_scal,dmp_exp)
    call  hessian(n,max_elem,emiss,xyz,iz,p,xva,xvb,pbc,lat,damp,base,dmp_scal,dmp_exp)
endif

! restore original cardinal numbers
iz=iz2
! set cart. frozen coords to zero
if(maxval(ifrez,n).eq.1.and.dograd) then
if(echo) write(*,'(a)') ' setting gradient of frozen cartesians to zero'
  do i=1,n
   if(ifrez(i).eq.1) gcp_g(1:3,i)=0.0d0
  enddo
endif

!output
if(echo) then
 write(*,*) ' '
 write(*,*) '|S|=',abs(sum(gcp_glat(1:3,1:3)))
 if(.not.verb) then
 write(*,'(a)')  'Write cell gradient to file gcp_cellgradient'
 open(unit=1,file='gcp_cellgradient')
 write(1,*) gcp_glat
 close(1)
 end if
 write(*,*) ' '
 end if
end if

!for scripting
open(unit=1,file='.CPC')
write(1,*) gcp_e
close(1)
open(unit=1,file='.CP')
write(1,'(F22.16)') gcp_e
close(1)


end subroutine gcp_call



!***************************************************
!* This subroutine calculates the gcp correction *
!***************************************************
subroutine gcp_egrad(n,max_elem,emiss,xyz,iz,p,gcp,g,grad,echo,xva,xvb,pbc,lat,damp,base,dmp_scal,dmp_exp)
implicit none
integer n,iz(n),max_elem,np
integer iat,jat
real*8  xyz(3,n), xyzjat(3)
real*8  g  (3,n)
real*8  gcp,tg,t0,t1
real*8  emiss(max_elem),dum22
real*8  p(4),thrR,thrE
real*8 xva(*),xvb(*)
real*8 va,vb
real*8 r,rscal,rscalexp,rscalexpm1,r0abij
real*8 tmp,ea,ecp,dum,tmpa,tmpb,dum2,tmpc,tmpd
real*8 sab,p1,p2,p3,p4,ene_old_num,ene_old_den
real*8 vec(3),gs(3),gab,ebas,gbas(3,n)
real*8 za(36),zb(36),lat(3,3)
logical echo,grad,damp,base,pbc
integer tau_max(3),a,b,c

real*8 rad(94)
data rad /&
0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0,&
0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0,&
2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0,&
1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0,&
2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0,&
1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0,&
2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0,&
1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0,&
1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0,&
1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,2.01D0,1.90D0,&
1.84D0,1.83D0,1.80D0,1.80D0 /

! cut-off radii for all element pairs
real*8 r0ab(max_elem,max_elem),autoang
!real*8 autoang
!for damping
real*8 dampval, grdfirst,grdsecond,ene_old,ene_dmp,grd_dmp,dmp_scal,dmp_exp
parameter (autoang =0.5291772083d0)
! Two threshold. thrR: distance cutoff thrE: numerical noise cutoff
thrR=60            ! 60 bohr
thrE=epsilon(1.d0) ! cut below machine precision rounding

if(echo) then
write(*,*) '  '
write(*,'('' cutoffs: '',F5.1)')
write(*,'(''   distance [bohr]'',F5.1)') thrR
write(*,'(''   noise    [a.u.]'',Es8.1)') thrE
write(*,'(''   SR-damping     '',L2)') damp
write(*,*) '  '
endif

g=0.0d0
ecp=0.0d0
gcp=0.0d0
dum=0.0d0
tg=0d0

p1=abs(p(1))
p2=abs(p(2))
p3=abs(p(3))
p4=abs(p(4))

!get r0
call setr0ab(max_elem,autoang,r0ab)

if(abs(p(1)-1.0d0).lt.1.d-6.and.abs(p(2)-1.315d0).lt.1.d-6)then  !SG dirty way to set the special r2SCAN-3c mode
   call setzet(p2,1.15d0,za,zb)
else
   call setzet(p2,1.00d0,za,zb)
endif

   gs=0
   if(echo) write(*,'(2x,a5,2x,a5,2x,a5,2x,a7,2x,a14,4x,a15)') &
        '#','ON','sites','Nvirt','Emiss','BSSE [kcal/mol]'
if(pbc) then
    !Determine supercell
    call criteria(thrR, lat, tau_max)
   ! Loop over all i atoms

   do iat=1,n
      va=xva(iat)
      ea=0.0d0
      np=0
      ! the BSSE due to atom jat, Loop over all j atoms in supercell
      do jat=1,n
      ! # of bf that are available from jat
      vb=xvb(jat)
      if(vb.lt.0.5) cycle
      do a=-tau_max(1),tau_max(1)
      do b=-tau_max(2),tau_max(2)
      do c=-tau_max(3),tau_max(3)
      !remove selfinteraction
      if((iat.eq.jat).and.(abs(a)+abs(b)+abs(c).eq.0)) cycle
         xyzjat(1)=xyz(1,jat)+a*lat(1,1)+b*lat(2,1)+c*lat(3,1)
         xyzjat(2)=xyz(2,jat)+a*lat(1,2)+b*lat(2,2)+c*lat(3,2)
         xyzjat(3)=xyz(3,jat)+a*lat(1,3)+b*lat(2,3)+c*lat(3,3)
!         write(*,*) 'a1 a2 a3',  lat(1,1),lat(1,2),lat(1,3)
!         write(*,*) 'b1 b2 b3',  lat(2,1),lat(2,2),lat(2,3)
!         write(*,*) 'c1 c2 c3',  lat(3,1),lat(3,2),lat(3,3)
!         stop
!         xyzjat(1)=xyz(1,jat)+a*lat(1,1)+b*lat(1,2)+c*lat(1,3)
!         xyzjat(2)=xyz(2,jat)+a*lat(2,1)+b*lat(2,2)+c*lat(2,3)
!         xyzjat(3)=xyz(3,jat)+a*lat(3,1)+b*lat(3,2)+c*lat(3,3)
         vec(1)=(xyz(1,iat)-xyzjat(1))
         vec(2)=(xyz(2,iat)-xyzjat(2))
         vec(3)=(xyz(3,iat)-xyzjat(3))
         r=sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
         ! distance cutoff
         if(r.gt.thrR) cycle
         ! calulate slater overlap sab
         call ssovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),sab)
         ! noise cutoff(sqrt(sab))
         if(sqrt(abs(sab)).lt.thrE) cycle
         ! evaluate gcp central expression
         ene_old_num=exp(-p3*r**p4)
         ene_old_den=sqrt(vb*Sab)
         ene_old=ene_old_num/ene_old_den
         ! noise cutoff(damp)
         if(abs(ene_old).lt.thrE) cycle

         if(damp) then
!D3 r0ab radii
            r0abij=r0ab(iz(iat),iz(jat))
            rscal=r/r0abij
!covalent radii
!          rscal=r/(rad(iz(iat))+rad(iz(jat)))
          rscalexp=rscal**dmp_exp
          dampval=(1d0-1d0/(1d0+dmp_scal*rscalexp))
          ene_dmp=ene_old*dampval
          ea=ea+emiss(iz(iat))*ene_dmp
         else
          ea=ea+emiss(iz(iat))*ene_old
         endif

         ! sites counter (i.e. # atoms contributing to the 'atomic bsse')
         np=np+1

         ! gradient for i,j pair
         if(grad)then
            call cpu_time(t0)
            call gsovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),gab)

            gs(1)=gab*vec(1)
            gs(2)=gab*vec(2)
            gs(3)=gab*vec(3)
            dum=exp(-p3*r**p4)*(-1d0/2d0)
            dum2=2d0*p3*p4*r**p4*sab/r
            dum22=r*sab**(3d0/2d0)
            tmpb=dum22*sqrt(vb)

            if(damp) then
              rscalexpm1=rscal**(dmp_exp-1)
              grd_dmp=dmp_scal*dmp_exp*rscalexpm1/r0abij
              grd_dmp=grd_dmp/((dmp_scal*rscalexp+1.0d0)**2)
            endif

            tmpa=dum2*vec(1)+gs(1)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(vec(1)/r)
            endif
            g(1,iat)=g(1,iat)+tmp*emiss(iz(iat))

            tmpa=dum2*vec(2)+gs(2)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(vec(2)/r)
            endif
            g(2,iat)=g(2,iat)+tmp*emiss(iz(iat))

            tmpa=dum2*vec(3)+gs(3)
            tmp=dum*tmpa/tmpb
            if(damp) then
              tmp=tmp*dampval+ene_old*grd_dmp*(vec(3)/r)
            endif
            g(3,iat)=g(3,iat)+tmp*emiss(iz(iat))

            if(va.lt.0.5) cycle
            if(damp) then
              ene_old_den=sqrt(va*Sab)
              ene_old=ene_old_num/ene_old_den
            endif

            tmpb=dum22*sqrt(va)

            tmpa=dum2*(-vec(1))-gs(1)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(-vec(1)/r)
            endif
            g(1,iat)=g(1,iat)-tmp*emiss(iz(jat))

            tmpa=dum2*(-vec(2))-gs(2)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(-vec(2)/r)
            endif
            g(2,iat)=g(2,iat)-tmp*emiss(iz(jat))

            tmpa=dum2*(-vec(3))-gs(3)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(-vec(3)/r)
            endif
            g(3,iat)=g(3,iat)-tmp*emiss(iz(jat))

            call cpu_time(t1)
            tg=tg+t1-t0
         endif
! end of j-loop
      enddo
!end supercell loop
   enddo
   enddo
   enddo
   if(echo) &
    write(*,'(2x,3(I5,2x),2x,F5.1,2x,F14.4,2x,F14.4,2x,a)')  &
             iat,iz(iat),np,va,emiss(iz(iat)), ea*627.5099*p(1)
    ecp=ecp+ea
! end of i-loop
   enddo
else
   ! Loop over all i atoms
   do iat=1,n
      va=xva(iat)
      ea=0.0d0

      np=0
      ! the BSSE due to atom jat, Loop over all j atoms
      do jat=1,n
         if(iat.eq.jat) cycle
         vec(1)=(xyz(1,iat)-xyz(1,jat))
         vec(2)=(xyz(2,iat)-xyz(2,jat))
         vec(3)=(xyz(3,iat)-xyz(3,jat))
         r=sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))

         ! # of bf that are available from jat
         vb=xvb(jat)
         if(vb.lt.0.5) cycle
         ! distance cutoff
         if(r.gt.thrR) cycle
         ! calulate slater overlap sab
         call ssovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),sab)
         ! noise cutoff(sqrt(sab))
         if(sqrt(abs(sab)).lt.thrE) cycle
         ! evaluate gcp central expression
         ene_old_num=exp(-p3*r**p4)
         ene_old_den=sqrt(vb*Sab)
         ene_old=ene_old_num/ene_old_den
         ! noise cutoff(damp)
         if(abs(ene_old).lt.thrE) cycle
         if(damp) then
!D3 r0ab radii
            r0abij=r0ab(iz(iat),iz(jat))
            rscal=r/r0abij
!covalent radii
!          rscal=r/(rad(iz(iat))+rad(iz(jat)))
          rscalexp=rscal**dmp_exp
          dampval=(1d0-1d0/(1d0+dmp_scal*rscalexp))
          ene_dmp=ene_old*dampval
          ea=ea+emiss(iz(iat))*ene_dmp
         else
          ea=ea+emiss(iz(iat))*ene_old
       endif

         ! sites counter (i.e. # atoms contributing to the 'atomic bsse')
         np=np+1

         ! gradient for i,j pair
         if(grad)then
            call cpu_time(t0)
            call gsovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),gab)

            gs(1)=gab*vec(1)
            gs(2)=gab*vec(2)
            gs(3)=gab*vec(3)
            dum=exp(-p3*r**p4)*(-1d0/2d0)
            dum2=2d0*p3*p4*r**p4*sab/r
            dum22=r*sab**(3d0/2d0)
            tmpb=dum22*sqrt(vb)

            if(damp) then
              rscalexpm1=rscal**(dmp_exp-1)
              grd_dmp=dmp_scal*dmp_exp*rscalexpm1/r0abij
              grd_dmp=grd_dmp/((dmp_scal*rscalexp+1.0d0)**2)
            endif

            tmpa=dum2*vec(1)+gs(1)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(vec(1)/r)
            endif
            g(1,iat)=g(1,iat)+tmp*emiss(iz(iat))

            tmpa=dum2*vec(2)+gs(2)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(vec(2)/r)
            endif
            g(2,iat)=g(2,iat)+tmp*emiss(iz(iat))

            tmpa=dum2*vec(3)+gs(3)
            tmp=dum*tmpa/tmpb
            if(damp) then
              tmp=tmp*dampval+ene_old*grd_dmp*(vec(3)/r)
            endif
            g(3,iat)=g(3,iat)+tmp*emiss(iz(iat))

            if(va.lt.0.5) cycle
            if(damp) then
              ene_old_den=sqrt(va*Sab)
              ene_old=ene_old_num/ene_old_den
            endif

            tmpb=dum22*sqrt(va)

            tmpa=dum2*(-vec(1))-gs(1)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(-vec(1)/r)
            endif
            g(1,iat)=g(1,iat)-tmp*emiss(iz(jat))

            tmpa=dum2*(-vec(2))-gs(2)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(-vec(2)/r)
            endif
            g(2,iat)=g(2,iat)-tmp*emiss(iz(jat))

            tmpa=dum2*(-vec(3))-gs(3)
            tmp=dum*tmpa/tmpb
            if(damp) then
               tmp=tmp*dampval+ene_old*grd_dmp*(-vec(3)/r)
            endif
            g(3,iat)=g(3,iat)-tmp*emiss(iz(jat))

            call cpu_time(t1)
            tg=tg+t1-t0
         endif
! end of j-loop
      enddo
      if(echo) &
       write(*,'(2x,3(I5,2x),2x,F5.1,2x,F14.4,2x,F14.4,2x,a)')  &
                iat,iz(iat),np,va,emiss(iz(iat)), ea*627.5099*p(1)

         ecp=ecp+ea
! end of i-loop
      enddo
endif
gcp=ecp*p(1)
g=g*p(1)


!Special HF-3c correction
if(base) then
 call basegrad(n,max_elem,iz,xyz,lat,pbc,0.7d0,0.03d0,ebas,gbas,echo)
 gcp=gcp+ebas
 if(grad) g=g+gbas
endif

if(grad.and.echo)  call prtim(6,tg,'t','gradient')

end subroutine gcp_egrad



!*********************************
!* subroutine calc stress tensor *
!*********************************
subroutine sgcp(n,max_elem,emiss,xyz,iz,p,gcp,g,grad,echo,xva,xvb,lat,sigma,base,damp,dmp_scal,dmp_exp)
implicit none
integer n,iz(n),max_elem,np
integer iat,jat
real*8  xyz(3,n), abc(3,n), xyzdum(3,n)
real*8  g(3,n)
real*8  gcp,tg,t0,t1
real*8  emiss(max_elem),dum22
real*8  p(4),thrR,thrE
real*8 xva(*),xvb(*)
real*8 va,vb
real*8 r,vec(3)
real*8 tmp,ea,ecp,dum,tmpa,tmpb,dum2
real*8 sab,p1,p2,p3,p4,dmp_scal,dmp_exp
real*8 gs(3),gab
real*8 za(36),zb(36), detlat
real*8 lat(3,3), xyzjat(3,n), sigma(3,3), delta, gcpl, gcpr, vol, latdum(3,3), gdum(3,n), edum
logical echo,grad,damp,base
integer a, b, c, tau_max(3), abcsum, selftest, i, j,k

delta=10E-6

!cell volume
latdum=0
detlat=0
detlat =detlat+lat(1,1)*lat(2,2)*lat(3,3)+lat(1,2)*lat(2,3)*lat(3,1)+lat(1,3)*lat(2,1)*lat(3,2)
detlat=detlat-lat(3,1)*lat(2,2)*lat(1,3)-lat(3,2)*lat(2,3)*lat(1,1)-lat(3,3)*lat(2,1)*lat(1,1)
vol=abs(detlat)


! Transform to direct system
! simultaniously shilft all atoms
latdum=lat
call xyz2abc(n,xyz,abc,lat)
do i=1,3
   do j=1,3
      latdum(i,j)=lat(i,j)-delta
      call abc2xyz(n,xyzdum,abc,latdum)
      call gcp_egrad(n,max_elem,emiss,xyzdum,iz,p,gcpl,gdum,.false.,.false.,xva,xvb,.true.,latdum,damp,base,dmp_scal,dmp_exp)
!       call e(n,max_elem,emiss,xyzdum,iz,p,gcpl,gdum,.false.,.false.,xva,xvb,.true.,latdum,damp,dmp_scal,dmp_exp)
      latdum(i,j)=lat(i,j)+delta
      call abc2xyz(n,xyzdum,abc,latdum)
!       call e(n,max_elem,emiss,xyzdum,iz,p,gcpr,gdum,.false.,.false.,xva,xvb,.true.,latdum,damp,dmp_scal,dmp_exp)
      call gcp_egrad(n,max_elem,emiss,xyzdum,iz,p,gcpr,gdum,.false.,.false.,xva,xvb,.true.,latdum,damp,base,dmp_scal,dmp_exp)
      ! symmetric numerical derivative
      sigma(i,j)=(gcpr-gcpl)/(2*delta)
      !sigma(i,j)=sigma(i,j)/vol
   end do
end do

end subroutine sgcp

!*****************
!* print timings *
!*****************
subroutine prtim(io,tt,is,string)
integer io
real*8 tt,t,sec
integer day,hour,min
character*(*) string
character*1 is

t=tt
day=idint(t/86400)
t=t-day*86400
hour=idint(t/3600)
t=t-hour*3600
min=idint(t/60)
t=t-60*min
sec=t

if(day.ne.0)then
   if(is.eq.'w')write(io,2)trim(string),day,hour,min,sec
   if(is.eq.'t')write(io,22)trim(string),day,hour,min,sec
   return
endif
if(hour.ne.0)then
   if(is.eq.'w')write(io,3)trim(string),hour,min,sec
   if(is.eq.'t')write(io,33)trim(string),hour,min,sec
   return
endif
if(min .ne.0)then
   if(is.eq.'w')write(io,4)trim(string),min,sec
   if(is.eq.'t')write(io,44)trim(string),min,sec
   return
endif
  if(is.eq.'w')write(io,5)trim(string),sec
  if(is.eq.'t')write(io,55)trim(string),sec
return

 1    format('======================================')
 2    format('wall-time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 3    format('wall-time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 4    format('wall-time for ',a,2x,i3,' m',f5.1,' s')
 5    format('wall-time for ',a,2x,f5.1,' s')

 22    format('cpu-time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 33    format('cpu-time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 44    format('cpu-time for ',a,2x,i3,' m',f5.1,' s')
 55    format('cpu-time for ',a,2x,f5.1,' s')

return
end

!*********************************
!* Load all necessary parameters *
!*********************************
subroutine setparam(emiss,nbas,p,method)
implicit none
integer, parameter :: mpar=86 ! maximal
integer, parameter :: apar=36 ! actual
integer nbas(mpar)
character*(*) method
real(8) emiss(mpar),p(*)
real(8) HFsv(apar),HFminis(apar),HF631gd(apar),HFsvp(apar),HFtz(apar),&
     HFvmb(apar),HFminisd(apar),oldHFsvp(apar),HFpobtz(apar),HFpobdzvp(apar),&
     HF2gcore(apar),HF2g(apar), HFdef1tzvp(apar),HFccdz(apar),HFaccdz(apar),&
     HFdzp(apar),HFhsv(apar),HFdz(apar),HFmsvp(apar),HFdef2mtzvp(apar),def2mtzvpp(apar) !SG
integer BASsv(apar),BASminis(apar),BAS631gd(apar),BAStz(apar),&
     BASsvp(apar),BASvmb(apar),BASminisd(apar),oldBASsvp(apar),BASpobtz(apar),BASpobdzvp(apar),&
     BAS2gcore(apar),BAS2g(apar),BASdef1tzvp(apar),BASccdz(apar),BASaccdz(apar),&
     BASdzp(apar),BAShsv(apar),BASdz(apar),BASmsvp(apar),BASdef2mtzvp(apar)
real(8) HFlanl2(10)
integer BASlanl2(10)

! ***************
! *  Emiss data *
! ***************
data HFsv/ & !H-Kr (no 3d)
0.009037,0.008843,&  ! He,He
0.204189,0.107747,0.049530,0.055482,0.072823,0.100847,0.134029,0.174222,&  ! Li-Ne
0.315616,0.261123,0.168568,0.152287,0.146909,0.168248,0.187882,0.211160,&  !Na -Ar
0.374252,0.460972,&  ! K-Ca
0.444886,0.404993,0.378406,0.373439,0.361245,0.360014,0.362928,0.243801,0.405299,0.396510,&   ! 3d-TM
0.362671,0.360457,0.363355,0.384170,0.399698,0.417307/ !Ga-Kr

data HFminis/ &
0.042400,0.028324,&
0.252661,0.197201,0.224237,0.279950,0.357906,0.479012,0.638518,0.832349, &
1.232920,1.343390,1.448280,1.613360,1.768140,1.992010,2.233110,2.493230, &
3.029640,3.389980,&  ! H-Ca
10*0,6*0/

data HF631GD/ &! H-Ca + Br (no 3d)
0.010083,0.008147,&
0.069260,0.030540,0.032736,0.021407,0.024248,0.036746,0.052733,0.075120,&
0.104255,0.070691,0.100260,0.072534,0.054099,0.056408,0.056025,0.057578,&
0.079198,0.161462,&
10*0.0, &
0.000000,0.000000,0.000000,0.000000,0.381049,0.000000/

data oldHFsvp / & ! Li,Na,Mg,K had wrong parameters
0.008107,0.008045,&
0.142957,0.028371,0.049369,0.055376,0.072785,0.100310,0.133273,0.173600,&
0.191109,0.222839,0.167188,0.149843,0.145396,0.164308,0.182990,0.205668,&
0.221189,0.299661,&
0.325995,0.305488,0.291723,0.293801,0.29179,0.296729,0.304603,0.242041,0.354186,0.350715,&
0.350021,0.345779,0.349532,0.367305,0.382008,0.399709/

data HFsvp /  & ! H-Kr
0.008107,0.008045,&
0.113583,0.028371,0.049369,0.055376,0.072785,0.100310,0.133273,0.173600,&
0.181140,0.125558,0.167188,0.149843,0.145396,0.164308,0.182990,0.205668,&
0.200956,0.299661, &
0.325995,0.305488,0.291723,0.293801,0.29179,0.296729,0.304603,0.242041,0.354186,0.350715,&
0.350021,0.345779,0.349532,0.367305,0.382008,0.399709/

data HFtz /&  ! H-Kr !def2-TZVP
0.007577,0.003312,&
0.086763,0.009962,0.013964,0.005997,0.004731,0.005687,0.006367,0.007511,&
0.077721,0.050003,0.068317,0.041830,0.025796,0.025512,0.023345,0.022734,&
0.097241,0.099167,&
0.219194,0.189098,0.164378,0.147238,0.137298,0.12751,0.118589,0.0318653,0.120985,0.0568313, &
0.090996,0.071820,0.063562,0.064241,0.061848,0.061021/

data HFdef2mtzvp /&  ! m def2-TZVP, no f for B-Ne
0.007930,0.003310,&
0.086760,0.009960,0.013960,0.006000,0.003760,0.004430,0.005380,0.006750,&
0.077720,0.050000,0.068320,0.041830,0.025800,0.025510,0.023340,0.022730,&
0.097240,0.099170,&
0.219190,0.189100,0.164380,0.147240,0.137300,0.127510,0.118590,0.031870,0.120990,0.056830,&
0.091000,0.071820,0.063560,0.064240,0.061850,0.061020/

data HFvmb/&
0.042400,0.028324,&
0.252661,0.197201,0.156009,0.164586,0.169273,0.214704,0.729138,0.336072,&
0.262329,0.231722,0.251169,0.287795,0.213590,0.250524,0.728579,0.260658, &
2*0,&
16*0/

data HFminisd/& !Al-Ar MINIS + Ahlrichs "P" funktions
0.0,0.0,&
8*0.0,&
2*0.0,1.446950,1.610980,1.766610,1.988230,2.228450,2.487960,&
2*0.0,&
16*0.0/

data HFlanl2/ & !  LANL2TZ+ vs LANL2DZ (ORCA), only Sc-Zn
0.102545,0.0719529,0.0491798,0.0362976,0.0266369,0.0235484,0.0171578,0.0438906,0.0100259,0.016208/

data HFpobtz /  & ! H-Kr, no RG
0.010077,0.000000,&
0.173239,0.101973,0.131181,0.032234,0.011630,0.008475,0.011673,0.000000,&
0.240653,0.661819,0.522306,0.14163,0.052456,0.184547,0.040837,0.000000,&
0.318136,0.564721,&
0.523194,0.767449,0.620122,0.390227,0.237571,0.247940,0.249589,0.117864,0.325725,0.197183,&
0.264489,0.180375,0.111262,0.147239,0.081747,0.000000/


data HF2gcore /  & ! only HCNOF yet
0.000539,0.000000,&
0.000000,0.000000,0.000000,0.173663,0.269952,0.364341,0.384923,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
0.000000,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000/

data HFdef1tzvp /&  !  org
0.007577,0.003312,&
0.136371,0.011163,0.017129,0.008140,0.005826,0.006777,0.007108,0.008132,&
0.134992,0.147417,0.085253,0.054238,0.033790,0.032862,0.029038,0.026555,&
0.141595,0.207980,&
0.223252,0.193038,0.167892,0.148726,0.140473,0.130220,0.121166,0.113839,0.121855,0.107138,&
0.105637,0.086639,0.075084,0.075089,0.070868,0.068706/

data def2mtzvpp /&    !SG
0.027000,0.000000,&
0.000000,0.000000,0.200000,0.020000,0.180000,0.080000,0.070000,0.065000,&
0.000000,0.000000,0.000000,0.200000,0.600000,0.600000,0.600000,0.300000,&
0.000000,0.000000,&
10*0.3,&
0.300000,0.300000,0.300000,0.300000,0.300000,0.000000/

data HF2g / & !no ne, ar ecp
0.0539181,0.161846,&
0.1581960,0.214318,0.808955,0.470398,0.724457,1.260960,2.014430,0.000000,&
0.3072290,0.265975,0.354139,0.307818,0.356962,0.461661,0.588346,0.000000,&
0.000000,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000/

data HFccdz / &
0.007907,0.008287,&
0.047380,0.014240,0.022133,0.014999,0.018148,0.028240,0.042261,0.061485,&
0.073185,0.056218,0.082660,0.052975,0.033874,0.034056,0.031433,0.030034,&
0.000000,0.078016,& !no k cc-pVDZ Basis
0.036885,0.038540,0.036474,0.036061,0.030289,0.027959,0.025177,0.022709,0.027386,0.015816,&
0.135176,0.115515,0.102761,0.102967,0.097891,0.097331/


data HFaccdz / & !for li,be,na,mg,k-zn energy below def2-QZVPD reference
0.001183,0.005948,&
0.000000,0.000000,0.005269,0.006380,0.011700,0.021199,0.034160,0.051481,&
0.000000,0.000000,0.016018,0.009268,0.010076,0.015153,0.016889,0.018563,&
0.000000,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
0.069963,0.065687,0.072944,0.077585,0.078777,0.080746/

data HFdzp / &
0.008107,0.008045,&
0.136751,0.016929,0.026729,0.021682,0.027391,0.040841,0.058747,0.082680,&
0.153286,0.162296,0.102704,0.073144,0.056217,0.061333,0.065045,0.071398,&
0.145642,0.212865,&
0.232821,0.204796,0.182933,0.169554,0.164701,0.160112,0.157723,0.158037,0.179104,0.169782,&
0.159396,0.140611,0.129645,0.132664,0.132121,0.134081/

data HFhsv / &
0.030224,0.028324,&
0.125379,0.064094,0.059751,0.079387,0.108929,0.167264,0.245786,0.347818,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
0.000000,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
0.000000,0.000000,0.000000,0.000000,0.000000,0.000000/

data HFdz / &
0.009037,0.008843,&
0.198254,0.000000,0.026921,0.021817,0.027458,0.041391,0.059495,0.083286,&
0.268608,0.202374,0.104146,0.075686,0.057826,0.065300,0.069912,0.076845,&
0.296046,0.370399,&
0.349482,0.302284,0.267639,0.244306,0.232237,0.221488,0.214153,0.032694,0.226865,0.213902,&
0.172296,0.155496,0.143646,0.149642,0.149871,0.151705/


data HFmsvp / &     !H-Kr modified Ahlrichs DZ, supplemented by def2-SV(P)
0.000000d0,0.000000d0,& !RG,H set to zero,  F adjusted empirically, Be corrected due to ROHF problems
0.107750d0,0.020000d0,0.026850d0,0.021740d0,0.027250d0,0.039930d0,0.030000d0,0.000000d0,&
0.153290d0,0.162300d0,0.102700d0,0.073140d0,0.056220d0,0.061330d0,0.065040d0,0.000000d0,&
0.200960d0,0.299660d0,&
0.325990d0,0.305490d0,0.291720d0,0.293800d0,0.291790d0,0.296730d0,0.304600d0,0.242040d0,0.354190d0,0.350720d0,&
0.350020d0,0.345780d0,0.349530d0,0.367310d0,0.382010d0,0.000000d0/


! *********************
! * nr. of basis fkt  *
! *********************
data BASsv/2*2,2*3,6*9,2*7,6*13,2*11,10*21,6*27/
data BASminis/2*1,2*2,6*5,2*6,6*9,2*10,16*0/
data BAS631gd/2,5,8*14,8*18,2*22,10*0,4*0,32,0/
data BASsvp/2*5,9,9,6*14,15,18,6*18,24,24,10*31,6*32/
data oldBASsvp/2*5,6,9,6*14,2*10,6*18,14,24,10*31,6*32/
data BAStz/       2*6,14,19,6*31,2*32,6*37,33,36,9*45,48,6*48/
!data BASdef2mtzvp/2*6,14,19,6*24,2*32,6*37,33,36,9*45,48,6*48/   !def2-TZVP, no f for B-Ne
data BASdef2mtzvp/2*3,8,11,3*19,24,2*19,2*14,6*22,18,28,10*31,6*36/   !def2-mTZVP
data BASvmb/2*1,2*2,6*4,2*1,6*4,2*0,16*0/ ! minimal basis set with ECPs
data BASminisd/2*0,2*0,6*0,2*0,6*14,2*0,16*0/
data BASlanl2/22, 22, 22, 22, 22, 22, 22, 22, 22, 18/ ! Sc-Zn LANL2DZ
data BASpobtz / & ! H-Kr no RG
6,0,&
7,7,18,18,18,18,18,0,&
19,19,22,22,22,22,22,0,&
23,23,&
40,40,40,40,40,40,40,40,40,40,&
41,41,41,41,41,0 /

data BASpobdzvp / & ! H-KR, no RG
5,0,&
6,6,14,14,14,14,14,0,&
15,18,18,18,18,18,18,0,&
19,19,&
31,31,31,31,31,31,31,31,31,31,&
32,32,32,32,32,0/

data BAS2gcore / & ! Only HCNOF yet
1,0,&
5,5,5,5,5,5,5,5,&
9,9,14,14,14,14,14,14,&
0,0,&
0,0,0,0,0,0,0,0,0,0,&
0,0,0,0,0,0/

data BAS2g / &
1,1,&
5,5,5,5,5,5,5,5,&
9,9,14,14,14,14,14,14,&
0,0,&
0,0,0,0,0,0,0,0,0,0,&
0,0,0,0,0,0/

data BASdef1tzvp / &
6,6,&
8,11,19,19,19,19,19,19,&
14,14,22,22,22,22,22,22,&
33,33,33,33,33,33,33,33,33,33,&
18,28,&
36,36,36,36,36,36/

data BASccdz / &
5,5,&
14,14,14,14,14,14,14,14,&
18,18,18,18,18,18,18,18,&
0,27,&
43,43,43,43,43,43,43,43,43,43,&
27,27,27,27,27,27/

data BASaccdz / &
9,9,&
23,23,23,23,23,23,23,23,&
27,27,27,27,27,27,27,27,&
0,0,&
59,59,59,59,59,59,59,59,59,59,&
36,36,36,36,36,36/

data BASdzp / &
5,5,&
7,10,15,15,15,15,15,15,&
15,15,23,23,23,23,23,23,&
26,36,&
41,41,41,41,41,41,41,41,41,41,&
41,41,41,41,41,41/

data BASdz / &
2,2,&
4,0,10,10,10,10,10,10,&
12,12,18,18,18,18,18,18,&
23,23,&
38,38,38,38,38,38,38,38,38,38,&
36,36,36,36,36,36/

data BASmsvp / &  ! modified Ahlrichs DZ, supplemented by def2-SV(P)
2,2,&
10,10,15,15,15,15,15,15,&
15,18,18,18,18,18,18,18,&
24,24,&
31,31,31,31,31,31,31,31,31,31,&
32,32,32,32,32,32/

! **************************************
! * load data into emiss() and nbas()  *
! **************************************

emiss=0d0
nbas =0d0

select case (method)
!*****************
!* Hartree-Fock  *
!*****************
  case ('hf/sv') ! RMS=0.3218975
     emiss(1:apar)=HFsv(1:apar)
     nbas(1:apar)=BASsv(1:apar)
     p(1)=0.1724d0
     p(2)=1.2804d0
     p(3)=0.8568d0
     p(4)=1.2342d0
  case ('hf/svp', 'hf/def2svp') ! RMS=0.4065  ! fit does not include Li,Na,Mg,K, so not change here
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     p(1)=0.2054d0
     p(2)=1.3157d0
     p(3)=0.8136d0
     p(4)=1.2572d0
  case ('hf/svp_old') ! RMS=0.4065
     emiss(1:apar)=oldHFsvp(1:apar)
     nbas(1:apar)=oldBASsvp(1:apar)
     p(1)=0.2054d0
     p(2)=1.3157d0
     p(3)=0.8136d0
     p(4)=1.2572d0
  case('hf/sv(p)', 'hf/def2sv(p)','hf/sv_p') !RMS=0.3502
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     p(1)=0.1373d0
     p(2)=1.4271d0
     p(3)=0.8141d0
     p(4)=1.2760d0
  case('hf/dzp') !RMS=0.4571
     emiss(1:apar)=HFdzp(1:apar)
     nbas(1:apar)=BASdzp(1:apar)
     p(1)=0.1443d0
     p(2)=1.4547d0
     p(3)=0.3711d0
     p(4)=1.6300d0
  case ('hf/631gd', 'hf/631gs') ! RMS= 0.40476
     emiss(1:apar)=HF631gd(1:apar)
     nbas(1:apar)=BAS631gd(1:apar)
     p(1)=0.2048d0
     p(2)=1.5652d0
     p(3)=0.9447d0
     p(4)=1.2100d0
  case ('hf/minis') ! RMS= 0.3040
    emiss(1:apar)=HFminis(1:apar)
    nbas(1:apar)=BASminis(1:apar)
     p(1)=0.1290d0
     p(2)=1.1526d0
     p(3)=1.1549d0
     p(4)=1.1763d0
  case ('hf/minix', 'hf3c')
!H-Mg: MINIS
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
!Al-Ar: MINIS+d (stan. TM exp.)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
!K-Zn: SV
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
!Ga-Kr: SVP
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
! Li,Be,Na,Mg MINIS + p-fkt
    emiss(3)=0.177871
    emiss(4)=0.171596
    nbas(3:4)=5
    emiss(11)=1.114110
    emiss(12)=1.271150
    nbas(11:12)=9
! fit param of HF/MINIS
     p(1)=0.1290d0
     p(2)=1.1526d0
     p(3)=1.1549d0
     p(4)=1.1763d0
  case ('hf/tz', 'hf/def2tzvp') !  RMS= 0.1150
     emiss(1:apar)=HFtz(1:apar)
     nbas(1:apar)=BAStz(1:apar)
     p(1)=0.3127d0
     p(2)=1.9914d0
     p(3)=1.0216d0
     p(4)=1.2833d0
  case ('hf/deftzvp', 'hf/tzvp') ! RMS=0.209
     emiss(1:apar)=HFdef1tzvp(1:apar)
     nbas(1:apar)=BASdef1tzvp(1:apar)
     p(1)=0.2600d0
     p(2)=2.2448d0
     p(3)=0.7998d0
     p(4)=1.4381d0
  case ('hf/ccdz', 'hf/ccpvdz') ! RMS=0.4968
     emiss(1:apar)=HFccdz(1:apar)
     nbas(1:apar)=BASccdz(1:apar)
     p(1)=0.4416d0
     p(2)=1.5185d0
     p(3)=0.6902d0
     p(4)=1.3713d0
case ('hf/accdz','hf/augccpvdz') !RMS=0.2222
     emiss(1:apar)=HFaccdz(1:apar)
     nbas(1:apar)=BASaccdz(1:apar)
     p(1)=0.0748d0
     p(2)=0.0663d0
     p(3)=0.3811d0
     p(4)=1.0155d0
  case ('hf/2g')
     emiss(1:apar)=HF2g(1:apar)
     nbas(1:apar)=BAS2g(1:apar)
     p(1)=0.2461d0
     p(2)=1.1616d0
     p(3)=0.7335d0
     p(4)=1.4709d0
  case ('hf/dz')  !RMS=0.3754
     emiss(1:apar)=HFdz(1:apar)
     nbas(1:apar)=BASdz(1:apar)
     p(1)=0.1059d0
     p(2)=1.4554d0
     p(3)=0.3711d0
     p(4)=1.6342d0

!************************
!* KS-DFT hybrid(B3LYP) *
!************************
  case ('dft/lanl','b3lyp/lanl')
    emiss(1:apar)=HF631gd(1:apar)
    nbas(1:apar)=BAS631gd(1:apar)
     p(1)=0.3405d0
     p(2)=1.6127d0
     p(3)=0.8589d0
     p(4)=1.2830d0
     emiss(21:30)=HFlanl2(1:10)
     nbas(21:30)=BASlanl2(1:10)
  case ('dft/sv','b3lyp/sv') ! RMS= 0.557
     emiss(1:apar)=HFsv(1:apar)
     nbas(1:apar)=BASsv(1:apar)
     p(1)=0.4048d0
     p(2)=1.1626d0
     p(3)=0.8652d0
     p(4)=1.2375d0
  case ('dft/sv(p)','b3lyp/sv(p)','b3lyp/def2sv(p)','pbe/def2sv(p)','dft/def2sv(p)','dft/sv_p') ! RMS= 0.57 ! def2-SV(P)
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     p(1)=0.2424d0
     p(2)=1.2371d0
     p(3)=0.6076d0
     p(4)=1.4078d0
  case ('dft/svx','b3lyp/svx') ! RMS=  0.56 ! def2-SV(P/h,c)  = SV at h,c
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     emiss(6)=HFsv(6)
     nbas(6)=BASsv(6)
     p(1)=0.1861d0
     p(2)=1.3200d0
     p(3)=0.6171d0
     p(4)=1.4019d0
  case ('dft/svp','b3lyp/svp', 'dft/def2svp', 'b3lyp/def2svp') ! RMS=0.6498
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     p(1)=0.2990d0
     p(2)=1.2605d0
     p(3)=0.6438d0
     p(4)=1.3694d0
  case ('dft/svp_old','b3lyp/svp_old') ! RMS=0.6498
     emiss(1:apar)=oldHFsvp(1:apar)
     nbas(1:apar)=oldBASsvp(1:apar)
     p(1)=0.2990d0
     p(2)=1.2605d0
     p(3)=0.6438d0
     p(4)=1.3694d0
  case('dft/dzp','b3lyp/dzp') !RMS=0.7184
     emiss(1:apar)=HFdzp(1:apar)
     nbas(1:apar)=BASdzp(1:apar)
     p(1)=0.2687d0
     p(2)=1.4634d0
     p(3)=0.3513d0
     p(4)=1.6880d0
  case ('dft/631gd','b3lyp/631gd', 'dft/631gs', 'b3lyp/631gs') ! RMS=  0.47856
    emiss(1:apar)=HF631gd(1:apar)
    nbas(1:apar)=BAS631gd(1:apar)
     p(1)=0.3405d0
     p(2)=1.6127d0
     p(3)=0.8589d0
     p(4)=1.2830d0
  case ('dft/minis','b3lyp/minis') ! RMS= 0.3400
    emiss(1:apar)=HFminis(1:apar)
    nbas(1:apar)=BASminis(1:apar)
     p(1)=0.2059d0
     p(2)=0.9722d0
     p(3)=1.1961d0
     p(4)=1.1456d0
  case('dft/minix','b3lyp/minix')
!H-Mg: MINIS
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
!Al-Ar: MINIS+d (stan. TM exp.)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
!K-Zn: SV
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
    !Ga-Kr: SVP
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
! Li,Be,Na,Mg MINIS + p-fkt
    emiss(3)=0.177871
    emiss(4)=0.171596
    nbas(3:4)=5
    emiss(11)=1.114110
    emiss(12)=1.271150
    nbas(11:12)=9
! fit param of DFT/MINIS
    p(1)=0.2059d0
    p(2)=0.9722d0
    p(3)=1.1961d0
    p(4)=1.1456d0
  case ('dft/tz','b3lyp/tz', 'dft/def2tzvp', 'b3lyp/def2tzvp') ! RMS=0.19648
     emiss(1:apar)=HFtz(1:apar)
     nbas(1:apar)=BAStz(1:apar)
     p(1)=0.2905d0
     p(2)=2.2495d0
     p(3)=0.8120d0
     p(4)=1.4412d0
  case ('dft/deftzvp','b3lyp/deftzvp', 'dft/tzvp', 'b3lyp/tzvp') ! RMS=0.1817
     emiss(1:apar)=HFdef1tzvp(1:apar)
     nbas(1:apar)=BASdef1tzvp(1:apar)
     p(1)=0.2393d0
     p(2)=2.2247d0
     p(3)=0.8185d0
     p(4)=1.4298d0
  case ('dft/ccdz','dft/ccpvdz','b3lyp/ccdz','b3lyp/ccpvdz') !RMS=0.7610
     emiss(1:apar)=HFccdz(1:apar)
     nbas(1:apar)=BASccdz(1:apar)
     p(1)=0.5383d0
     p(2)=1.6482d0
     p(3)=0.6230d0
     p(4)=1.4523d0
  case ('dft/accdz','dft/augccpvdz','b3lyp/accdz','b3lyp/augccpvdz') !RMS=0.1840
     emiss(1:apar)=HFaccdz(1:apar)
     nbas(1:apar)=BASaccdz(1:apar)
     p(1)=0.1465d0
     p(2)=0.0500d0
     p(3)=0.6003d0
     p(4)=0.8761d0
  case ('dft/pobtz','b3lyp/pobtz') ! RMS=
     emiss(1:apar)=HFpobtz(1:apar)
     nbas(1:apar)=BASpobtz(1:apar)
     p(1)=0.1300d0
     p(2)=1.3743d0
     p(3)=0.4792d0
     p(4)=1.3962d0
  case('dft/dz','b3lyp/dz')
     emiss(1:apar)=HFdz(1:apar)
     nbas(1:apar)=BASdz(1:apar)
     p(1)=0.2687d0
     p(2)=1.4634d0
     p(3)=0.3513d0
     p(4)=1.6880d0


!*****************
!* special cases *
!*****************
  case ('blyp/minis','gga/minis') ! RMS= 0.3462
  emiss(1:apar)=HFminis(1:apar)
  nbas(1:apar)=BASminis(1:apar)
  p(1)=0.1566d0
  p(2)=1.0271d0
  p(3)=1.0732d0
  p(4)=1.1968d0
  case ('tpss/minis') ! RMS=
  emiss(1:apar)=HFminis(1:apar)
  nbas(1:apar)=BASminis(1:apar)
  p(1)=0.22982d0
  p(2)=1.35401d0
  p(3)=1.47633d0
  p(4)=1.11300d0
  case ('tpss/svp') ! RMS=  0.618
  emiss(1:apar)=HFsvp(1:apar)
  nbas(1:apar)=BASsvp(1:apar)
  p(1)=0.6647d0
  p(2)=1.3306d0
  p(3)=1.0792d0
  p(4)=1.1651d0
  case ('gga/svp','blyp/svp') ! RMS=
  emiss(1:apar)=HFsvp(1:apar)
  nbas(1:apar)=BASsvp(1:apar)
  p(1)=0.6823d0
  p(2)=1.2491d0
  p(3)=0.8225d0
  p(4)=1.2811d0
  case ('blyp/sv','gga/sv') ! RMS = 0.6652
  emiss(1:apar)=HFsv(1:apar)
  nbas(1:apar)=BASsv(1:apar)
  p(1)=0.2727d0
  p(2)=1.4022d0
  p(3)=0.8055d0
  p(4)=1.3000d0
  case ('blyp/tz','gga/tz') !RMS = 0.21408
  emiss(1:apar)=HFtz(1:apar)
  nbas(1:apar)=BAStz(1:apar)
  p(1)=0.1182d0
  p(2)=1.0631d0
  p(3)=1.0510d0
  p(4)=1.1287d0
  case ('pw6b95/minis')  ! RMS = 0.3279929
  emiss(1:apar)=HFminis(1:apar)
  nbas(1:apar)=BASminis(1:apar)
  p(1)=0.21054d0
  p(2)=1.25458d0
  p(3)=1.35003d0
  p(4)=1.14061d0
  case ('pw6b95/svp')  ! RMS = 0.58312
  emiss(1:apar)=HFsvp(1:apar)
  nbas(1:apar)=BASsvp(1:apar)
  p(1)=0.3098d0
  p(2)=1.2373d0
  p(3)=0.6896d0
  p(4)=1.3347d0
  case ('pbeh3c', 'pbeh3c/msvp')
  emiss(1:apar)=HFmsvp(1:apar)
  emiss(19:apar)=HFdzp(19:apar)
  emiss(36)=0.0d0
  nbas(1:apar)=BASmsvp(1:apar)
  p(1)=1.00000d0
  p(2)=1.32492d0
  p(3)=0.27649d0
  p(4)=1.95600d0
  case ('hse3c', 'hse3c/msvp')
  emiss(1:apar)=HFmsvp(1:apar)
  emiss(19:apar)=HFdzp(19:apar)
  emiss(36)=0.0d0
  nbas(1:apar)=BASmsvp(1:apar)
  p(1)=1.00000d0
  p(2)=1.32378d0
  p(3)=0.28314d0
  p(4)=1.94527d0
  case('b3pbe3c')
  emiss(1:apar)=HFdef2mtzvp(1:apar)
  nbas(1:apar)=BASdef2mtzvp(1:apar)
  p(1)=1.0000d0
  p(2)=2.98561d0
  p(3)=0.3011d0
  p(4)=2.4405d0

!***********************
!* load just emiss/nbf *
!**********************
case('sv')
     emiss(1:apar)=HFsv(1:apar)
     nbas(1:apar)=BASsv(1:apar)
case('sv(p)','def2sv(p)')
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=0.009037d0
     nbas(1)=2
case ('svx') ! RMS=  ! def2-SV(P/h,c)  = SV at h,c
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     emiss(6)=HFsv(6)
     nbas(6)=BASsv(6)
case('svp')
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
case('minis')
     emiss(1:apar)=HFminis(1:apar)
     nbas(1:apar)=BASminis(1:apar)
case('631gd')
     emiss(1:apar)=HF631gd(1:apar)
     nbas(1:apar)=BAS631gd(1:apar)
case('tz')
     emiss(1:apar)=HFtz(1:apar)
     nbas(1:apar)=BAStz(1:apar)
case('deftzvp')
     emiss(1:apar)=HFdef1tzvp(1:apar)
     nbas(1:apar)=BASdef1tzvp(1:apar)
case('ccdz')
     emiss(1:apar)=HFccdz(1:apar)
     nbas(1:apar)=BASccdz(1:apar)
case('accdz')
     emiss(1:apar)=HFaccdz(1:apar)
     nbas(1:apar)=BASaccdz(1:apar)
case('pobtz')
   emiss(1:apar)=HFpobtz(1:apar)
   nbas(1:apar)=BASpobtz(1:apar)
case('pobdzvp')
   emiss(1:apar)=HFpobdzvp(1:apar)
   nbas(1:apar)=BASpobdzvp(1:apar)
case ('minix')
!H-Mg: MINIS
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
!Al-Ar: MINIS+d (stan. TM exp.)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
!K-Zn: SV
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
!Ga-Kr: SVP
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
! Li,Be,Na,Mg MINIS + p-fkt
    emiss(3)=0.177871
    emiss(4)=0.171596
    nbas(3:4)=5
    emiss(11)=1.114110
    emiss(12)=1.271150
    nbas(11:12)=9
case('gcore')
   emiss(1:apar)=HF2gcore(1:apar)
   nbas(1:apar)=BAS2gcore(1:apar)
case('twog')
   emiss(1:apar)=HF2g(1:apar)
   nbas(1:apar)=BAS2g(1:apar)
case('fitg')
   emiss(1:apar)=HF2g(1:apar)
   nbas(1:apar)=BAS2g(1:apar)
case('dzp')
   emiss(1:apar)=HFdzp(1:apar)
   nbas(1:apar)=BASdzp(1:apar)
case('hsv')
   emiss(1:apar)=HFhsv(1:apar)
   nbas(1:apar)=BAShsv(1:apar)
case('dz')
   emiss(1:apar)=HFdz(1:apar)
   nbas(1:apar)=BASdz(1:apar)
case('msvp')
  emiss(1:apar)=HFmsvp(1:apar)
  nbas(1:apar)=BASmsvp(1:apar)
case('def2mtzvp')
  emiss(1:apar)=HFdef2mtzvp(1:apar)
  nbas(1:apar)=BASdef2mtzvp(1:apar)
case('def2mtzvpp')  !SG
  emiss(1:apar)=def2mtzvpp(1:apar)
  nbas(1:apar)=1
  p(1)=1.0000d0
  p(2)=1.3150d0
  p(3)=0.9410d0
  p(4)=1.4636d0
case default
  write(*,'(3a)')  '** ',trim(method),' **'
  error stop 'not implemented'
end select

end



!******************************************************************************
!* calculates the s-type overlap integral over 1s, 2s and 3s slater functions
!* added support for 3s functions
!* ovl = overlap integral
!* za  = slater exponent atom A
!* zb  = slater exponent atom B
!* R   = distance between atom A and B
!* Inspired by mopac7.0
!******************************************************************************
subroutine ssovl(r,iat,jat,iz,xza,xzb,ovl)
implicit none
integer ii,shell(72)
logical debug
real(8) za,zb,R,ovl,ax,bx,norm,R05
integer na,nb
real(8) Bxx0,Bxx1,Bxx2,xx,Bxx4,Bxx6
real(8) Bxx3,Bxx5
data shell/                 &
!          h, he
          1,1               &
!         li-ne
          ,2,2,2,2,2,2,2,2, &
!         na-ar
          3,3,3,3,3,3,3,3,  &
! 4s,5s will be treated as 3s
!         k-rn , no f-elements
          54*3/
! ...
real(kind=8) xza,xzb
integer iat,jat,iz(*)

       za=xza
       zb=xzb
       na=iz(iat)
       nb=iz(jat)
debug=.false.
!debug=.true.

! ii selects kind of ovl by multiplying the shell
! kind    <1s|1s>  <2s|1s>  <2s|2s>  <1s|3s>  <2s|3s>  <3s|3s>
! case:      1        2        4       3        6         9
!
ii=shell(na)*shell(nb)
if(debug) write(*,*) 'shell', ii

R05=R*0.5
ax=(za+zb)*R05
bx=(zb-za)*R05

! same elements
if(za.eq.zb.OR.abs(za-zb).lt.0.1) then
  select case (ii)
   case (1)
    ovl=0.25d0*sqrt((za*zb*R*R)**3)*(A2(ax)*Bint(bx,0)-Bint(bx,2)*A0(ax))
   case (2)
    ovl = SQRT(1.D00/3.D00)
    if(shell(na).lt.shell(nb)) then
    ! <1s|2s>
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00
      ovl=ovl*norm*(A3(ax)*Bint(bx,0)-Bint(bx,3)*A0(ax)+A2(ax)*Bint(bx,1)-Bint(bx,2)*A1(ax))
     else
    ! switch za/zb to get <2s|1s>
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00
      ovl=ovl*norm*(A3(ax)*Bint(bx,0)-Bint(bx,3)*A0(ax)+A2(ax)*Bint(bx,1)-Bint(bx,2)*A1(ax))
    endif
   case (4)
    norm=SQRT((ZA*ZB)**5)*(R**5)*0.0625d0
    ovl=norm* (A4(ax)*Bint(bx,0)+Bint(bx,4)*A0(ax)-2.0d0*A2(ax)*Bint(bx,2))*(1d0/3d0)
   case(3)
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*Bint(bx,0)-Bint(bx,4)*A0(ax)+2.d0*(A3(ax)*Bint(bx,1)-Bint(bx,3)*A1(ax)))/sqrt(3.d0)
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*Bint(bx,0)-Bint(bx,4)*A0(ax)+2.d0*(A3(ax)*Bint(bx,1)-Bint(bx,3)*A1(ax)))/sqrt(3.d0)
    endif
   case(6)
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*Bint(bx,0)+A4(ax)*Bint(bx,1) &
         & -2d0*(A3(ax)*Bint(bx,2)+A2(ax)*Bint(bx,3)) &
         & +A1(ax)*Bint(bx,4)+A0(ax)*Bint(bx,5))/3.d0
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*Bint(bx,0)+A4(ax)*Bint(bx,1) &
         & -2d0*(A3(ax)*Bint(bx,2)+A2(ax)*Bint(bx,3)) &
         & +A1(ax)*Bint(bx,4)+A0(ax)*Bint(bx,5))/3.d0
    endif
   case(9)
      norm=sqrt((ZA*ZB*R*R)**7)/480.d0
      ovl=norm*(A6(ax)*Bint(bx,0)-3.d0*(A4(ax)*Bint(bx,2) &
         & -A2(ax)*Bint(bx,4))-A0(ax)*Bint(bx,6))/3.D00
   end select
else ! different elements
   select case (ii)
   case (1)
      norm=0.25d0*sqrt((za*zb*R*R)**3)
      ovl=(A2(ax)*B0(bx)-B2(bx)*A0(ax))*norm
   case (2)
      ovl = SQRT(1.D00/3.D00)
    if(shell(na).lt.shell(nb)) then
    ! <1s|2s>
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00
      ovl=ovl*norm*(A3(ax)*B0(bx)-B3(bx)*A0(ax)+A2(ax)*B1(bx)-B2(bx)*A1(ax))
     else
    ! switch za/zb to get <2s|1s>
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00
      ovl=ovl*norm*(A3(ax)*B0(bx)-B3(bx)*A0(ax)+A2(ax)*B1(bx)-B2(bx)*A1(ax))
    endif
   case (4) ! <2s|2s>
      norm=SQRT((ZA*ZB)**5)*(R**5)*0.0625D00
      ovl=norm* (A4(ax)*B0(bx)+B4(bx)*A0(ax)-2.0D00*A2(ax)*B2(bx))*(1d0/3d0)
   case(3)  ! <1s|3s> + <3s|1s>
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*B0(bx)-B4(bx)*A0(ax)+2.d0*(A3(ax)*B1(bx)-B3(bx)*A1(ax)))/sqrt(3.d0)
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*B0(bx)-B4(bx)*A0(ax)+2.d0*(A3(ax)*B1(bx)-B3(bx)*A1(ax)))/sqrt(3.d0)
    endif
   case(6)  ! <2s|3s> + <3s|2s>
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*B0(bx)+A4(ax)*B1(bx)-2d0*(A3(ax)*B2(bx)+A2(ax)*B3(bx))+A1(ax)*B4(bx)+A0(ax)*B5(bx))/3.d0
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*B0(bx)+A4(ax)*B1(bx)-2d0*(A3(ax)*B2(bx)+A2(ax)*B3(bx))+A1(ax)*B4(bx)+A0(ax)*B5(bx))/3.d0
    endif
    case(9) ! <3s|3>
      norm=sqrt((ZA*ZB*R*R)**7)/1440.d0
!      ovl=norm*(A6(ax)*B0(bx)-3.d0*(A4(ax)*B2(bx)-A2(ax)*B4(bx))-A0(ax)*Bint(bx,6))
      ovl=norm*(A6(ax)*B0(bx)-3.d0*(A4(ax)*B2(bx)-A2(ax)*B4(bx))-A0(ax)*B6(bx))
   end select
endif
return
end subroutine ssovl


!****************************************
!* A(x) auxiliary integrals             *
!* Quantenchemie - Ein Lehrgang Vol 5   *
!* p. 570  eq. 11.4.14                  *
!****************************************

real(8) pure function A0(x)
! Hilfsintegral A_0
implicit none
real(8), intent(in) :: x
A0=exp(-x)/x
return
end function

real(8) pure function A1(x)
! Hilfsintegral A_1
implicit none
real(8), intent(in) :: x
A1=((1+x)*exp(-x))/(x**2)
return
end function


real(8) pure function A2(x)
! Hilfsintegral A_2
implicit none
real(8), intent(in) :: x
A2=((2d0+2d0*x+x**2)*exp(-x))/x**3
return
end function


real(8) pure function A3(x)
! Hilfsintegral A_3
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4
x2=x*x
x3=x2*x
x4=x3*x
xx=(6d0+6d0*x+3d0*x2+x3)
A3=(xx*exp(-x))/x4
return
end function


real(8) pure function A4(x)
! Hilfsintegral A_4
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4,x5
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
xx=(24d0+24d0*x+12d0*x2+4d0*x3+x4)
A4=(xx*exp(-x))/x5
return
end function

real(8) pure function A5(x)
! Hilfsintegral A_5
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4,x5,x6
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
xx=(120d0+120d0*x+60d0*x2+20d0*x3+5d0*x4+x5)
A5=(xx*exp(-x))/x6
return
end function

real(8) pure function A6(x)
! Hilfsintegral A_6
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4,x5,x6,x7
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
x7=x6*x
xx=(720d0+720d0*x+360d0*x2+120d0*x3+30d0*x4+6d0*x5+x6)
A6=(xx*exp(-x))/x7
return
end function


!**************************************
!* B(x) auxiliary integrals           *
!* Quantenchemie - Ein Lehrgang Vol 5 *
!* p. 570  eq. 11.4.14b               *
!**************************************


real(8) pure function B0(x)
! Hilfsintegral B_0
implicit none
real(8), intent(in) :: x
B0=(exp(x)-exp(-x))/x
return
end function

real(8) pure function B1(x)
! Hilfsintegral B_1
implicit none
real(8), intent(in) :: x
real(8) x2,x3
x2=x*x
x3=x2*x
B1=((1d0-x)*exp(x)-(1d0+x)*exp(-x))/x2
return
end function

real(8) pure function B2(x)
! Hilfsintegral B_2
implicit none
real(8), intent(in) :: x
real(8) x2,x3
x2=x*x
x3=x2*x
B2=(((2d0-2*x+x2)*exp(x)) - ((2d0+2d0*x+x2)*exp(-x)))/x3
return
end function

real(8) pure function B3(x)
! Hilfsintegral B_3
implicit none
real(8), intent(in) :: x
real(8) xx,yy
real(8) x2,x3,x4
x2=x*x
x3=x2*x
x4=x3*x
xx=(6d0-6d0*x+3d0*x2-x3)*exp(x)/x4
yy=(6d0+6d0*x+3d0*x2+x3)*exp(-x)/x4
B3=xx-yy
return
end function


real(8) pure function B4(x)
! Hilfsintegral B_4
implicit none
real(8), intent(in) :: x
real(8) xx,yy
real(8) x2,x3,x4,x5
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
xx=(24d0-24d0*x+12d0*x2-4d0*x3+x4)*exp(x)/x5
yy=(24d0+24d0*x+12d0*x2+4d0*x3+x4)*exp(-x)/x5
B4=xx-yy
return
end function

real(8) pure function B5(x)
! Hilfsintegral B_5
implicit none
real(8), intent(in) :: x
real(8) xx,yy
real(8) x2,x3,x4,x5,x6
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
xx=(120d0-120*x+60*x2-20*x3+5*x4-x5)*exp(x)/x6
yy=(120d0+120*x+60*x2+20*x3+5*x4+x5)*exp(-x)/x6
B5=xx-yy
return
end function

real(8) function B6(x)
! Hilfsintegral B_6
implicit none
real(8), intent(in) :: x
real(8) x2,x3,x4,x5,x6,x7,yy,xx
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
x7=x6*x
xx=(720d0 - 720d0*x+ 360d0*x2 - 120d0*x3 + 30d0*x4 - 6d0*x5 + x6)*exp(x)/x7
yy=(720d0 + 720d0*x + 360d0*x2 + 120d0*x3 + 30d0*x4 + 6d0*x5 + x6)*exp(-x)/x7
B6=xx-yy
return
end function


real*8 function bint(x,k)
! calculates B_k(x)
! general summation formula
! 'infinite' sum is numerically unstable. 12 terms seem
! accurate enough
implicit none
real(8), intent(in) :: x
real(8) xx,yy
integer, intent(in) :: k
integer i
bint=0

if(abs(x).lt.1e-6) then
do i=0,k
   bint=(1.d0+(-1d0)**i)/(dble(i)+1.d0)
end do
return
endif

do i=0,12
xx=1d0-((-1d0)**(k+i+1))
yy=dble(fact(i))*dble((k+i+1))
bint=bint+xx/yy*(-x)**i
enddo


end function bint


! faculty function
integer(8) function fact(N)
implicit none
integer j,n
fact=1
do j=2,n
  fact=fact*j
enddo
return

end




subroutine gsovl(r,iat,jat,iz,xza,xzb,g)
! GRADIENT
! calculates the s-type overlap integral over 1s,2s,3s slater functions
! ovl = overlap integral
! za  = slater exponent atom A
! zb  = slater exponent atom B
! R   = distance between atom A and B
implicit none
integer ii,shell(72)
logical debug
real(8) ax,bx,R05,za,zb,R
integer na,nb
data shell/                 &
!          h, he
          1,1               &
!         li-ne
          ,2,2,2,2,2,2,2,2, &
!         na-ar
          3,3,3,3,3,3,3,3,  &
! 4s,5s will be treated as 3s
!         k-rn , no f-elements
          54*3/
! ...
real*8 g,Fa,Fb
!--------------------- set exponents ---------------------------------------
real(kind=8) xza,xzb
real(kind=8) xx
integer iat,jat,iz(*)
logical lsame

       za=xza
       zb=xzb
       na=iz(iat)
       nb=iz(jat)
!----------------------------------------------------------------------------

debug=.false.
!debug=.true.

! ii selects kind of ovl by multiplying the shell
! kind    <1s|1s>  <2s|1s>  <2s|2s>  <1s|3s>  <2s|3s>  <3s|3s>
! case:      1        2        4       3        6         9
!
ii=shell(na)*shell(nb)
if(debug) write(*,*) 'gshell', ii
R05=R*0.5
ax=(za+zb)*R05
Fa=(za+zb)
bx=(zb-za)*R05
Fb=(zb-za)
lsame=.false.
!
! same elements
if(za.eq.zb.OR.abs(za-zb).lt.0.1) then
lsame=.true.
! call arguments: gtype(exponents,argumentDeriv.,distance,gradient,(Switch shell),sameElement)
  select case (ii)
   case (1)
     call g1s1s(za,zb,Fa,Fb,R,g,lsame)
   case (2)
    if(shell(na).lt.shell(nb)) then
      call g2s1s(za,zb,Fa,Fb,R,g,.false.,lsame)
     else
      xx=za
      za=zb
      zb=xx
      call g2s1s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case (4)
     call g2s2s(za,zb,Fa,Fb,R,g,lsame)
   case(3)
    if(shell(na).lt.shell(nb)) then
    call g1s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call g1s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case(6)
    if(shell(na).lt.shell(nb)) then
    call g2s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call g2s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case(9)
    call g3s3s(za,zb,Fa,Fb,R,g,lsame)
   end select
else ! different elements
lsame=.false.
   select case (ii)
   case (1)
     call g1s1s(za,zb,Fa,Fb,R,g,lsame)
   return
   case (2)  ! <1s|2s>
    if(shell(na).lt.shell(nb)) then
      call g2s1s(za,zb,Fa,Fb,R,g,.false.,lsame)
     else
      xx=za
      za=zb
      zb=xx
      call g2s1s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case (4) ! <2s|2s>
      call g2s2s(za,zb,Fa,Fb,R,g,lsame)
   case(3)  ! <1s|3s> + <3s|1s>
    if(shell(na).lt.shell(nb)) then
    call g1s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call g1s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case(6)  ! <2s|3s> + <3s|2s>
    if(shell(na).lt.shell(nb)) then
    call g2s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call g2s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
    case(9) ! <3s|3>
    call g3s3s(za,zb,Fa,Fb,R,g,lsame)
   end select
endif

return
end subroutine gsovl


!-------------------------------------------------------------
! Maple was used to find the analy. derivatives of
! the slater integrals (expressions over A,B aux. integrals)
! Optimized fortran code by maple with some human-corrections
!-------------------------------------------------------------
subroutine g1s1s(za,zb,Fa,Fb,R,g,sameElement)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit real(8) (t)
real(8) za,zb,Fa,Fb
real(8) g,R
logical sameElement

if(sameElement) then
  t1 = za ** 2
  t3 = zb ** 2
  t5 = t1 * za * t3 * zb
  t6 = R ** 2
  t7 = t6 ** 2
  t10 = Fa * R
  t14 = exp(-0.5D0 * t10)
  t17 = sqrt(t5 * t7 * t6)
  g = -(1d0/3d0) * t5 * t7 / Fa * (0.2D1 + t10) * t14 / t17
  return
else
  t1 = za ** 2
  t3 = zb ** 2
  t5 = t1 * za * t3 * zb
  t6 = Fb ** 2
  t7 = Fb * R
  t8 = 0.5D0 * t7
  t9 = exp(t8)
  t12 = exp(-t8)
  t15 = t6 * Fa
  t22 = Fa ** 2
  t23 = t22 * t9
  t27 = t22 * t12
  t31 = t6 * Fb
  t32 = R * t31
  t37 = t22 * Fa
  t38 = R * t37
  t43 = R ** 2
  t44 = t43 * t31
  t51 = t43 * t37
  t56 = 0.4D1 * t6 * t9 - 0.4D1 * t6 * t12 + 0.2D1 * t15 * R * t9 -          &
  0.2D1 * t15 * R * t12 - 0.4D1 * t23 + 0.2D1 * t23 * t7 + 0.4D1 * t27       &
  + 0.2D1 * t27 * t7 - 0.2D1 * t32 * t9 - 0.2D1 * t32 * t12 - 0.2D1 * t38 *  &
  t9 + 0.2D1 * t38 * t12 - 0.1D1 * t44 * Fa * t9 - 0.1D1                     &
  * t44 * Fa * t12 + t51 * t9 * Fb + t51 * t12 * Fb
  t61 = exp(-0.5D0 * Fa * R)
  t62 = t43 ** 2
  t65 = sqrt(t5 * t62 * t43)
  g = -0.2D1 * t5 * R * t56 * t61 / t65 / t31 / t37
  return
endif


end subroutine g1s1s


subroutine g2s1s(za,zb,Fa,Fb,R,g,switch,lsame)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical switch
logical lsame
norm=(1d0/24d0)*sqrt(za**3*zb**5*3d0)

if(switch) then
Fb=-Fb
endif

if(lsame) then

      t1 = Fa * R
      t3 = exp(-0.5000000000D0 * t1)
      t6 = Fa ** 2
      t7 = R ** 2
      g = -0.1000000000D-8 * R * t3 * (0.5333333380D10 + 0.2666666670D10 &
      * t1 + 0.1333333333D10 * t6 * t7) / t6
      g=g*norm

else

      t3 = exp(-0.5000000000D0 * Fa * R)
      t4 = Fa ** 2
      t5 = t4 * Fa
      t6 = Fb * R
      t7 = 0.5000000000D0 * t6
      t8 = exp(t7)
      t9 = t5 * t8
      t11 = Fb ** 2
      t12 = t11 * Fa
      t15 = exp(-t7)
      t18 = t4 ** 2
      t19 = R * t18
      t22 = t11 ** 2
      t29 = Fb * t4
      t36 = R ** 2
      t37 = t36 * t18
      t44 = t36 * R
      t48 = -0.12D2 * t9 + 0.4D1 * t12 * t8 - 0.4D1 * t12 * t15 &
      - 0.6D1 * t19 * t8 - 0.6D1 * t22 * t8 * R - 0.6D1 * t22 * t15 &
      * R + 0.4D1 * t29 * t15 - 0.4D1 * t29 * t8 + 0.6D1 * t19 * t15 &
      + 0.2D1 * t37 * t8 * Fb + 0.4D1 * t37 * t15 * Fb + t44 * t18 * t15 * t11
      t49 = t5 * t15
      t51 = t11 * Fb
      t58 = t51 * Fa
      t59 = R * t8
      t76 = t36 * t15
      t79 = t22 * Fa
      t87 = 0.12D2 * t49 - 0.12D2 * t51 * t15 - 0.1D1 * t22 * t4 * t15 * t44 &
      + 0.4D1 * t58 * t59 - 0.8D1 * t58 * R * t15 + 0.4D1 * t9 * t6 + 0.8D1 * &
      t49 * t6 + 0.2D1 * t49 * t11 * t36 + 0.4D1 * t11 * t4 * t59 - 0.2D1 * t51 &
      * t4 * t76 - 0.2D1 * t79 * t36 * t8 - 0.4D1 * t79 * t76 + 0.12D2 * t51 * t8
      g = -0.16D2 * t3 * (t48 + t87) / t36 / t22 / t18
      g=g*norm
endif

return
end subroutine g2s1s

subroutine g2s2s(za,zb,Fa,Fb,R,g,SameElement)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical SameElement

norm=1d0/(16d0*3d0)*SQRT((ZA*ZB)**5)

if(SameElement) then

      t2 = R ** 2
      t5 = Fa ** 2
      t9 = t5 * Fa
      t10 = t2 ** 2
      t16 = exp(-Fa * R / 0.2D1)
      g = (-0.4266666666D2 * R - 0.2133333333D2 * Fa * t2 - 0.2133333333D1 &
          * t5 * t2 * R - 0.1066666666D1 * t9 * t10) * t16 / t9
      g=g*norm


return
else
      t1 = R ** 2
      t3 = 0.3840000000D3 * t1 * Fb
      t4 = t1 * R
      t5 = Fb ** 2
      t7 = 0.6400000000D2 * t4 * t5
      t8 = 0.7680000000D3 * R
      t10 = Fa ** 2
      t11 = t10 ** 2
      t12 = t11 * Fa
      t14 = Fb * R
      t15 = 0.768000000D3 * t14
      t17 = 0.1280000000D3 * t5 * t1
      t21 = 0.256000000D3 * t5 * R
      t22 = t5 * Fb
      t24 = 0.1280000000D3 * t22 * t1
      t26 = t10 * Fa
      t28 = t5 ** 2
      t30 = 0.1280000000D3 * t1 * t28
      t32 = 0.256000000D3 * t22 * R
      t33 = 0.512000000D3 * t5
      t34 = t28 * Fb
      t36 = 0.6400000000D2 * t4 * t34
      t40 = 0.768000000D3 * t28 * R
      t42 = 0.3840000000D3 * t1 * t34
      t45 = 0.1536000000D4 * t28
      t47 = 0.7680000000D3 * t34 * R
      t51 = exp(-0.5D0 * Fa * R)
      t53 = 0.5D0 * t14
      t54 = exp(-t53)
      t68 = exp(t53)
      g = (((t3 + t7 + t8) * t12 + (0.1536000000D4 + t15 + t17) * t11 + &
      (-t21 - t24) * t26 + (t30 - t32 - t33 + t36) * t10 + (t40 + t42) *&
      Fa + t45 + t47) * t51 * t54 + ((t3 - t8 - t7) * t12 + (-0.1536000000D4 &
      + t15 - t17) * t11 + (-t24 + t21) * t26 + (-t30 + t33 - t32 &
      + t36) * t10 + (-t40 + t42) * Fa + t47 - t45) * t51 * t68) / t1 /  &
      t12 / t34


      g=g*norm


return
endif

end subroutine g2s2s






subroutine g1s3s(za,zb,Fa,Fb,R,g,switch,lsame)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical switch
logical lsame

if(switch) Fb=-Fb

norm=SQRT((ZA**3)*(ZB**7)/7.5d00)/(16d0*sqrt(3d0))

if(lsame) then

  t1 = Fa * R
  t3 = exp(-0.5000000000D0 * t1)
  t4 = R ** 2
  g = -0.1600000000D1 * t3 * t4 * R * (0.2D1 + t1) / Fa
  g=g*norm
else

      t3 = exp(-0.5000D0 * Fa * R)
      t4 = Fb ** 2
      t5 = t4 ** 2
      t6 = t5 * Fb
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = Fb * R
      t10 = 0.50D0 * t9
      t11 = exp(t10)
      t15 = exp(-t10)
      t16 = t8 * t15
      t19 = Fa ** 2
      t21 = t8 * R
      t22 = t21 * t15
      t25 = t19 * Fa
      t27 = t8 ** 2
      t31 = t19 ** 2
      t32 = t31 * Fa
      t33 = t8 * t32
      t45 = t4 * Fb
      t48 = t31 * t15
      t55 = t4 * t25
      t56 = t11 * R
      t59 = t15 * R
      t62 = t5 * Fa
      t73 = -0.6D1 * t7 * t8 * t11 - 0.18D2 * t7 * t16 - 0.6D1 * t6 * t19 &
      * t22 - 0.1D1 * t6 * t25 * t27 * t15 + 0.6D1 * t33 * t11 * Fb + 0.18D2 &
      * t33 * t15 * Fb + 0.6D1 * t21 * t32 * t15 * t4 + t27 * t32* t15 * t45 &
      + 0.2D1 * t48 * t45 * t21 + 0.12D2 * t48 * t4 * t8 + 0.12D2 * t55 * t56 &
      + 0.12D2 * t55 * t59 + 0.12D2 * t62 * t56 - 0.36D2 * t62 * t59 - 0.12D2 &
      * t5 * t19 * t16 - 0.2D1 * t5 * t25 * t22
      t74 = t31 * t11
      t79 = t45 * t19
      t92 = R * t32
      t95 = t45 * Fa
      t100 = Fb * t25
      t111 = 0.12D2 * t74 * t9 + 0.36D2 * t48 * t9 + 0.12D2 * t79 * t56 - 0.12D2  &
      * t79 * t59 + 0.48D2 * t5 * t11 - 0.24D2 * t6 * t11 * R - 0.24D2 * t6 * t15 &
      * R - 0.24D2 * t92 * t11 + 0.24D2 * t95 * t11 - 0.24D2 * t95 * t15 + 0.24D2 &
      * t100 * t15 + 0.24D2 * t92 * t15 - 0.24D2 * t100 * t11 - 0.48D2 * t5 * t15 &
      - 0.48D2 * t74 + 0.48D2 * t48
      g = -0.32D2 * t3 * (t73 + t111) / t8 / t6 / t32

      g=g*norm
endif

return
end subroutine g1s3s


subroutine g2s3s(za,zb,Fa,Fb,R,g,switch,lsame)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical switch
logical lsame
norm=sqrt((za**5)*(zb**7)/7.5D00)/96.d0

if(switch) Fb=-Fb

if(lsame) then
      t1 = Fa * R
      t3 = exp(-0.5000000000D0 * t1)
      t6 = Fa ** 2
      t7 = R ** 2
      t14 = t6 ** 2
      t15 = t7 ** 2
      g = -0.2000000000D-8 * R * t3 * (0.1280000000D12 + 0.6400000000D11 &
      * t1 + 0.1280000000D11 * t6 * t7 + 0.1066666670D10 * t6 * Fa * t7 &
      * R + 0.533333333D9 * t14 * t15) / t14
      g=g*norm
else

      t3 = exp(-0.5D0 * Fa * R)
      t4 = Fb ** 2
      t5 = t4 ** 2
      t6 = Fa ** 2
      t7 = t6 * Fa
      t8 = t5 * t7
      t9 = R ** 2
      t11 = 0.50D0 * Fb * R
      t12 = exp(t11)
      t13 = t9 * t12
      t16 = t6 ** 2
      t17 = t16 * Fa
      t18 = exp(-t11)
      t21 = t5 * Fb
      t28 = t9 * t18
      t32 = t9 * R
      t33 = t32 * t18
      t36 = t5 * t4
      t38 = t9 ** 2
      t39 = t38 * t18
      t41 = t21 * Fa
      t42 = R * t12
      t45 = t16 * t6
      t46 = t4 * Fb
      t49 = t46 * t16
      t52 = -0.6D1 * t8 * t13 + 0.120D3 * t17 * t18 + 0.120D3 * t21 * t18 &
       - 0.120D3 * t17 * t12 - 0.120D3 * t21 * t12 - 0.6D1 * t8 * t28 - 0.2D1 &
      * t5 * t16 * t33 + t36 * t7 * t39 - 0.48D2 * t41 * t42 + t45 * t46 * t39 - 0.6D1 * t49 * t13
      t54 = R * t18
      t60 = t46 * t6
      t63 = Fb * t16
      t66 = t5 * Fa
      t69 = t4 * t7
      t72 = t36 * t6
      t75 = t32 * t12
      t78 = Fb * t9
      t84 = Fb * t17
      t87 = -0.24D2 * t46 * t7 * t54 - 0.24D2 * t5 * t6 * t42 + 0.24D2 *&
       t60 * t12 + 0.24D2 * t63 * t18 - 0.24D2 * t66 * t12 - 0.24D2 * t69 &
      * t18 + 0.9D1 * t72 * t33 + 0.3D1 * t72 * t75 + 0.24D2 * t78 * t45 &
      * t12 - 0.6D1 * t49 * t28 + 0.48D2 * t84 * t42
      t102 = t21 * t6
      t105 = t4 * t17
      t113 = t45 * t4
      t118 = 0.72D2 * t84 * t54 + 0.72D2 * t41 * t54 + 0.36D2 * t78 * t45 &
      * t18 + 0.2D1 * t46 * t17 * t33 + 0.24D2 * t4 * t16 * t42 - 0.6D1   &
      * t102 * t13 - 0.6D1 * t105 * t13 + 0.18D2 * t105 * t28 + 0.2D1 &
      * t21 * t7 * t33 - 0.3D1 * t113 * t75 + 0.9D1 * t113 * t33
      t121 = t36 * Fa
      t130 = R * t45
      t145 = 0.18D2 * t102 * t28 + 0.24D2 * t121 * t13 + 0.36D2 * t121 * &
       t28 - 0.24D2 * t60 * t18 - 0.24D2 * t63 * t12 + 0.60D2 * t130 * t18 &
      + 0.60D2 * t36 * t18 * R + 0.24D2 * t69 * t12 + 0.60D2 * t36 * t12 * &
      R - 0.60D2 * t130 * t12 + 0.24D2 * t66 * t18
      g = 0.128D3 * t3 * (t52 + t87 + t118 + t145) / t9 / t36 / t45

  g=g*norm
endif

return
end subroutine g2s3s


subroutine g3s3s(za,zb,Fa,Fb,R,g,SameElement)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical SameElement

norm=sqrt((ZA*ZB)**7)/1440.d0

if(SameElement) then

      t1 = Fa * R
      t3 = exp(-0.5000000000D0 * t1)
      t5 = Fa ** 2
      t6 = t5 ** 2
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = t8 ** 2
      g = -0.2000000000D-8 * t3 * R * (0.457142857D9 * t7 * t9 &
      * R + 0.7680000000D12 * t1 + 0.1536000000D12 * t5 * t8 &
      + 0.1280000000D11 * t5 * Fa * t8 * R + 0.914285715D9 * t6 * t9 + 0.1536000000D13) / t7


      g=g*norm
return
else


      t3 = exp(-0.5000000000D0 * Fa * R)
      t4 = Fa ** 2
      t5 = t4 ** 2
      t6 = t5 * t4
      t7 = Fb * R
      t8 = 0.5000000000D0 * t7
      t9 = exp(-t8)
      t10 = t6 * t9
      t13 = Fb ** 2
      t14 = t13 * Fb
      t15 = t13 ** 2
      t16 = t15 * t14
      t17 = R ** 2
      t18 = t17 * R
      t19 = t16 * t18
      t23 = exp(t8)
      t24 = t6 * t23
      t27 = t5 * Fa
      t28 = t27 * t13
      t29 = R * t23
      t32 = t6 * t13
      t33 = t17 * t9
      t36 = t15 * Fb
      t37 = t4 * t36
      t38 = t9 * R
      t43 = t17 * t23
      t46 = t4 * Fa
      t47 = t5 * t46
      t48 = t47 * t18
      t52 = t47 * t17
      t65 = 0.120D3 * t10 * t7 - 0.12D2 * t19 * t4 * t9 + 0.120D3 &
      * t24 * t7 + 0.24D2 * t28 * t29 + 0.24D2 * t32 * t33 + 0.24D2 * t37 &
      * t38 - 0.24D2 * t28 * t38 - 0.24D2 * t32 * t43 - 0.12D2 * t48 * t13 &
      * t23 + 0.60D2 * t52 * t23 * Fb + 0.12D2 * t48 * t13 * t9 + 0.60D2 &
      * t52 * t9 * Fb - 0.12D2 * t19 * t4 * t23
      t66 = t17 ** 2
      t67 = t16 * t66
      t74 = t27 * t14
      t77 = t6 * t14
      t78 = t18 * t23
      t81 = t46 * t15
      t86 = t27 * t15
      t89 = t5 * t36
      t90 = t18 * t9
      t97 = t46 * t36
      t104 = -0.1D1 * t67 * t46 * t9 - 0.1D1 * t67 * t46 * t23 - 0.12D2 &
      * t74 * t43 + 0.2D1 * t77 * t78 - 0.24D2 * t81 * t29 + 0.24D2 * t81 &
      * t38 + 0.2D1 * t86 * t78 + 0.2D1 * t89 * t90 - 0.2D1 * t86 * t90 &
      + 0.24D2 * t37 * t29 + 0.12D2 * t97 * t33 + 0.2D1 * t89 * t78 - 0.12D2 * t74 * t33
      t108 = t5 * t14
      t111 = t15 * t13
      t112 = t111 * t4
      t117 = t111 * t46
      t122 = t111 * Fa
      t129 = t4 * t15
      t132 = t47 * R
      t139 = 0.2D1 * t77 * t90 - 0.24D2 * t108 * t38 + 0.24D2 * t112 * t43 &
      - 0.24D2 * t112 * t33 + 0.2D1 * t117 * t78 - 0.2D1 * t117 * t90 + 0.120D3 &
      * t122 * t29 - 0.120D3 * t122 * t38 + 0.12D2 * t97 * t43 - 0.48D2 * t129 &
      * t23 + 0.120D3 * t132 * t9 - 0.120D3 * t132 * t23 + 0.240D3 * t111 * t23
      t140 = t47 * t66
      t145 = t16 * R
      t150 = t16 * t17
      t160 = t5 * t13
      t170 = t140 * t14 * t23 + t140 * t14 * t9 - 0.120D3 * t145 * t9 - 0.24D2 &
      * t108 * t29 - 0.60D2 * t150 * Fa * t23 - 0.240D3 * t111 * t9 - 0.240D3 &
      * t24 + 0.240D3 * t10 + 0.48D2 * t129 * t9 - 0.48D2 * t160 * t9 + 0.48D2 &
      * t160 * t23 - 0.120D3 * t145 * t23 - 0.60D2 * t150 * Fa * t9
      g = -0.768D3 * t3 * (t65 + t104 + t139 + t170) / t17 / t47 / t16

      g=g*norm
return
endif

end subroutine g3s3s

subroutine done(aa,io)
! normal program termination
implicit none
integer io
character(*)aa

 write(io,'(a)') ' '
 write(io,'(5x,''normal termination of '',a)') adjustl(trim(aa))
 error stop
end subroutine done

!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
subroutine tmolrd(maxat,xyz,iat,ifrez,nat,infile,echo)
implicit none
integer maxat
character*2 cc,ff
character*80  atmp
character*(*) infile
real(kind=8) xyz(3,maxat),xx(5)
real(kind=8) bohr
integer iat(maxat),nat,i,nn,j,ifrez(maxat),istat,iff
logical da,echo,gaus,tm

bohr=0.52917726d0
i=0
ifrez=0
iff=0

inquire(file=infile,exist=da)
select case (da)
case (.true.)
      if(echo) write(*,'('' reading...'',$)')

open(unit=3,file=infile)
! test for tmol or xyz file

! start anew
do
 read(3,'(a)',end=125) atmp
 if(index(atmp,'$coord').ne.0) then

 ! count number of atoms
 do while (da)
  read(3,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit
   i=i+1
  enddo

 100 continue
 nat=i

 rewind(unit=3)

 ! read TMOL file
 read(3,*) atmp ! $coord
 do j=1,nat
    read(3,'(a)') atmp ! $coord
    backspace(3)
!   if(index(atmp,' F ').ne.0) then
!   read(3,*) xyz(1,j),xyz(2,j),xyz(3,j),cc,ff
!    ifrez(j)=1
!   iff=iff+1
!   else ! default
     read(3,*) xyz(1,j),xyz(2,j),xyz(3,j),cc
!  endif
   call elem(cc,iat(j))
   xyz(1:3,j)=xyz(1:3,j)
  enddo
 if(echo) write(*,*) ' Turbomole file [bohr] :  ', trim(infile)

 close(3)
 exit
 elseif (index(atmp,'G A U S S I A N').ne.0) then
   do
    read(3,*,end=124) atmp ! $coord
     if (index(atmp,'???').ne.0) then
      ! get position Z
      !
     endif
   enddo
124 continue
   ! go to position Z
   ! read in

 exit
 else ! hopefully an xyz file !
 rewind(3)
       read(3,'(a)',end=101) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do
            nat=nat+1
            read(3,'(a)',end=123) atmp
           enddo
          else
            nat=idint(xx(1))
           read(3,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(3,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,iat(i))
            xyz(1:3,i)=xx(1:3)*1d0/bohr
       enddo
      if(echo) write(*,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
101  exit
endif
enddo
125 close(3)

if(maxval(ifrez,nat).eq.1) then
 if(echo) then
  write(*,'(a,x,I4,x,a)') '  found ',iff,' frozen cart. coordinates'
  if(iff.lt.50) then ! dont spam the output to much ...
   write(*,'(a,$)') '  atom nr: '
   do i=1,nat
     if(ifrez(i).eq.1) write(*,'(1x,I2,$)') i
   enddo
   write(*,'(a)') ' '
  endif
 endif
endif

case (.false.)
  write(*,*) ' no input file <',trim(infile) ,'> found !! '
end select
end subroutine

!********************************
!* convert a word to lower case *
!********************************
      subroutine lower_case(word)
      character (len=*) , intent(in out) :: word
      integer :: i,ic,nlen
      nlen = len(word)
      do i=1,nlen
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
      end do
      end subroutine lower_case


!********************************************
!* split a string s into n separate words w *
!********************************************
subroutine charsplit(s,n,w)
implicit none
integer i,n,k
character*80, intent(in) :: s
character*80, intent(out) :: w(n)
character*80   a,aa

aa=adjustl(s)
do i=1,n
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
enddo
return
end subroutine


!***************************************************
!* split a string s iand the return the x'ths word *
!***************************************************
subroutine charXsplit(s,wx,x)
implicit none
integer i,k,x
character*80, intent(in) :: s
character*80, intent(out) ::wx
character*80   w(20)
character*80   a,aa

aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
  if(i.gt.50) error stop '!* string split error: subroutine charXsplit *!'
enddo
wx=w(x)
return
end subroutine



subroutine help(h)
implicit none
logical h

if(.not.h) return

print*,'                                  '
print*,' *    ----   HELP    ----    *    '
print*,'                                   '
print*,' gcp <coordinates> -level <method>   '
print*,'                                   '
print*,' example:                                 '
print*,' -> gcp water.xyz -level dft/631gd     '
print*,' -> gcp POSCAR -pbc -level hf/svp '
print*,'                                  '
print*,'                                  '
print*,' Valid input coordinates are      '
print*,'   - TURBOMOLE files (in bohr)    '
print*,'   - XMOL files (in angstrom)     '
print*,'   - VASP POSCAR files (in angstrom)'
print*,'   - CRYSTAL STRUC.INCOOR files (in bohr)'
print*,'                                  '
print*,' available parametrisation for <method>        '
print*,'      HF/MINIS       DFT/MINIS'
print*,'      HF/MINIX       DFT/MINIX'
print*,'      HF/SV          DFT/SV '
print*,'      HF/def2-SV(P)  DFT/def2-SV(P)'
print*,'      HF/def2-SVP    DFT/def2-SVP'
!print*,'      HF/SVP_OLD     DFT/SVP_old'
print*,'      HF/DZP         DFT/DZP'
print*,'      HF/def-TZVP    DFT/def-TZVP'
print*,'      HF/def2-TZVP   DFT/def2-TZVP'
print*,'      HF/631Gd       DFT/631Gd'
print*,'      HF/def2-TZVP   DFT/def2-TZVP'
print*,'      HF/cc-pVDZ     DFT/cc-pVDZ'
print*,'      HF/aug-cc-pVDZ DFT/aug-cc-pVDZ'
print*,'        -            DFT/SV(P/h,c)'
print*,'        -            DFT/LANL'
print*,'                                  '
print*,' note: <method> does not need to be          '
print*,'       capitelized                           '
print*,'                                  '
print*,' command line options:            '
print*,'    -h                       this print out   '
print*,'    -pbc                     use periodic version   '
print*,'    -level <string>    specify <method>   '
print*,'      if <method> = file or left empty a parameter   '
print*,'      file from ~/.gcppar.$HOSTNAME will be read in    '
print*,'    -l                       same as -level   '
print*,'    -grad                    request gradient '
print*,'    -stress                  request cell gradient '
print*,'    -noprint                 surpress output   '
print*,'    -parfile                 print <gcp.param> '
print*,'    -local                   use  .gcppar in work-dir '
print*,'    -hess                    request hessian  '
print*,'    -test                    stop after parameter setup  '
print*,'    -v                       verbose output: print gradient to stdout  '
print*,'                             instead of gcp_gradient     '
print*,'    -vasp                    define input as VASP POSCAR file '
print*,'    -crystal                 define input as CRYSTAL STRUC.INCOOR file '
print*,'                                  '
print*,'                                  '
print*,' *    ----   ****    ----    *    '

call done('help',6)
end subroutine help

subroutine setzet(eta,etaspec,za,zb)
!*************************
!* set slater exponents  *
!*************************
implicit none
integer i
real(8) za(36),zb(36)
real(8) ZS(36),ZP(36),ZD(36),eta,etaspec
data ZS /1.2000,1.6469,0.6534,1.0365,1.3990,1.7210,2.0348,2.2399,2.5644,2.8812,&
0.8675,1.1935,1.5143,1.7580,1.9860,2.1362,2.3617,2.5796,0.9362,1.2112,&
1.2870,1.3416,1.3570,1.3804,1.4761,1.5465,1.5650,1.5532,1.5781,1.7778,&
2.0675,2.2702,2.4546,2.5680,2.7523,2.9299/
data ZP /0.0000,0.0000,0.5305,0.8994,1.2685,1.6105,1.9398,2.0477,2.4022,2.7421,&
0.6148,0.8809,1.1660,1.4337,1.6755,1.7721,2.0176,2.2501,0.6914,0.9329,&
0.9828,1.0104,0.9947,0.9784,1.0641,1.1114,1.1001,1.0594,1.0527,1.2448,&
1.5073,1.7680,1.9819,2.0548,2.2652,2.4617/
data ZD /0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,&
0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,&
2.4341,2.6439,2.7809,2.9775,3.2208,3.4537,3.6023,3.7017,3.8962,2.0477,&
2.4022,2.7421,0.6148,0.8809,1.1660,1.4337/

  do i=1,36
select case (i)
  case(:2)
    za(i)=ZS(i)
  case(3:20,31:)
    za(i)=( ZS(i)+ZP(i) )/2d0
  case(21:30)
    za(i)=( ZS(i)+ZP(i)+ZD(i) )/3d0
end select
  enddo

  za(11:36)=za(11:36)*etaspec !SG r2scan-3c/def2-mtzvpp change

  za=za*eta
  zb=za

 return
end subroutine setzet



!C     *****************************************************************

      SUBROUTINE READL(A1,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) A1
      DIMENSION X(*)
      I=0
      IS=1
  10  I=I+1
      X(I)=READAA(A1,IS,IB,IE)
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END


      FUNCTION READAA(A,ISTART,IEND,IEND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 READAA
      CHARACTER*(*) A
      NINE=ICHAR('9')
      IZERO=ICHAR('0')
      MINUS=ICHAR('-')
      IDOT=ICHAR('.')
      ND=ICHAR('D')
      NE=ICHAR('E')
      IBL=ICHAR(' ')
      IEND=0
      IEND2=0
      IDIG=0
      C1=0
      C2=0
      ONE=1.D0
      X = 1.D0
      NL=LEN(A)
      DO 10 J=ISTART,NL-1
         N=ICHAR(A(J:J))
         M=ICHAR(A(J+1:J+1))
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20

   10 CONTINUE
      READAA=0.D0
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I:I))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C1=C1*10+N-IZERO
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
            ONE=-1.D0
         ELSEIF(N.EQ.IDOT) THEN
            GOTO 40
         ELSE
            GOTO 60
         ENDIF
   30 CONTINUE
   40 CONTINUE
      IDIG=0
      DO 50 II=I+1,NL
         N=ICHAR(A(II:II))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C2=C2*10+N-IZERO
            X = X /10
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
            X=-X
         ELSE
            GOTO 60
         ENDIF
   50 CONTINUE
!C
!C PUT THE PIECES TOGETHER
!C
   60 CONTINUE
      READAA= ONE * ( C1 + C2 * X)
      DO J=IEND,NL
         N=ICHAR(A(J:J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
      ENDDO
      IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN

   57 C1=0.0D0
      ONE=1.0D0
      DO I=J+1,NL
         N=ICHAR(A(I:I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
         IF(N.EQ.MINUS)ONE=-1.0D0
      ENDDO
   61 CONTINUE
   70 READAA=READAA*10**(ONE*C1)
      RETURN
      END

!C     *****************************************************************


 SUBROUTINE ELEM(KEY1, NAT)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 CHARACTER*(*) KEY1
 CHARACTER*2 ELEMNT(94),E

 DATA ELEMNT/'h ','he',                                      &
 'li','be','b ','c ','n ','o ','f ','ne',                    &
 'na','mg','al','si','p ','s ','cl','ar',                    &
 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',     &
 'zn','ga','ge','as','se','br','kr',                         &
 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',     &
 'cd','in','sn','sb','te','i ','xe',                         &
 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
 'au','hg','tl','pb','bi','po','at','rn',                    &
 'fr','ra','ac','th','pa','u ','np','pu'/

 nat=0
 e='  '
 k=1
 DO J=1,len(key1)
    if (k.gt.2)exit
    N=ICHAR(key1(J:J))
    if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
       e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
       k=k+1
    endif
    if(n.ge.ichar('a') .and. n.le.ichar('z') )then
       e(k:k)=key1(j:j)
       k=k+1
    endif
 enddo

 DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

FUNCTION ESYM(I)
CHARACTER*2 ESYM
CHARACTER*2 ELEMNT(94)
DATA ELEMNT/'h ','he',                                           &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu'/
  ESYM=ELEMNT(I)
  RETURN
  END




! printout for gradient
subroutine wregrad(maxat,nat,xyz,iat,ifrez,edisp,gin,echo)
implicit none
integer nat,iat(nat),maxat
integer ifrez(nat)
real*8 edisp
real*8,  intent(in) :: gin(3,nat)
real*8  g(3,nat)
real*8 xyz(3,nat)

integer i,j,nn,nl
character*128 a,a1
real*8 xx(10),gsum,x,y,z
real*8 gr(3,nat)
logical ex,echo

g=gin

open(unit=142,file='gcp_gradient')
! write(*,*) 'Ggcp'
do i=1,nat
   write(142,'(F22.19,x,F22.19,x,F22.19)')g(1:3,i)
enddo
close(142)
return
end

!*******************************************************
!* For Turbomole:                                      *
!* add gcp energy in file <energy>                     *
!* and g to file <gradient>                            *
!*******************************************************
subroutine wregrad_tm(maxat,nat,xyz,iat,ifrez,edisp,gin,echo)
implicit none
integer nat,iat(nat),maxat
integer ifrez(nat)
real*8 edisp
real*8,  intent(in) :: gin(3,nat)
real*8  g(3,nat)
real*8 xyz(3,nat)
integer i,j,nn,nl
character*128 a,a1
real*8 xx(10),gsum,x,y,z
real*8 gr(3,nat)
logical ex,echo

g=gin

inquire(file='gradient',exist=ex)
if(.not.ex) then
 if(echo) then
    write(*,*) 'no gradient file found to add G(gcp)!'
    write(*,*) 'hence written to file gcp_gradient'
 endif
open(unit=142,file='gcp_gradient')
do i=1,nat
   write(142,'(3D22.14)')g(1:3,i)
enddo
close(142)
    return
endif
! write file gradient
      j=0
      open(unit=42,file='gradient')
201   read(42,'(a)',end=301)a1
      j=j+1
      if(index(a1,'cycle').ne.0)nl=j
      goto 201
301   continue

      if(nl.lt.2)then
         write(*,*) 'illegal gradient file to add Gdisp!'
         return
      endif

      rewind 42
      do i=1,nl
      read(42,'(a)')a1
      enddo
      do i=1,nat
         read(42,*)x,y,z
         xyz(1,i)=x
         xyz(2,i)=y
         xyz(3,i)=z
      enddo
      gsum=0
      do i=1,nat
         read(42,*)gr(1,i),gr(2,i),gr(3,i)
         g(1:3,i)=g(1:3,i)+gr(1:3,i)
      enddo
      gsum=sqrt(sum(g(1:3,1:nat)**2))

      rewind 42
      open(unit=43,file='gradient.tmp')
      j=0
401   read(42,'(a)',end=501)a1
      j=j+1
      if(j.lt.nl)then
         write(43,'(a)')trim(a1)
      else
         call readl(a1,xx,nn)
         j=idint(xx(1))
         write(43,'(''  cycle = '',i6,4x,''SCF energy ='',F18.11,3x,  &
                   &''|dE/dxyz| ='',F10.6)')j,xx(2)+edisp,gsum
         do i=1,nat
          if(ifrez(i).eq.1) then
            write(43,'(3(F20.14,2x),4x,2a2)')xyz(1,i),xyz(2,i),xyz(3,i), &
                                           esym(iat(i)),' f'
          else
            write(43,'(3(F20.14,2x),4x,a2)')xyz(1,i),xyz(2,i),xyz(3,i), &
                                           esym(iat(i))
          endif
         enddo
         do i=1,nat
            write(43,'(3D22.13)')g(1,i),g(2,i),g(3,i)
         enddo
         a='$end'
         write(43,'(a)')trim(a)
         goto 501
      endif
      goto 401
501   continue
      close(42)
      close(43)

      call system('mv gradient.tmp gradient')

! write file energy
      j=1
      open(unit=42,file='energy')
      open(unit=43,file='energy.tmp')
 50   read(42,'(a)',end=100)a
      call readl(a,xx,nn)
      if(nn.gt.3)nl=j
      j=j+1
      goto 50
100   continue

      rewind 42
      j=0
  20  read(42,'(a)',end=200)a
      j=j+1
      if(j.lt.nl)then
         write(43,'(a)')trim(a)
         call readl(a,xx,nn)
      else
         call readl(a,xx,nn)
         xx(2)=xx(2)+edisp
         write(43,'(i6,4F20.11)')idint(xx(1)),xx(2:nn)
         a='$end'
         write(43,'(a)')trim(a)
         goto 200
      endif
      goto 20
 200  continue
close(42)
close(43)

call system('mv energy.tmp energy')

end


! subroutine hessian(n,max_elem,emiss,xyz,iz,p,gcp,xva,xvb,DoHess,damp,base,dmp_scal,dmp_exp)
subroutine hessian(n,max_elem,emiss,xyz,iz,p,xva,xvb,pbc,lat,damp,base,dmp_scal,dmp_exp)
implicit none
integer iz(*),n,max_elem
integer ii,jj,hi,hj,i,j
real*8  xyz(3,n),xva(*),xvb(*)
real*8  p(4),gcp,step
real*8  emiss(max_elem),dmp_scal,dmp_exp
real*8  t1,t0,xx,thrE,lat(3,3)
real*8, allocatable :: h(:,:), hess(:,:), gl(:,:), gr(:,:)
logical da,damp,base,pbc
da=.false.

thrE=epsilon(1d0)
step=0.005d0

write(*,'(a)') 'Doing Hessian numerically ..'
  call cpu_time(t0)
  allocate(hess(3*n,3*n),gl(3,n),gr(3,n))
  hess=0.0d0
  do i=1,n
   do j=1,3
    hi=(i-1)*3+j
    xyz(j,i)=xyz(j,i)+step
    call gcp_egrad(n,max_elem,emiss, &
            xyz,iz,p,gcp,gr,.true.,.false.,xva,xvb,pbc,lat,damp,base,dmp_scal,dmp_exp)
    xyz(j,i)=xyz(j,i)-2d0*step
    call gcp_egrad(n,max_elem,emiss, &
            xyz,iz,p,gcp,gl,.true.,.false.,xva,xvb,pbc,lat,damp,base,dmp_scal,dmp_exp)
    xyz(j,i)=xyz(j,i)+step

  do ii=1,n
   do jj=1,3
    hj=(ii-1)*3+jj
    xx=(gr(jj,ii)-gl(jj,ii))/(2d0*step)
    if(abs(xx).gt.thrE)  hess(hi,hj)=xx
   enddo ! jj-loop
  enddo  ! ii-loop

 enddo ! j-loop
 enddo  ! i-loop


!  deallocate(gl,gr)
  allocate(h(3*n,3*n))
  h=0d0

! symmetrize
  write(*,'(a)')  'Symmetrizing Hessian ..'
  do i=1,n*3
     do j=1,n*3
        h(j,i)=(hess(i,j)+hess(j,i))*0.5d0
        h(i,j)=h(j,i)
     enddo
  enddo

  write(*,'(a)')  'gCP Hessian correction written into .. ''gcp_hessian'' '
call wrhess(3*n,h,'gcp_hessian')

deallocate(h,hess)

  call cpu_time(t1)
  call prtim(6,t1-t0,'t','hessian')

end subroutine hessian



subroutine rdhess(nat3,h,fname)
implicit none
integer nat3
real*8 h(nat3,nat3)
character*(*) fname
integer iunit,i,j,mincol,maxcol
character*5 adum
character*80 a80

write(*,*)
write(*,*) 'Reading <',trim(fname),'>'
iunit=11
open(unit=iunit,file=fname)

do
read(iunit,'(a)')a80
  if(index(a80,'$hessian').ne.0)then
   do  i=1,nat3
      maxcol = 0
 200   mincol = maxcol + 1
       maxcol = min(maxcol+5,nat3)
       read(iunit,'(a5   ,5f15.10)')adum ,(h(i,j),j=mincol,maxcol)
       if (maxcol.lt.nat3) goto 200
   enddo
 close(iunit,status='keep')
 exit
 endif
enddo

  300 return
end

subroutine wrhess(nat3,h,fname)
implicit none
integer nat3
real*8 h(nat3,nat3)
character*(*) fname
integer iunit,i,j,mincol,maxcol
character*5 adum
character*80 a80

adum='   '
iunit=11
open(unit=iunit,file=fname)
a80='$hessian'
write(iunit,'(a)')a80
do i=1,nat3
  maxcol = 0
  do while (maxcol.lt.nat3)
    mincol = maxcol + 1
    maxcol = min(maxcol+5,nat3)
    write(iunit,'(a5,5f15.10)')adum,(h(i,j),j=mincol,maxcol)
  enddo !if (maxcol.lt.nat3) goto 200
enddo
write(iunit,'(''$end'')')
close(iunit)
end

character(80) function getlevel(method)
implicit none
character(*) method
character(80) bas
integer i

!do while (index(method,'/').ne.0)
  i=scan(method,'/')
  bas=trim(method((i+1):))
!enddo
getlevel=' '
select case(bas)
case('minis')
 getlevel='basis: MINIS'
case('minix')
 getlevel='basis: MINIX'
case('sv')
 getlevel='basis: SV (Ahlrichs)'
case('sv(p)', 'def2-sv(p)')
 getlevel='basis: def2-SV(P)'
case('svp')
 getlevel='basis: def2-SVP'
case('svx')
 getlevel='basis: def2-SV(P/h,c)'
case('svp_old')
 getlevel='basis: def2-SVP (old parameters)'
case('631gd')
 getlevel='basis: 6-31G(d) (Turbomole)'
case('lanl')
 getlevel='basis: 6-31G(d)+LANL2DZ(Sc-Zn) (Turbomole)'
case('tz')
 getlevel='basis: def2-TZVP'
case('def1tzvp')
 getlevel='basis: def1-TZVP'
case('pbeh3c')
 getlevel='basis: def2-mSVP'
end select

end function


subroutine basegrad(n,max_elem,iz,coord,lat,pbc,rscal,qscal,e,g,echo)
implicit none
integer n,max_elem,maxat,iz(*)
real*8 coord(3,*),e,g(3,*),lat(3,3),xyzjat(3)
real*8 rscal,qscal
real*8 fi,fj,ff,rf,r,expt
!c cut-off radii for all element pairs
real*8 r0ab(max_elem,max_elem),autoang
parameter (autoang =0.5291772083d0)
logical echo,pbc
real*8 r0,thrR,vec(3)
integer i,j,tau_max(3),a,b,c

!threshold
thrR=30            ! 30 bohr
g(1:3,1:n)=0
e=0

call setr0ab(max_elem,autoang,r0ab)

!c paramter of the method are rscal and qscal
!c     rscal=0.7
!c     qscal=0.03

if(pbc) then
!Determine supercell
  call criteria(thrR, lat, tau_max)
  ! Loop over all i atoms
  do i=1,n
  ! the BSSE due to atom jat, Loop over all j atoms in supercell
  do j=1,i
  if(iz(i).lt.1.or.iz(i).gt.18) cycle
  if(iz(j).lt.1.or.iz(j).gt.18) cycle
  fi=float(iz(i))
  fj=float(iz(j))
  if (i.eq.j) then
     ff=-0.5d0*(fi*fj)**1.5d0
  else
     ff=-(fi*fj)**1.5d0
  end if
  r0=rscal*r0ab(iz(i),iz(j))**0.75d0
  do a=-tau_max(1),tau_max(1)
  do b=-tau_max(2),tau_max(2)
  do c=-tau_max(3),tau_max(3)
  !remove selfinteraction
  if((i.eq.j).and.(abs(a)+abs(b)+abs(c).eq.0)) cycle
     xyzjat(1)=coord(1,j)+a*lat(1,1)+b*lat(2,1)+c*lat(3,1)
     xyzjat(2)=coord(2,j)+a*lat(1,2)+b*lat(2,2)+c*lat(3,2)
     xyzjat(3)=coord(3,j)+a*lat(1,3)+b*lat(2,3)+c*lat(3,3)
     vec(1)=(coord(1,i)-xyzjat(1))
     vec(2)=(coord(2,i)-xyzjat(2))
     vec(3)=(coord(3,i)-xyzjat(3))
     r=sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
     if(r.gt.thrR) cycle
     expt=exp(-r0*r)
     e=e+ff*expt
     rf=qscal/r
     g(1,i)=g(1,i)-ff*r0*vec(1)*expt*rf
     g(1,j)=g(1,j)+ff*r0*vec(1)*expt*rf
     g(2,i)=g(2,i)-ff*r0*vec(2)*expt*rf
     g(2,j)=g(2,j)+ff*r0*vec(2)*expt*rf
     g(3,i)=g(3,i)-ff*r0*vec(3)*expt*rf
     g(3,j)=g(3,j)+ff*r0*vec(3)*expt*rf
  enddo
  enddo
  enddo
  enddo
enddo
e=e*qscal
else
do i=1,n-1
 do j=i+1,n
  if(iz(i).lt.1.or.iz(i).gt.18) cycle
  if(iz(j).lt.1.or.iz(j).gt.18) cycle
  vec(1)=coord(1,i)-coord(1,j)
  vec(2)=coord(2,i)-coord(2,j)
  vec(3)=coord(3,i)-coord(3,j)
  r=sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
  if(r.gt.thrR) cycle
  r0=rscal*r0ab(iz(i),iz(j))**0.75d0
  fi=float(iz(i))
  fj=float(iz(j))
  ff=-(fi*fj)**1.5d0
  expt=exp(-r0*r)
  e=e+ff*expt
  rf=qscal/r
  g(1,i)=g(1,i)-ff*r0*vec(1)*expt*rf
  g(1,j)=g(1,j)+ff*r0*vec(1)*expt*rf
  g(2,i)=g(2,i)-ff*r0*vec(2)*expt*rf
  g(2,j)=g(2,j)+ff*r0*vec(2)*expt*rf
  g(3,i)=g(3,i)-ff*r0*vec(3)*expt*rf
  g(3,j)=g(3,j)+ff*r0*vec(3)*expt*rf
  enddo
enddo
e=e*qscal
endif

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    rdatomnumbervasp
!c    reads the number of atoms from vasp file POSCAR
subroutine rdatomnumbervasp(infile,n,echo)
use gcp_strings
implicit none
!c input
character*20 infile    !inputfile
logical echo
!c output
integer n              !atomnumber
!c intern
integer nargs,i,k,ios
logical da,old
character*80 line, args(90)
old=.false.
n=0
inquire(file=infile,exist=da)
if(da.neqv..true.) then
   error stop 'no input file'
end if
open(unit=3,file=infile)
read(3,'(a)') line !comment or atom species
read(3,'(a)') line !scale factor
read(3,'(a)') line !cell matrix
read(3,'(a)') line
read(3,'(a)') line
read(3,'(a)') line !atom species or #atoms per specie
call parse(line,' ',args,nargs)
do i=1,nargs
   call value(args(i),k,ios)
   if(ios.eq.0) then
      old=.true.
   endif
   n=n+k
end do
if(.not.old) then
   read(3,'(a)') line !#atoms per species
   n=0
   call parse(line,' ',args,nargs)
   do i=1,nargs
      call value(args(i),k,ios)
      n=n+k
   end do
endif
close(3)
end subroutine rdatomnumbervasp

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    rdcoordvasp
!c    read geometry from vasp file POSCAR
!c    reads fixed coordinates (selective dynamics)
subroutine rdcoordvasp(xyz,lat,iat,nat,infile,echo)
use gcp_strings
implicit none
!input
character*(*) infile    !input file
integer nat             !atombumber
logical echo
!output
real*8 xyz(3,nat)       !cartesian coordinates in au
integer sdxyz(3,nat)    !atom gradient will be multiplied with sdxyz
                        !sdxyz=1: move atom  sdxyz=0: fix atom
real*8 lat(3,3)         !cellmatrix in au
integer iat(nat)        !element ordinal number
!intern
character*2 cc
character*80  atmp
character(100) :: comment
character*80 line,args(90),args2(90), speciex(94)
real(kind=8) xx(5), scalar, abc(3,nat)
integer i,nn,j, ios, nargs, nargs2, k, natspeciex(94), nspecies,l
real(kind=8) bohr
logical da,sd,old
bohr=0.52917726d0
sd=.false.
old=.false.
sdxyz=1
i=0
inquire(file=infile,exist=da)
if(da.eqv..true.) then
   if(echo.eqv..true.) then
      write(*,'('' Vasp structure file read '')')
   end if
end if

!c test for input version
open(unit=3,file=infile)
read(3,'(a)') line   !1
read(3,'(a)') line   !2
read(3,'(a)') line   !3
read(3,'(a)') line   !4
read(3,'(a)') line   !5
read(3,'(a)') line   !6
call parse(line,' ',args,nargs)
call value(args(1),k,ios)
if(ios.eq.0) then
   old=.true.
endif
close(3)

open(unit=3,file=infile)
read(3,'(a)') comment   !comment line or atom specie
read(3,'(a)') line      !scale factor
call parse(line,' ',args,nargs)
call value(args(1),scalar,ios)
!* Read cellmatrix in A, convert to a.u.
do k=1,3
   read(3,'(a)') line
   call parse(line,' ',args,nargs)
   do i=1,3
      call value(args(i),lat(k,i),ios)
      lat(k,i)=scalar*lat(k,i)/bohr
   end do
end do
!* Read species and numbers
!c (nat redundant, already as input)
nat=0
call readline(3,line,ios)
call parse(line, ' ',args, nargs)
if(.not.old) then
   call readline(3,line,ios)
   call parse(line, ' ',args2, nargs2)
else
   args2=args
   nargs2=nargs
   call parse(comment, ' ',args, nargs)
endif
if(nargs.eq.nargs2) then
   nspecies=nargs
   do i=1,nspecies
      speciex(i)=args(i)
      call value(args2(i), natspeciex(i),ios)
      nat=nat+natspeciex(i)
   end do
end if
!* Read Coordinate system
call readline(3,line,ios)
call parse(line,' ', args,nargs)
comment=lowercase(args(1))
if(comment(1:1).eq.'s') then
   sd=.true.
   call readline(3,line,ios)
   call parse(line,' ', args,nargs)
   comment=lowercase(args(1))
end if
if(comment(1:1).eq.'d') then
!   write(*,*) 'Direct coordinates found, convert to cartesian'
!* Read coordinates
   i=1
   do k=1,nspecies
      do l=1,natspeciex(k)
         call readline(3,line,ios)
         call parse(line,' ', args,nargs)
         do j=1,3
            call value(args(j),abc(j,i),ios)
         end do
         if(sd) then
            do j=1,3
               args(j+3)=lowercase(args(j+3))
               if(args(j+3).eq.'f')then
                  sdxyz(j,i)=0
               end if
            end do
         end if
         do j=1,3
            xyz(j,i)=lat(1,j)*abc(j,i)+lat(2,j)*abc(2,i)+lat(3,j)*abc(3,i)
         end do
         call elem(speciex(k),iat(i))
         xyz(1:3,i)=xyz(1:3,i)
         i=i+1
      end do
   end do
elseif(comment(1:1).eq.'c') then
   i=1
   do k=1,nspecies
      do l=1,natspeciex(k)
         call readline(3,line,ios)
         call parse(line,' ', args,nargs)
         do j=1,3
            call value(args(j),xyz(j,i),ios)
            xyz(j,i)=scalar*xyz(j,i)/bohr
         end do
         !test for fixed atoms
         if(sd) then
            do j=1,3
               args(j+3)=lowercase(args(j+3))
               if(args(j+3).eq.'f')then
                  sdxyz(j,i)=0
               end if
            end do
         end if
         call elem(speciex(k),iat(i))
         xyz(1:3,i)=xyz(1:3,i)
         i=i+1
      end do
   end do
end if
close(3)
end subroutine rdcoordvasp


SUBROUTINE criteria(r_cutoff,lat_trans,tau_max)

REAL*8 :: r_cutoff
REAL*8 :: lat(3,3)
REAL*8 :: lat_trans(3,3)
integer :: tau_max(3)
REAL*8 :: norm1(3),norm2(3),norm3(3)
REAL*8 :: cos10,cos21,cos32
lat = transpose(lat_trans)

!c find normal to the plane...
call crossprodukt(lat(:,2),lat(:,3),norm1)
call crossprodukt(lat(:,3),lat(:,1),norm2)
call crossprodukt(lat(:,1),lat(:,2),norm3)
!c ...normalize it...
norm1=norm1/vectornorm(norm1)
norm2=norm2/vectornorm(norm2)
norm3=norm3/vectornorm(norm3)
!c cos angles between normals and lattice vectors
cos10=SUM(norm1*lat(:,1))
cos21=SUM(norm2*lat(:,2))
cos32=SUM(norm3*lat(:,3))
tau_max(1)= ceiling(abs(r_cutoff/cos10))
tau_max(2)= ceiling(abs(r_cutoff/cos21))
tau_max(3)= ceiling(abs(r_cutoff/cos32))
!write(*,'(3f8.4)')criteria(1),criteria(2),criteria(3)
END SUBROUTINE CRITERIA


FUNCTION vectornorm(VECT)
REAL*8 :: VECT(3)
REAL*8 :: SVECT(3)
REAL*8 :: vectornorm
SVECT=VECT*VECT
vectornorm=SUM(SVECT)
vectornorm=vectornorm**(0.5)
END FUNCTION vectornorm

SUBROUTINE crossprodukt(A,B,C)
IMPLICIT NONE
REAL*8 :: A(3),B(3)
REAL*8 :: X,Y,Z
REAL*8 :: C(3)
X=A(2)*B(3)-B(2)*A(3)
Y=A(3)*B(1)-B(3)*A(1)
Z=A(1)*B(2)-B(1)*A(2)
C=(/X,Y,Z/)
END SUBROUTINE crossprodukt

subroutine abc2xyz(n,xyz,abc,lat)
implicit none
integer :: n,i,j
real*8 :: xyz(3,n), abc(3,n), lat(3,3),latt(3,3)
xyz=0
!transpose matrix to get
! a1 b1 c1
! a2 b2 c2
! a3 b3 c3
do i=1,3
   do j=1,3
      latt(i,j)=lat(j,i)
   end do
end do
do i=1,n
   do j=1,3
      xyz(j,i)=latt(j,1)*abc(1,i)+latt(j,2)*abc(2,i)+latt(j,3)*abc(3,i)
   end do
end do
end subroutine abc2xyz

subroutine xyz2abc(n,xyz,abc,lat)
implicit none
integer :: n,i,j
real*8 :: xyz(3,n), abc(3,n), lat(3,3), ilat(3,3), detlat, cut,latt(3,3)
ilat=0
abc=0
!transpose matrix to get
! a1 b1 c1
! a2 b2 c2
! a3 b3 c3
do i=1,3
   do j=1,3
      latt(i,j)=lat(j,i)
   end do
end do


! Calculate inverse cell metrix
detlat=0
detlat =detlat+latt(1,1)*latt(2,2)*latt(3,3)+latt(1,2)*latt(2,3)*latt(3,1)+latt(1,3)*latt(2,1)*latt(3,2)
detlat=detlat-latt(3,1)*latt(2,2)*latt(1,3)-latt(3,2)*latt(2,3)*latt(1,1)-latt(3,3)*latt(2,1)*latt(1,1)

if(detlat.eq.0) then
   write(*,*) ' Error: singular cell matrix'
   error stop 'Singular cellmatrix'
end if
ilat(1,1)=latt(2,2)*latt(3,3)-latt(3,2)*latt(2,3)
ilat(1,2)=latt(1,3)*latt(3,2)-latt(3,3)*latt(1,2)
ilat(1,3)=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)
ilat(2,1)=latt(2,3)*latt(3,1)-latt(3,3)*latt(2,1)
ilat(2,2)=latt(1,1)*latt(3,3)-latt(3,1)*latt(1,3)
ilat(2,3)=latt(1,3)*latt(2,1)-latt(2,3)*latt(1,1)
ilat(3,1)=latt(2,1)*latt(3,2)-latt(3,1)*latt(2,2)
ilat(3,2)=latt(1,2)*latt(3,1)-latt(3,2)*latt(1,1)
ilat(3,3)=latt(1,1)*latt(2,2)-latt(2,1)*latt(1,2)
do i=1,3
   do j=1,3
      ilat(i,j)=ilat(i,j)/detlat
   end do
end do

do i=1,n
   do j=1,3
      abc(j,i)=ilat(j,1)*xyz(1,i)+ilat(j,2)*xyz(2,i)+ilat(j,3)*xyz(3,i)
      !cut=int(abc(j,i))
      !abc(j,i)=abc(j,i)-cut
      !if(abc(j,i).lt.0) then
      !   abc(j,i)=abc(j,i)+1
      !end if
   end do
end do

end subroutine xyz2abc

subroutine set_criteria_gcp(r_cut, lat, tau_max)
implicit none
real*8 :: r_cut, lat(3,3), latlen(3)
integer :: i, tau_max(3)
do i=1,3
   latlen(i) = lat(i,1)+lat(i,2)+lat(i,3)
   tau_max(i) = int(ceiling(r_cut/latlen(i)))
end do
end subroutine set_criteria_gcp


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C set cut-off radii
!C in parts due to INTEL compiler bug
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine setr0ab(max_elem,autoang,r)
implicit none
integer max_elem,i,j,k
real*8 r(max_elem,max_elem),autoang
real*8 r0ab(4465)

      r0ab(   1:  70)=(/ &
        2.1823,  1.8547,  1.7347,  2.9086,  2.5732,  3.4956,  2.3550 &
     ,  2.5095,  2.9802,  3.0982,  2.5141,  2.3917,  2.9977,  2.9484 &
     ,  3.2160,  2.4492,  2.2527,  3.1933,  3.0214,  2.9531,  2.9103 &
     ,  2.3667,  2.1328,  2.8784,  2.7660,  2.7776,  2.7063,  2.6225 &
     ,  2.1768,  2.0625,  2.6395,  2.6648,  2.6482,  2.5697,  2.4846 &
     ,  2.4817,  2.0646,  1.9891,  2.5086,  2.6908,  2.6233,  2.4770 &
     ,  2.3885,  2.3511,  2.2996,  1.9892,  1.9251,  2.4190,  2.5473 &
     ,  2.4994,  2.4091,  2.3176,  2.2571,  2.1946,  2.1374,  2.9898 &
     ,  2.6397,  3.6031,  3.1219,  3.7620,  3.2485,  2.9357,  2.7093 &
     ,  2.5781,  2.4839,  3.7082,  2.5129,  2.7321,  3.1052,  3.2962 &
     /)
      r0ab(  71: 140)=(/ &
        3.1331,  3.2000,  2.9586,  3.0822,  2.8582,  2.7120,  3.2570 &
     ,  3.4839,  2.8766,  2.7427,  3.2776,  3.2363,  3.5929,  3.2826 &
     ,  3.0911,  2.9369,  2.9030,  2.7789,  3.3921,  3.3970,  4.0106 &
     ,  2.8884,  2.6605,  3.7513,  3.1613,  3.3605,  3.3325,  3.0991 &
     ,  2.9297,  2.8674,  2.7571,  3.8129,  3.3266,  3.7105,  3.7917 &
     ,  2.8304,  2.5538,  3.3932,  3.1193,  3.1866,  3.1245,  3.0465 &
     ,  2.8727,  2.7664,  2.6926,  3.4608,  3.2984,  3.5142,  3.5418 &
     ,  3.5017,  2.6190,  2.4797,  3.1331,  3.0540,  3.0651,  2.9879 &
     ,  2.9054,  2.8805,  2.7330,  2.6331,  3.2096,  3.5668,  3.3684 &
     ,  3.3686,  3.3180,  3.3107,  2.4757,  2.4019,  2.9789,  3.1468 &
     /)
     r0ab( 141: 210)=(/ &
        2.9768,  2.8848,  2.7952,  2.7457,  2.6881,  2.5728,  3.0574 &
     ,  3.3264,  3.3562,  3.2529,  3.1916,  3.1523,  3.1046,  2.3725 &
     ,  2.3289,  2.8760,  2.9804,  2.9093,  2.8040,  2.7071,  2.6386 &
     ,  2.5720,  2.5139,  2.9517,  3.1606,  3.2085,  3.1692,  3.0982 &
     ,  3.0352,  2.9730,  2.9148,  3.2147,  2.8315,  3.8724,  3.4621 &
     ,  3.8823,  3.3760,  3.0746,  2.8817,  2.7552,  2.6605,  3.9740 &
     ,  3.6192,  3.6569,  3.9586,  3.6188,  3.3917,  3.2479,  3.1434 &
     ,  4.2411,  2.7597,  3.0588,  3.3474,  3.6214,  3.4353,  3.4729 &
     ,  3.2487,  3.3200,  3.0914,  2.9403,  3.4972,  3.7993,  3.6773 &
     ,  3.8678,  3.5808,  3.8243,  3.5826,  3.4156,  3.8765,  4.1035 &
     /)
     r0ab( 211: 280)=(/ &
        2.7361,  2.9765,  3.2475,  3.5004,  3.4185,  3.4378,  3.2084 &
     ,  3.2787,  3.0604,  2.9187,  3.4037,  3.6759,  3.6586,  3.8327 &
     ,  3.5372,  3.7665,  3.5310,  3.3700,  3.7788,  3.9804,  3.8903 &
     ,  2.6832,  2.9060,  3.2613,  3.4359,  3.3538,  3.3860,  3.1550 &
     ,  3.2300,  3.0133,  2.8736,  3.4024,  3.6142,  3.5979,  3.5295 &
     ,  3.4834,  3.7140,  3.4782,  3.3170,  3.7434,  3.9623,  3.8181 &
     ,  3.7642,  2.6379,  2.8494,  3.1840,  3.4225,  3.2771,  3.3401 &
     ,  3.1072,  3.1885,  2.9714,  2.8319,  3.3315,  3.5979,  3.5256 &
     ,  3.4980,  3.4376,  3.6714,  3.4346,  3.2723,  3.6859,  3.8985 &
     ,  3.7918,  3.7372,  3.7211,  2.9230,  2.6223,  3.4161,  2.8999 &
     /)
     r0ab( 281: 350)=(/ &
        3.0557,  3.3308,  3.0555,  2.8508,  2.7385,  2.6640,  3.5263 &
     ,  3.0277,  3.2990,  3.7721,  3.5017,  3.2751,  3.1368,  3.0435 &
     ,  3.7873,  3.2858,  3.2140,  3.1727,  3.2178,  3.4414,  2.5490 &
     ,  2.7623,  3.0991,  3.3252,  3.1836,  3.2428,  3.0259,  3.1225 &
     ,  2.9032,  2.7621,  3.2490,  3.5110,  3.4429,  3.3845,  3.3574 &
     ,  3.6045,  3.3658,  3.2013,  3.6110,  3.8241,  3.7090,  3.6496 &
     ,  3.6333,  3.0896,  3.5462,  2.4926,  2.7136,  3.0693,  3.2699 &
     ,  3.1272,  3.1893,  2.9658,  3.0972,  2.8778,  2.7358,  3.2206 &
     ,  3.4566,  3.3896,  3.3257,  3.2946,  3.5693,  3.3312,  3.1670 &
     ,  3.5805,  3.7711,  3.6536,  3.5927,  3.5775,  3.0411,  3.4885 &
     /)
     r0ab( 351: 420)=(/ &
        3.4421,  2.4667,  2.6709,  3.0575,  3.2357,  3.0908,  3.1537 &
     ,  2.9235,  3.0669,  2.8476,  2.7054,  3.2064,  3.4519,  3.3593 &
     ,  3.2921,  3.2577,  3.2161,  3.2982,  3.1339,  3.5606,  3.7582 &
     ,  3.6432,  3.5833,  3.5691,  3.0161,  3.4812,  3.4339,  3.4327 &
     ,  2.4515,  2.6338,  3.0511,  3.2229,  3.0630,  3.1265,  2.8909 &
     ,  3.0253,  2.8184,  2.6764,  3.1968,  3.4114,  3.3492,  3.2691 &
     ,  3.2320,  3.1786,  3.2680,  3.1036,  3.5453,  3.7259,  3.6090 &
     ,  3.5473,  3.5327,  3.0018,  3.4413,  3.3907,  3.3593,  3.3462 &
     ,  2.4413,  2.6006,  3.0540,  3.1987,  3.0490,  3.1058,  2.8643 &
     ,  2.9948,  2.7908,  2.6491,  3.1950,  3.3922,  3.3316,  3.2585 &
     /)
     r0ab( 421: 490)=(/ &
        3.2136,  3.1516,  3.2364,  3.0752,  3.5368,  3.7117,  3.5941 &
     ,  3.5313,  3.5164,  2.9962,  3.4225,  3.3699,  3.3370,  3.3234 &
     ,  3.3008,  2.4318,  2.5729,  3.0416,  3.1639,  3.0196,  3.0843 &
     ,  2.8413,  2.7436,  2.7608,  2.6271,  3.1811,  3.3591,  3.3045 &
     ,  3.2349,  3.1942,  3.1291,  3.2111,  3.0534,  3.5189,  3.6809 &
     ,  3.5635,  3.5001,  3.4854,  2.9857,  3.3897,  3.3363,  3.3027 &
     ,  3.2890,  3.2655,  3.2309,  2.8502,  2.6934,  3.2467,  3.1921 &
     ,  3.5663,  3.2541,  3.0571,  2.9048,  2.8657,  2.7438,  3.3547 &
     ,  3.3510,  3.9837,  3.6871,  3.4862,  3.3389,  3.2413,  3.1708 &
     ,  3.6096,  3.6280,  3.6860,  3.5568,  3.4836,  3.2868,  3.3994 &
     /)
     r0ab( 491: 560)=(/ &
        3.3476,  3.3170,  3.2950,  3.2874,  3.2606,  3.9579,  2.9226 &
     ,  2.6838,  3.7867,  3.1732,  3.3872,  3.3643,  3.1267,  2.9541 &
     ,  2.8505,  2.7781,  3.8475,  3.3336,  3.7359,  3.8266,  3.5733 &
     ,  3.3959,  3.2775,  3.1915,  3.9878,  3.8816,  3.5810,  3.5364 &
     ,  3.5060,  3.8097,  3.3925,  3.3348,  3.3019,  3.2796,  3.2662 &
     ,  3.2464,  3.7136,  3.8619,  2.9140,  2.6271,  3.4771,  3.1774 &
     ,  3.2560,  3.1970,  3.1207,  2.9406,  2.8322,  2.7571,  3.5455 &
     ,  3.3514,  3.5837,  3.6177,  3.5816,  3.3902,  3.2604,  3.1652 &
     ,  3.7037,  3.6283,  3.5858,  3.5330,  3.4884,  3.5789,  3.4094 &
     ,  3.3473,  3.3118,  3.2876,  3.2707,  3.2521,  3.5570,  3.6496 &
     /)
     r0ab( 561: 630)=(/ &
        3.6625,  2.7300,  2.5870,  3.2471,  3.1487,  3.1667,  3.0914 &
     ,  3.0107,  2.9812,  2.8300,  2.7284,  3.3259,  3.3182,  3.4707 &
     ,  3.4748,  3.4279,  3.4182,  3.2547,  3.1353,  3.5116,  3.9432 &
     ,  3.8828,  3.8303,  3.7880,  3.3760,  3.7218,  3.3408,  3.3059 &
     ,  3.2698,  3.2446,  3.2229,  3.4422,  3.5023,  3.5009,  3.5268 &
     ,  2.6026,  2.5355,  3.1129,  3.2863,  3.1029,  3.0108,  2.9227 &
     ,  2.8694,  2.8109,  2.6929,  3.1958,  3.4670,  3.4018,  3.3805 &
     ,  3.3218,  3.2815,  3.2346,  3.0994,  3.3937,  3.7266,  3.6697 &
     ,  3.6164,  3.5730,  3.2522,  3.5051,  3.4686,  3.4355,  3.4084 &
     ,  3.3748,  3.3496,  3.3692,  3.4052,  3.3910,  3.3849,  3.3662 &
     /)
     r0ab( 631: 700)=(/ &
        2.5087,  2.4814,  3.0239,  3.1312,  3.0535,  2.9457,  2.8496 &
     ,  2.7780,  2.7828,  2.6532,  3.1063,  3.3143,  3.3549,  3.3120 &
     ,  3.2421,  3.1787,  3.1176,  3.0613,  3.3082,  3.5755,  3.5222 &
     ,  3.4678,  3.4231,  3.1684,  3.3528,  3.3162,  3.2827,  3.2527 &
     ,  3.2308,  3.2029,  3.3173,  3.3343,  3.3092,  3.2795,  3.2452 &
     ,  3.2096,  3.2893,  2.8991,  4.0388,  3.6100,  3.9388,  3.4475 &
     ,  3.1590,  2.9812,  2.8586,  2.7683,  4.1428,  3.7911,  3.8225 &
     ,  4.0372,  3.7059,  3.4935,  3.3529,  3.2492,  4.4352,  4.0826 &
     ,  3.9733,  3.9254,  3.8646,  3.9315,  3.7837,  3.7465,  3.7211 &
     ,  3.7012,  3.6893,  3.6676,  3.7736,  4.0660,  3.7926,  3.6158 &
     /)
     r0ab( 701: 770)=(/ &
        3.5017,  3.4166,  4.6176,  2.8786,  3.1658,  3.5823,  3.7689 &
     ,  3.5762,  3.5789,  3.3552,  3.4004,  3.1722,  3.0212,  3.7241 &
     ,  3.9604,  3.8500,  3.9844,  3.7035,  3.9161,  3.6751,  3.5075 &
     ,  4.1151,  4.2877,  4.1579,  4.1247,  4.0617,  3.4874,  3.9848 &
     ,  3.9280,  3.9079,  3.8751,  3.8604,  3.8277,  3.8002,  3.9981 &
     ,  3.7544,  4.0371,  3.8225,  3.6718,  4.3092,  4.4764,  2.8997 &
     ,  3.0953,  3.4524,  3.6107,  3.6062,  3.5783,  3.3463,  3.3855 &
     ,  3.1746,  3.0381,  3.6019,  3.7938,  3.8697,  3.9781,  3.6877 &
     ,  3.8736,  3.6451,  3.4890,  3.9858,  4.1179,  4.0430,  3.9563 &
     ,  3.9182,  3.4002,  3.8310,  3.7716,  3.7543,  3.7203,  3.7053 &
     /)
     r0ab( 771: 840)=(/ &
        3.6742,  3.8318,  3.7631,  3.7392,  3.9892,  3.7832,  3.6406 &
     ,  4.1701,  4.3016,  4.2196,  2.8535,  3.0167,  3.3978,  3.5363 &
     ,  3.5393,  3.5301,  3.2960,  3.3352,  3.1287,  2.9967,  3.6659 &
     ,  3.7239,  3.8070,  3.7165,  3.6368,  3.8162,  3.5885,  3.4336 &
     ,  3.9829,  4.0529,  3.9584,  3.9025,  3.8607,  3.3673,  3.7658 &
     ,  3.7035,  3.6866,  3.6504,  3.6339,  3.6024,  3.7708,  3.7283 &
     ,  3.6896,  3.9315,  3.7250,  3.5819,  4.1457,  4.2280,  4.1130 &
     ,  4.0597,  3.0905,  2.7998,  3.6448,  3.0739,  3.2996,  3.5262 &
     ,  3.2559,  3.0518,  2.9394,  2.8658,  3.7514,  3.2295,  3.5643 &
     ,  3.7808,  3.6931,  3.4723,  3.3357,  3.2429,  4.0280,  3.5589 &
     /)
     r0ab( 841: 910)=(/ &
        3.4636,  3.4994,  3.4309,  3.6177,  3.2946,  3.2376,  3.2050 &
     ,  3.1847,  3.1715,  3.1599,  3.5555,  3.8111,  3.7693,  3.5718 &
     ,  3.4498,  3.3662,  4.1608,  3.7417,  3.6536,  3.6154,  3.8596 &
     ,  3.0301,  2.7312,  3.5821,  3.0473,  3.2137,  3.4679,  3.1975 &
     ,  2.9969,  2.8847,  2.8110,  3.6931,  3.2076,  3.4943,  3.5956 &
     ,  3.6379,  3.4190,  3.2808,  3.1860,  3.9850,  3.5105,  3.4330 &
     ,  3.3797,  3.4155,  3.6033,  3.2737,  3.2145,  3.1807,  3.1596 &
     ,  3.1461,  3.1337,  3.4812,  3.6251,  3.7152,  3.5201,  3.3966 &
     ,  3.3107,  4.1128,  3.6899,  3.6082,  3.5604,  3.7834,  3.7543 &
     ,  2.9189,  2.6777,  3.4925,  2.9648,  3.1216,  3.2940,  3.0975 &
     /)
     r0ab( 911: 980)=(/ &
        2.9757,  2.8493,  2.7638,  3.6085,  3.1214,  3.4006,  3.4793 &
     ,  3.5147,  3.3806,  3.2356,  3.1335,  3.9144,  3.4183,  3.3369 &
     ,  3.2803,  3.2679,  3.4871,  3.1714,  3.1521,  3.1101,  3.0843 &
     ,  3.0670,  3.0539,  3.3890,  3.5086,  3.5895,  3.4783,  3.3484 &
     ,  3.2559,  4.0422,  3.5967,  3.5113,  3.4576,  3.6594,  3.6313 &
     ,  3.5690,  2.8578,  2.6334,  3.4673,  2.9245,  3.0732,  3.2435 &
     ,  3.0338,  2.9462,  2.8143,  2.7240,  3.5832,  3.0789,  3.3617 &
     ,  3.4246,  3.4505,  3.3443,  3.1964,  3.0913,  3.8921,  3.3713 &
     ,  3.2873,  3.2281,  3.2165,  3.4386,  3.1164,  3.1220,  3.0761 &
     ,  3.0480,  3.0295,  3.0155,  3.3495,  3.4543,  3.5260,  3.4413 &
     /)
     r0ab( 981:1050)=(/ &
        3.3085,  3.2134,  4.0170,  3.5464,  3.4587,  3.4006,  3.6027 &
     ,  3.5730,  3.4945,  3.4623,  2.8240,  2.5960,  3.4635,  2.9032 &
     ,  3.0431,  3.2115,  2.9892,  2.9148,  2.7801,  2.6873,  3.5776 &
     ,  3.0568,  3.3433,  3.3949,  3.4132,  3.3116,  3.1616,  3.0548 &
     ,  3.8859,  3.3719,  3.2917,  3.2345,  3.2274,  3.4171,  3.1293 &
     ,  3.0567,  3.0565,  3.0274,  3.0087,  2.9939,  3.3293,  3.4249 &
     ,  3.4902,  3.4091,  3.2744,  3.1776,  4.0078,  3.5374,  3.4537 &
     ,  3.3956,  3.5747,  3.5430,  3.4522,  3.4160,  3.3975,  2.8004 &
     ,  2.5621,  3.4617,  2.9154,  3.0203,  3.1875,  2.9548,  2.8038 &
     ,  2.7472,  2.6530,  3.5736,  3.0584,  3.3304,  3.3748,  3.3871 &
     /)
     r0ab(1051:1120)=(/ &
        3.2028,  3.1296,  3.0214,  3.8796,  3.3337,  3.2492,  3.1883 &
     ,  3.1802,  3.4050,  3.0756,  3.0478,  3.0322,  3.0323,  3.0163 &
     ,  3.0019,  3.3145,  3.4050,  3.4656,  3.3021,  3.2433,  3.1453 &
     ,  3.9991,  3.5017,  3.4141,  3.3520,  3.5583,  3.5251,  3.4243 &
     ,  3.3851,  3.3662,  3.3525,  2.7846,  2.5324,  3.4652,  2.8759 &
     ,  3.0051,  3.1692,  2.9273,  2.7615,  2.7164,  2.6212,  3.5744 &
     ,  3.0275,  3.3249,  3.3627,  3.3686,  3.1669,  3.0584,  2.9915 &
     ,  3.8773,  3.3099,  3.2231,  3.1600,  3.1520,  3.4023,  3.0426 &
     ,  3.0099,  2.9920,  2.9809,  2.9800,  2.9646,  3.3068,  3.3930 &
     ,  3.4486,  3.2682,  3.1729,  3.1168,  3.9952,  3.4796,  3.3901 &
     /)
     r0ab(1121:1190)=(/ &
        3.3255,  3.5530,  3.5183,  3.4097,  3.3683,  3.3492,  3.3360 &
     ,  3.3308,  2.5424,  2.6601,  3.2555,  3.2807,  3.1384,  3.1737 &
     ,  2.9397,  2.8429,  2.8492,  2.7225,  3.3875,  3.4910,  3.4520 &
     ,  3.3608,  3.3036,  3.2345,  3.2999,  3.1487,  3.7409,  3.8392 &
     ,  3.7148,  3.6439,  3.6182,  3.1753,  3.5210,  3.4639,  3.4265 &
     ,  3.4075,  3.3828,  3.3474,  3.4071,  3.3754,  3.3646,  3.3308 &
     ,  3.4393,  3.2993,  3.8768,  3.9891,  3.8310,  3.7483,  3.3417 &
     ,  3.3019,  3.2250,  3.1832,  3.1578,  3.1564,  3.1224,  3.4620 &
     ,  2.9743,  2.8058,  3.4830,  3.3474,  3.6863,  3.3617,  3.1608 &
     ,  3.0069,  2.9640,  2.8427,  3.5885,  3.5219,  4.1314,  3.8120 &
     /)
     r0ab(1191:1260)=(/ &
        3.6015,  3.4502,  3.3498,  3.2777,  3.8635,  3.8232,  3.8486 &
     ,  3.7215,  3.6487,  3.4724,  3.5627,  3.5087,  3.4757,  3.4517 &
     ,  3.4423,  3.4139,  4.1028,  3.8388,  3.6745,  3.5562,  3.4806 &
     ,  3.4272,  4.0182,  3.9991,  4.0007,  3.9282,  3.7238,  3.6498 &
     ,  3.5605,  3.5211,  3.5009,  3.4859,  3.4785,  3.5621,  4.2623 &
     ,  3.0775,  2.8275,  4.0181,  3.3385,  3.5379,  3.5036,  3.2589 &
     ,  3.0804,  3.0094,  2.9003,  4.0869,  3.5088,  3.9105,  3.9833 &
     ,  3.7176,  3.5323,  3.4102,  3.3227,  4.2702,  4.0888,  3.7560 &
     ,  3.7687,  3.6681,  3.6405,  3.5569,  3.4990,  3.4659,  3.4433 &
     ,  3.4330,  3.4092,  3.8867,  4.0190,  3.7961,  3.6412,  3.5405 &
     /)
     r0ab(1261:1330)=(/ &
        3.4681,  4.3538,  4.2136,  3.9381,  3.8912,  3.9681,  3.7909 &
     ,  3.6774,  3.6262,  3.5999,  3.5823,  3.5727,  3.5419,  4.0245 &
     ,  4.1874,  3.0893,  2.7917,  3.7262,  3.3518,  3.4241,  3.5433 &
     ,  3.2773,  3.0890,  2.9775,  2.9010,  3.8048,  3.5362,  3.7746 &
     ,  3.7911,  3.7511,  3.5495,  3.4149,  3.3177,  4.0129,  3.8370 &
     ,  3.7739,  3.7125,  3.7152,  3.7701,  3.5813,  3.5187,  3.4835 &
     ,  3.4595,  3.4439,  3.4242,  3.7476,  3.8239,  3.8346,  3.6627 &
     ,  3.5479,  3.4639,  4.1026,  3.9733,  3.9292,  3.8667,  3.9513 &
     ,  3.8959,  3.7698,  3.7089,  3.6765,  3.6548,  3.6409,  3.5398 &
     ,  3.8759,  3.9804,  4.0150,  2.9091,  2.7638,  3.5066,  3.3377 &
     /)
     r0ab(1331:1400)=(/ &
        3.3481,  3.2633,  3.1810,  3.1428,  2.9872,  2.8837,  3.5929 &
     ,  3.5183,  3.6729,  3.6596,  3.6082,  3.5927,  3.4224,  3.2997 &
     ,  3.8190,  4.1865,  4.1114,  4.0540,  3.6325,  3.5697,  3.5561 &
     ,  3.5259,  3.4901,  3.4552,  3.4315,  3.4091,  3.6438,  3.6879 &
     ,  3.6832,  3.7043,  3.5557,  3.4466,  3.9203,  4.2919,  4.2196 &
     ,  4.1542,  3.7573,  3.7039,  3.6546,  3.6151,  3.5293,  3.4849 &
     ,  3.4552,  3.5192,  3.7673,  3.8359,  3.8525,  3.8901,  2.7806 &
     ,  2.7209,  3.3812,  3.4958,  3.2913,  3.1888,  3.0990,  3.0394 &
     ,  2.9789,  2.8582,  3.4716,  3.6883,  3.6105,  3.5704,  3.5059 &
     ,  3.4619,  3.4138,  3.2742,  3.7080,  3.9773,  3.9010,  3.8409 &
     /)
     r0ab(1401:1470)=(/ &
        3.7944,  3.4465,  3.7235,  3.6808,  3.6453,  3.6168,  3.5844 &
     ,  3.5576,  3.5772,  3.5959,  3.5768,  3.5678,  3.5486,  3.4228 &
     ,  3.8107,  4.0866,  4.0169,  3.9476,  3.6358,  3.5800,  3.5260 &
     ,  3.4838,  3.4501,  3.4204,  3.3553,  3.6487,  3.6973,  3.7398 &
     ,  3.7405,  3.7459,  3.7380,  2.6848,  2.6740,  3.2925,  3.3386 &
     ,  3.2473,  3.1284,  3.0301,  2.9531,  2.9602,  2.8272,  3.3830 &
     ,  3.5358,  3.5672,  3.5049,  3.4284,  3.3621,  3.3001,  3.2451 &
     ,  3.6209,  3.8299,  3.7543,  3.6920,  3.6436,  3.3598,  3.5701 &
     ,  3.5266,  3.4904,  3.4590,  3.4364,  3.4077,  3.5287,  3.5280 &
     ,  3.4969,  3.4650,  3.4304,  3.3963,  3.7229,  3.9402,  3.8753 &
     /)
     r0ab(1471:1540)=(/ &
        3.8035,  3.5499,  3.4913,  3.4319,  3.3873,  3.3520,  3.3209 &
     ,  3.2948,  3.5052,  3.6465,  3.6696,  3.6577,  3.6388,  3.6142 &
     ,  3.5889,  3.3968,  3.0122,  4.2241,  3.7887,  4.0049,  3.5384 &
     ,  3.2698,  3.1083,  2.9917,  2.9057,  4.3340,  3.9900,  4.6588 &
     ,  4.1278,  3.8125,  3.6189,  3.4851,  3.3859,  4.6531,  4.3134 &
     ,  4.2258,  4.1309,  4.0692,  4.0944,  3.9850,  3.9416,  3.9112 &
     ,  3.8873,  3.8736,  3.8473,  4.6027,  4.1538,  3.8994,  3.7419 &
     ,  3.6356,  3.5548,  4.8353,  4.5413,  4.3891,  4.3416,  4.3243 &
     ,  4.2753,  4.2053,  4.1790,  4.1685,  4.1585,  4.1536,  4.0579 &
     ,  4.1980,  4.4564,  4.2192,  4.0528,  3.9489,  3.8642,  5.0567 &
     /)
     r0ab(1541:1610)=(/ &
        3.0630,  3.3271,  4.0432,  4.0046,  4.1555,  3.7426,  3.5130 &
     ,  3.5174,  3.2884,  3.1378,  4.1894,  4.2321,  4.1725,  4.1833 &
     ,  3.8929,  4.0544,  3.8118,  3.6414,  4.6373,  4.6268,  4.4750 &
     ,  4.4134,  4.3458,  3.8582,  4.2583,  4.1898,  4.1562,  4.1191 &
     ,  4.1069,  4.0639,  4.1257,  4.1974,  3.9532,  4.1794,  3.9660 &
     ,  3.8130,  4.8160,  4.8272,  4.6294,  4.5840,  4.0770,  4.0088 &
     ,  3.9103,  3.8536,  3.8324,  3.7995,  3.7826,  4.2294,  4.3380 &
     ,  4.4352,  4.1933,  4.4580,  4.2554,  4.1072,  5.0454,  5.1814 &
     ,  3.0632,  3.2662,  3.6432,  3.8088,  3.7910,  3.7381,  3.5093 &
     ,  3.5155,  3.3047,  3.1681,  3.7871,  3.9924,  4.0637,  4.1382 &
     /)
     r0ab(1611:1680)=(/ &
        3.8591,  4.0164,  3.7878,  3.6316,  4.1741,  4.3166,  4.2395 &
     ,  4.1831,  4.1107,  3.5857,  4.0270,  3.9676,  3.9463,  3.9150 &
     ,  3.9021,  3.8708,  4.0240,  4.1551,  3.9108,  4.1337,  3.9289 &
     ,  3.7873,  4.3666,  4.5080,  4.4232,  4.3155,  3.8461,  3.8007 &
     ,  3.6991,  3.6447,  3.6308,  3.5959,  3.5749,  4.0359,  4.3124 &
     ,  4.3539,  4.1122,  4.3772,  4.1785,  4.0386,  4.7004,  4.8604 &
     ,  4.6261,  2.9455,  3.2470,  3.6108,  3.8522,  3.6625,  3.6598 &
     ,  3.4411,  3.4660,  3.2415,  3.0944,  3.7514,  4.0397,  3.9231 &
     ,  4.0561,  3.7860,  3.9845,  3.7454,  3.5802,  4.1366,  4.3581 &
     ,  4.2351,  4.2011,  4.1402,  3.5381,  4.0653,  4.0093,  3.9883 &
     /)
     r0ab(1681:1750)=(/ &
        3.9570,  3.9429,  3.9112,  3.8728,  4.0682,  3.8351,  4.1054 &
     ,  3.8928,  3.7445,  4.3415,  4.5497,  4.3833,  4.3122,  3.8051 &
     ,  3.7583,  3.6622,  3.6108,  3.5971,  3.5628,  3.5408,  4.0780 &
     ,  4.0727,  4.2836,  4.0553,  4.3647,  4.1622,  4.0178,  4.5802 &
     ,  4.9125,  4.5861,  4.6201,  2.9244,  3.2241,  3.5848,  3.8293 &
     ,  3.6395,  3.6400,  3.4204,  3.4499,  3.2253,  3.0779,  3.7257 &
     ,  4.0170,  3.9003,  4.0372,  3.7653,  3.9672,  3.7283,  3.5630 &
     ,  4.1092,  4.3347,  4.2117,  4.1793,  4.1179,  3.5139,  4.0426 &
     ,  3.9867,  3.9661,  3.9345,  3.9200,  3.8883,  3.8498,  4.0496 &
     ,  3.8145,  4.0881,  3.8756,  3.7271,  4.3128,  4.5242,  4.3578 &
     /)
     r0ab(1751:1820)=(/ &
        4.2870,  3.7796,  3.7318,  3.6364,  3.5854,  3.5726,  3.5378 &
     ,  3.5155,  4.0527,  4.0478,  4.2630,  4.0322,  4.3449,  4.1421 &
     ,  3.9975,  4.5499,  4.8825,  4.5601,  4.5950,  4.5702,  2.9046 &
     ,  3.2044,  3.5621,  3.8078,  3.6185,  3.6220,  3.4019,  3.4359 &
     ,  3.2110,  3.0635,  3.7037,  3.9958,  3.8792,  4.0194,  3.7460 &
     ,  3.9517,  3.7128,  3.5474,  4.0872,  4.3138,  4.1906,  4.1593 &
     ,  4.0973,  3.4919,  4.0216,  3.9657,  3.9454,  3.9134,  3.8986 &
     ,  3.8669,  3.8289,  4.0323,  3.7954,  4.0725,  3.8598,  3.7113 &
     ,  4.2896,  4.5021,  4.3325,  4.2645,  3.7571,  3.7083,  3.6136 &
     ,  3.5628,  3.5507,  3.5155,  3.4929,  4.0297,  4.0234,  4.2442 &
     /)
     r0ab(1821:1890)=(/ &
        4.0112,  4.3274,  4.1240,  3.9793,  4.5257,  4.8568,  4.5353 &
     ,  4.5733,  4.5485,  4.5271,  2.8878,  3.1890,  3.5412,  3.7908 &
     ,  3.5974,  3.6078,  3.3871,  3.4243,  3.1992,  3.0513,  3.6831 &
     ,  3.9784,  3.8579,  4.0049,  3.7304,  3.9392,  3.7002,  3.5347 &
     ,  4.0657,  4.2955,  4.1705,  4.1424,  4.0800,  3.4717,  4.0043 &
     ,  3.9485,  3.9286,  3.8965,  3.8815,  3.8500,  3.8073,  4.0180 &
     ,  3.7796,  4.0598,  3.8470,  3.6983,  4.2678,  4.4830,  4.3132 &
     ,  4.2444,  3.7370,  3.6876,  3.5935,  3.5428,  3.5314,  3.4958 &
     ,  3.4730,  4.0117,  4.0043,  4.2287,  3.9939,  4.3134,  4.1096 &
     ,  3.9646,  4.5032,  4.8356,  4.5156,  4.5544,  4.5297,  4.5083 &
     /)
     r0ab(1891:1960)=(/ &
        4.4896,  2.8709,  3.1737,  3.5199,  3.7734,  3.5802,  3.5934 &
     ,  3.3724,  3.4128,  3.1877,  3.0396,  3.6624,  3.9608,  3.8397 &
     ,  3.9893,  3.7145,  3.9266,  3.6877,  3.5222,  4.0448,  4.2771 &
     ,  4.1523,  4.1247,  4.0626,  3.4530,  3.9866,  3.9310,  3.9115 &
     ,  3.8792,  3.8641,  3.8326,  3.7892,  4.0025,  3.7636,  4.0471 &
     ,  3.8343,  3.6854,  4.2464,  4.4635,  4.2939,  4.2252,  3.7169 &
     ,  3.6675,  3.5739,  3.5235,  3.5126,  3.4768,  3.4537,  3.9932 &
     ,  3.9854,  4.2123,  3.9765,  4.2992,  4.0951,  3.9500,  4.4811 &
     ,  4.8135,  4.4959,  4.5351,  4.5105,  4.4891,  4.4705,  4.4515 &
     ,  2.8568,  3.1608,  3.5050,  3.7598,  3.5665,  3.5803,  3.3601 &
     /)
     r0ab(1961:2030)=(/ &
        3.4031,  3.1779,  3.0296,  3.6479,  3.9471,  3.8262,  3.9773 &
     ,  3.7015,  3.9162,  3.6771,  3.5115,  4.0306,  4.2634,  4.1385 &
     ,  4.1116,  4.0489,  3.4366,  3.9732,  3.9176,  3.8983,  3.8659 &
     ,  3.8507,  3.8191,  3.7757,  3.9907,  3.7506,  4.0365,  3.8235 &
     ,  3.6745,  4.2314,  4.4490,  4.2792,  4.2105,  3.7003,  3.6510 &
     ,  3.5578,  3.5075,  3.4971,  3.4609,  3.4377,  3.9788,  3.9712 &
     ,  4.1997,  3.9624,  4.2877,  4.0831,  3.9378,  4.4655,  4.7974 &
     ,  4.4813,  4.5209,  4.4964,  4.4750,  4.4565,  4.4375,  4.4234 &
     ,  2.6798,  3.0151,  3.2586,  3.5292,  3.5391,  3.4902,  3.2887 &
     ,  3.3322,  3.1228,  2.9888,  3.4012,  3.7145,  3.7830,  3.6665 &
     /)
     r0ab(2031:2100)=(/ &
        3.5898,  3.8077,  3.5810,  3.4265,  3.7726,  4.0307,  3.9763 &
     ,  3.8890,  3.8489,  3.2706,  3.7595,  3.6984,  3.6772,  3.6428 &
     ,  3.6243,  3.5951,  3.7497,  3.6775,  3.6364,  3.9203,  3.7157 &
     ,  3.5746,  3.9494,  4.2076,  4.1563,  4.0508,  3.5329,  3.4780 &
     ,  3.3731,  3.3126,  3.2846,  3.2426,  3.2135,  3.7491,  3.9006 &
     ,  3.8332,  3.8029,  4.1436,  3.9407,  3.7998,  4.1663,  4.5309 &
     ,  4.3481,  4.2911,  4.2671,  4.2415,  4.2230,  4.2047,  4.1908 &
     ,  4.1243,  2.5189,  2.9703,  3.3063,  3.6235,  3.4517,  3.3989 &
     ,  3.2107,  3.2434,  3.0094,  2.8580,  3.4253,  3.8157,  3.7258 &
     ,  3.6132,  3.5297,  3.7566,  3.5095,  3.3368,  3.7890,  4.1298 &
     /)
     r0ab(2101:2170)=(/ &
        4.0190,  3.9573,  3.9237,  3.2677,  3.8480,  3.8157,  3.7656 &
     ,  3.7317,  3.7126,  3.6814,  3.6793,  3.6218,  3.5788,  3.8763 &
     ,  3.6572,  3.5022,  3.9737,  4.3255,  4.1828,  4.1158,  3.5078 &
     ,  3.4595,  3.3600,  3.3088,  3.2575,  3.2164,  3.1856,  3.8522 &
     ,  3.8665,  3.8075,  3.7772,  4.1391,  3.9296,  3.7772,  4.2134 &
     ,  4.7308,  4.3787,  4.3894,  4.3649,  4.3441,  4.3257,  4.3073 &
     ,  4.2941,  4.1252,  4.2427,  3.0481,  2.9584,  3.6919,  3.5990 &
     ,  3.8881,  3.4209,  3.1606,  3.1938,  2.9975,  2.8646,  3.8138 &
     ,  3.7935,  3.7081,  3.9155,  3.5910,  3.4808,  3.4886,  3.3397 &
     ,  4.1336,  4.1122,  3.9888,  3.9543,  3.8917,  3.5894,  3.8131 &
     /)
     r0ab(2171:2240)=(/ &
        3.7635,  3.7419,  3.7071,  3.6880,  3.6574,  3.6546,  3.9375 &
     ,  3.6579,  3.5870,  3.6361,  3.5039,  4.3149,  4.2978,  4.1321 &
     ,  4.1298,  3.8164,  3.7680,  3.7154,  3.6858,  3.6709,  3.6666 &
     ,  3.6517,  3.8174,  3.8608,  4.1805,  3.9102,  3.8394,  3.8968 &
     ,  3.7673,  4.5274,  4.6682,  4.3344,  4.3639,  4.3384,  4.3162 &
     ,  4.2972,  4.2779,  4.2636,  4.0253,  4.1168,  4.1541,  2.8136 &
     ,  3.0951,  3.4635,  3.6875,  3.4987,  3.5183,  3.2937,  3.3580 &
     ,  3.1325,  2.9832,  3.6078,  3.8757,  3.7616,  3.9222,  3.6370 &
     ,  3.8647,  3.6256,  3.4595,  3.9874,  4.1938,  4.0679,  4.0430 &
     ,  3.9781,  3.3886,  3.9008,  3.8463,  3.8288,  3.7950,  3.7790 &
     /)
     r0ab(2241:2310)=(/ &
        3.7472,  3.7117,  3.9371,  3.6873,  3.9846,  3.7709,  3.6210 &
     ,  4.1812,  4.3750,  4.2044,  4.1340,  3.6459,  3.5929,  3.5036 &
     ,  3.4577,  3.4528,  3.4146,  3.3904,  3.9014,  3.9031,  4.1443 &
     ,  3.8961,  4.2295,  4.0227,  3.8763,  4.4086,  4.7097,  4.4064 &
     ,  4.4488,  4.4243,  4.4029,  4.3842,  4.3655,  4.3514,  4.1162 &
     ,  4.2205,  4.1953,  4.2794,  2.8032,  3.0805,  3.4519,  3.6700 &
     ,  3.4827,  3.5050,  3.2799,  3.3482,  3.1233,  2.9747,  3.5971 &
     ,  3.8586,  3.7461,  3.9100,  3.6228,  3.8535,  3.6147,  3.4490 &
     ,  3.9764,  4.1773,  4.0511,  4.0270,  3.9614,  3.3754,  3.8836 &
     ,  3.8291,  3.8121,  3.7780,  3.7619,  3.7300,  3.6965,  3.9253 &
     /)
     r0ab(2311:2380)=(/ &
        3.6734,  3.9733,  3.7597,  3.6099,  4.1683,  4.3572,  4.1862 &
     ,  4.1153,  3.6312,  3.5772,  3.4881,  3.4429,  3.4395,  3.4009 &
     ,  3.3766,  3.8827,  3.8868,  4.1316,  3.8807,  4.2164,  4.0092 &
     ,  3.8627,  4.3936,  4.6871,  4.3882,  4.4316,  4.4073,  4.3858 &
     ,  4.3672,  4.3485,  4.3344,  4.0984,  4.2036,  4.1791,  4.2622 &
     ,  4.2450,  2.7967,  3.0689,  3.4445,  3.6581,  3.4717,  3.4951 &
     ,  3.2694,  3.3397,  3.1147,  2.9661,  3.5898,  3.8468,  3.7358 &
     ,  3.9014,  3.6129,  3.8443,  3.6054,  3.4396,  3.9683,  4.1656 &
     ,  4.0394,  4.0158,  3.9498,  3.3677,  3.8718,  3.8164,  3.8005 &
     ,  3.7662,  3.7500,  3.7181,  3.6863,  3.9170,  3.6637,  3.9641 &
     /)
     r0ab(2381:2450)=(/ &
        3.7503,  3.6004,  4.1590,  4.3448,  4.1739,  4.1029,  3.6224 &
     ,  3.5677,  3.4785,  3.4314,  3.4313,  3.3923,  3.3680,  3.8698 &
     ,  3.8758,  4.1229,  3.8704,  4.2063,  3.9987,  3.8519,  4.3832 &
     ,  4.6728,  4.3759,  4.4195,  4.3952,  4.3737,  4.3551,  4.3364 &
     ,  4.3223,  4.0861,  4.1911,  4.1676,  4.2501,  4.2329,  4.2208 &
     ,  2.7897,  3.0636,  3.4344,  3.6480,  3.4626,  3.4892,  3.2626 &
     ,  3.3344,  3.1088,  2.9597,  3.5804,  3.8359,  3.7251,  3.8940 &
     ,  3.6047,  3.8375,  3.5990,  3.4329,  3.9597,  4.1542,  4.0278 &
     ,  4.0048,  3.9390,  3.3571,  3.8608,  3.8056,  3.7899,  3.7560 &
     ,  3.7400,  3.7081,  3.6758,  3.9095,  3.6552,  3.9572,  3.7436 &
     /)
     r0ab(2451:2520)=(/ &
        3.5933,  4.1508,  4.3337,  4.1624,  4.0916,  3.6126,  3.5582 &
     ,  3.4684,  3.4212,  3.4207,  3.3829,  3.3586,  3.8604,  3.8658 &
     ,  4.1156,  3.8620,  4.1994,  3.9917,  3.8446,  4.3750,  4.6617 &
     ,  4.3644,  4.4083,  4.3840,  4.3625,  4.3439,  4.3253,  4.3112 &
     ,  4.0745,  4.1807,  4.1578,  4.2390,  4.2218,  4.2097,  4.1986 &
     ,  2.8395,  3.0081,  3.3171,  3.4878,  3.5360,  3.5145,  3.2809 &
     ,  3.3307,  3.1260,  2.9940,  3.4741,  3.6675,  3.7832,  3.6787 &
     ,  3.6156,  3.8041,  3.5813,  3.4301,  3.8480,  3.9849,  3.9314 &
     ,  3.8405,  3.8029,  3.2962,  3.7104,  3.6515,  3.6378,  3.6020 &
     ,  3.5849,  3.5550,  3.7494,  3.6893,  3.6666,  3.9170,  3.7150 &
     /)
     r0ab(2521:2590)=(/ &
        3.5760,  4.0268,  4.1596,  4.1107,  3.9995,  3.5574,  3.5103 &
     ,  3.4163,  3.3655,  3.3677,  3.3243,  3.2975,  3.7071,  3.9047 &
     ,  3.8514,  3.8422,  3.8022,  3.9323,  3.7932,  4.2343,  4.4583 &
     ,  4.3115,  4.2457,  4.2213,  4.1945,  4.1756,  4.1569,  4.1424 &
     ,  4.0620,  4.0494,  3.9953,  4.0694,  4.0516,  4.0396,  4.0280 &
     ,  4.0130,  2.9007,  2.9674,  3.8174,  3.5856,  3.6486,  3.5339 &
     ,  3.2832,  3.3154,  3.1144,  2.9866,  3.9618,  3.8430,  3.9980 &
     ,  3.8134,  3.6652,  3.7985,  3.5756,  3.4207,  4.4061,  4.2817 &
     ,  4.1477,  4.0616,  3.9979,  3.6492,  3.8833,  3.8027,  3.7660 &
     ,  3.7183,  3.6954,  3.6525,  3.9669,  3.8371,  3.7325,  3.9160 &
     /)
     r0ab(2591:2660)=(/ &
        3.7156,  3.5714,  4.6036,  4.4620,  4.3092,  4.2122,  3.8478 &
     ,  3.7572,  3.6597,  3.5969,  3.5575,  3.5386,  3.5153,  3.7818 &
     ,  4.1335,  4.0153,  3.9177,  3.8603,  3.9365,  3.7906,  4.7936 &
     ,  4.7410,  4.5461,  4.5662,  4.5340,  4.5059,  4.4832,  4.4604 &
     ,  4.4429,  4.2346,  4.4204,  4.3119,  4.3450,  4.3193,  4.3035 &
     ,  4.2933,  4.1582,  4.2450,  2.8559,  2.9050,  3.8325,  3.5442 &
     ,  3.5077,  3.4905,  3.2396,  3.2720,  3.0726,  2.9467,  3.9644 &
     ,  3.8050,  3.8981,  3.7762,  3.6216,  3.7531,  3.5297,  3.3742 &
     ,  4.3814,  4.2818,  4.1026,  4.0294,  3.9640,  3.6208,  3.8464 &
     ,  3.7648,  3.7281,  3.6790,  3.6542,  3.6117,  3.8650,  3.8010 &
     /)
     r0ab(2661:2730)=(/ &
        3.6894,  3.8713,  3.6699,  3.5244,  4.5151,  4.4517,  4.2538 &
     ,  4.1483,  3.8641,  3.7244,  3.6243,  3.5589,  3.5172,  3.4973 &
     ,  3.4715,  3.7340,  4.0316,  3.9958,  3.8687,  3.8115,  3.8862 &
     ,  3.7379,  4.7091,  4.7156,  4.5199,  4.5542,  4.5230,  4.4959 &
     ,  4.4750,  4.4529,  4.4361,  4.1774,  4.3774,  4.2963,  4.3406 &
     ,  4.3159,  4.3006,  4.2910,  4.1008,  4.1568,  4.0980,  2.8110 &
     ,  2.8520,  3.7480,  3.5105,  3.4346,  3.3461,  3.1971,  3.2326 &
     ,  3.0329,  2.9070,  3.8823,  3.7928,  3.8264,  3.7006,  3.5797 &
     ,  3.7141,  3.4894,  3.3326,  4.3048,  4.2217,  4.0786,  3.9900 &
     ,  3.9357,  3.6331,  3.8333,  3.7317,  3.6957,  3.6460,  3.6197 &
     /)
     r0ab(2731:2800)=(/ &
        3.5779,  3.7909,  3.7257,  3.6476,  3.5729,  3.6304,  3.4834 &
     ,  4.4368,  4.3921,  4.2207,  4.1133,  3.8067,  3.7421,  3.6140 &
     ,  3.5491,  3.5077,  3.4887,  3.4623,  3.6956,  3.9568,  3.8976 &
     ,  3.8240,  3.7684,  3.8451,  3.6949,  4.6318,  4.6559,  4.4533 &
     ,  4.4956,  4.4641,  4.4366,  4.4155,  4.3936,  4.3764,  4.1302 &
     ,  4.3398,  4.2283,  4.2796,  4.2547,  4.2391,  4.2296,  4.0699 &
     ,  4.1083,  4.0319,  3.9855,  2.7676,  2.8078,  3.6725,  3.4804 &
     ,  3.3775,  3.2411,  3.1581,  3.1983,  2.9973,  2.8705,  3.8070 &
     ,  3.7392,  3.7668,  3.6263,  3.5402,  3.6807,  3.4545,  3.2962 &
     ,  4.2283,  4.1698,  4.0240,  3.9341,  3.8711,  3.5489,  3.7798 &
     /)
     r0ab(2801:2870)=(/ &
        3.7000,  3.6654,  3.6154,  3.5882,  3.5472,  3.7289,  3.6510 &
     ,  3.6078,  3.5355,  3.5963,  3.4480,  4.3587,  4.3390,  4.1635 &
     ,  4.0536,  3.7193,  3.6529,  3.5512,  3.4837,  3.4400,  3.4191 &
     ,  3.3891,  3.6622,  3.8934,  3.8235,  3.7823,  3.7292,  3.8106 &
     ,  3.6589,  4.5535,  4.6013,  4.3961,  4.4423,  4.4109,  4.3835 &
     ,  4.3625,  4.3407,  4.3237,  4.0863,  4.2835,  4.1675,  4.2272 &
     ,  4.2025,  4.1869,  4.1774,  4.0126,  4.0460,  3.9815,  3.9340 &
     ,  3.8955,  2.6912,  2.7604,  3.6037,  3.4194,  3.3094,  3.1710 &
     ,  3.0862,  3.1789,  2.9738,  2.8427,  3.7378,  3.6742,  3.6928 &
     ,  3.5512,  3.4614,  3.4087,  3.4201,  3.2607,  4.1527,  4.0977 &
     /)
     r0ab(2871:2940)=(/ &
        3.9523,  3.8628,  3.8002,  3.4759,  3.7102,  3.6466,  3.6106 &
     ,  3.5580,  3.5282,  3.4878,  3.6547,  3.5763,  3.5289,  3.5086 &
     ,  3.5593,  3.4099,  4.2788,  4.2624,  4.0873,  3.9770,  3.6407 &
     ,  3.5743,  3.5178,  3.4753,  3.3931,  3.3694,  3.3339,  3.6002 &
     ,  3.8164,  3.7478,  3.7028,  3.6952,  3.7669,  3.6137,  4.4698 &
     ,  4.5488,  4.3168,  4.3646,  4.3338,  4.3067,  4.2860,  4.2645 &
     ,  4.2478,  4.0067,  4.2349,  4.0958,  4.1543,  4.1302,  4.1141 &
     ,  4.1048,  3.9410,  3.9595,  3.8941,  3.8465,  3.8089,  3.7490 &
     ,  2.7895,  2.5849,  3.6484,  3.0162,  3.1267,  3.2125,  3.0043 &
     ,  2.9572,  2.8197,  2.7261,  3.7701,  3.2446,  3.5239,  3.4696 &
     /)
     r0ab(2941:3010)=(/ &
        3.4261,  3.3508,  3.1968,  3.0848,  4.1496,  3.6598,  3.5111 &
     ,  3.4199,  3.3809,  3.5382,  3.2572,  3.2100,  3.1917,  3.1519 &
     ,  3.1198,  3.1005,  3.5071,  3.5086,  3.5073,  3.4509,  3.3120 &
     ,  3.2082,  4.2611,  3.8117,  3.6988,  3.5646,  3.6925,  3.6295 &
     ,  3.5383,  3.4910,  3.4625,  3.4233,  3.4007,  3.2329,  3.6723 &
     ,  3.6845,  3.6876,  3.6197,  3.4799,  3.3737,  4.4341,  4.0525 &
     ,  3.9011,  3.8945,  3.8635,  3.8368,  3.8153,  3.7936,  3.7758 &
     ,  3.4944,  3.4873,  3.9040,  3.7110,  3.6922,  3.6799,  3.6724 &
     ,  3.5622,  3.6081,  3.5426,  3.4922,  3.4498,  3.3984,  3.4456 &
     ,  2.7522,  2.5524,  3.5742,  2.9508,  3.0751,  3.0158,  2.9644 &
     /)
     r0ab(3011:3080)=(/ &
        2.8338,  2.7891,  2.6933,  3.6926,  3.1814,  3.4528,  3.4186 &
     ,  3.3836,  3.2213,  3.1626,  3.0507,  4.0548,  3.5312,  3.4244 &
     ,  3.3409,  3.2810,  3.4782,  3.1905,  3.1494,  3.1221,  3.1128 &
     ,  3.0853,  3.0384,  3.4366,  3.4562,  3.4638,  3.3211,  3.2762 &
     ,  3.1730,  4.1632,  3.6825,  3.5822,  3.4870,  3.6325,  3.5740 &
     ,  3.4733,  3.4247,  3.3969,  3.3764,  3.3525,  3.1984,  3.5989 &
     ,  3.6299,  3.6433,  3.4937,  3.4417,  3.3365,  4.3304,  3.9242 &
     ,  3.7793,  3.7623,  3.7327,  3.7071,  3.6860,  3.6650,  3.6476 &
     ,  3.3849,  3.3534,  3.8216,  3.5870,  3.5695,  3.5584,  3.5508 &
     ,  3.4856,  3.5523,  3.4934,  3.4464,  3.4055,  3.3551,  3.3888 &
     /)
     r0ab(3081:3150)=(/ &
        3.3525,  2.7202,  2.5183,  3.4947,  2.8731,  3.0198,  3.1457 &
     ,  2.9276,  2.7826,  2.7574,  2.6606,  3.6090,  3.0581,  3.3747 &
     ,  3.3677,  3.3450,  3.1651,  3.1259,  3.0147,  3.9498,  3.3857 &
     ,  3.2917,  3.2154,  3.1604,  3.4174,  3.0735,  3.0342,  3.0096 &
     ,  3.0136,  2.9855,  2.9680,  3.3604,  3.4037,  3.4243,  3.2633 &
     ,  3.1810,  3.1351,  4.0557,  3.5368,  3.4526,  3.3699,  3.5707 &
     ,  3.5184,  3.4085,  3.3595,  3.3333,  3.3143,  3.3041,  3.1094 &
     ,  3.5193,  3.5745,  3.6025,  3.4338,  3.3448,  3.2952,  4.2158 &
     ,  3.7802,  3.6431,  3.6129,  3.5853,  3.5610,  3.5406,  3.5204 &
     ,  3.5036,  3.2679,  3.2162,  3.7068,  3.4483,  3.4323,  3.4221 &
     /)
     r0ab(3151:3220)=(/ &
        3.4138,  3.3652,  3.4576,  3.4053,  3.3618,  3.3224,  3.2711 &
     ,  3.3326,  3.2950,  3.2564,  2.5315,  2.6104,  3.2734,  3.2299 &
     ,  3.1090,  2.9942,  2.9159,  2.8324,  2.8350,  2.7216,  3.3994 &
     ,  3.4475,  3.4354,  3.3438,  3.2807,  3.2169,  3.2677,  3.1296 &
     ,  3.7493,  3.8075,  3.6846,  3.6104,  3.5577,  3.2052,  3.4803 &
     ,  3.4236,  3.3845,  3.3640,  3.3365,  3.3010,  3.3938,  3.3624 &
     ,  3.3440,  3.3132,  3.4035,  3.2754,  3.8701,  3.9523,  3.8018 &
     ,  3.7149,  3.3673,  3.3199,  3.2483,  3.2069,  3.1793,  3.1558 &
     ,  3.1395,  3.4097,  3.5410,  3.5228,  3.5116,  3.4921,  3.4781 &
     ,  3.4690,  4.0420,  4.1759,  4.0078,  4.0450,  4.0189,  3.9952 &
     /)
     r0ab(3221:3290)=(/ &
        3.9770,  3.9583,  3.9434,  3.7217,  3.8228,  3.7826,  3.8640 &
     ,  3.8446,  3.8314,  3.8225,  3.6817,  3.7068,  3.6555,  3.6159 &
     ,  3.5831,  3.5257,  3.2133,  3.1689,  3.1196,  3.3599,  2.9852 &
     ,  2.7881,  3.5284,  3.3493,  3.6958,  3.3642,  3.1568,  3.0055 &
     ,  2.9558,  2.8393,  3.6287,  3.5283,  4.1511,  3.8259,  3.6066 &
     ,  3.4527,  3.3480,  3.2713,  3.9037,  3.8361,  3.8579,  3.7311 &
     ,  3.6575,  3.5176,  3.5693,  3.5157,  3.4814,  3.4559,  3.4445 &
     ,  3.4160,  4.1231,  3.8543,  3.6816,  3.5602,  3.4798,  3.4208 &
     ,  4.0542,  4.0139,  4.0165,  3.9412,  3.7698,  3.6915,  3.6043 &
     ,  3.5639,  3.5416,  3.5247,  3.5153,  3.5654,  4.2862,  4.0437 &
     /)
     r0ab(3291:3360)=(/ &
        3.8871,  3.7741,  3.6985,  3.6413,  4.2345,  4.3663,  4.3257 &
     ,  4.0869,  4.0612,  4.0364,  4.0170,  3.9978,  3.9834,  3.9137 &
     ,  3.8825,  3.8758,  3.9143,  3.8976,  3.8864,  3.8768,  3.9190 &
     ,  4.1613,  4.0566,  3.9784,  3.9116,  3.8326,  3.7122,  3.6378 &
     ,  3.5576,  3.5457,  4.3127,  3.1160,  2.8482,  4.0739,  3.3599 &
     ,  3.5698,  3.5366,  3.2854,  3.1039,  2.9953,  2.9192,  4.1432 &
     ,  3.5320,  3.9478,  4.0231,  3.7509,  3.5604,  3.4340,  3.3426 &
     ,  4.3328,  3.8288,  3.7822,  3.7909,  3.6907,  3.6864,  3.5793 &
     ,  3.5221,  3.4883,  3.4649,  3.4514,  3.4301,  3.9256,  4.0596 &
     ,  3.8307,  3.6702,  3.5651,  3.4884,  4.4182,  4.2516,  3.9687 &
     /)
     r0ab(3361:3430)=(/ &
        3.9186,  3.9485,  3.8370,  3.7255,  3.6744,  3.6476,  3.6295 &
     ,  3.6193,  3.5659,  4.0663,  4.2309,  4.0183,  3.8680,  3.7672 &
     ,  3.6923,  4.5240,  4.4834,  4.1570,  4.3204,  4.2993,  4.2804 &
     ,  4.2647,  4.2481,  4.2354,  3.8626,  3.8448,  4.2267,  4.1799 &
     ,  4.1670,  3.8738,  3.8643,  3.8796,  4.0575,  4.0354,  3.9365 &
     ,  3.8611,  3.7847,  3.7388,  3.6826,  3.6251,  3.5492,  4.0889 &
     ,  4.2764,  3.1416,  2.8325,  3.7735,  3.3787,  3.4632,  3.5923 &
     ,  3.3214,  3.1285,  3.0147,  2.9366,  3.8527,  3.5602,  3.8131 &
     ,  3.8349,  3.7995,  3.5919,  3.4539,  3.3540,  4.0654,  3.8603 &
     ,  3.7972,  3.7358,  3.7392,  3.8157,  3.6055,  3.5438,  3.5089 &
     /)
     r0ab(3431:3500)=(/ &
        3.4853,  3.4698,  3.4508,  3.7882,  3.8682,  3.8837,  3.7055 &
     ,  3.5870,  3.5000,  4.1573,  4.0005,  3.9568,  3.8936,  3.9990 &
     ,  3.9433,  3.8172,  3.7566,  3.7246,  3.7033,  3.6900,  3.5697 &
     ,  3.9183,  4.0262,  4.0659,  3.8969,  3.7809,  3.6949,  4.2765 &
     ,  4.2312,  4.1401,  4.0815,  4.0580,  4.0369,  4.0194,  4.0017 &
     ,  3.9874,  3.8312,  3.8120,  3.9454,  3.9210,  3.9055,  3.8951 &
     ,  3.8866,  3.8689,  3.9603,  3.9109,  3.9122,  3.8233,  3.7438 &
     ,  3.7436,  3.6981,  3.6555,  3.5452,  3.9327,  4.0658,  4.1175 &
     ,  2.9664,  2.8209,  3.5547,  3.3796,  3.3985,  3.3164,  3.2364 &
     ,  3.1956,  3.0370,  2.9313,  3.6425,  3.5565,  3.7209,  3.7108 &
     /)
     r0ab(3501:3570)=(/ &
        3.6639,  3.6484,  3.4745,  3.3492,  3.8755,  4.2457,  3.7758 &
     ,  3.7161,  3.6693,  3.6155,  3.5941,  3.5643,  3.5292,  3.4950 &
     ,  3.4720,  3.4503,  3.6936,  3.7392,  3.7388,  3.7602,  3.6078 &
     ,  3.4960,  3.9800,  4.3518,  4.2802,  3.8580,  3.8056,  3.7527 &
     ,  3.7019,  3.6615,  3.5768,  3.5330,  3.5038,  3.5639,  3.8192 &
     ,  3.8883,  3.9092,  3.9478,  3.7995,  3.6896,  4.1165,  4.5232 &
     ,  4.4357,  4.4226,  4.4031,  4.3860,  4.3721,  4.3580,  4.3466 &
     ,  4.2036,  4.2037,  3.8867,  4.2895,  4.2766,  4.2662,  4.2598 &
     ,  3.8408,  3.9169,  3.8681,  3.8250,  3.7855,  3.7501,  3.6753 &
     ,  3.5499,  3.4872,  3.5401,  3.8288,  3.9217,  3.9538,  4.0054 &
     /)
     r0ab(3571:3640)=(/ &
        2.8388,  2.7890,  3.4329,  3.5593,  3.3488,  3.2486,  3.1615 &
     ,  3.1000,  3.0394,  2.9165,  3.5267,  3.7479,  3.6650,  3.6263 &
     ,  3.5658,  3.5224,  3.4762,  3.3342,  3.7738,  4.0333,  3.9568 &
     ,  3.8975,  3.8521,  3.4929,  3.7830,  3.7409,  3.7062,  3.6786 &
     ,  3.6471,  3.6208,  3.6337,  3.6519,  3.6363,  3.6278,  3.6110 &
     ,  3.4825,  3.8795,  4.1448,  4.0736,  4.0045,  3.6843,  3.6291 &
     ,  3.5741,  3.5312,  3.4974,  3.4472,  3.4034,  3.7131,  3.7557 &
     ,  3.7966,  3.8005,  3.8068,  3.8015,  3.6747,  4.0222,  4.3207 &
     ,  4.2347,  4.2191,  4.1990,  4.1811,  4.1666,  4.1521,  4.1401 &
     ,  3.9970,  3.9943,  3.9592,  4.0800,  4.0664,  4.0559,  4.0488 &
     /)
     r0ab(3641:3710)=(/ &
        3.9882,  4.0035,  3.9539,  3.9138,  3.8798,  3.8355,  3.5359 &
     ,  3.4954,  3.3962,  3.5339,  3.7595,  3.8250,  3.8408,  3.8600 &
     ,  3.8644,  2.7412,  2.7489,  3.3374,  3.3950,  3.3076,  3.1910 &
     ,  3.0961,  3.0175,  3.0280,  2.8929,  3.4328,  3.5883,  3.6227 &
     ,  3.5616,  3.4894,  3.4241,  3.3641,  3.3120,  3.6815,  3.8789 &
     ,  3.8031,  3.7413,  3.6939,  3.4010,  3.6225,  3.5797,  3.5443 &
     ,  3.5139,  3.4923,  3.4642,  3.5860,  3.5849,  3.5570,  3.5257 &
     ,  3.4936,  3.4628,  3.7874,  3.9916,  3.9249,  3.8530,  3.5932 &
     ,  3.5355,  3.4757,  3.4306,  3.3953,  3.3646,  3.3390,  3.5637 &
     ,  3.7053,  3.7266,  3.7177,  3.6996,  3.6775,  3.6558,  3.9331 &
     /)
     r0ab(3711:3780)=(/ &
        4.1655,  4.0879,  4.0681,  4.0479,  4.0299,  4.0152,  4.0006 &
     ,  3.9883,  3.8500,  3.8359,  3.8249,  3.9269,  3.9133,  3.9025 &
     ,  3.8948,  3.8422,  3.8509,  3.7990,  3.7570,  3.7219,  3.6762 &
     ,  3.4260,  3.3866,  3.3425,  3.5294,  3.7022,  3.7497,  3.7542 &
     ,  3.7494,  3.7370,  3.7216,  3.4155,  3.0522,  4.2541,  3.8218 &
     ,  4.0438,  3.5875,  3.3286,  3.1682,  3.0566,  2.9746,  4.3627 &
     ,  4.0249,  4.6947,  4.1718,  3.8639,  3.6735,  3.5435,  3.4479 &
     ,  4.6806,  4.3485,  4.2668,  4.1690,  4.1061,  4.1245,  4.0206 &
     ,  3.9765,  3.9458,  3.9217,  3.9075,  3.8813,  3.9947,  4.1989 &
     ,  3.9507,  3.7960,  3.6925,  3.6150,  4.8535,  4.5642,  4.4134 &
     /)
     r0ab(3781:3850)=(/ &
        4.3688,  4.3396,  4.2879,  4.2166,  4.1888,  4.1768,  4.1660 &
     ,  4.1608,  4.0745,  4.2289,  4.4863,  4.2513,  4.0897,  3.9876 &
     ,  3.9061,  5.0690,  5.0446,  4.6186,  4.6078,  4.5780,  4.5538 &
     ,  4.5319,  4.5101,  4.4945,  4.1912,  4.2315,  4.5534,  4.4373 &
     ,  4.4224,  4.4120,  4.4040,  4.2634,  4.7770,  4.6890,  4.6107 &
     ,  4.5331,  4.4496,  4.4082,  4.3095,  4.2023,  4.0501,  4.2595 &
     ,  4.5497,  4.3056,  4.1506,  4.0574,  3.9725,  5.0796,  3.0548 &
     ,  3.3206,  3.8132,  3.9720,  3.7675,  3.7351,  3.5167,  3.5274 &
     ,  3.3085,  3.1653,  3.9500,  4.1730,  4.0613,  4.1493,  3.8823 &
     ,  4.0537,  3.8200,  3.6582,  4.3422,  4.5111,  4.3795,  4.3362 &
     /)
     r0ab(3851:3920)=(/ &
        4.2751,  3.7103,  4.1973,  4.1385,  4.1129,  4.0800,  4.0647 &
     ,  4.0308,  4.0096,  4.1619,  3.9360,  4.1766,  3.9705,  3.8262 &
     ,  4.5348,  4.7025,  4.5268,  4.5076,  3.9562,  3.9065,  3.8119 &
     ,  3.7605,  3.7447,  3.7119,  3.6916,  4.1950,  4.2110,  4.3843 &
     ,  4.1631,  4.4427,  4.2463,  4.1054,  4.7693,  5.0649,  4.7365 &
     ,  4.7761,  4.7498,  4.7272,  4.7076,  4.6877,  4.6730,  4.4274 &
     ,  4.5473,  4.5169,  4.5975,  4.5793,  4.5667,  4.5559,  4.3804 &
     ,  4.6920,  4.6731,  4.6142,  4.5600,  4.4801,  4.0149,  3.8856 &
     ,  3.7407,  4.1545,  4.2253,  4.4229,  4.1923,  4.5022,  4.3059 &
     ,  4.1591,  4.7883,  4.9294,  3.3850,  3.4208,  3.7004,  3.8800 &
     /)
     r0ab(3921:3990)=(/ &
        3.9886,  3.9040,  3.6719,  3.6547,  3.4625,  3.3370,  3.8394 &
     ,  4.0335,  4.2373,  4.3023,  4.0306,  4.1408,  3.9297,  3.7857 &
     ,  4.1907,  4.3230,  4.2664,  4.2173,  4.1482,  3.6823,  4.0711 &
     ,  4.0180,  4.0017,  3.9747,  3.9634,  3.9383,  4.1993,  4.3205 &
     ,  4.0821,  4.2547,  4.0659,  3.9359,  4.3952,  4.5176,  4.3888 &
     ,  4.3607,  3.9583,  3.9280,  3.8390,  3.7971,  3.7955,  3.7674 &
     ,  3.7521,  4.1062,  4.3633,  4.2991,  4.2767,  4.4857,  4.3039 &
     ,  4.1762,  4.6197,  4.8654,  4.6633,  4.5878,  4.5640,  4.5422 &
     ,  4.5231,  4.5042,  4.4901,  4.3282,  4.3978,  4.3483,  4.4202 &
     ,  4.4039,  4.3926,  4.3807,  4.2649,  4.6135,  4.5605,  4.5232 &
     /)
     r0ab(3991:4060)=(/ &
        4.4676,  4.3948,  4.0989,  3.9864,  3.8596,  4.0942,  4.2720 &
     ,  4.3270,  4.3022,  4.5410,  4.3576,  4.2235,  4.6545,  4.7447 &
     ,  4.7043,  3.0942,  3.2075,  3.5152,  3.6659,  3.8289,  3.7459 &
     ,  3.5156,  3.5197,  3.3290,  3.2069,  3.6702,  3.8448,  4.0340 &
     ,  3.9509,  3.8585,  3.9894,  3.7787,  3.6365,  4.1425,  4.1618 &
     ,  4.0940,  4.0466,  3.9941,  3.5426,  3.8952,  3.8327,  3.8126 &
     ,  3.7796,  3.7635,  3.7356,  4.0047,  3.9655,  3.9116,  4.1010 &
     ,  3.9102,  3.7800,  4.2964,  4.3330,  4.2622,  4.2254,  3.8195 &
     ,  3.7560,  3.6513,  3.5941,  3.5810,  3.5420,  3.5178,  3.8861 &
     ,  4.1459,  4.1147,  4.0772,  4.3120,  4.1207,  3.9900,  4.4733 &
     /)
     r0ab(4061:4130)=(/ &
        4.6157,  4.4580,  4.4194,  4.3954,  4.3739,  4.3531,  4.3343 &
     ,  4.3196,  4.2140,  4.2339,  4.1738,  4.2458,  4.2278,  4.2158 &
     ,  4.2039,  4.1658,  4.3595,  4.2857,  4.2444,  4.1855,  4.1122 &
     ,  3.7839,  3.6879,  3.5816,  3.8633,  4.1585,  4.1402,  4.1036 &
     ,  4.3694,  4.1735,  4.0368,  4.5095,  4.5538,  4.5240,  4.4252 &
     ,  3.0187,  3.1918,  3.5127,  3.6875,  3.7404,  3.6943,  3.4702 &
     ,  3.4888,  3.2914,  3.1643,  3.6669,  3.8724,  3.9940,  4.0816 &
     ,  3.8054,  3.9661,  3.7492,  3.6024,  4.0428,  4.1951,  4.1466 &
     ,  4.0515,  4.0075,  3.5020,  3.9158,  3.8546,  3.8342,  3.8008 &
     ,  3.7845,  3.7549,  3.9602,  3.8872,  3.8564,  4.0793,  3.8835 &
     /)
     r0ab(4131:4200)=(/ &
        3.7495,  4.2213,  4.3704,  4.3300,  4.2121,  3.7643,  3.7130 &
     ,  3.6144,  3.5599,  3.5474,  3.5093,  3.4853,  3.9075,  4.1115 &
     ,  4.0473,  4.0318,  4.2999,  4.1050,  3.9710,  4.4320,  4.6706 &
     ,  4.5273,  4.4581,  4.4332,  4.4064,  4.3873,  4.3684,  4.3537 &
     ,  4.2728,  4.2549,  4.2032,  4.2794,  4.2613,  4.2491,  4.2375 &
     ,  4.2322,  4.3665,  4.3061,  4.2714,  4.2155,  4.1416,  3.7660 &
     ,  3.6628,  3.5476,  3.8790,  4.1233,  4.0738,  4.0575,  4.3575 &
     ,  4.1586,  4.0183,  4.4593,  4.5927,  4.4865,  4.3813,  4.4594 &
     ,  2.9875,  3.1674,  3.4971,  3.6715,  3.7114,  3.6692,  3.4446 &
     ,  3.4676,  3.2685,  3.1405,  3.6546,  3.8579,  3.9637,  4.0581 &
     /)
     r0ab(4201:4270)=(/ &
        3.7796,  3.9463,  3.7275,  3.5792,  4.0295,  4.1824,  4.1247 &
     ,  4.0357,  3.9926,  3.4827,  3.9007,  3.8392,  3.8191,  3.7851 &
     ,  3.7687,  3.7387,  3.9290,  3.8606,  3.8306,  4.0601,  3.8625 &
     ,  3.7269,  4.2062,  4.3566,  4.3022,  4.1929,  3.7401,  3.6888 &
     ,  3.5900,  3.5350,  3.5226,  3.4838,  3.4594,  3.8888,  4.0813 &
     ,  4.0209,  4.0059,  4.2810,  4.0843,  3.9486,  4.4162,  4.6542 &
     ,  4.5005,  4.4444,  4.4196,  4.3933,  4.3741,  4.3552,  4.3406 &
     ,  4.2484,  4.2413,  4.1907,  4.2656,  4.2474,  4.2352,  4.2236 &
     ,  4.2068,  4.3410,  4.2817,  4.2479,  4.1921,  4.1182,  3.7346 &
     ,  3.6314,  3.5168,  3.8582,  4.0927,  4.0469,  4.0313,  4.3391 &
     /)
     r0ab(4271:4340)=(/ &
        4.1381,  3.9962,  4.4429,  4.5787,  4.4731,  4.3588,  4.4270 &
     ,  4.3957,  2.9659,  3.1442,  3.4795,  3.6503,  3.6814,  3.6476 &
     ,  3.4222,  3.4491,  3.2494,  3.1209,  3.6324,  3.8375,  3.9397 &
     ,  3.8311,  3.7581,  3.9274,  3.7085,  3.5598,  4.0080,  4.1641 &
     ,  4.1057,  4.0158,  3.9726,  3.4667,  3.8802,  3.8188,  3.7989 &
     ,  3.7644,  3.7474,  3.7173,  3.9049,  3.8424,  3.8095,  4.0412 &
     ,  3.8436,  3.7077,  4.1837,  4.3366,  4.2816,  4.1686,  3.7293 &
     ,  3.6709,  3.5700,  3.5153,  3.5039,  3.4684,  3.4437,  3.8663 &
     ,  4.0575,  4.0020,  3.9842,  4.2612,  4.0643,  3.9285,  4.3928 &
     ,  4.6308,  4.4799,  4.4244,  4.3996,  4.3737,  4.3547,  4.3358 &
     /)
     r0ab(4341:4410)=(/ &
        4.3212,  4.2275,  4.2216,  4.1676,  4.2465,  4.2283,  4.2161 &
     ,  4.2045,  4.1841,  4.3135,  4.2562,  4.2226,  4.1667,  4.0932 &
     ,  3.7134,  3.6109,  3.4962,  3.8352,  4.0688,  4.0281,  4.0099 &
     ,  4.3199,  4.1188,  3.9768,  4.4192,  4.5577,  4.4516,  4.3365 &
     ,  4.4058,  4.3745,  4.3539,  2.8763,  3.1294,  3.5598,  3.7465 &
     ,  3.5659,  3.5816,  3.3599,  3.4024,  3.1877,  3.0484,  3.7009 &
     ,  3.9451,  3.8465,  3.9873,  3.7079,  3.9083,  3.6756,  3.5150 &
     ,  4.0829,  4.2780,  4.1511,  4.1260,  4.0571,  3.4865,  3.9744 &
     ,  3.9150,  3.8930,  3.8578,  3.8402,  3.8073,  3.7977,  4.0036 &
     ,  3.7604,  4.0288,  3.8210,  3.6757,  4.2646,  4.4558,  4.2862 &
     /)
     r0ab(4411:4465)=(/ &
        4.2122,  3.7088,  3.6729,  3.5800,  3.5276,  3.5165,  3.4783 &
     ,  3.4539,  3.9553,  3.9818,  4.2040,  3.9604,  4.2718,  4.0689 &
     ,  3.9253,  4.4869,  4.7792,  4.4918,  4.5342,  4.5090,  4.4868 &
     ,  4.4680,  4.4486,  4.4341,  4.2023,  4.3122,  4.2710,  4.3587 &
     ,  4.3407,  4.3281,  4.3174,  4.1499,  4.3940,  4.3895,  4.3260 &
     ,  4.2725,  4.1961,  3.7361,  3.6193,  3.4916,  3.9115,  3.9914 &
     ,  3.9809,  3.9866,  4.3329,  4.1276,  3.9782,  4.5097,  4.6769 &
     ,  4.5158,  4.3291,  4.3609,  4.3462,  4.3265,  4.4341          &
     /)

      k=0
      do i=1,max_elem
         do j=1,i
            k=k+1
            r(i,j)=r0ab(k)/autoang
            r(j,i)=r0ab(k)/autoang
         enddo
      enddo

      end subroutine setr0ab


! short-range bond length correction
! modified form derived from HF-3c SRB potential
! requires TRUE ordinal numbers in array iz
! the empirical parameters are qscal (prefactor) and rscal (radii scaling)
! SG, Nov. 2016
subroutine srb_egrad2(xyz,iz,Hlat,n,energy,g,grad,rscal,qscal,echo,pbc)
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
parameter (autokcal=627.509541d0,max_elem=94,max_para=36)
DIMENSION XYZ(3,N),G(3,N),Hlat(3,3),IZ(N),ITAU_MAX(3),pp(max_elem),co(n)
REAL*8,DIMENSION(:,:),ALLOCATABLE :: r0ab
integer version,iz
logical echo,grad,pbc
character*2  esym
!AUTOANG=PAR(32)
parameter (AUTOANG =0.5291772083d0)

! Two threshold. thrR: distance cutoff thrE: numerical noise cutoff
thrR=30.0d0            ! X bohr
thrE=epsilon(1.d0)
energy=0.0d0
g=0.0d0

!Determine supercell
Itau_max=0
if(pbc) then
call criteria(thrR,Hlat,Itau_max)
endif
!get R0AB radii
allocate(R0AB(MAX_ELEM,MAX_ELEM))
!CALL CRYALLOC(R0AB,MAX_ELEM,MAX_ELEM,'SRB_EGRAD','R0AB')
call setr0ab(max_elem,AUTOANG,r0ab)
!call setr0ab(max_elem,r0ab)

!Parameter for B97-3c
!paramter of the method are rscal and qscal
!if(echo) then
!   write(IOUT,'(2x,a18)')'Parameter for SRB: '
!   write(IOUT,'(2x,a6,f10.4)')'rscal ',rscal
!   write(IOUT,'(2x,a6,f10.4)')'qscal ',qscal
!endif
ethr=0.d0

!loop over all i atoms
if (pbc) then
do iat=1,n
! all elements
!   if(iz(iat).lt.1.or.iz(iat).gt.18) cycle
   ! the SRB  due to atom jat, Loop over all j atoms in supercell
   ITAU1=Itau_max(1)
   ITAU2=Itau_max(2)
   ITAU3=Itau_max(3)
   do i=-Itau1,Itau1
      do j=-Itau2,Itau2
         do k=-Itau3,Itau3
            do jat=1,n
!               if(iz(jat).lt.1.or.iz(jat).gt.18) cycle
               iselftest=0
               ! Test for equal atoms, remove selfinteraction
               if(iat.eq.jat) then
                  iselftest=iselftest+1
                  if(i.eq.0)iselftest=iselftest+1
                  if(j.eq.0)iselftest=iselftest+1
                  if(k.eq.0)iselftest=iselftest+1
               end if
               if(iselftest.eq.4)cycle
               dx=xyz(1,iat)-xyz(1,jat)+i*Hlat(1,1)+j*Hlat(2,1)+k*Hlat(3,1)
               dy=xyz(2,iat)-xyz(2,jat)+i*Hlat(1,2)+j*Hlat(2,2)+k*Hlat(3,2)
               dz=xyz(3,iat)-xyz(3,jat)+i*Hlat(1,3)+j*Hlat(2,3)+k*Hlat(3,3)
               r=sqrt(dx*dx+dy*dy+dz*dz)
               ! distance cutoff
               if(r.gt.thrR) cycle
! Do SRB for B97-3c
               r0abij=r0ab(iz(iat),iz(jat))
               r0=rscal/r0ab(iz(iat),iz(jat))
               fi=real(iz(iat))
               fj=real(iz(jat))
               ff=-(fi*fj)**0.5d0
               ener_dum=qscal*ff*exp(-r0*r)
               !factor 1/2 from double counting
               ener_dum=ener_dum*0.5d0
               if(abs(ener_dum).lt.Ethr) cycle
               co(iat)=co(iat)+ener_dum
               energy=energy+ener_dum
               if(grad) then
                  rf=qscal/r
                  g(1,iat)=g(1,iat)-ff*r0*dx*exp(-r0*r)*rf
                  g(2,iat)=g(2,iat)-ff*r0*dy*exp(-r0*r)*rf
                  g(3,iat)=g(3,iat)-ff*r0*dz*exp(-r0*r)*rf
               endif
            enddo !jat
         enddo !k
      enddo !j
   enddo !i
enddo !iat
else
do iat=1,n
! all elements
!   if(iz(iat).lt.1.or.iz(iat).gt.18) cycle
   ! the SRB  due to atom jat, Loop over all j atoms in supercell
   do jat=1,n
      !               if(iz(jat).lt.1.or.iz(jat).gt.18) cycle
      ! Test for equal atoms, remove selfinteraction
      if(iat.eq.jat) then
         cycle
      end if
      dx=xyz(1,iat)-xyz(1,jat)
      dy=xyz(2,iat)-xyz(2,jat)
      dz=xyz(3,iat)-xyz(3,jat)
      r=sqrt(dx*dx+dy*dy+dz*dz)
      ! distance cutoff
      if(r.gt.thrR) cycle
      ! Do SRB for B97-3c
      r0abij=r0ab(iz(iat),iz(jat))
      r0=rscal/r0ab(iz(iat),iz(jat))
      fi=real(iz(iat))
      fj=real(iz(jat))
      ff=-(fi*fj)**0.5d0
      ener_dum=qscal*ff*exp(-r0*r)
      !factor 1/2 from double counting
      ener_dum=ener_dum*0.5d0
      if(abs(ener_dum).lt.Ethr) cycle
      co(iat)=co(iat)+ener_dum
      energy=energy+ener_dum
      if(grad) then
         rf=qscal/r
         g(1,iat)=g(1,iat)-ff*r0*dx*exp(-r0*r)*rf
         g(2,iat)=g(2,iat)-ff*r0*dy*exp(-r0*r)*rf
         g(3,iat)=g(3,iat)-ff*r0*dz*exp(-r0*r)*rf
      endif
   enddo !jat
enddo !iat
end if
if(echo)then
   write(IOUT,'(/2x,a5,2x,a5,4x,a15)') &
        '#','ON','SRB [kcal/mol]'
   do i=1,n
      write(IOUT,'(2x,2(i5,2x),F9.3)') i,iz(i),co(i)*AUTOKCAL
   enddo
!   write(IOUT,*)'** SRB correction **'
!   write(IOUT,'(2x,a7,F18.10,'' / (a.u.) || '',x,F11.4,'' / (kcal/mol)'')')'Esrb:  ',energy,energy*AUTOKCAL
!   if(grad)write(IOUT,*)'|G|=',sum(abs(g(1:3,1:n)))
endif

end subroutine srb_egrad2
end module gcp
