! This file is part of mctc-gcp.
! SPDX-Identifier: LGPL-3.0-or-later
!
! mctc-gcp is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! mctc-gcp is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with mctc-gcp.  If not, see <https://www.gnu.org/licenses/>.

!> Geometrical counterpoise correction.
!>
!> The correction itself is implemented in s-dftd3, this module only provides
!> the legacy interface of the standalone gCP program on top of it.
module gcp
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   use dftd3_cutoff, only : realspace_cutoff
   use dftd3_gcp, only : gcp_param, get_gcp_param, get_geometric_counterpoise
   use dftd3_output, only : ascii_gcp_param, turbomole_gradient
   use gcp_version, only : get_gcp_version
   use mctc_env, only : wp
   use mctc_io, only : structure_type, new
   use mctc_io_math, only : matinv_3x3
   implicit none
   private

   public :: gcp_call, wregrad_tm, get_gcp_version

   !> Step size of the numerical Hessian in Bohr
   real(wp), parameter :: hess_step = 0.005_wp

contains

!> Evaluate the geometrical counterpoise correction, interface for
!> orca/turbomole/crystal
subroutine gcp_call(n, xyz, lat, iz, gcp_e, gcp_g, gcp_glat, dograd, dohess, &
      & pbc, method, echo, parfile)

   !> Number of atoms
   integer, intent(in) :: n

   !> Cartesian coordinates in Bohr
   real(wp), intent(in) :: xyz(3, n)

   !> Lattice vectors in Bohr, stored row-wise
   real(wp), intent(in) :: lat(3, 3)

   !> Atomic numbers
   integer, intent(in) :: iz(n)

   !> Counterpoise energy
   real(wp), intent(out) :: gcp_e

   !> Counterpoise gradient
   real(wp), intent(out) :: gcp_g(3, n)

   !> Counterpoise lattice gradient
   real(wp), intent(out) :: gcp_glat(3, 3)

   !> Evaluate gradient and lattice gradient
   logical, intent(in) :: dograd

   !> Evaluate the numerical Hessian
   logical, intent(in) :: dohess

   !> Periodic system
   logical, intent(in) :: pbc

   !> Level of theory, either <method>/<basis> or a composite method
   character(len=*), intent(in) :: method

   !> Print information
   logical, intent(in) :: echo

   !> Write a parameter file, no longer supported
   logical, intent(in) :: parfile

   type(structure_type) :: mol
   type(gcp_param) :: param
   real(wp) :: sigma(3, 3)
   integer :: unit

   call new_gcp_structure(mol, n, xyz, lat, iz, pbc)
   call new_gcp_param(param, mol, method)

   if (parfile) write(error_unit, '(a)') &
      & "[Warn] Parameter files are not supported anymore"
   if (echo) call ascii_gcp_param(output_unit, mol, param, method)

   gcp_e = 0.0_wp
   gcp_g(:, :) = 0.0_wp
   gcp_glat(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (dograd) then
      call get_geometric_counterpoise(mol, param, realspace_cutoff(), gcp_e, &
         & gcp_g, sigma)
      if (pbc) then
         gcp_glat(:, :) = transpose(matmul(sigma, matinv_3x3(lat)))
         open(file="gcp_cellgradient", newunit=unit)
         write(unit, *) gcp_glat
         close(unit)
      end if
   else
      call get_geometric_counterpoise(mol, param, realspace_cutoff(), gcp_e)
   end if

   if (dohess) call numerical_hessian(mol, param, echo)

   ! for scripting
   open(file=".CPC", newunit=unit)
   write(unit, *) gcp_e
   close(unit)
   open(file=".CP", newunit=unit)
   write(unit, '(f22.16)') gcp_e
   close(unit)

end subroutine gcp_call


!> Add the counterpoise gradient to a Turbomole gradient file
subroutine wregrad_tm(maxat, nat, xyz, iat, ifrez, edisp, gin, echo)

   !> Maximum number of atoms, unused
   integer, intent(in) :: maxat

   !> Number of atoms
   integer, intent(in) :: nat

   !> Cartesian coordinates in Bohr
   real(wp), intent(in) :: xyz(3, nat)

   !> Atomic numbers
   integer, intent(in) :: iat(nat)

   !> Frozen atoms, unused
   integer, intent(in) :: ifrez(nat)

   !> Counterpoise energy
   real(wp), intent(in) :: edisp

   !> Counterpoise gradient
   real(wp), intent(in) :: gin(3, nat)

   !> Print information
   logical, intent(in) :: echo

   type(structure_type) :: mol
   integer :: stat

   call new(mol, iat, xyz)
   call turbomole_gradient(mol, "gradient", edisp, gin, stat)
   if (echo) then
      if (stat == 0) then
         write(output_unit, '(a)') "[Info] gCP gradient added to Turbomole gradient file"
      else
         write(output_unit, '(a)') "[Warn] Could not add to Turbomole gradient file"
      end if
   end if

end subroutine wregrad_tm


!> Create the molecular structure data from the legacy plain arrays
subroutine new_gcp_structure(mol, n, xyz, lat, iz, pbc)
   type(structure_type), intent(out) :: mol
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(in) :: lat(3, 3)
   integer, intent(in) :: iz(n)
   logical, intent(in) :: pbc

   if (pbc) then
      call new(mol, iz, xyz, lattice=transpose(lat), periodic=spread(.true., 1, 3))
   else
      call new(mol, iz, xyz)
   end if
end subroutine new_gcp_structure


!> Obtain the counterpoise parameters for a level of theory string
subroutine new_gcp_param(param, mol, input)
   type(gcp_param), intent(out) :: param
   type(structure_type), intent(in) :: mol
   character(len=*), intent(in) :: input

   character(len=:), allocatable :: method, basis

   call split_level(input, method, basis)
   if (allocated(basis)) then
      call get_gcp_param(param, mol, method, basis)
   else
      call get_gcp_param(param, mol, method)
   end if

   if (.not.allocated(param%emiss) .and. .not.param%srb) then
      write(error_unit, '(a)') "[Fatal] No gCP parameters available for '"//trim(input)//"'"
      error stop
   end if
end subroutine new_gcp_param


!> Split a level of theory string into method and basis set, normalizing
!> case and hyphenation, and resolve the aliases of the legacy implementation
subroutine split_level(input, method, basis)
   character(len=*), intent(in) :: input
   character(len=:), allocatable, intent(out) :: method, basis

   character(len=len(input)) :: level
   integer :: ic, is, nc

   nc = 0
   do ic = 1, len_trim(input)
      if (input(ic:ic) == "-") cycle
      nc = nc + 1
      is = index("ABCDEFGHIJKLMNOPQRSTUVWXYZ", input(ic:ic))
      if (is > 0) then
         level(nc:nc) = achar(is + iachar("a") - 1)
      else
         level(nc:nc) = input(ic:ic)
      end if
   end do

   is = index(level(:nc), "/")
   if (is > 0) then
      method = level(:is-1)
      basis = level(is+1:nc)
      select case(basis)
      case("def2svp")
         basis = "svp"
      case("tzvp")
         basis = "deftzvp"
      end select
      ! composite methods imply their own basis set
      if (index(method, "3c") > 0) deallocate(basis)
   else
      method = level(:nc)
      select case(method)
      case("mtzvpp", "def2mtzvpp")
         method = "r2scan3c"
      case("b3pbe3c")
         basis = "def2mtzvp"
      end select
   end if

   ! generic (hybrid) functionals are parametrized as B3LYP in s-dftd3
   select case(method)
   case("dft", "pbe")
      method = "b3lyp"
   end select
end subroutine split_level


!> Numerical Hessian of the counterpoise correction, written to gcp_hessian
subroutine numerical_hessian(mol, param, echo)
   type(structure_type), intent(inout) :: mol
   type(gcp_param), intent(in) :: param
   logical, intent(in) :: echo

   real(wp) :: energy, sigma(3, 3)
   real(wp), allocatable :: hessian(:, :), gr(:, :), gl(:, :)
   integer :: iat, ic, ii, unit, ndim

   ndim = 3 * mol%nat
   allocate(hessian(ndim, ndim), source=0.0_wp)
   allocate(gr(3, mol%nat), gl(3, mol%nat))

   if (echo) write(output_unit, '(a)') "[Info] Doing Hessian numerically ..."
   do iat = 1, mol%nat
      do ic = 1, 3
         ii = 3*(iat - 1) + ic
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + hess_step
         energy = 0.0_wp
         gr(:, :) = 0.0_wp
         sigma(:, :) = 0.0_wp
         call get_geometric_counterpoise(mol, param, realspace_cutoff(), energy, gr, sigma)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*hess_step
         energy = 0.0_wp
         gl(:, :) = 0.0_wp
         sigma(:, :) = 0.0_wp
         call get_geometric_counterpoise(mol, param, realspace_cutoff(), energy, gl, sigma)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + hess_step
         hessian(:, ii) = reshape(gr - gl, [ndim]) / (2*hess_step)
      end do
   end do
   hessian(:, :) = 0.5_wp * (hessian + transpose(hessian))

   open(file="gcp_hessian", newunit=unit)
   write(unit, '(a)') "$hessian"
   do ii = 1, ndim
      write(unit, '(5x,5f15.10)') hessian(ii, :)
   end do
   write(unit, '(a)') "$end"
   close(unit)
   if (echo) write(output_unit, '(a)') "[Info] gCP Hessian written to 'gcp_hessian'"

end subroutine numerical_hessian

end module gcp
