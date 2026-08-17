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

   public :: gcp_call, wregrad_tm, get_gcp_version, new_gcp_param

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
   logical :: found

   call read_legacy_param(param, mol, input, found)
   if (found) return

   call split_level(input, method, basis)
   if (allocated(basis)) then
      call get_gcp_param(param, mol, method, basis)
   else
      call get_gcp_param(param, mol, method)
   end if
   call fix_msvp_virtuals(param, method, basis)

   if (.not.allocated(param%emiss) .and. .not.param%srb) then
      write(error_unit, '(a)') "[Fatal] No gCP parameters available for '"//trim(input)//"'"
      error stop
   end if
end subroutine new_gcp_param


subroutine read_legacy_param(param, mol, input, found)
   type(gcp_param), intent(out) :: param
   type(structure_type), intent(in) :: mol
   character(len=*), intent(in) :: input
   logical, intent(out) :: found

   character(len=:), allocatable :: file
   character(len=256) :: home, home2
   logical :: exists

   found = .false.

   if (len_trim(input) > 0) then
      inquire(file=trim(input), exist=exists)
      if (exists) then
         call read_legacy_param_file(param, mol, trim(input), trim(input), found)
         if (found) return
      end if
   end if

   home = ''
   call get_environment_variable("HOME", home)
   if (len_trim(home) == 0) then
      call get_environment_variable("USERPROFILE", home)
   end if
   if (len_trim(home) == 0) then
      call get_environment_variable("HOMEDRIVE", home)
      if (len_trim(home) > 0) then
         home2 = ''
         call get_environment_variable("HOMEPATH", home2)
         if (len_trim(home2) > 0) then
            home = trim(home)//trim(home2)
         end if
      end if
   end if
   if (len_trim(home) > 0) then
      file = trim(home)//"/.gcppar"
      call read_legacy_param_file(param, mol, file, trim(input), found)
      if (found) return
   end if

   inquire(file=".gcppar", exist=exists)
   if (exists) then
      call read_legacy_param_file(param, mol, ".gcppar", trim(input), found)
   end if
end subroutine read_legacy_param


subroutine read_legacy_param_file(param, mol, file, target, found)
   type(gcp_param), intent(out) :: param
   type(structure_type), intent(in) :: mol
   character(len=*), intent(in) :: file, target
   logical, intent(out) :: found

   integer :: unit, stat
   logical :: exists, use_first, exact_match
   character(len=256) :: line
   character(len=32) :: key, first_key
   character(len=:), allocatable :: basis, method
   real(wp) :: sigma, eta, alpha, beta, first_sigma, first_eta, first_alpha, first_beta
   real(wp), allocatable :: values(:)

   found = .false.
   exact_match = .false.
   use_first = .false.
   key = ''
   first_key = ''
   sigma = 0.0_wp
   eta = 0.0_wp
   alpha = 0.0_wp
   beta = 0.0_wp
   first_sigma = 0.0_wp
   first_eta = 0.0_wp
   first_alpha = 0.0_wp
   first_beta = 0.0_wp

   inquire(file=trim(file), exist=exists)
   if (.not.exists) return

   open(file=trim(file), newunit=unit, status='old', action='read')
   do
      read(unit, '(a)', iostat=stat) line
      if (stat /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle
      read(line, *, iostat=stat) key, sigma, eta, alpha, beta
      if (stat /= 0) cycle
      if (.not.use_first .and. (trim(target) == 'file' .or. trim(target) == trim(file))) then
         first_key = key
         first_sigma = sigma
         first_eta = eta
         first_alpha = alpha
         first_beta = beta
         use_first = .true.
      end if
      if (trim(key) == trim(target)) then
         exact_match = .true.
         exit
      end if
   end do
   close(unit)
   if (.not.exact_match .and. .not.use_first) return

   if (exact_match) then
      key = key
   else
      key = first_key
      sigma = first_sigma
      eta = first_eta
      alpha = first_alpha
      beta = first_beta
   end if

   call split_level(trim(key), method, basis)
   if (.not.allocated(basis) .or. len_trim(basis) == 0) then
      basis = trim(key)
   end if
   if (len_trim(basis) == 0) basis = 'sv'
   call get_gcp_param(param, mol, method='hf', basis=basis)
   param%sigma = sigma
   param%alpha = alpha
   param%beta = beta
   if (allocated(param%slater) .and. eta > 0.0_wp) then
      values = legacy_slater_exp(param%zeff)
      param%slater = eta * values
   end if

   found = .true.
end subroutine read_legacy_param_file


pure function legacy_slater_exp(zeff) result(values)
   integer, intent(in) :: zeff(:)
   real(wp) :: values(size(zeff))

   real(wp), parameter :: slater_s(*) = [ &
      & 1.2000_wp, 1.6469_wp, 0.6534_wp, 1.0365_wp, 1.3990_wp, 1.7210_wp, 1.9398_wp, 2.2399_wp, 2.5644_wp, 2.8812_wp, &
      & 0.8675_wp, 1.1935_wp, 1.5143_wp, 1.7580_wp, 1.9860_wp, 2.1362_wp, 2.3617_wp, 2.5796_wp, 0.9362_wp, 1.2112_wp, &
      & 1.2870_wp, 1.3416_wp, 1.3570_wp, 1.3804_wp, 1.4761_wp, 1.5465_wp, 1.5650_wp, 1.5532_wp, 1.5781_wp, 1.7778_wp, &
      & 2.0675_wp, 2.2702_wp, 2.4546_wp, 2.5680_wp, 2.7523_wp, 2.9299_wp]
   real(wp), parameter :: slater_p(*) = [ &
      & 0.0000_wp, 0.0000_wp, 0.5305_wp, 0.8994_wp, 1.2685_wp, 1.6105_wp, 1.9398_wp, 2.0477_wp, 2.4022_wp, 2.7421_wp, &
      & 0.6148_wp, 0.8809_wp, 1.1660_wp, 1.4337_wp, 1.6755_wp, 1.7721_wp, 2.0176_wp, 2.2501_wp, 0.6914_wp, 0.9329_wp, &
      & 0.9828_wp, 1.0104_wp, 0.9947_wp, 0.9784_wp, 1.0641_wp, 1.1114_wp, 1.1001_wp, 1.0594_wp, 1.0527_wp, 1.2448_wp, &
      & 1.5073_wp, 1.7680_wp, 1.9819_wp, 2.0548_wp, 2.2652_wp, 2.4617_wp]
   real(wp), parameter :: slater_d(*) = [ &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 2.4341_wp, 2.6439_wp, 2.7809_wp, 2.9775_wp, 3.2208_wp, 3.4537_wp, 3.6023_wp, 3.7017_wp, 3.8962_wp, 2.0477_wp, &
      & 2.4022_wp, 2.7421_wp, 0.6148_wp, 0.8809_wp, 1.1660_wp, 1.4337_wp]
   integer :: ip

   values(:) = 0.0_wp
   do ip = 1, size(zeff)
      if (zeff(ip) < 1 .or. zeff(ip) > 36) cycle
      if (zeff(ip) <= 2) then
         values(ip) = slater_s(zeff(ip))
      else if (zeff(ip) <= 20) then
         values(ip) = 0.5_wp * (slater_s(zeff(ip)) + slater_p(zeff(ip)))
      else if (zeff(ip) <= 30) then
         values(ip) = (slater_s(zeff(ip)) + 2.0_wp * slater_d(zeff(ip))) / 3.0_wp
      else
         values(ip) = slater_s(zeff(ip))
      end if
   end do
end function legacy_slater_exp


! s-dftd3 1.5.0 still counts 10 instead of 9 basis functions for Li and
! Be in def2-mSVP, drop this once the fix is released upstream
subroutine fix_msvp_virtuals(param, method, basis)
   type(gcp_param), intent(inout) :: param
   character(len=*), intent(in) :: method
   character(len=*), intent(in), optional :: basis

   logical :: msvp

   msvp = method == "pbeh3c" .or. method == "hse3c"
   if (present(basis)) msvp = msvp .or. basis == "msvp" .or. basis == "def2msvp"
   if (.not.msvp .or. .not.allocated(param%xv)) return

   where(param%zeff == 3 .or. param%zeff == 4) param%xv = param%xv - 1.0_wp
end subroutine fix_msvp_virtuals


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
