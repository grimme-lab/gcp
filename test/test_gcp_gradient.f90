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

module test_gcp_gradient
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mctc_io_math, only : matinv_3x3
   use dftd3_gcp, only : gcp_param
   use mstore, only : get_structure
   use gcp
   implicit none
   private

   public :: collect_gcp_gradient

   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_gcp_gradient(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("HF/DZ", test_hf_dz), &
      & new_unittest("DFT/SV", test_dft_sv), &
      & new_unittest("HF-3c", test_hf3c), &
      & new_unittest("B97-3c", test_b973c), &
      & new_unittest("HSE-3c", test_hse3c), &
      & new_unittest("Custom parameter file", test_paramfile), &
      & new_unittest("HF-3c (lattice)", test_hf3c_latt) &
      & ]

end subroutine collect_gcp_gradient


subroutine test_paramfile(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   integer :: unit
   character(len=256) :: home

   call get_environment_variable('HOME', home)
   if (len_trim(home) == 0) home = '/tmp'
   home = trim(home)//'/gcp-paramfile-test'
   call execute_command_line('rm -rf '//trim(home))
   call execute_command_line('mkdir -p '//trim(home))
   open(newunit=unit, file=trim(home)//'/legacy.gcppar', status='replace')
   write(unit, '(a)') 'sv 0.2 0.4 0.8 0.9'
   close(unit)

   call get_structure(mol, 'X23', 'CO2')
   call new_gcp_param(param, mol, trim(home)//'/legacy.gcppar')

   if (abs(param%sigma - 0.2_wp) > 1.0e-12_wp .or. &
      & abs(param%alpha - 0.8_wp) > 1.0e-12_wp .or. &
      & abs(param%beta - 0.9_wp) > 1.0e-12_wp) then
      call test_failed(error, 'Legacy parameter file values were not applied')
      call execute_command_line('rm -rf '//trim(home))
      return
   end if

   call execute_command_line('rm -rf '//trim(home))
end subroutine test_paramfile


subroutine test_numgrad(error, mol, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Method name
   character(len=*), intent(in) :: method

   integer :: iat, ic, mat
   real(wp) :: energy, er, el
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp) :: gradlatt(3, 3), lattice(3, 3)
   logical :: pbc
   character(len=20) :: method_str
   logical, parameter :: dohess = .false.
   logical, parameter :: echo = .false.
   logical, parameter :: parfile = .false.
   real(wp), parameter :: step = 1.0e-6_wp

   method_str = method
   pbc = any(mol%periodic)
   lattice = 0.0_wp
   if (allocated(mol%lattice)) lattice = transpose(mol%lattice)

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))

   if (pbc) then
      mat = min(mol%nat, 3)
   else
      mat = mol%nat
   end if
   do iat = 1, mat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call gcp_call(mol%nat, mol%xyz, lattice, mol%num(mol%id), &
            & er, gradient, gradlatt, .false., dohess, pbc, method_str, &
            & echo, parfile)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call gcp_call(mol%nat, mol%xyz, lattice, mol%num(mol%id), &
            & el, gradient, gradlatt, .false., dohess, pbc, method_str, &
            & echo, parfile)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call gcp_call(mol%nat, mol%xyz, lattice, mol%num(mol%id), &
      & energy, gradient, gradlatt, .true., dohess, pbc, method_str, echo, parfile)

   if (any(abs(gradient(:, :mat)-numgrad(:, :mat)) > thr)) then
      call test_failed(error, "Numerical and analytical gradient do not match")
      do iat = 1, mat
         print'(3es14.5,3x,a)', gradient(:, iat)-numgrad(:, iat), mol%sym(mol%id(iat))
      end do
   end if

end subroutine test_numgrad


subroutine test_hf_dz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, "hf/dz")

end subroutine test_hf_dz


subroutine test_dft_sv(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numgrad(error, mol, "dft/sv")

end subroutine test_dft_sv


subroutine test_hse3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_numgrad(error, mol, "hse3c")

end subroutine test_hse3c


subroutine test_hf3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ammonia")
   call test_numgrad(error, mol, "hf3c")

end subroutine test_hf3c


subroutine test_b973c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, "b973c")

end subroutine test_b973c


!> Compare the lattice gradient against numerical differentiation of the
!> energy with respect to the cell parameters at fixed fractional coordinates
subroutine test_numlatt(error, mol, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Method name
   character(len=*), intent(in) :: method

   integer :: ic, jc
   real(wp) :: energy, er, el
   real(wp) :: gradlatt(3, 3), numlatt(3, 3), lattice(3, 3), trial(3, 3)
   real(wp), allocatable :: gradient(:, :), frac(:, :), xyz(:, :)
   character(len=20) :: method_str
   logical, parameter :: dohess = .false.
   logical, parameter :: echo = .false.
   logical, parameter :: parfile = .false.
   real(wp), parameter :: step = 1.0e-5_wp, thr2 = 1.0e-6_wp

   method_str = method
   lattice = transpose(mol%lattice)
   frac = matmul(matinv_3x3(mol%lattice), mol%xyz)
   allocate(gradient(3, mol%nat), xyz(3, mol%nat))

   do ic = 1, 3
      do jc = 1, 3
         trial(:, :) = lattice
         trial(ic, jc) = lattice(ic, jc) + step
         xyz(:, :) = matmul(transpose(trial), frac)
         call gcp_call(mol%nat, xyz, trial, mol%num(mol%id), &
            & er, gradient, gradlatt, .false., dohess, .true., method_str, echo, parfile)
         trial(ic, jc) = lattice(ic, jc) - step
         xyz(:, :) = matmul(transpose(trial), frac)
         call gcp_call(mol%nat, xyz, trial, mol%num(mol%id), &
            & el, gradient, gradlatt, .false., dohess, .true., method_str, echo, parfile)
         numlatt(ic, jc) = 0.5_wp*(er - el)/step
      end do
   end do

   call gcp_call(mol%nat, mol%xyz, lattice, mol%num(mol%id), &
      & energy, gradient, gradlatt, .true., dohess, .true., method_str, echo, parfile)

   if (any(abs(gradlatt - numlatt) > thr2)) then
      call test_failed(error, "Numerical and analytical lattice gradient do not match")
      do ic = 1, 3
         print'(3es14.5)', gradlatt(ic, :) - numlatt(ic, :)
      end do
   end if

end subroutine test_numlatt


subroutine test_hf3c_latt(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "formamide")
   call test_numlatt(error, mol, "hf3c")

end subroutine test_hf3c_latt


end module test_gcp_gradient
