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

module test_gcp_3c_methods
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use gcp
   implicit none
   private

   public :: collect_gcp_3c_methods

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_gcp_3c_methods(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("HF-3c", test_hf3c), &
      & new_unittest("HF-3c (PBC)", test_hf3c_pbc), &
      & new_unittest("PBEh-3c", test_pbeh3c), &
      & new_unittest("PBEh-3c (PBC)", test_pbeh3c_pbc), &
      & new_unittest("HSE-3c", test_hse3c), &
      & new_unittest("HSE-3c (PBC)", test_hse3c_pbc), &
      & new_unittest("B97-3c", test_b973c), &
      & new_unittest("B97-3c (PBC)", test_b973c_pbc), &
      & new_unittest("r2SCAN-3c", test_r2scan3c), &
      & new_unittest("r2SCAN-3c (PBC)", test_r2scan3c_pbc) &
      & ]

end subroutine collect_gcp_3c_methods


subroutine test_generic(error, mol, method, energy_ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Method name
   character(len=*), intent(in) :: method

   !> Reference energy
   real(wp), intent(in) :: energy_ref

   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :)
   real(wp) :: gradlatt(3, 3)
   logical :: pbc
   character(len=20) :: method_str
   logical, parameter :: dograd = .false.
   logical, parameter :: dohess = .false.
   logical, parameter :: echo = .false.
   logical, parameter :: parfile = .false.

   method_str = method
   pbc = any(mol%periodic)

   allocate(gradient(3, mol%nat))

   call gcp_call(mol%nat, mol%xyz, mol%lattice, mol%num(mol%id), &
      & energy, gradient, gradlatt, dograd, dohess, pbc, method_str, echo, parfile)

   call check(error, energy_ref, energy, thr=thr)
   if (allocated(error)) then
      print*, energy
   end if

end subroutine test_generic


subroutine test_hf3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, "hf3c", -3.0475695153528881E-2_wp)

end subroutine test_hf3c


subroutine test_pbeh3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, "pbeh3c", 2.0602039298887861E-2_wp)

end subroutine test_pbeh3c


subroutine test_hse3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_generic(error, mol, "hse3c", 1.9885888051458127E-2_wp)

end subroutine test_hse3c


subroutine test_r2scan3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_generic(error, mol, "r2scan3c", 1.1501835545305577E-2_wp)

end subroutine test_r2scan3c


subroutine test_b973c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_generic(error, mol, "b973c", -3.6973735977079793E-2_wp)

end subroutine test_b973c


subroutine test_hf3c_pbc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "acetic")
   call test_generic(error, mol, "hf3c", -0.21073394289693193_wp)

end subroutine test_hf3c_pbc


subroutine test_pbeh3c_pbc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "adaman")
   call test_generic(error, mol, "pbeh3c", 9.6863188972415457E-2_wp)

end subroutine test_pbeh3c_pbc


subroutine test_hse3c_pbc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ammonia")
   call test_generic(error, mol, "hse3c", 2.0735966360832053E-2_wp)

end subroutine test_hse3c_pbc


subroutine test_r2scan3c_pbc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "benzene")
   call test_generic(error, mol, "mtzvpp", 1.1008015279278560E-2_wp)

end subroutine test_r2scan3c_pbc


subroutine test_b973c_pbc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_generic(error, mol, "b973c", -4.7346715088819775E-2_wp)

end subroutine test_b973c_pbc


end module test_gcp_3c_methods
