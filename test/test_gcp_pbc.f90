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

module test_gcp_pbc
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use gcp
   implicit none
   private

   public :: collect_gcp_pbc

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_gcp_pbc(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("anthracene", test_anthracene), &
      & new_unittest("ethcar", test_ethcar), &
      & new_unittest("formamide", test_formamide), &
      & new_unittest("hexdio", test_hexdio), &
      & new_unittest("naph", test_naph), &
      & new_unittest("oxacb", test_oxacb), &
      & new_unittest("trioxane", test_trioxane) &
      & ]

end subroutine collect_gcp_pbc


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


subroutine test_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "anthracene")
   call test_generic(error, mol, "hf3c", -0.59161990455181912_wp)

end subroutine test_anthracene


subroutine test_ethcar(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ethcar")
   call test_generic(error, mol, "hf3c", -0.11814220512782711_wp)

end subroutine test_ethcar


subroutine test_formamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "formamide")
   call test_generic(error, mol, "hf3c", -0.15500632261441180_wp)

end subroutine test_formamide


subroutine test_hexdio(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "hexdio")
   call test_generic(error, mol, "hf3c", -4.4642369521389225E-002_wp)

end subroutine test_hexdio


subroutine test_naph(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "naph")
   call test_generic(error, mol, "hf3c", -1.0399600000679481_wp)

end subroutine test_naph


subroutine test_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, "hf3c", -0.15716399742412498_wp)

end subroutine test_oxacb


subroutine test_trioxane(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "trioxane")
   call test_generic(error, mol, "b973c", -0.36904665894177763_wp)

end subroutine test_trioxane


end module test_gcp_pbc
