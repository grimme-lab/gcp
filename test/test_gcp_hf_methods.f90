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

module test_gcp_hf_methods
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use gcp
   implicit none
   private

   public :: collect_gcp_hf_methods

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_gcp_hf_methods(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("HF/MINIS", test_hf_minis), &
      & new_unittest("HF/MINIX", test_hf_minix), &
      & new_unittest("HF/SV", test_hf_sv), &
      & new_unittest("HF/def2-SV(P)", test_hf_def2sv_p), &
      & new_unittest("HF/def2-SVP", test_hf_def2svp), &
      & new_unittest("HF/DZP", test_hf_dzp), &
      & new_unittest("HF/def-TZVP", test_hf_deftzvp), &
      & new_unittest("HF/def2-TZVP", test_hf_def2tzvp), &
      & new_unittest("HF/cc-pVDZ", test_hf_ccpvdz), &
      & new_unittest("HF/aug-cc-pVDZ", test_hf_augccpvdz) &
      & ]

end subroutine collect_gcp_hf_methods


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


subroutine test_hf_minis(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_generic(error, mol, "hf/minis", 1.5336656998809515E-2_wp)

end subroutine test_hf_minis


subroutine test_hf_minix(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_generic(error, mol, "hf/minix", 1.8559832674618637E-2_wp)

end subroutine test_hf_minix


subroutine test_hf_sv(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_generic(error, mol, "hf/sv", 6.4701104644333107E-3_wp)

end subroutine test_hf_sv


subroutine test_hf_def2sv_p(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_generic(error, mol, "hf/def2-sv(p)", 1.5009548340274959E-2_wp)

end subroutine test_hf_def2sv_p


subroutine test_hf_def2svp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "08")
   call test_generic(error, mol, "hf/def2-svp", 1.8249078757662872E-2_wp)

end subroutine test_hf_def2svp


subroutine test_hf_dzp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "08")
   call test_generic(error, mol, "hf/dzp", 2.4002048454141937E-2_wp)

end subroutine test_hf_dzp


subroutine test_hf_deftzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "09")
   call test_generic(error, mol, "hf/def-tzvp", 4.9799363418100567E-3_wp)

end subroutine test_hf_deftzvp


subroutine test_hf_def2tzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "10")
   call test_generic(error, mol, "hf/def2-tzvp", 2.1614549442281130E-3_wp)

end subroutine test_hf_def2tzvp


subroutine test_hf_ccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "11")
   call test_generic(error, mol, "hf/cc-pvdz", 9.5654476866715681E-3_wp)

end subroutine test_hf_ccpvdz


subroutine test_hf_augccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "12")
   call test_generic(error, mol, "hf/aug-cc-pvdz", 2.6289439007561133E-3_wp)

end subroutine test_hf_augccpvdz


end module test_gcp_hf_methods
