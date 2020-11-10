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

module test_gcp_dft_methods
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use gcp
   implicit none
   private

   public :: collect_gcp_dft_methods

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_gcp_dft_methods(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("DFT/MINIS", test_dft_minis), &
      & new_unittest("DFT/MINIX", test_dft_minix), &
      & new_unittest("DFT/SV", test_dft_sv), &
      & new_unittest("DFT/def2-SV(P)", test_dft_def2sv_p), &
      & new_unittest("DFT/def2-SVP", test_dft_def2svp), &
      & new_unittest("DFT/DZP", test_dft_dzp), &
      & new_unittest("DFT/def-TZVP", test_dft_deftzvp), &
      & new_unittest("DFT/def2-TZVP", test_dft_def2tzvp), &
      & new_unittest("DFT/cc-pVDZ", test_dft_ccpvdz), &
      & new_unittest("DFT/aug-cc-pVDZ", test_dft_augccpvdz) &
      & ]

end subroutine collect_gcp_dft_methods


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


subroutine test_dft_minis(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "13")
   call test_generic(error, mol, "dft/minis", 2.8831779044393737E-2_wp)

end subroutine test_dft_minis


subroutine test_dft_minix(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "14")
   call test_generic(error, mol, "dft/minix", 1.8191237759356119E-2_wp)

end subroutine test_dft_minix


subroutine test_dft_sv(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_generic(error, mol, "dft/sv", 8.9178105933847233E-3_wp)

end subroutine test_dft_sv


subroutine test_dft_def2sv_p(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_generic(error, mol, "dft/def2-sv(p)", 9.0695660381885229E-3_wp)

end subroutine test_dft_def2sv_p


subroutine test_dft_def2svp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_generic(error, mol, "dft/def2-svp", 2.4834643199826180E-2_wp)

end subroutine test_dft_def2svp


subroutine test_dft_dzp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "18")
   call test_generic(error, mol, "dft/dzp", 2.2195075119243783E-2_wp)

end subroutine test_dft_dzp


subroutine test_dft_deftzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "19")
   call test_generic(error, mol, "dft/def-tzvp", 3.6258064909241677E-3_wp)

end subroutine test_dft_deftzvp


subroutine test_dft_def2tzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "20")
   call test_generic(error, mol, "dft/def2-tzvp", 3.9044561047188592E-3_wp)

end subroutine test_dft_def2tzvp


subroutine test_dft_ccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "21")
   call test_generic(error, mol, "dft/cc-pvdz", 9.9649294926933245E-3_wp)

end subroutine test_dft_ccpvdz


subroutine test_dft_augccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "22")
   call test_generic(error, mol, "dft/aug-cc-pvdz", 2.9849348931724447E-3_wp)

end subroutine test_dft_augccpvdz


end module test_gcp_dft_methods
