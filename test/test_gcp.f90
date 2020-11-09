! This file is part of mctc-gcp.
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

module test_gcp
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io_structure, only : structure_type
   use testsuite_structure, only : get_structure
   use gcp
   implicit none
   private

   public :: collect_gcp

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_gcp(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("HF-3c", test_hf3c), &
      & new_unittest("PBEh-3c", test_pbeh3c), &
      & new_unittest("HSE-3c", test_hse3c), &
      & new_unittest("HF/MINIS", test_hf_minis), &
      & new_unittest("HF/MINIX", test_hf_minix), &
      & new_unittest("HF/SV", test_hf_sv), &
      & new_unittest("HF/def2-SV(P)", test_hf_def2sv_p), &
      & new_unittest("HF/def2-SVP", test_hf_def2svp), &
      & new_unittest("HF/DZP", test_hf_dzp), &
      & new_unittest("HF/def-TZVP", test_hf_deftzvp), &
      & new_unittest("HF/def2-TZVP", test_hf_def2tzvp), &
      & new_unittest("HF/cc-pVDZ", test_hf_ccpvdz), &
      & new_unittest("HF/aug-cc-pVDZ", test_hf_augccpvdz), &
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

end subroutine collect_gcp


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

   call get_structure(mol, "mindless01")
   call test_generic(error, mol, "hf3c", -0.10185940933506053_wp)

end subroutine test_hf3c


subroutine test_pbeh3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless02")
   call test_generic(error, mol, "pbeh3c", 2.0602039298887861E-2_wp)

end subroutine test_pbeh3c


subroutine test_hse3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless03")
   call test_generic(error, mol, "hse3c", 1.9885888051458127E-2_wp)

end subroutine test_hse3c


subroutine test_hf_minis(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless04")
   call test_generic(error, mol, "hf/minis", 1.5336656998809515E-2_wp)

end subroutine test_hf_minis


subroutine test_hf_minix(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless05")
   call test_generic(error, mol, "hf/minix", 1.8559832674618637E-2_wp)

end subroutine test_hf_minix


subroutine test_hf_sv(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless06")
   call test_generic(error, mol, "hf/sv", 6.4701104644333107E-3_wp)

end subroutine test_hf_sv


subroutine test_hf_def2sv_p(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless07")
   call test_generic(error, mol, "hf/def2-sv(p)", 1.5009548340274959E-2_wp)

end subroutine test_hf_def2sv_p


subroutine test_hf_def2svp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless08")
   call test_generic(error, mol, "hf/def2-svp", 1.8249078757662872E-2_wp)

end subroutine test_hf_def2svp


subroutine test_hf_dzp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless08")
   call test_generic(error, mol, "hf/dzp", 2.4002048454141937E-2_wp)

end subroutine test_hf_dzp


subroutine test_hf_deftzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless09")
   call test_generic(error, mol, "hf/def-tzvp", 4.9799363418100567E-3_wp)

end subroutine test_hf_deftzvp


subroutine test_hf_def2tzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless10")
   call test_generic(error, mol, "hf/def2-tzvp", 2.1614549442281130E-3_wp)

end subroutine test_hf_def2tzvp


subroutine test_hf_ccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless11")
   call test_generic(error, mol, "hf/cc-pvdz", 9.5654476866715681E-3_wp)

end subroutine test_hf_ccpvdz


subroutine test_hf_augccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless12")
   call test_generic(error, mol, "hf/aug-cc-pvdz", 2.6289439007561133E-3_wp)

end subroutine test_hf_augccpvdz


subroutine test_dft_minis(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless13")
   call test_generic(error, mol, "dft/minis", 2.8831779044393737E-2_wp)

end subroutine test_dft_minis


subroutine test_dft_minix(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless14")
   call test_generic(error, mol, "dft/minix", 1.8191237759356119E-2_wp)

end subroutine test_dft_minix


subroutine test_dft_sv(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless15")
   call test_generic(error, mol, "dft/sv", 8.9178105933847233E-3_wp)

end subroutine test_dft_sv


subroutine test_dft_def2sv_p(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless16")
   call test_generic(error, mol, "dft/def2-sv(p)", 9.0695660381885229E-3_wp)

end subroutine test_dft_def2sv_p


subroutine test_dft_def2svp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless17")
   call test_generic(error, mol, "dft/def2-svp", 2.4834643199826180E-2_wp)

end subroutine test_dft_def2svp


subroutine test_dft_dzp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless18")
   call test_generic(error, mol, "dft/dzp", 2.2195075119243783E-2_wp)

end subroutine test_dft_dzp


subroutine test_dft_deftzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless19")
   call test_generic(error, mol, "dft/def-tzvp", 3.6258064909241677E-3_wp)

end subroutine test_dft_deftzvp


subroutine test_dft_def2tzvp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless20")
   call test_generic(error, mol, "dft/def2-tzvp", 3.9044561047188592E-3_wp)

end subroutine test_dft_def2tzvp


subroutine test_dft_ccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless21")
   call test_generic(error, mol, "dft/cc-pvdz", 9.9649294926933245E-3_wp)

end subroutine test_dft_ccpvdz


subroutine test_dft_augccpvdz(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "mindless22")
   call test_generic(error, mol, "dft/aug-cc-pvdz", 2.9849348931724447E-3_wp)

end subroutine test_dft_augccpvdz


end module test_gcp
