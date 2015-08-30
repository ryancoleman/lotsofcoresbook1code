!-------------------------------------------------------------------------------
! Program  : LL-2D-STD-MPI
! Revision : 1.0 (2008-06-15)
! Author   : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Driver for the parallel D2Q9 Lee-Lin multiphase LBM.
!> @details
!! Driver for the standard implementation of the Lee-Lin multiphase LBM
!! using D2Q9 discretization and periodic boundary conditions. For details:
!!
!! Journal of Computational Physics 206: 16-47, 2005.
!!
!! Parallel implementation using MPI. Convergence test and post-processing
!!  functions designed for single drop/bubble cases.
!!
!! The average velocity, mass conservation factor, effective radius of the drop,
!! pressure difference between the inside and the outside of the drop and the
!! error with respect to the analytical value given by Laplace's equation are
!! written to file "stats.out"

!-------------------------------------------------------------------------------
! Copyright 2008 Carlos Rosales Fernandez, David S. Whyte and IHPC (A*STAR).
!
! This file is part of MP-LABS.
!
! MP-LABS is free software: you can redistribute it and/or modify it under the
! terms of the GNU GPL version 3 or (at your option) any later version.
!
! MP-LABS is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! MP-LABS, in the file COPYING.txt. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------

 PROGRAM main

! Common Variables
 USE Domain,      ONLY : iStep, MaxStep, RelaxStep, tDump, tStat
 USE FluidParams, ONLY : eps, pConv
 USE MPIParams,   ONLY : nprocs, proc
 USE MPI
 IMPLICIT NONE

! Local quantities
 INTEGER :: MPI_ERR

! Initialize MPI Environment
 CALL MPI_INIT(MPI_ERR)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, MPI_ERR)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc, MPI_ERR)

! Read input parameters, broadcast and reate virtual CPU grid
 CALL Parameters
 CALL VGrid

! Allocate memory for common arrays and initialize
 CALL MemAlloc(1)
 CALL Init

! Main Loop
 DO iStep = 1, MaxStep

   CALL PreStream
   CALL MPI_UpdateFG
   CALL MPI_UpdateRho
   CALL Hydrodynamics
   CALL PostStream

! Save data if required
   IF ( MOD(iStep,tStat) == 0 ) CALL Stats
   IF ( MOD(iStep,tDump) == 0 ) CALL VtkPlane

! Test pressure convergence
   IF ( eps <= pConv ) EXIT

 END DO
 RelaxStep = iStep - 1

! Save final data
 CALL FinalDump
 CALL VtkPlane

! Free common arrays, finalize MPI communications and end program
 CALL MemAlloc(2)
 CALL MPI_FINALIZE(MPI_ERR)

 END PROGRAM
