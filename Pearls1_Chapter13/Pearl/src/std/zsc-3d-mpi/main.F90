!-------------------------------------------------------------------------------
! Program  : ZSC-3D-STD-MPI
! Revision : 1.0 (2008-06-15)
! Author   : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!            Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Driver for the parallel D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM
!> @details
!! Driver for the standard parallel implementation of the Zheng-Shu-Chew
!! multiphase LBM using D3Q7/D3Q19 discretization and periodic boundary
!! conditions, including the gravitational force. For details:
!!
!! Journal of Computational Physics 218: 353-371, 2006.
!!
!! Parallel Implementation using MPI. Convergence test and post-processing
!! functions designed for single drop/bubble cases.
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
 USE Domain,      ONLY : iStep, MaxStep, now, nxt, RelaxStep, tDump, tStat
 USE FluidParams, ONLY : eps, pConv
 USE MPIParams,   ONLY : nprocs, proc
 USE MPI
 IMPLICIT NONE

! Llocal variables
 INTEGER :: MPI_ERR


!--------- Initialize MPI Environment ------------------------------------------
 CALL MPI_INIT(MPI_ERR)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, MPI_ERR)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc, MPI_ERR) 

! Read input parameters, broadcast, and create virtual CPU grid
 CALL Parameters
 CALL VGrid

! Allocate memory for common arrays and initialize
 CALL MemAlloc(1)
 CALL Init

! Interface relaxation loop
 DO iStep = 1, RelaxStep

   CALL CollisionRelax
   CALL PostCollision
   CALL Stream
   CALL PostStream

   now = 1 - now
   nxt = 1 - nxt

   IF( MOD(iStep,tStat) == 0 ) CALL RelaxStats
   IF( MOD(iStep,tDump) == 0 ) CALL VtkDump

   IF( eps < pConv ) EXIT

 END DO
 RelaxStep = iStep

! Save relaxed configuration data
 CALL Stats
 CALL VtkDump

! Main Iteration Loop
 DO iStep = RelaxStep, Maxstep

   CALL Collision
   CALL PostCollision
   CALL Stream
   CALL PostStream

   now = 1 - now
   nxt = 1 - nxt

   IF( MOD(iStep,tStat) == 0 ) CALL Stats
   IF( MOD(iStep,tDump) == 0 ) CALL VtkDump

 END DO

! Save final data
 CALL FinalDump
 CALL VtkDump

! Free memory, finalize MPI and end program
 CALL MemAlloc(2)
 CALL MPI_FINALIZE(MPI_ERR)

 END PROGRAM
