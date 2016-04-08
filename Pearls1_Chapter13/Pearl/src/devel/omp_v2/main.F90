!-------------------------------------------------------------------------------
! Program    : LBS3D-DEVEL-OMPv2
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Driver for the OMP_V2 development implementation of the Zheng-Shu-Chew
! multiphase LBM using D3Q7/D3Q19 discretization and periodic boundary
! conditions, including the gravitational force. For details:
!
! Journal of Computational Physics 218: 353-371, 2006.
!
! The average velocity, mass conservation factor, effective radius of the drop,
! pressure difference between the inside and the outside of the drop and the
! error with respect to the analytical value given by Laplace's equation are
! written to file "stats.out"
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
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
 USE Timers
 USE Domain,      ONLY : iStep, MaxStep, RelaxStep, tSave, tStat, &
                         STAGE, CompletedRelaxSteps, NG, now, nxt
 USE FluidParams, ONLY : eps, pConv, phiStar, grav, gravity
 USE OMP_LIB
 IMPLICIT NONE

! Read input parameters, broadcast, and create virtual CPU grid
 tStart = omp_get_wtime()
 CALL Parameters
 now = 0
 nxt = 1

! Allocate memory for common arrays and initialize
 CALL MemAlloc(1)
 CALL Init

! Save initialized data
 CALL RelaxStats
#ifdef LBMIO
 CALL VtkSave
#endif

! Interface relaxation loop
 grav = 0.D0
 tCurrent = omp_get_wtime()
 tSetup = tCurrent - tStart
 DO iStep = 1, RelaxStep

   CALL Collision
   CALL PostCollision
   CALL Stream
   CALL PostStream

   now = 1 - now
   nxt = 1 - nxt

#ifdef LBMIO
   IF( MOD(iStep,tStat) == 0 ) CALL RelaxStats
   IF( MOD(iStep,tSave) == 0 ) CALL VtkSave
#endif

   IF( eps < pConv ) EXIT

 END DO
 tCurrent = omp_get_wtime()
 tRelax = tCurrent - tSetup - tStart
 CompletedRelaxSteps = iStep

! Save relaxed configuration data
 STAGE = 1
#ifdef LBMIO
 CALL Stats
 CALL VtkSave
 tCurrent = omp_get_wtime()
 tIO = tCurrent - tRelax - tSetup - tStart
#else
 tIO = 0.D0
#endif

! Main Iteration Loop
 grav = 2.D0*phiStar*gravity
 DO iStep = 1, Maxstep

   CALL Collision
   CALL PostCollision
   CALL Stream
   CALL PostStream

   now = 1 - now
   nxt = 1 - nxt

#ifdef LBMIO
   IF( MOD(iStep,tStat) == 0 ) CALL Stats
   IF( MOD(iStep,tSave) == 0 ) CALL VtkSave
#endif

 END DO
 tCurrent = omp_get_wtime()
! nIO is the approximate number of times we save the full domain data to file
! In most sane implementations this dominates the IO cost.
 nIO   = (CompletedRelaxSteps + MaxStep) / tSave
 tEvol = tCurrent - tIO*nIO - tRelax - tSetup - tStart
 tRun  = tCurrent - tStart

! Save final data
 CALL FinalSave
#ifdef LBMIO
 CALL VtkSave
#endif

! Free memory, finalize MPI and end program
 CALL MemAlloc(2)

 END PROGRAM
