!-------------------------------------------------------------------------------
! Program  : ZSC-2D-DGR
! Revision : 1.0 (2008-06-15)
! Author   : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Driver for the dual grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM
!> @details
!! Driver for the dual grid implementation of the Zheng-Shu-Chew multiphase LBM
!! using D2Q5/D2Q9 discretization and periodic boundary conditions. For details:
!!
!! Journal of Computational Physics 218: 353-371, 2006.
!!
!! Serial Implementation. Convergence test and post-processing functions
!! designed for single drop/bubble cases.
!!
!! In this dual grid implementation the order parameter grid is twice larger
!! than the coarser momentum grid. This provides improved numerical
!! differentials at the interface and leads to better mass conservation and
!! greater stability. All differential terms are calculated in the order
!! parameter grid and then copied over to the momentum grid (the grids overlap).
!! The velocity and pressure are calculated in the momentum grid and
!! interpolated in the order parameter grid. Two collision-stream steps are
!! are taken in the order parameter grid for every step taken in the momentum
!! grid to ensure synchronization.
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
 IMPLICIT NONE


! Read input parameters, allocate memory for common arrays and initialize
 CALL Parameters
 CALL MemAlloc(1)
 CALL Init

! Main iteration loop
 DO iStep = 1, MaxStep

! Full Step for momentum
   CALL CollisionG
   CALL Update

! Fractional Step 1 for order parameter
   CALL CollisionF
   CALL Stream
   CALL Differentials

! Fractional Step 2 for order parameter
   CALL CollisionF
   CALL Stream
   CALL Differentials

! Update momentum mesh index
   now = 1 - now
   nxt = 1 - nxt

! Save data if required
   IF ( MOD(istep,tStat) == 0 ) CALL Stats
   IF ( MOD(istep,tDump) == 0 ) CALL VtkPlane

! Test pressure convergence
   IF ( eps < pConv ) EXIT

 END DO
 RelaxStep = iStep - 1

! Save final data
 CALL VtkPlane
 CALL FinalDump

! Free memory
 CALL MemAlloc(2)


END PROGRAM
