!-------------------------------------------------------------------------------
! Subroutine : Parameters
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Read input parameters and define simulation constants
!> @details
!! Read input parameters from files "properties.in" and "discrete.in" and define
!! constants for the simulation in the parallel D2Q9 Lee-Lin multiphase LBM.

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

 SUBROUTINE Parameters

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : inv12, inv6, MaxStep, tDump, tStat, xmax, xmin, ymax, ymin
 USE FluidParams
 USE LBMParams, ONLY : Cs_sq
 USE MPIParams, ONLY : master, mpi_xdim, mpi_ydim, proc
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: i
 INTEGER :: IO_ERR, MPI_ERR
 INTEGER, DIMENSION (1:10) :: intParam
 REAL(KIND = DBL), DIMENSION (1:7) :: dblParam
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bubblesVec

! Read parameter data from file "properties.in" in master node
 IF (proc == master) THEN
   OPEN(UNIT = 10, FILE = "properties.in", STATUS = "OLD", ACTION = "READ", &
        IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     READ(10,*)
     READ(10,*) Maxstep
     READ(10,*)
     READ(10,*) tStat, tDump
     READ(10,*)
     READ(10,*) xmin, xmax, ymin, ymax
     READ(10,*)
     READ(10,*) rhoL, rhoH
     READ(10,*)
     READ(10,*) tauL, tauH
     READ(10,*)
     READ(10,*) IntWidth, sigma
     READ(10,*)
     READ(10,*) pConv
     READ(10,*)
     READ(10,*) mpi_xdim, mpi_ydim
     CLOSE(UNIT = 10)
   ELSE
     CALL MPI_FINALIZE(MPI_ERR)
     STOP "Error: Unable to open input file 'properties.in'."
   END IF

   intParam(1) = MaxStep
   intParam(2) = tStat
   intParam(3) = tDump
   intParam(4) = xmin
   intParam(5) = xmax
   intParam(6) = ymin
   intParam(7) = ymax
   intParam(8) = mpi_xdim
   intParam(9) = mpi_ydim

   dblParam(1) = rhoL
   dblParam(2) = rhoH
   dblParam(3) = tauL
   dblParam(4) = tauH
   dblParam(5) = IntWidth
   dblParam(6) = sigma
   dblParam(7) = pConv


! Read bubble/drop positions from file "discrete.in"
   OPEN(UNIT = 10, FILE = "discrete.in", STATUS = "OLD", ACTION = "READ", &
      IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     READ(10,*)
     READ(10,*) nBubbles
     READ(10,*)
     ALLOCATE( bubbles(1:nBubbles,1:3) )
     ALLOCATE( bubblesVec(1:nBubbles*3) )
     DO i = 1, nBubbles
       READ(10,*) bubblesVec(i), bubblesVec(i+nBubbles), bubblesVec(i+2*nBubbles)
     END DO
     CLOSE(UNIT = 10)
   ELSE
     CALL MPI_FINALIZE(MPI_ERR)
     STOP "Error: Unable to open input file 'discrete.in'."
   END IF

   intParam(10) = nBubbles
 END IF

! Broadcast values of parameters from the master to all other processors
 CALL MPI_BCAST(intParam, 10, MPI_INTEGER, master, MPI_COMM_WORLD, MPI_ERR)
 CALL MPI_BCAST(dblParam,  7, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, MPI_ERR)

! Assign parameter data to locally named variables in slave nodes
 IF(proc /= master) THEN
   MaxStep  = intParam(1)
   tStat    = intParam(2)
   tDump    = intParam(3)
   xmin     = intParam(4)
   xmax     = intParam(5)
   ymin     = intParam(6)
   ymax     = intParam(7)
   mpi_xdim = intParam(8)
   mpi_ydim = intParam(9)
   nBubbles = intParam(10)

   rhoL     = dblParam(1)
   rhoH     = dblParam(2)
   tauL     = dblParam(3)
   tauH     = dblParam(4)
   IntWidth = dblParam(5)
   sigma    = dblParam(6)
   pConv    = dblParam(7)

   ALLOCATE( bubbles(1:nBubbles,1:3) )
   ALLOCATE( bubblesVec(1:nBubbles*3) )
 END IF

! Broadcast bubble/drop information from master to all other processors
 CALL MPI_BCAST(bubblesVec, nBubbles*3, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, MPI_ERR)

! Save bubble/drop data in local array bubbles(xc,yc,Rc)
 DO i = 1, nBubbles
   bubbles(i,1) = bubblesVec(i)            ! xc(i)
   bubbles(i,2) = bubblesVec(i+nBubbles)   ! yc(i)
   bubbles(i,3) = bubblesVec(i+2*nBubbles) ! Rc(i)
 END DO
 DEALLOCATE( bubblesVec )

! Constants related to the interface tension and the free energy
 beta     = 12.D0*sigma/(IntWidth*( (rhoH - rhoL)**4 ))
 beta4    = 4.D0*beta
 kappa    = 1.5D0*sigma*IntWidth/( (rhoH - rhoL)**2 )
 kappa_6  = kappa*inv6
 kappa_12 = kappa*inv12
 rhostar  = 0.5D0*(rhoL + rhoH)
 tauRhoStar = (tauH - tauL)/(rhoH - rhoL)

! Dynamic viscosity of both liquids
 muL = tauL*rhoL*Cs_sq
 muH = tauH*rhoH*Cs_sq

! Redefine MaxStep (iterate between 0 and MaxStep-1)
 MaxStep = MaxStep - 1

 RETURN
 END SUBROUTINE Parameters
