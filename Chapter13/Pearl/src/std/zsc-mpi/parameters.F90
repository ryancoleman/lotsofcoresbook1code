!-------------------------------------------------------------------------------
! Subroutine : Parameters
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!              Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Read input parameters and define simulation constants
!> @details
!! Read input parameters from files "properties.in" and "discrete.in" and define
!! constants for the simulation in the parallel D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM.
!!
!! Data files are read from a single processor (master), broadcasted to all
!! processors, and saved into local variables.

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
 USE NTypes, ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams, ONLY : master, mpi_xdim, mpi_ydim, proc
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: i
 INTEGER :: IO_ERR, MPI_ERR
 INTEGER, DIMENSION (1:10) :: intParam
 REAL(KIND = DBL), DIMENSION (1:8) :: dblParam
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bubblesVec


! Read parameter data from file in master node
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
     READ(10,*) tauRho, tauPhi
     READ(10,*)
     READ(10,*) IntWidth, sigma, Gamma
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
   dblParam(3) = tauRho
   dblParam(4) = tauPhi
   dblParam(5) = IntWidth
   dblParam(6) = sigma
   dblParam(7) = Gamma
   dblParam(8) = pConv

!  Read bubble positions
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
 CALL MPI_BCAST(dblParam,  8, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, MPI_ERR)

! Assign parameter data to local-named variables in slave nodes
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
   tauRho   = dblParam(3)
   tauPhi   = dblParam(4)
   IntWidth = dblParam(5)
   sigma    = dblParam(6)
   Gamma    = dblParam(7)
   pConv    = dblParam(8)

   ALLOCATE( bubbles(1:nBubbles,1:3) )
   ALLOCATE( bubblesVec(1:nBubbles*3) )
 END IF

!  Broadcast bubble(s) information from master to all other processors
 CALL MPI_BCAST(bubblesVec, nBubbles*3, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, MPI_ERR)

!  Save in local array bubbles(xc,yc,Rc)
 DO i = 1, nBubbles
   bubbles(i,1) = bubblesVec(i)            ! xc(i)
   bubbles(i,2) = bubblesVec(i+nBubbles)   ! yc(i)
   bubbles(i,3) = bubblesVec(i+2*nBubbles) ! Rc(i)
 END DO
 DEALLOCATE( bubblesVec )

! Fluid properties
 phiStar    = 0.5D0*(rhoH - rhoL)
 phiStar2   = phistar*phistar
 phiStar4   = phistar2*phistar2
 eta        = 1.D0/(tauphi + 0.5D0)
 eta2       = 1.D0 - eta
 invEta2    = 0.5D0/eta
 invTauRho  = 1.D0/taurho
 invTauRho2 = 1.D0 - 0.5D0*invtaurho
 invTauPhi  = 1.D0/tauphi

! Chemical potential parameters
 alpha  = 0.75D0*sigma/(IntWidth*phiStar4)
 alpha4 = alpha*4.D0
 kappa  = 0.50D0*alpha*(IntWidth*phiStar)*(IntWidth*phiStar)

! Modified distribution weights for use in Collision
 Eg0T = Eg0C*invTauRho2
 Eg1T = Eg1C*invTauRho2
 Eg2T = Eg2C*invTauRho2

! Redefine MaxStep (iterate between 0 and Maxstep-1)
 MaxStep = MaxStep - 1

 RETURN
 END SUBROUTINE Parameters
