!-------------------------------------------------------------------------------
! Subroutine : Parameters
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Read input parameters from files "properties.in" and "discrete.in" and define
! constants for the simulation in the hybrid MPI+OMP D3Q7/D3Q19 Zheng-Shu-Chew
! multiphase LBM.
!
! DISTRO determines if the bubbles are generated randomly using the parameters
! (xb,yb,zb) or following given positions:
! DISTRO == 1 -> Random bubble distribution
! DISTRO != 2 -> Given bubble positions
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

 SUBROUTINE Parameters

! Common variables
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams, ONLY : master, mpi_xdim, mpi_ydim, mpi_zdim, proc
 USE MPI
 IMPLICIT NONE

! Local variables
 INTEGER :: i, j, k, xb, yb, zb, Rb, DISTRO, IO_ERR, MPI_ERR
 INTEGER, DIMENSION (1:14) :: intParam
 REAL(KIND = DBL) :: dx, dy, dz, xdmax, ydmax, zdmax, x, y, z, Lx, Ly, Lz
 REAL(KIND = DBL), DIMENSION (1:9) :: dblParam
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bubblesVec


! Read parameter data from file in master node
 IF (proc == master) THEN
   OPEN(UNIT = 10, FILE = "properties.in", STATUS = "OLD", ACTION = "READ", &
        IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     READ(10,*)
     READ(10,*) Maxstep, RelaxStep
     READ(10,*)
     READ(10,*) tStat, tSave
     READ(10,*)
     READ(10,*) xmin, xmax, ymin, ymax, zmin, zmax
     READ(10,*)
     READ(10,*) rhoL, rhoH
     READ(10,*)
     READ(10,*) tauRho, tauPhi
     READ(10,*)
     READ(10,*) IntWidth, sigma, Gamma
     READ(10,*)
     READ(10,*) Eo, pConv
     READ(10,*)
     READ(10,*) mpi_xdim, mpi_ydim, mpi_zdim
     CLOSE(UNIT = 10)
   ELSE
     STOP "Error: Unable to open input file 'properties.in'."
   END IF

   intParam(1)  = MaxStep
   intParam(2)  = RelaxStep
   intParam(3)  = tStat
   intParam(4)  = tSave
   intParam(5)  = xmin
   intParam(6)  = xmax
   intParam(7)  = ymin
   intParam(8)  = ymax
   intParam(9) = zmin
   intParam(10) = zmax
   intParam(11) = mpi_xdim
   intParam(12) = mpi_ydim
   intParam(13) = mpi_zdim

   dblParam(1) = rhoL
   dblParam(2) = rhoH
   dblParam(3) = tauRho
   dblParam(4) = tauPhi
   dblParam(5) = Gamma
   dblParam(6) = IntWidth
   dblParam(7) = sigma
   dblParam(8) = Eo
   dblParam(9) = pConv

! Read bubble positions
   OPEN(UNIT = 10, FILE = "discrete.in", STATUS = "OLD", ACTION = "READ", &
      IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     READ(10,*)
     READ(10,*) DISTRO
     IF(DISTRO == 1) THEN
       READ(10,*)
       READ(10,*) xb, yb, zb, Rb
       CLOSE(UNIT = 10)
       nBubbles = xb*yb*zb
       Lx = DBLE(xmax)/DBLE(xb)
       Ly = DBLE(ymax)/DBLE(yb)
       Lz = DBLE(zmax)/DBLE(zb)
       xdmax = 0.5D0*Lx - DBLE(Rb)
       ydmax = 0.5D0*Ly - DBLE(Rb)
       zdmax = 0.5D0*Lz - DBLE(Rb)
       CALL RANDOM_SEED
       ALLOCATE( bubbles(1:nBubbles,1:4) )
       ALLOCATE( bubblesVec(1:nBubbles*4) )
       OPEN(UNIT = 10, FILE = "discrete.out", STATUS = "NEW", IOSTAT = IO_ERR)
       IF ( IO_ERR == 0 ) THEN
         WRITE(10,*)'Lx = ',Lx,' Ly = ',Ly,' Lz = ',Lz
         DO i = 1, xb
           DO j = 1, yb
             DO k = 1, zb
               CALL RANDOM_NUMBER(dx)
               CALL RANDOM_NUMBER(dy)
               CALL RANDOM_NUMBER(dz)
               x = (i - 0.5D0)*Lx + xdmax*(1.0D0 - 2.0D0*dx)
               y = (j - 0.5D0)*Ly + ydmax*(1.0D0 - 2.0D0*dy)
               z = (k - 0.5D0)*Lz + zdmax*(1.0D0 - 2.0D0*dz)
               bubblesVec(i)              = x
               bubblesVec(i + nBubbles)   = y
               bubblesVec(i + 2*nBubbles) = z
               bubblesVec(i + 3*nBubbles) = Rb
               WRITE(10,*)dx, dy, dz, x, y, z, Rb
             END DO
           END DO
         END DO
         CLOSE(UNIT = 10)
       ELSE
         STOP "Error: Unable to open output file 'discrete.out'."
       END IF
     ELSE
       READ(10,*)
       READ(10,*) nBubbles
       READ(10,*)
       ALLOCATE( bubbles(1:nBubbles,1:4) )
       ALLOCATE( bubblesVec(1:nBubbles*4) )
       DO i = 1, nBubbles
         READ(10,*) bubblesVec(i), bubblesVec(i+nBubbles), bubblesVec(i+2*nBubbles), bubblesVec(i+3*nBubbles)
       END DO
     END IF
   CLOSE(UNIT=10)
   ELSE
     STOP "Error: Unable to open input file 'discrete.in'."
   END IF
   intParam(14) = nBubbles
 END IF

! Broadcast values of parameters from the master to all other processors
 CALL MPI_BCAST(intParam, 14, MPI_INTEGER, master, MPI_COMM_WORLD, MPI_ERR)
 CALL MPI_BCAST(dblParam,  9, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, MPI_ERR)

! Assign parameter data to local-named variables in slave nodes
 IF (proc /= master) THEN
   MaxStep   = intParam(1)
   RelaxStep = intParam(2)
   tStat     = intParam(3)
   tSave     = intParam(4)
   xmin      = intParam(5)
   xmax      = intParam(6)
   ymin      = intParam(7)
   ymax      = intParam(8)
   zmin      = intParam(9)
   zmax      = intParam(10)
   mpi_xdim  = intParam(11)
   mpi_ydim  = intParam(12)
   mpi_zdim  = intParam(13)
   nBubbles  = intParam(14)

   rhoL     = dblParam(1)
   rhoH     = dblParam(2)
   taurho   = dblParam(3)
   tauphi   = dblParam(4)
   Gamma    = dblParam(5)
   IntWidth = dblParam(6)
   sigma    = dblParam(7)
   Eo       = dblParam(8)
   pConv    = dblParam(9)

   ALLOCATE( bubbles(1:nBubbles,1:4) )
   ALLOCATE( bubblesVec(1:nBubbles*4) )
 END IF


! Broadcast bubble(s) information from master to all other processors
 CALL MPI_BCAST( bubblesVec, nBubbles*4, MPI_DOUBLE_PRECISION, master, &
                 MPI_COMM_WORLD, MPI_ERR )

! Save in local array bubbles(xc,yc,zc,Rc)
 DO i = 1, nBubbles
   bubbles(i,1) = bubblesVec(i)              ! xc(i)
   bubbles(i,2) = bubblesVec(i + nBubbles)   ! yc(i)
   bubbles(i,3) = bubblesVec(i + 2*nBubbles) ! zc(i)
   bubbles(i,4) = bubblesVec(i + 3*nBubbles) ! Rc(i)
 END DO
 DEALLOCATE( bubblesVec )

! Fluid properties
 phiStar       = 0.5D0*(rhoH - rhoL)
 phiStar2      = phiStar*phiStar
 invTauRho     = 1.D0/tauRho
 invTauRhoOne  = 1.D0 - invTauRho
 invTauRhoHalf = 1.D0 - 0.5D0*invTauRho

 eta        = 1.D0/(tauPhi + 0.5D0)
 eta2       = 1.D0 - eta
 invEta2    = 0.5D0/eta
 invTauPhi  = 1.D0/tauPhi
 invTauPhi1 = 1.D0 - invTauPhi

! Chemical Potential Stuff
 alpha  = 0.75D0*sigma/(IntWidth*phistar**4.D0)
 alpha4 = alpha*4.D0
 kappa  = (IntWidth*phistar)**2.D0*alpha/2.D0

! Rc = bubbles(1,4) Use the first bubble as the largest in the system
 gravity = 0.25D0*Eo*sigma/( (rhoH - rhoL)*bubbles(1,4)*bubbles(1,4) )
 grav    = 0.D0

! Modified LBM parameters
 Eg0n  = invTauRho*Eg0
 Eg1n  = invTauRho*Eg1
 Eg2n  = invTauRho*Eg2

 EgC0n = invTauRhoHalf*Eg0C
 EgC1n = invTauRhoHalf*Eg1C
 EgC2n = invTauRhoHalf*Eg2C

 RETURN
 END SUBROUTINE Parameters
