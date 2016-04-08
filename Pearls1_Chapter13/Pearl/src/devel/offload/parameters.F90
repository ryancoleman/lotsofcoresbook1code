!-------------------------------------------------------------------------------
! Subroutine : Parameters
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Read input parameters from files "properties.in" and "discrete.in" and define
! constants for the simulation in the offload D3Q7/D3Q19 Zheng-Shu-Chew
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
 USE NTypes, ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local variables
 INTEGER :: i, j, k, xb, yb, zb, Rb, DISTRO, IO_ERR
 REAL(KIND = DBL) :: dx, dy, dz, xdmax, ydmax, zdmax, x, y, z, Lx, Ly, Lz
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bubblesVec


! Read parameter data from file in master node
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
     READ(10,*) MIC_FRACTION
     CLOSE(UNIT = 10)
   ELSE
     STOP "Error: Unable to open input file 'properties.in'."
   END IF

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
               bubbles(i,1) = x
               bubbles(i,2) = y
               bubbles(i,3) = z
               bubbles(i,4) = Rb
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
       DO i = 1, nBubbles
         READ(10,*) bubbles(i,1), bubbles(i,2), bubbles(i,3), bubbles(i,4)
       END DO
     END IF
   CLOSE(UNIT=10)
   ELSE
     STOP "Error: Unable to open input file 'discrete.in'."
   END IF

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

! Middle planes (x = xo), (y = yo), (z = zo)
 xo = INT(0.5D0*(xmax + xmin))
 yo = INT(0.5D0*(ymax + ymin))
 zo = INT(0.5D0*(zmax + zmin))

 xl = xmin; xlg = xmin - 1
 xu = xmax; xug = xmax + 1
 yl = ymin; ylg = ymin - 1
 yu = ymax; yug = ymax + 1
 zl = zmin; zlg = zmin - 1
 zu = zmax; zug = zmax + 1

! Set domain limits
 NX = xmax
 NY = ymax
 NZ = zmax

! Set ghost limits
 NXG = NX + 2
 NYG = NY + 2
 NZG = NZ + 2
 NG  = NXG*NYG*NZG

! Set array limits in MIC coprocessor
! Note I am assigning 50% of work to CPU and 50% to MIC
! This should be a variable.
 NY_MIC  = NY*MIC_FRACTION
 NYG_MIC = NY_MIC+2
 NG_MIC  = NXG*NYG_MIC*NZG

 NXYG_MIC = NXG*NYG_MIC

 RETURN
 END SUBROUTINE Parameters
