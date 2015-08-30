!-------------------------------------------------------------------------------
! Subroutine : Stats
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Calculate intermediate simulation results and save to file 'stats.out'
!> @details
!! Calculate average velocity, mass conservation factor and effective radius of
!! the drop, and write them to file "stats.out" for the parallel D2Q5/D2Q9
!! Zheng-Shu-Chew multiphase LBM.
!!
!! The effective radius is calculated assuming the drop is a perfect circle with
!! area given by Vol = (4/3)*Pi*R*R*R.
!!
!! Parallel implementation using MPI. Variables with the "Local" ending refer to
!! quantities calculated locally in the current processor, vproc. Variables with
!! the ending "Global" refer to qunatities calculated on the complete domain.
!! The parameters stored in the MPI excahnge buffers are defined below.
!!
!! @param data(1) : Volume
!! @param data(2) : Ux
!! @param data(3) : Uy
!! @param data(4) : Uz

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

 SUBROUTINE Stats

! Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain,      ONLY : invInitVol, invPi, iStep, now, xl, xu, yl, yu, zl, zu
 USE FluidParams, ONLY : phi
 USE LBMParams,   ONLY : g
 USE MPIParams,   ONLY : master, MPI_COMM_VGRID, vproc
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k
 INTEGER :: IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: ux, uy, uz, rhon, invRhon
 REAL(KIND = DBL) :: invVol, Ref, Vef, Vol
 REAL(KIND = DBL), DIMENSION(1:4) :: dataLocal, dataGlobal

! Initialize
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0

! Loop through all nodes inside the drop for the current processor, vproc
 DO k = zl, zu
   DO j = yl, yu
     DO i = xl, xu

       rhon    = SUM( g(i,j,k,:,now) )
       invRhon = 1.D0/rhon

! Calculate the velocity
       ux = ( g(i,j,k, 1,now) - g(i,j,k, 2,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) + g(i,j,k, 9,now) - g(i,j,k,10,now) &
          +   g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) &
          -   g(i,j,k,14,now) )*invRhon

       uy = ( g(i,j,k, 3,now) - g(i,j,k, 4,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) - g(i,j,k, 9,now) + g(i,j,k,10,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) &
          -   g(i,j,k,18,now) )*invRhon

       uz = ( g(i,j,k, 5,now) - g(i,j,k, 6,now) + g(i,j,k,11,now) &
          -   g(i,j,k,12,now) - g(i,j,k,13,now) + g(i,j,k,14,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) &
          +   g(i,j,k,18,now) )*invRhon

! Calculate the accumulated quantities
       IF ( phi(i,j,k) >= 0.D0 ) THEN
         dataLocal(1) = dataLocal(1) + 1.D0
         dataLocal(2) = dataLocal(2) + ux
         dataLocal(3) = dataLocal(3) + uy
         dataLocal(4) = dataLocal(4) + uz
       END IF

     END DO
   END DO
 END DO

! Gather global information
 CALL MPI_ALLREDUCE(dataLocal, dataGlobal, 4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_VGRID, MPI_ERR)
 Vol    = dataGlobal(1)
 invVol = 1.D0/Vol

! Define average velocity of the drop and effective radius and volume
 ux  = dataGlobal(2)*invVol
 uy  = dataGlobal(3)*invVol
 uz  = dataGlobal(4)*invVol
 Ref = ( Vol*invPi*0.75D0 )**(1.D0/3.D0)
 Vef = Vol*invInitVol

! Save velocity and effective radius and volume data from the master node only
 IF( vproc == master ) THEN
   OPEN(UNIT = 10, FILE = "stats.out", STATUS = "UNKNOWN", POSITION = "APPEND",&
        IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     WRITE(10,'(I9,7ES19.9)')iStep,ux,uy,uz,Vef,Ref
     CLOSE(UNIT = 10)
   ELSE
     CALL MemAlloc(2)
     CALL MPI_FINALIZE(MPI_ERR)
     STOP "Error: Unable to open output file 'stats.out'."
   END IF
 END IF

 RETURN
 END SUBROUTINE Stats
