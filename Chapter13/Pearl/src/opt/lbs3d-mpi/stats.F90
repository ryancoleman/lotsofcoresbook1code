!-------------------------------------------------------------------------------
! Subroutine : Stats
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Calculate average velocity, mass conservation factor and effective radius of
! the drop, and write them to file "stats.out" for the hybrid D3Q7/D3Q19
! Zheng-Shu-Chew multiphase LBM.
!
! The effective radius is calculated assuming the drop is a perfect circle with
! area given by Vol = (4/3)*Pi*R*R*R.
!
! Parallel implementation using MPI. Variables with the "Local" ending refer to
! quantities calculated locally in the current processor, vproc. Variables with
! the ending "Global" refer to qunatities calculated on the complete domain.
! The parameters stored in the MPI excahnge buffers are defined below.
!
! data(1) : Volume
! data(2) : Ux
! data(3) : Uy
! data(4) : Uz
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

 SUBROUTINE Stats

! Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams
 USE MPIParams
 USE MPI
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m, IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: ux, uy, uz, rhon, invRhon
 REAL(KIND = DBL) :: invVol, Ref, Vef, Vol
 REAL(KIND = DBL), DIMENSION(1:4) :: dataLocal, dataGlobal

! Initialize
 Vol = 0.D0
 ux  = 0.D0
 uy  = 0.D0
 uz  = 0.D0
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0

! Loop through all nodes inside the drop for the current processor, vproc
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX

       m = i + NXG*( j + NYG*k ) + now

! Calculate the accumulated quantities
       IF ( phi(m-now) >= 0.D0 ) THEN
         rhon    = g0(m-now)  + g1(m)  + g2(m)  + g3(m)  + g4(m)  + g5(m)  + g6(m)  & 
                 + g7(m)  + g8(m)  + g9(m)  + g10(m) + g11(m) + g12(m) + g13(m) & 
                 + g14(m) + g15(m) + g16(m) + g17(m) + g18(m)
         invRhon = 1.D0/rhon

         dataLocal(1) = dataLocal(1) + 1.D0

         dataLocal(2) = dataLocal(2) + ( g1(m)  - g2(m)  + g7(m) &
                      - g8(m) + g9(m) - g10(m) + g11(m)          &
                      - g12(m) + g13(m) - g14(m) )*invRhon
         dataLocal(3) = dataLocal(3) + ( g3(m)  - g4(m)  + g7(m) &
                      - g8(m) - g9(m) + g10(m) + g15(m)          &
                      - g16(m) + g17(m) - g18(m) )*invRhon
         dataLocal(4) = dataLocal(4) + ( g5(m)  - g6(m)  + g11(m) &
                      - g12(m) - g13(m) + g14(m)                  &
                      + g15(m) - g16(m) - g17(m) + g18(m) )*invRhon
       END IF

     END DO
   END DO
 END DO

! Gather global information
 CALL MPI_ALLREDUCE( dataLocal, dataGlobal, 4, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_VGRID, MPI_ERR )
 Vol    = dataGlobal(1)
 invVol = 1.D0/Vol

! Define average velocity of the drop and effective radius and volume
 ux  = dataGlobal(2)*invVol
 uy  = dataGlobal(3)*invVol
 uz  = dataGlobal(4)*invVol
 Ref = ( Vol*invPi*0.75D0 )**(1.D0/3.D0)
 Vef = Vol*invInitVol

! Save velocity and effective radius and volume data from the master node only
 IF ( proc == master ) THEN
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
