!-------------------------------------------------------------------------------
! Subroutine : VtkSave
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Save order parameter, pressure, and velocity data in VTK format for the
! serial D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM. Filenames are
! 'EVOL_YYY_XXXXXXX.VTK' where YYY is the processor number and XXXXXXX is the
! time step. These files are understood by Paraview.
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

 SUBROUTINE VtkSave

! Common Variables
 USE Ntypes,      ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, k, m, mg, IO_ERR
 REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
 REAL(KIND = DBL) :: gradPhiSq
 CHARACTER(20) :: filename


! Set output files formats and names
 IF( STAGE == 0 ) THEN
    WRITE( filename, '(A,I3.3,A,I7.7,A)')"INIT_",vproc,"_",iStep,".vtk"
 ELSE
    WRITE( filename, '(A,I3.3,A,I7.7,A)')"EVOL_",vproc,"_",iStep,".vtk"
 END IF  
   

 OPEN(UNIT = 12, FILE = filename, STATUS = "NEW", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(12,'(A)')"# vtk DataFile Version 2.0"
   WRITE(12,'(A)')"MP-LABS v1.3"
   WRITE(12,'(A)')"ASCII"
   WRITE(12,*)
   WRITE(12,'(A)')"DATASET STRUCTURED_POINTS"
   WRITE(12,'(A,3(1x,I4))')"DIMENSIONS",NX,NY,NZ
   WRITE(12,'(A,3(1x,I4))')"ORIGIN",xmin,ymin,zmin
   WRITE(12,'(A,3(1x,I1))')"SPACING",1,1,1
   WRITE(12,*)
   WRITE(12,'(A,1x,I12)')"POINT_DATA",NX*NY*NZ
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Phi double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO k = 1, NZ
      DO j = 1, NY
         DO i = 1, NX
            WRITE(UNIT = 12,FMT = '(1(ES21.11E3))')phi( i + NXG*( j + NYG*k ) )
         END DO
      END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Pressure double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO k = 1, NZ
      DO j = 1, NY
         DO i = 1, NX

            m  = i + NXG*( j + NYG*k )
            mg = i + NXG*( j + NYG*k ) + now

            rhon  = g0(mg)  + g1(mg)  + g2(mg)  + g3(mg)  + g4(mg)  + g5(mg)  + g6(mg)  & 
                  + g7(mg)  + g8(mg)  + g9(mg)  + g10(mg) + g11(mg) + g12(mg) + g13(mg) & 
                  + g14(mg) + g15(mg) + g16(mg) + g17(mg) + g18(mg)
            phin  = phi(m)
            phin2 = phin*phin

            gradPhiSq = gradPhiX(m)*gradPhiX(m) &
                      + gradPhiY(m)*gradPhiY(m) &
                      + gradPhiZ(m)*gradPhiZ(m)

! Calculate the pressure
            pressure = alpha*( phin2*( 3.D0*phin2 - 2.D0*phiStar2 ) - phiStar4 ) &
                     - kappa*( phin*lapPhi(m) - 0.5D0*gradPhiSq ) + Cs_sq*rhon

            WRITE(UNIT = 12,FMT = '(1(ES21.11E3))')pressure
         END DO
      END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"VECTORS Velocity double"
   DO k = 1, NZ
      DO j = 1, NY
         DO i = 1, NX

            m = i + NXG*( j + NYG*k ) + now

            rhon  = g0(m)  + g1(m)  + g2(m)  + g3(m)  + g4(m)  + g5(m)  + g6(m)  & 
                  + g7(m)  + g8(m)  + g9(m)  + g10(m) + g11(m) + g12(m) + g13(m) & 
                  + g14(m) + g15(m) + g16(m) + g17(m) + g18(m)
           invRhon = 1.D0/rhon

           ux = ( g1(m) - g2(m) + g7(m) - g8(m) + g9(m) - g10(m) + g11(m) &
              -   g12(m) + g13(m) - g14(m) )*invRhon
           uy = ( g3(m) - g4(m) + g7(m) - g8(m) - g9(m) + g10(m) + g15(m) &
              -   g16(m) + g17(m) - g18(m) )*invRhon
           uz = ( g5(m) - g6(m) + g11(m) - g12(m) - g13(m) + g14(m) &
              +   g15(m) - g16(m) - g17(m) + g18(m) )*invRhon

            WRITE(UNIT = 12,FMT = '(3(ES21.11E3))')ux,uy,uz
         END DO
      END DO
   END DO
   CLOSE(UNIT = 12)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output vtk file."
 END IF

 RETURN
 END SUBROUTINE VtkSave

