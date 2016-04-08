!-------------------------------------------------------------------------------
! Subroutine : VtkPlane
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Save simulation data in VTK format
!> @details
!! Save order parameter, pressure, and velocity data in VTK format for the
!! parallel D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM. Filenames are
!! 'DATA_YYY_XXXXXXX.VTK' where YYY is the processor number and XXXXXXX is the
!! time step. These files are understood by Paraview.

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

 SUBROUTINE VtkPlane

! Common Variables
 USE Domain,      ONLY : iStep, NX, NY, xl, xu, yl, yu
 USE FluidParams, ONLY : p, phi, u
 USE MPI
 IMPLICIT NONE

! Required functions
 CHARACTER(12) :: filer

! Local Variables
 INTEGER :: i, j
 INTEGER :: IO_ERR, MPI_ERR
 CHARACTER(14) :: filer1
 CHARACTER(16) :: fmtd1, fmtd3
 CHARACTER(20) :: filer2

! Set output files formats and names
 fmtd1 = '(1(ES21.11E3))'
 fmtd3 = '(3(ES21.11E3))'
 filer1 = filer()
 filer2 = "DATA"//filer1(1:12)//".vtk"

 OPEN(UNIT = 12, FILE = filer2, STATUS = "NEW", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(12,'(A)')"# vtk DataFile Version 2.0"
   WRITE(12,'(A)')"MP-LABS v1.2"
   WRITE(12,'(A)')"ASCII"
   WRITE(12,*)
   WRITE(12,'(A)')"DATASET STRUCTURED_POINTS"
   WRITE(12,*)"DIMENSIONS",NX,NY,1
   WRITE(12,*)"ORIGIN",xl,yl,1
   WRITE(12,*)"SPACING",1,1,1
   WRITE(12,*)
   WRITE(12,*)"POINT_DATA",NX*NY
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Phi double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO j = yl, yu
       DO i = xl, xu
          WRITE(UNIT = 12,FMT = fmtd1)phi(i,j)
       END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Pressure double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO j = yl, yu
       DO i = xl, xu
          WRITE(UNIT = 12,FMT = fmtd1)p(i,j)
       END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"VECTORS Velocity double"
   DO j = yl, yu
       DO i = xl, xu
         WRITE(UNIT = 12,FMT = fmtd3)u(i,j,1),u(i,j,2),0.D0
      END DO
   END DO
   CLOSE(UNIT = 12)
 ELSE
   CALL MemAlloc(2)
   CALL MPI_FINALIZE(MPI_ERR)
   STOP "Error: Unable to open output vtk file."
 END IF

 RETURN
 END SUBROUTINE VtkPlane

!-------------------------------------------------------------------------------
! Function : filer
! Revision : 1.0 (2008-06-15)
! Author   : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @brief Set vtk file name as a function of processor number and time step.
 FUNCTION filer()

! Common variables
 USE Domain,    ONLY : iStep
 USE MPIParams, ONLY : vproc
 IMPLICIT NONE

! Function type
 CHARACTER(12) :: filer

!-------- Files name depending on processor (vproc) and timestep (iStep) -------
 IF (vproc < 10) THEN
   IF (iStep < 10) THEN
     WRITE(filer,9001)vproc,iStep
   ELSE IF (iStep < 100) THEN
     WRITE(filer,9002)vproc,iStep
   ELSE IF (iStep < 1000) THEN
     WRITE(filer,9003)vproc,iStep
   ELSE IF (iStep < 10000) THEN
     WRITE(filer,9004)vproc,iStep
   ELSE IF (iStep < 100000) THEN
     WRITE(filer,9005)vproc,iStep
   ELSE IF (iStep < 1000000) THEN
     WRITE(filer,9006)vproc,iStep
   ELSE IF (iStep < 10000000) THEN
     WRITE(filer,9007)vproc,iStep
   END IF
 ELSE IF (vproc < 100) THEN
   IF (iStep < 10) THEN
     WRITE(filer,9011)vproc,iStep
   ELSE IF (iStep < 100) THEN
     WRITE(filer,9012)vproc,iStep
   ELSE IF (iStep < 1000) THEN
     WRITE(filer,9013)vproc,iStep
   ELSE IF (iStep < 10000) THEN
     WRITE(filer,9014)vproc,iStep
   ELSE IF (iStep < 100000) THEN
     WRITE(filer,9015)vproc,iStep
   ELSE IF (iStep < 1000000) THEN
     WRITE(filer,9016)vproc,iStep
   ELSE IF (iStep < 10000000) THEN
     WRITE(filer,9017)vproc,iStep
   END IF
  ELSE
   IF (iStep < 10) THEN
     WRITE(filer,9111)vproc,iStep
   ELSE IF (iStep < 100) THEN
     WRITE(filer,9112)vproc,iStep
   ELSE IF (iStep < 1000) THEN
     WRITE(filer,9113)vproc,iStep
   ELSE IF (iStep < 10000) THEN
     WRITE(filer,9114)vproc,iStep
   ELSE IF (iStep < 100000) THEN
     WRITE(filer,9115)vproc,iStep
   ELSE IF (iStep < 1000000) THEN
     WRITE(filer,9116)vproc,iStep
   ELSE IF (iStep < 10000000) THEN
     WRITE(filer,9117)vproc,iStep
   END IF
 END IF

!-------- File name format depending on processor (vproc) and timestep (iStep) -
 9001 FORMAT('_00',i1,'_000000',i1)
 9002 FORMAT('_00',i1,'_00000',i2)
 9003 FORMAT('_00',i1,'_0000',i3)
 9004 FORMAT('_00',i1,'_000',i4)
 9005 FORMAT('_00',i1,'_00',i5)
 9006 FORMAT('_00',i1,'_0',i6)
 9007 FORMAT('_00',i1,'_',i7)
 9011 FORMAT('_0',i2,'_000000',i1)
 9012 FORMAT('_0',i2,'_00000',i2)
 9013 FORMAT('_0',i2,'_0000',i3)
 9014 FORMAT('_0',i2,'_000',i4)
 9015 FORMAT('_0',i2,'_00',i5)
 9016 FORMAT('_0',i2,'_0',i6)
 9017 FORMAT('_0',i2,'_',i7)
 9111 FORMAT('_',i3,'_000000',i1)
 9112 FORMAT('_',i3,'_00000',i2)
 9113 FORMAT('_',i3,'_0000',i3)
 9114 FORMAT('_',i3,'_000',i4)
 9115 FORMAT('_',i3,'_00',i5)
 9116 FORMAT('_',i3,'_0',i6)
 9117 FORMAT('_',i3,'_',i7)


 RETURN
 END FUNCTION filer
