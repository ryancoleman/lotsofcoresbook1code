!-------------------------------------------------------------------------------
! Subroutine : PackMIC_F
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Copy exchange values to MIC buffer for transfer to host CPU.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
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

 SUBROUTINE PackMIC_F
!DIR$ ATTRIBUTES OFFLOAD:MIC :: PackMIC_f

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, k, offset

 offset = NX*NZ

!$OMP PARALLEL

!-------------------------------------------------------------------------------
! We copy the f3 from the bottom layer in the domain (owned by the mic) to the 
! exchange buffer
!-------------------------------------------------------------------------------
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
!DIR$ IVDEP
     DO i = 1, NX
       f_buff_mic( i + NX*(k-1) ) = f3_mic( i + NXG*( 1 + NYG_MIC*k ) + now_mic )
     END DO
 END DO

!-------------------------------------------------------------------------------
! We copy the f4 from the top layer in the domain (owned by the mic) to the 
! exchange buffer
!-------------------------------------------------------------------------------
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
!DIR$ IVDEP
     DO i = 1, NX
       f_buff_mic( i + NX*(k-1) + offset ) = f4_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + now_mic )
     END DO
 END DO
!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PackMIC_F
