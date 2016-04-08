!-------------------------------------------------------------------------------
! Subroutine : PackCPU
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Copy exchange values to CPU buffer for transfer to MIC.
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

 SUBROUTINE PackCPU

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, id_buf, k, id_cpu, offset

 offset = NX*NZ

!$OMP PARALLEL

!-------------------------------------------------------------------------------
! Because the MIC computes all point in the [1:NY_MIC] range and the CPU in the 
! [NY_MIC+1:NY] range we need to exchange all the downward pointing values on 
! the bottommost ghost layer of the CPU (NY_MIC)
!-------------------------------------------------------------------------------
!$OMP DO PRIVATE(i,id_buf,id_cpu)
 DO k = 1, NZ
!DIR$ IVDEP
     DO i = 1, NX
       id_buf = i + NX*(k-1)
       id_cpu = i + NXG*( NY_MIC + NYG*k ) + nxt

       buff_cpu(id_buf           ) = f4(id_cpu)
       buff_cpu(id_buf +   offset) = g4(id_cpu)
       buff_cpu(id_buf + 2*offset) = g8(id_cpu)
       buff_cpu(id_buf + 3*offset) = g9(id_cpu)
       buff_cpu(id_buf + 4*offset) = g16(id_cpu)
       buff_cpu(id_buf + 5*offset) = g18(id_cpu)
     END DO

!-------------------------------------------------------------------------------
! In the multiphase case, because we are using periodic boundary conditions,
! we also have to exchange the outward pointing top ghosts in the CPU (NY+1)
!-------------------------------------------------------------------------------
!DIR$ IVDEP
     DO i = 1, NX
       id_buf = i + NX*(k-1)
       id_cpu = i + NXG*( NY+1 + NYG*k ) + nxt

       buff_cpu(id_buf + 6*offset)  = f3(id_cpu)
       buff_cpu(id_buf + 7*offset)  = g3(id_cpu)
       buff_cpu(id_buf + 8*offset)  = g7(id_cpu)
       buff_cpu(id_buf + 9*offset)  = g10(id_cpu)
       buff_cpu(id_buf + 10*offset) = g15(id_cpu)
       buff_cpu(id_buf + 11*offset) = g17(id_cpu)
     END DO

 END DO


!-------------------------------------------------------------------
! And, of course, we also need the edge exchanges along X and Z ...
!-------------------------------------------------------------------

!$OMP DO PRIVATE(i)
 DO i = 1, NX
   buff_cpu_edge_x(i       ) = g15( i + NXG*( NY+1   + NYG*(NZ+1) ) + nxt )
   buff_cpu_edge_x(i +   NX) = g16( i + NXG*( NY_MIC              ) + nxt )
   buff_cpu_edge_x(i + 2*NX) = g17( i + NXG*( NY+1                ) + nxt )
   buff_cpu_edge_x(i + 3*NX) = g18( i + NXG*( NY_MIC + NYG*(NZ+1) ) + nxt )
 END DO

 !$OMP DO PRIVATE(k)
 DO k = 1, NZ
   buff_cpu_edge_z(k       ) = g7( NX+1 + NXG*( NY+1   + NYG*k ) + nxt )
   buff_cpu_edge_z(k +   NZ) = g8(        NXG*( NY_MIC + NYG*k ) + nxt )
   buff_cpu_edge_z(k + 2*NZ) = g9( NX+1 + NXG*( NY_MIC + NYG*k ) + nxt )
   buff_cpu_edge_z(k + 3*NZ) = g10(       NXG*( NY+1   + NYG*k ) + nxt )
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PackCPU
