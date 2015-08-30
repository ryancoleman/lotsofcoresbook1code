!-------------------------------------------------------------------------------
! Subroutine : PackMIC
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

 SUBROUTINE PackMIC
!DIR$ ATTRIBUTES OFFLOAD:MIC :: PackMIC 

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, id_buf, k, id_mic, offset

 offset = NX*NZ

!$OMP PARALLEL

!-------------------------------------------------------------------------------
! Because the MIC computes all point in the [1:NY_MIC] range and the CPU in the 
! [NY_MIC+1:NY] range we need to exchange all the upward pointing values on the
! topmost ghost layer of the MIC (NY_MIC+1) and the downward pointing values on 
! the bottom ghost layer of the MIC (0).
!-------------------------------------------------------------------------------
!$OMP DO PRIVATE(i,id_buf,id_mic)
 DO k = 1, NZ

!DIR$ IVDEP
     DO i = 1, NX
       id_buf = i + NX*(k-1)
       id_mic = i + NXG*( NYG_MIC*k ) + nxt_mic

       buff_mic(id_buf           )  = f4_mic(id_mic)
       buff_mic(id_buf +   offset)  = g4_mic(id_mic)
       buff_mic(id_buf + 2*offset)  = g8_mic(id_mic)
       buff_mic(id_buf + 3*offset)  = g9_mic(id_mic)
       buff_mic(id_buf + 4*offset) = g16_mic(id_mic)
       buff_mic(id_buf + 5*offset) = g18_mic(id_mic)

     END DO

!DIR$ IVDEP
     DO i = 1, NX
       id_buf = i + NX*(k-1)
       id_mic = i + NXG*( NY_MIC+1 + NYG_MIC*k ) + nxt_mic

       buff_mic(id_buf +  6*offset) = f3_mic(id_mic)
       buff_mic(id_buf +  7*offset) = g3_mic(id_mic)
       buff_mic(id_buf +  8*offset) = g7_mic(id_mic)
       buff_mic(id_buf +  9*offset) = g10_mic(id_mic)
       buff_mic(id_buf + 10*offset) = g15_mic(id_mic)
       buff_mic(id_buf + 11*offset) = g17_mic(id_mic)

     END DO

 END DO

!-------------------------------------------------------------------
! And, of course, we also need the edge exchanges along X and Z ...
!-------------------------------------------------------------------

!$OMP DO PRIVATE(i)
 DO i = 1, NX
   buff_mic_edge_x(i       ) = g15_mic( i + NXG*( NY_MIC+1 + NYG_MIC*(NZ+1) ) + nxt_mic )
   buff_mic_edge_x(i +   NX) = g16_mic( i                                     + nxt_mic )
   buff_mic_edge_x(i + 2*NX) = g17_mic( i + NXG*( NY_MIC+1                  ) + nxt_mic )
   buff_mic_edge_x(i + 3*NX) = g18_mic( i + NXG*(            NYG_MIC*(NZ+1) ) + nxt_mic )
 END DO

 !$OMP DO PRIVATE(k)
 DO k = 1, NZ
   buff_mic_edge_z(k       ) = g7_mic( NX+1 + NXG*( NY_MIC+1 + NYG_MIC*k ) + nxt_mic )
   buff_mic_edge_z(k +   NZ) = g8_mic(        NXG*(            NYG_MIC*k ) + nxt_mic )
   buff_mic_edge_z(k + 2*NZ) = g9_mic( NX+1 + NXG*(            NYG_MIC*k ) + nxt_mic )
   buff_mic_edge_z(k + 3*NZ) = g10_mic(       NXG*( NY_MIC+1 + NYG_MIC*k ) + nxt_mic )
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PackMIC
