!-------------------------------------------------------------------------------
! Subroutine : MemAlloc
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
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

 SUBROUTINE MemAlloc(FLAG)

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi, lapPhi, gradPhiX, gradPhiY, gradPhiZ
 USE LBMParams
 USE MPIParams
 IMPLICIT NONE

 !   Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN
! Near neighbors and order parameter arrays
 ALLOCATE( phi( 0:NG ) )
! ALLOCATE( nb(xmin:xmax,ymin:ymax,zmin:zmax,1:6) )
 ALLOCATE( lapPhi( 0:NG ) )
 ALLOCATE( gradPhiX( 0:NG ) )
 ALLOCATE( gradPhiY( 0:NG ) )
 ALLOCATE( gradPhiZ( 0:NG ) )

! Distribution function arrays
 ALLOCATE( f0( 0:NG ) )
 ALLOCATE( f1( 0:2*NG ) )
 ALLOCATE( f2( 0:2*NG ) )
 ALLOCATE( f3( 0:2*NG ) )
 ALLOCATE( f4( 0:2*NG ) )
 ALLOCATE( f5( 0:2*NG ) )
 ALLOCATE( f6( 0:2*NG ) )

 ALLOCATE( g0( 0:NG ) )
 ALLOCATE( g1( 0:2*NG ) )
 ALLOCATE( g2( 0:2*NG ) )
 ALLOCATE( g3( 0:2*NG ) )
 ALLOCATE( g4( 0:2*NG ) )
 ALLOCATE( g5( 0:2*NG ) )
 ALLOCATE( g6( 0:2*NG ) )
 ALLOCATE( g7( 0:2*NG ) )
 ALLOCATE( g8( 0:2*NG ) )
 ALLOCATE( g9( 0:2*NG ) )
 ALLOCATE( g10( 0:2*NG ) )
 ALLOCATE( g11( 0:2*NG ) )
 ALLOCATE( g12( 0:2*NG ) )
 ALLOCATE( g13( 0:2*NG ) )
 ALLOCATE( g14( 0:2*NG ) )
 ALLOCATE( g15( 0:2*NG ) )
 ALLOCATE( g16( 0:2*NG ) )
 ALLOCATE( g17( 0:2*NG ) )
 ALLOCATE( g18( 0:2*NG ) )

! Allocation for data exchange in the x direction
   xsize  = NY*NZ
   xsize6 = NY*NZ*6
   
   ALLOCATE( f_west_snd(1:xsize) )
   ALLOCATE( f_east_snd(1:xsize) )
   ALLOCATE( f_west_rcv(1:xsize) )
   ALLOCATE( f_east_rcv(1:xsize) )

   ALLOCATE( west_snd(1:xsize6) )
   ALLOCATE( east_snd(1:xsize6) )
   ALLOCATE( west_rcv(1:xsize6) )
   ALLOCATE( east_rcv(1:xsize6) )

! Allocation for data exchange in the y direction
   ysize  = NX*NZ
   ysize6 = NX*NZ*6

   ALLOCATE( f_south_snd(1:ysize) )
   ALLOCATE( f_north_snd(1:ysize) )
   ALLOCATE( f_south_rcv(1:ysize) )
   ALLOCATE( f_north_rcv(1:ysize) )

   ALLOCATE( south_snd(1:ysize6) )
   ALLOCATE( north_snd(1:ysize6) )
   ALLOCATE( south_rcv(1:ysize6) )
   ALLOCATE( north_rcv(1:ysize6) )

! Allocation for data exchanged in the y direction
   zsize  = NX*NY
   zsize6 = NX*NY*6

   ALLOCATE( f_bot_snd(1:zsize) )
   ALLOCATE( f_bot_rcv(1:zsize) )
   ALLOCATE( f_top_snd(1:zsize) )
   ALLOCATE( f_top_rcv(1:zsize) )

   ALLOCATE( bot_snd(1:zsize6) )
   ALLOCATE( bot_rcv(1:zsize6) )
   ALLOCATE( top_snd(1:zsize6) )
   ALLOCATE( top_rcv(1:zsize6) )

! Allocation for data exchange in edges along the x direction
   xedge2 = NX*2
   
   ALLOCATE( tn_snd(1:xedge2) )
   ALLOCATE( ts_snd(1:xedge2) )
   ALLOCATE( bn_snd(1:xedge2) )
   ALLOCATE( bs_snd(1:xedge2) )
   ALLOCATE( tn_rcv(1:xedge2) )
   ALLOCATE( ts_rcv(1:xedge2) )
   ALLOCATE( bn_rcv(1:xedge2) )
   ALLOCATE( bs_rcv(1:xedge2) )

! Allocation for data exchange in edges along the y direction
   yedge2 = NY*2

   ALLOCATE( te_snd(1:yedge2) )
   ALLOCATE( tw_snd(1:yedge2) )
   ALLOCATE( be_snd(1:yedge2) )
   ALLOCATE( bw_snd(1:yedge2) )
   ALLOCATE( te_rcv(1:yedge2) )
   ALLOCATE( tw_rcv(1:yedge2) )
   ALLOCATE( be_rcv(1:yedge2) )
   ALLOCATE( bw_rcv(1:yedge2) )

! Allocation for data exchange in edges along the z direction
   zedge2 = NZ*2

   ALLOCATE( ne_snd(1:zedge2) )
   ALLOCATE( nw_snd(1:zedge2) )
   ALLOCATE( se_snd(1:zedge2) )
   ALLOCATE( sw_snd(1:zedge2) )
   ALLOCATE( ne_rcv(1:zedge2) )
   ALLOCATE( nw_rcv(1:zedge2) )
   ALLOCATE( se_rcv(1:zedge2) )
   ALLOCATE( sw_rcv(1:zedge2) )

 ELSE

   DEALLOCATE( phi )
   DEALLOCATE( lapPhi )
   DEALLOCATE( gradPhiX )
   DEALLOCATE( gradPhiY )
   DEALLOCATE( gradPhiZ )
   DEALLOCATE( f0 )
   DEALLOCATE( f1 )
   DEALLOCATE( f2 )
   DEALLOCATE( f3 )
   DEALLOCATE( f4 )
   DEALLOCATE( f5 )
   DEALLOCATE( f6 )
   DEALLOCATE( g0 )
   DEALLOCATE( g1 )
   DEALLOCATE( g2 )
   DEALLOCATE( g3 )
   DEALLOCATE( g4 )
   DEALLOCATE( g5 )
   DEALLOCATE( g6 )
   DEALLOCATE( g7 )
   DEALLOCATE( g8 )
   DEALLOCATE( g9 )
   DEALLOCATE( g10 )
   DEALLOCATE( g11 )
   DEALLOCATE( g12 )
   DEALLOCATE( g13 )
   DEALLOCATE( g14 )
   DEALLOCATE( g15 )
   DEALLOCATE( g16 )
   DEALLOCATE( g17 )
   DEALLOCATE( g18 )
   DEALLOCATE( f_west_snd )
   DEALLOCATE( f_west_rcv )
   DEALLOCATE( f_east_snd )
   DEALLOCATE( f_east_rcv )
   DEALLOCATE( west_snd )
   DEALLOCATE( east_snd )
   DEALLOCATE( west_rcv )
   DEALLOCATE( east_rcv )
   DEALLOCATE( f_south_snd )
   DEALLOCATE( f_south_rcv )
   DEALLOCATE( f_north_snd )
   DEALLOCATE( f_north_rcv )
   DEALLOCATE( south_snd )
   DEALLOCATE( north_snd )
   DEALLOCATE( south_rcv )
   DEALLOCATE( north_rcv )
   DEALLOCATE( f_bot_snd )
   DEALLOCATE( f_bot_rcv )
   DEALLOCATE( f_top_snd )
   DEALLOCATE( f_top_rcv )
   DEALLOCATE( bot_snd )
   DEALLOCATE( bot_rcv )
   DEALLOCATE( top_snd )
   DEALLOCATE( top_rcv )
   DEALLOCATE( tn_snd )
   DEALLOCATE( ts_snd )
   DEALLOCATE( te_snd )
   DEALLOCATE( tw_snd )
   DEALLOCATE( bn_snd )
   DEALLOCATE( bs_snd )
   DEALLOCATE( be_snd )
   DEALLOCATE( bw_snd )
   DEALLOCATE( ne_snd )
   DEALLOCATE( nw_snd )
   DEALLOCATE( se_snd )
   DEALLOCATE( sw_snd ) 
   DEALLOCATE( tn_rcv )
   DEALLOCATE( ts_rcv )
   DEALLOCATE( te_rcv )
   DEALLOCATE( tw_rcv )
   DEALLOCATE( bn_rcv )
   DEALLOCATE( bs_rcv )
   DEALLOCATE( be_rcv )
   DEALLOCATE( bw_rcv )
   DEALLOCATE( ne_rcv )
   DEALLOCATE( nw_rcv )
   DEALLOCATE( se_rcv )
   DEALLOCATE( sw_rcv )

 END IF

 RETURN
 END SUBROUTINE MemAlloc
