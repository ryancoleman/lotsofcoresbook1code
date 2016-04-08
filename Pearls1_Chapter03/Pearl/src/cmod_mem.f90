!Copyright (c) 2014, Per Berg and Jacob Weismann Poulsen, DMI
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met: 
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer. 
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution. 
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!The views and conclusions contained in the software and documentation are those
!of the authors and should not be interpreted as representing official policies, 
!either expressed or implied, of the FreeBSD Project.

module cmod_mem
  !=============================================================================
  !
  ! Module for memory allocation in cmod, especially for nested variables.
  ! Assuming shape of arrays is important.
  ! Names of predefined type defs follow the notation:
  !   cmTD   
  ! where cm is short for "Cmod Memory", T is a type identifier
  !   T  = i for integer(4)
  !   T  = r for real(8) 
  !   T  = l for logical 
  ! and D identifies the dimension (shape) of the array
  !   D  = 1 for 1D arrays
  !   D  = 2 for 2D arrays
  !   D  = 3 for 3D arrays
  ! 
  ! Example of usage: 
  ! Define a 3D integer(4) index array for each of N model domains 
  ! when N is specified at run-time, ie it is not known during compilation.
  ! Use the index arrays in subsequent calls of a subroutine requiring four
  ! arguments: the index array and its three dimensions. 
  ! Free memory when done.
  !   
  !   use cmod_mem
  !   interface
  !     subroutine sub( ix, n1, n2, n3 )
  !       integer(4), intent(in)    :: n1, n2, n3
  !       integer(4), intent(inout) :: ix(:,:,:)
  !     end subroutine sub
  !   end interface
  !   type (cmi3), allocatable :: idx(:)
  !   integer(4) :: N, ia
  !   integer(4), allocatable :: dims(:,:)
  !      ... statements getting N
  !   allocate( dims(N,3) )
  !      ... statements getting dimensions dims(:,1:3) in each area
  !   allocate( idx(N) )
  !   do ia=1,N
  !     allocate( idx(ia)%p( 1:dims(ia,1), 1:dims(ia,2), 1:dims(ia,3) ) )
  !   enddo
  !   do ia=1,N
  !     call sub( idx(ia), dims(ia,1), dims(ia,2), dims(ia,3) )
  !   enddo
  !   deallocate( idx )
  !=============================================================================
  implicit none
  private

  type cmi1
    integer(4), allocatable :: p(:)  
  end type

  type cmi2
    integer(4), allocatable :: p(:,:)  
  end type

  type cmi3
    integer(4), allocatable :: p(:,:,:)  
  end type

  type cmr1
    real(8), allocatable :: p(:)  
  end type

  type cmr2
    real(8), allocatable :: p(:,:)  
  end type

  type cmr3
    real(8), allocatable :: p(:,:,:)  
  end type

  type cml1
    logical, allocatable :: p(:)  
  end type

  type cml2
    logical, allocatable :: p(:,:)  
  end type

  type cml3
    logical, allocatable :: p(:,:,:)  
  end type

  public  :: cmi1, cmi2, cmi3, cmr1, cmr2, cmr3, cml1, cml2, cml3

end module cmod_mem
