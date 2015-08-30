        subroutine assert(isok, mesg, ival )
        implicit none
        logical isok
        character*(*) mesg
        integer ival

        if (.not.isok) then
                write(*,*) mesg,ival
                stop '** assertion error '
        endif
        return
        end
