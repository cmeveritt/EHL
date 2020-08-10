! Sobroutine to break if the out file is to large. If the out file is too large it indicates that something has gone wrong
    subroutine Stop_to_large_out(t)
    implicit none
    integer     out_size
    integer     t
    
    inquire(UNIT=4, SIZE=out_size)                            ! Read the size of the out file
    if( out_size .GT. 1e8)then
        WRITE(4,*)'Breaking due to to big out file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        WRITE(4,*)'At t = ', t
        OPEN(20,FILE='STOP_Too_large_out_file.DAT',STATUS='UNKNOWN')
        stop 'Out greater than 1e8' 
    endif
    
    return
    end
    