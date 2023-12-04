module timings
    implicit none
    save
    private
    public tic,toc,startClock,stopClock

    real(kind(0e0)) proc_time
    integer wall_time, rate
    
contains

    subroutine tic(startTime)
        ! tic() starts the cpu timer. The subroutines should record the processor time at execution of this command in a module variable. 
        ! Use toc()/ toc(elapsedTime) to print/get the elapsed time since the last call of tic().
        ! tic(startTime) returns the value of the processor time at execution of this command in the real(kind(0.e0)) variable startTime. 
        ! Calling this command shall not change the module variable associated with calling tic() (ie. without argument).
        real(kind(0e0)), intent(out), optional :: startTime

        if (present(startTime)) then
            call cpu_time(startTime)
        else
            call cpu_time(proc_time)
        endif
    end subroutine tic

    subroutine toc(elapsedTime, startTime)
	    ! toc() prints the elapsed cpu time in seconds since the most recent call of tic() (ie. tic called without output argument).
        ! toc(elapsedTime) returns the elapsed cpu time in seconds since the most recent call of tic() (ie. tic called without output argument) in the real(kind(0.e0)) variable elapsedTime.
        ! if toc() or toc(elapsedTime) is called without first calling tic(), a meaningless result may be returned.
        ! toc(startTime=startTime) prints the elapsed cpu time in seconds since the call of the tic command corresponding to startTime.
        ! toc(elapsedTime, startTime) returns the elapsed cpu time in seconds since the call of the tic command corresponding to startTime in the real(kind(0.e0)) variable elapsedTime.
        real(kind(0e0)), intent(out), optional :: elapsedTime
        real(kind(0e0)), intent(out), optional :: startTime
        real(kind(0e0))                        :: tmp

        if (.not. present(elapsedTime) .and. .not. present(startTime)) then
            tmp = proc_time
            call cpu_time(proc_time)
            print *, 'Elapsed CPU time : ', proc_time - tmp
        else if (.not. present(startTime)) then
            tmp = proc_time
            call cpu_time(proc_time)
            elapsedTime = proc_time - tmp
        else if (.not. present(elapsedTime)) then
            call cpu_time(proc_time)
            print *, 'Elapsed CPU time : ', proc_time - startTime
        else
            call cpu_time(proc_time)
            elapsedTime = proc_time - startTime
        endif
    end subroutine toc

    subroutine startClock(startTime)
        ! startClock() starts the wall clock timer. The subroutine should record the value of a real-time clock at the execution of this command in a module variable. Use stopClock()/stopClock(elapsedTime) to print/get the elapsed time.
        ! startClock(startTime) stores the value of a real-time clock at execution of this command in the integer(kind(0)) variable startTime. Calling this command shall not change the module variable associated with calling startClock() (ie.startClock called without argument).
        integer, intent(out), optional :: startTime

        if (present(startTime)) then
            call system_clock(startTime, rate)
        else
            call system_clock(wall_time, rate)
        endif
    end subroutine startClock
    
    subroutine stopClock(elapsedTime, startTime)
	 ! stopClock() prints the elapsed wall clock time in milliseconds since the most recent call of startClock() (ie. startClock called without output argument).
        ! stopClock(elapsedTime) returns the elapsed wall clock time in milliseconds since the most recent call of startClock() (ie. startClock called without output argument) in the integer(kind(0)) variable  elapsedTime.
        ! if stopClock() or stopClock(elapsedTime) is called without first calling startClock() a meaningless result may be returned.
        ! stopClock(startTime=startTime) prints the elapsed wall clock time in milliseconds since the call of the startClock command corresponding to startTime.
        ! stopClock(elapsedTime, startTime) returns the elapsed wall clock time in milliseconds since the call of the startClock command corresponding to startTime in the integer(kind(0)) variable elapsedTime.
        real(kind(0e0)), intent(out), optional :: elapsedTime
        integer, intent(out), optional :: startTime
        integer                        :: tmp

        if (.not. present(elapsedTime) .and. .not. present(startTime)) then
            tmp = wall_time
            call system_clock(wall_time)
            print *, 'Elapsed real time : ', real(wall_time - tmp)/real(rate)*10**3
        else if (.not. present(startTime)) then
            tmp = wall_time
            call system_clock(wall_time)
            elapsedTime =  real(wall_time - tmp)/real(rate)*10**3
        else if (.not. present(elapsedTime)) then
            call system_clock(wall_time)
            print *, 'Elapsed real time : ', real(wall_time - startTime)/real(rate)*10**3
        else
            call system_clock(wall_time)
            elapsedTime = real(wall_time - startTime)/real(rate)*10**3
        endif
    end subroutine stopClock
end
