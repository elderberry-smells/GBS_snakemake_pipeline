#ifndef __TIMER_HH
#define __TIMER_HH
//
//  File      : timer.hh
//

#include <sys/time.h>
#include <sys/resource.h>

#ifdef SUNOS
#include <sys/lwp.h>
#else
//#define _GNU_SOURCE
#endif
    
#include <unistd.h>

#include <stdio.h>

// type for different timer-measurment
typedef enum { REAL_TIME, CPU_TIME, USER_TIME, SYSTEM_TIME, THREAD_CPU }  timetype_t;

//
// the timer class
//
class TTimer 
{
public:
    
protected:
    // start/stop times
    double          _start, _stop;

    // temp structs
    struct timeval  _timeval_data;
#ifdef SUNOS
    struct lwpinfo  _lwpinfo_data;
#else
    struct rusage   _rusage_data;
#endif
    
    // what kind of time we should stop
    timetype_t      _type;
    
public:
    TTimer ( timetype_t ttype = CPU_TIME ) : _start(0), _stop(0), _type(ttype) {}
    
	// sets first/second timer (cpu- or real-time by _real-field)
    TTimer & start    () { _stop = _start = system_time(); return *this; }
    TTimer & stop     () { _stop =          system_time(); return *this; }

	// returns time between start and end in seconds
	double diff () const { return _stop - _start; }

    // get time of system (usertime or real)
    double system_time ();

    // set type of time to measure and return *this (for output)
    TTimer & real () { _type = REAL_TIME; return *this; }
    TTimer & cpu  () { _type = CPU_TIME;  return *this; }
    TTimer & user () { _type = USER_TIME; return *this; }
    TTimer & sys  () { _type = SYSTEM_TIME;  return *this; }
    TTimer & thread  () { _type = THREAD_CPU;  return *this; }
    
};

#endif  // __TIMER_HH
