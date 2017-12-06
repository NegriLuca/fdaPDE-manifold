#include <iostream>
#include <iomanip>
#include <ctime>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h> 

class timer {
public:
  void start() {
    clock_gettime(CLOCK_REALTIME, &begin);
  }

  timespec stop() {
    timespec end;
    clock_gettime(CLOCK_REALTIME, &end);
    
    timespec difference = diff(begin, end);
    Rprintf("It took %u.%09us\n", difference.tv_sec, difference.tv_nsec);
//    std::cout << "It took " << difference.tv_sec << "."
//              << std::setfill('0') << std::setw(9) << difference.tv_nsec
//              << "s" << std::endl;
    return difference;
  }
private:
  timespec diff(timespec start, timespec end)
  {
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
      temp.tv_sec = end.tv_sec-start.tv_sec-1;
      temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
      temp.tv_sec = end.tv_sec-start.tv_sec;
      temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
  }

private:
  timespec begin;
};
