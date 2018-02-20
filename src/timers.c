#include <map>
#include <string>
#include <cstdio>
#include <sys/time.h>
#include "timers.h"

class timer
{
public:
  double t_start;
  double t_total;
  double t_end;
  bool   exists;
};

std::map < std::string, timer > alltimers;

double mtime()
{
  struct timeval tmp;
  gettimeofday( &tmp, (struct timezone *)0 );
  return tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
}

extern void timer_start(std::string name)
{
  name.append("\0");
  if (alltimers.find(name) == alltimers.end())
  {
    timer new_timer;
    new_timer.t_total = 0.0;
    alltimers[name]   = new_timer;
  }
  alltimers[name].t_start = mtime();
}

extern void timer_pause(std::string name)
{
  name.append("\0");
  alltimers[name].t_total += (mtime() - alltimers[name].t_start);
  alltimers[name].t_start  = mtime();
}

extern void timer_stop(std::string name)
{
  name.append("\0");
  alltimers[name].t_total = mtime() - alltimers[name].t_start;
}

extern void print_timers()
{
  double total_time = alltimers["Total"].t_total;
  for (std::map<std::string, timer>::iterator it = alltimers.begin();
       it != alltimers.end(); ++it)
  {
    double percent = 100.0 * (it->second.t_total / total_time);
    printf("Timer %-20s: %12.6f seconds. (%5.2f%% of total)\n",
           it->first.c_str(), it->second.t_total, percent);
  }
}
