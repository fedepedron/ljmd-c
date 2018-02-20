#include <map>
#include <string>
#include <cstdio>
#include <sys/time.h>

extern void timer_start(std::string name);
extern void timer_pause(std::string name);
extern void timer_stop(std::string name);
extern void print_timers();
