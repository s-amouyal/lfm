#ifndef TIMING_H
#define TIMING_H

#include <unistd.h>

/* Measuring CPU time */

//***************************************************************************************************
class cpu_timing{

	public:
		uint64_t second;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline void calibrate_rdtsc(){

    		uint64_t t0, t1;

    		t0 = get_rdtsc();
    		sleep(1);
    		t1 = get_rdtsc();

    		second = t1 - t0;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline uint64_t get_rdtsc() {

    		//uint32_t low, high;
    		//asm volatile ("rdtsc" : "=a" (low), "=d" (high));
    		//return ((uint64_t)high << 32) | (uint64_t)low;
    		return 1;
		}
};

#endif
