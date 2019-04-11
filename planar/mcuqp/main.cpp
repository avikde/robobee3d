#include "stdio.h"
#include <string.h>
#include <sys/types.h> // for useconds_t
#include <unistd.h>
#include <Arduino.h>
#include <variant.h>

#include "workspace.h"
#include "osqp.h"

const uint8_t LEDB = PD8, LEDG = PD10, LEDR = PD15;

// libc_init_array
extern void (*__preinit_array_start[])(void) WEAK;
extern void (*__preinit_array_end[])(void) WEAK;
extern void (*__init_array_start[])(void) WEAK;
extern void (*__init_array_end[])(void) WEAK;
extern void (*__fini_array_start[])(void) WEAK;
extern void (*__fini_array_end[])(void) WEAK;

extern "C" int _write(int file, char *data, int len)
{
	// The special files number 0, 1, 2 are for stdin, stdout, stderr (typically)
	switch (file)
	{
	case STDIN_FILENO:
	case STDOUT_FILENO:
	case STDERR_FILENO:
		return Serial1.write((const uint8_t *)data, len);
		// TODO figure out how to use this for log file writing (user params) / open / close etc.
	default:
		return -1;
	}
}
/*
 read
 Read a character to a file. `libc' subroutines will use this system routine for input from all files, including stdin
 Returns -1 on error or blocks until the number of characters have been read.
 */
extern "C" int _read(int file, char *ptr, int len)
{
	static int avail;
	// The special files number 0, 1, 2 are for stdin, stdout, stderr (typically)
	switch (file)
	{
	case STDIN_FILENO:
	case STDOUT_FILENO:
	case STDERR_FILENO:
		avail = Serial1.available();
		// must fit in output buffer
		if (avail > len)
			avail = len;
		// Note this is a primitive Arduino-style read. This could be replaced with readLatestDMA
		for (int i = 0; i < avail; ++i)
		{
			ptr[i] = (char)Serial1.read();
		}
		return avail;
	default:
		return -1;
	}
}

extern "C" int usleep(useconds_t us)
{
	//	FIXME after scheduler has been started, use vTaskDelay
	// if (rtosSchedulerStarted)
	// {
	// }
	// else
	// {
	if (us < 1000)
		delayMicroseconds(us);
	else
	{
		//    trace_printf("aa = %u", aa);
		delayMicroseconds(us % (uint32_t)1000);
		delay(us / 1000);
	}
	// }
	return 0;
}

// ---

int main(int argc, char **argv)
{
	// This is basically libc_init_array -- handles global constructors
	unsigned int count;
	unsigned int i;
	count = __preinit_array_end - __preinit_array_start;
	for (i = 0; i < count; i++)
		__preinit_array_start[i]();
	count = __init_array_end - __init_array_start;
	for (i = 0; i < count; i++)
		__init_array_start[i]();
	count = __fini_array_end - __fini_array_start;
	for (i = 0; i < count; i++)
		__fini_array_start[i]();

	// FIXME needs to move up
	// /* volatile uint32_t regVal =  */ rebootInit();

	//4 bits for preemp priority 0 bit for sub priority
	NVIC_PriorityGroupConfig(NVIC_PriorityGroup_4);
	//  koduino init
	// From variant: init peripherals
	variantInit();
	// Timing: micros(), millis(), delay() etc.
	systemClockInit();
	// External interrupts
	interrupts();

	// // Solve Problem
	// osqp_solve(&workspace);

	// // Print status
	// printf("Status:                %s\n", (&workspace)->info->status);
	// printf("Number of iterations:  %d\n", (int)((&workspace)->info->iter));
	// printf("Objective value:       %.4e\n", (&workspace)->info->obj_val);
	// printf("Primal residual:       %.4e\n", (&workspace)->info->pri_res);
	// printf("Dual residual:         %.4e\n", (&workspace)->info->dua_res);

	pinMode(LEDB, OUTPUT);
	volatile int test = 0;

	while (1)
	{
		delay(1000);
		digitalWrite(LEDB, TOGGLE);
		test++;
	}
	return 0;
}