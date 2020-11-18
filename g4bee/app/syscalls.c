/**
 * @file syscalls.c
 * @author Avik De (avikde@gmail.com)
 * @brief Stuff needed to use C++
 * @version 0.1
 * @date 2019-06-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "main.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

char *__env[1] = { 0 };
char **environ = __env;

/*----------------------------------------------------------------------------
 *        Exported variables
 *----------------------------------------------------------------------------*/

// #undef errno
// extern int errno ;
// extern int  _end ;
/*----------------------------------------------------------------------------
 *        Exported functions
 *----------------------------------------------------------------------------*/
extern void _exit(int status);
extern void _kill(int pid, int sig);
extern int _getpid(void);

/*
 link
 Establish a new name for an existing file. Minimal implementation:
 */
extern int link(const char __attribute__((unused)) * cOld, const char __attribute__((unused)) * cNew)
{
	return 0;
}

/*
 sbrk
 Increase program data space.
 Malloc and related functions depend on this
 */
extern caddr_t _sbrk(int incr) {
	extern char _end; // Defined by the linker
	static char *heap_end;
	char *prev_heap_end;
	if (heap_end == 0) {
		heap_end = &_end;
	}
	prev_heap_end = heap_end;
	char * stack = (char*) __get_MSP();
	if (heap_end + incr > stack) {
		// _write(1, (char *)"Heap and stack collision\n", 25);
		return (caddr_t) -1;
	}
	heap_end += incr;
	return (caddr_t) prev_heap_end;
}

extern int _close(int __attribute__((unused)) file)
{
	return 0;
}

/*
 fstat
 Status of an open file. For consistency with other minimal implementations in these examples,
 all files are regarded as character special devices.
 The `sys/stat.h' header file required is distributed in the `include' subdirectory for this C library.
 */
extern int _fstat(int __attribute__((unused)) file, struct stat __attribute__((unused)) * st)
{
	// st->st_mode = S_IFCHR ;

	return 0;
}

/*
 isatty
 Query whether output stream is a terminal. For consistency with the other minimal implementations,
 */
extern int _isatty(int __attribute__((unused)) file)
{
	return 1;
}

/*
 lseek
 Set position in a file. Minimal implementation:
 */
extern int _lseek(int __attribute__((unused)) file, int __attribute__((unused)) ptr, int __attribute__((unused)) dir)
{
	return 0;
}

extern void _exit(int __attribute__((unused)) status)
{
	for (;;)
		;
}

/*
 kill
 Send a signal. Minimal implementation:
 */
extern void _kill(int __attribute__((unused)) pid, int __attribute__((unused)) sig)
{
	return;
}

/*
 getpid
 Process-ID; this is sometimes used to generate strings unlikely to conflict with other processes. Minimal implementation, for a system without processes:
 */
extern int _getpid(void) {
	return -1;
}

#ifdef USE_FULL_ASSERT
void assert_failed(uint8_t* file, uint32_t line) {
	trace_printf("%s:%d\r\n", file, line);
	while (1)
		;
}
#endif

#pragma GCC diagnostic pop
