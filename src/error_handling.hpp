#ifndef stack_traces_hpp
#define stack_traces_hpp

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <signal.h>

extern char const * program_name;

void handle_sqlite_error(void * p_arg, int error_code, const char * msg);

void register_signal_handler();
void print_stack_trace();
int  divide_by_zero();
void cause_segfault();
void stack_overflow();
void infinite_loop();
void illegal_instruction();

#endif // #define stack_traces_hpp
