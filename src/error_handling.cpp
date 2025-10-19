/***
    Parts adapted from https://gist.github.com/jvranish/4441299
    by Job Vranish, downloaded on 2017-10-07.
    Original code is in the public domain. 
***/

/* compile with:
 on linux:   gcc -g stack_traces.c
 on OS X:    gcc -g -fno-pie stack_traces.c
 on windows: gcc -g stack_traces.c -limagehlp
 */

#include "error_handling.hpp"

#include <signal.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdbool.h>
#include <errno.h>
#include <err.h>
#include <execinfo.h>

#ifdef __APPLE__
#include <dlfcn.h>
#endif

char const * program_name;


void handle_sqlite_error(void * p_arg, int error_code, const char * msg) {
    fprintf(stderr, "SQLite error occurred! Aborting.\n");
    fprintf(stderr, "(%d) %s\n", error_code, msg);
    print_stack_trace();
    exit(1);
}

/* Resolve symbol name and source location given the path to the executable 
 and an address */
int addr2line(char const * const program_name, void const * const addr)
{
    char addr2line_cmd[512] = {0};
    
    /* have addr2line map the address to the relent line in the code */
#ifdef __APPLE__
    /* apple does things differently... */
    
    Dl_info info;
    if(dladdr(addr, &info)) {
        sprintf(addr2line_cmd,"atos -o %.256s -l %p %p", info.dli_fname, info.dli_fbase, addr);
        //printf("%s\n", addr2line_cmd);
    }
    else {
        
    }
     
#else
    sprintf(addr2line_cmd,"addr2line -f -p -e %.256s %p", program_name, addr); 
#endif
    
    return system(addr2line_cmd);
}



#define MAX_STACK_FRAMES 64
static void *stack_traces[MAX_STACK_FRAMES];
void print_stack_trace()
{
    int i, trace_size = 0;
    char **messages = (char **)NULL;
    
    trace_size = backtrace(stack_traces, MAX_STACK_FRAMES);
    messages = backtrace_symbols(stack_traces, trace_size);
    
    /* skip the first couple stack frames (as they are this function and
     our handler) and also skip the last frame as it's (always?) junk. */
    // for (i = 3; i < (trace_size - 1); ++i)
    for (i = 0; i < trace_size; ++i) // we'll use this for now so you can see what's going on
    {
        //printf("line %d\n", i);
        //printf("%p\n", stack_traces[i]);
        //printf("%s\n", messages[i]);
        if(addr2line(program_name, stack_traces[i]) != 0)
        {
            printf("  error determining line # for: %s\n", messages[i]);
        }
        
    }
    if (messages) { free(messages); } 
}

void signal_handler(int sig, siginfo_t *siginfo, void *context)
{
    (void)context;
    switch(sig)
    {
        case SIGSEGV:
            fputs("Caught SIGSEGV: Segmentation Fault\n", stderr);
            break;
        case SIGINT:
            fputs("Caught SIGINT: Interactive attention signal, (usually ctrl+c)\n", stderr);
            break;
        case SIGFPE:
            switch(siginfo->si_code)
        {
            case FPE_INTDIV:
                fputs("Caught SIGFPE: (integer divide by zero)\n", stderr);
                break;
            case FPE_INTOVF:
                fputs("Caught SIGFPE: (integer overflow)\n", stderr);
                break;
            case FPE_FLTDIV:
                fputs("Caught SIGFPE: (floating-point divide by zero)\n", stderr);
                break;
            case FPE_FLTOVF:
                fputs("Caught SIGFPE: (floating-point overflow)\n", stderr);
                break;
            case FPE_FLTUND:
                fputs("Caught SIGFPE: (floating-point underflow)\n", stderr);
                break;
            case FPE_FLTRES:
                fputs("Caught SIGFPE: (floating-point inexact result)\n", stderr);
                break;
            case FPE_FLTINV:
                fputs("Caught SIGFPE: (floating-point invalid operation)\n", stderr);
                break;
            case FPE_FLTSUB:
                fputs("Caught SIGFPE: (subscript out of range)\n", stderr);
                break;
            default:
                fputs("Caught SIGFPE: Arithmetic Exception\n", stderr);
                break;
        }
        case SIGILL:
            switch(siginfo->si_code)
        {
            case ILL_ILLOPC:
                fputs("Caught SIGILL: (illegal opcode)\n", stderr);
                break;
            case ILL_ILLOPN:
                fputs("Caught SIGILL: (illegal operand)\n", stderr);
                break;
            case ILL_ILLADR:
                fputs("Caught SIGILL: (illegal addressing mode)\n", stderr);
                break;
            case ILL_ILLTRP:
                fputs("Caught SIGILL: (illegal trap)\n", stderr);
                break;
            case ILL_PRVOPC:
                fputs("Caught SIGILL: (privileged opcode)\n", stderr);
                break;
            case ILL_PRVREG:
                fputs("Caught SIGILL: (privileged register)\n", stderr);
                break;
            case ILL_COPROC:
                fputs("Caught SIGILL: (coprocessor error)\n", stderr);
                break;
            case ILL_BADSTK:
                fputs("Caught SIGILL: (internal stack error)\n", stderr);
                break;
            default:
                fputs("Caught SIGILL: Illegal Instruction\n", stderr);
                break;
        }
            break;
        case SIGTERM:
            fputs("Caught SIGTERM: a termination request was sent to the program\n", stderr);
            break;
        case SIGABRT:
            fputs("Caught SIGABRT: usually caused by an abort() or assert()\n", stderr);
            break;
        default:
            break;
    }
    print_stack_trace();
    _Exit(1);
}

/* static uint8_t alternate_stack[SIGSTKSZ]; */
void register_signal_handler()
{
    /* setup alternate stack */
    {
        stack_t ss = {};
        /* malloc is usually used here, I'm not 100% sure my static allocation
         is valid but it seems to work just fine. */
        /* ss.ss_sp = (void*)alternate_stack;
        ss.ss_size = SIGSTKSZ;
        ss.ss_flags = 0; */
        ss.ss_sp = malloc(SIGSTKSZ);
        ss.ss_size = SIGSTKSZ;
        ss.ss_flags = 0;
        
        if (sigaltstack(&ss, NULL) != 0) { err(1, "sigaltstack"); }
    }
    
    /* register our signal handlers */
    {
        struct sigaction sig_action = {};
        sig_action.sa_sigaction = signal_handler;
        sigemptyset(&sig_action.sa_mask);
        
#ifdef __APPLE__
        /* for some reason we backtrace() doesn't work on osx
         when we use an alternate stack */
        sig_action.sa_flags = SA_SIGINFO;
#else
        sig_action.sa_flags = SA_SIGINFO | SA_ONSTACK;
#endif
        
        if (sigaction(SIGSEGV, &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGFPE,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGINT,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGILL,  &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGTERM, &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGABRT, &sig_action, NULL) != 0) { err(1, "sigaction"); }
    }
}

int  divide_by_zero();
void cause_segfault();
void stack_overflow();
void infinite_loop();
void illegal_instruction();

int divide_by_zero()
{
    int a = 1;
    int b = 0; 
    return a / b;
}

void cause_segfault()
{
    int * p = (int*)0x12345678;
    *p = 0;
}

void stack_overflow()
{
    int foo[1000];
    (void)foo;
    stack_overflow();
}

/* break out with ctrl+c to test SIGINT handling */
void infinite_loop()
{
    while(1) {};
}

void illegal_instruction()
{
    /* I couldn't find an easy way to cause this one, so I'm cheating */
    raise(SIGILL);
}

