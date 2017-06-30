/************************************************************
 *
 * options -- an interface for getopt_long
 *
 * Copyright (C) 199 Sebastian Will
 *
 * wills@informatik.uni-muenchen.de
 *
 ************************************************************/

/************************************************************
 * 
 * interface to getopt
 *  
 *  - easier to use
 *  - may be restricted to certain cases
 *  - hopefully useable for many cases
 *
 ************************************************************/
 
//============================================================================
//
// Changed by Wu Aiping (wuaiping@moon.ibp.ac.cn), Aug, 2007
//
//============================================================================

#ifndef OPTIONS_H
#define OPTIONS_H

#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <assert.h>

using namespace std;


#ifndef bool
#define bool int
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

/* argument types */
#define O_NO_ARG 0
#define O_ARG_STRING 0
#define O_ARG_INT 1
#define O_ARG_FLOAT 2
#define O_ARG_DOUBLE 3

/* definition of a option */
typedef struct {
  char *longname; /* long option name */
  char shortname; /* short option char */
  bool *flag; /* pointer to flag that indicates if option given */
  int arg_type; /* type of argument */
  void *argument; /* pointer to variable that should hold argument, 0 indicates no arg */
  char *deflt; /* pointer to default argument, if arg optional. otherwise 0 */
  char *argname; /* optional name for an argument (shown in usage string)*/
  char *description; /*optionl description (shown in help) */
}
option_def;

/* longname==0 and shortname==0 is not allowed for regular options definition */

/*
Example for option_def array

bool opt_help;
.
.
.
int optVal_size;
int default_size=1000;

struct option_def my_options[] = {
  {"help",'h',&opt_help,O_NO_ARG,0,0,0,"This help"},
  {"num",'n',&opt_num,O_ARG_INT,&optVal_num,0,0,"Some arbitrary number"},
  {"output",'o',0,O_ARG_STRING,&outputfile,0,
      "output-file", "File for output"}, // mandatory
  {"size",'s',0,O_ARG_INT,&optVal_size,"100","size","Size of problem"},
  {0,0,0,O_ARG_STRING,&inputfile,0,"input-file","File for input"},
  {0,0,0,0,0,0,0,0}
}; 

the last entry shows how to define non-option command line
arguments. If there is more than one such definition, the order of
those argument definitions is important. This is different to the
option argument definitions.

*/

/************************************************************/
/* prototypes */

/* process the options */
bool process_options(int argc, char *argv[], option_def options[]);

/* print a usage string */
void print_usage(char *progname, option_def options[]);

/* print a longer help */
void print_help(char *progname, option_def options[]);

/* print the options which are read*/
void print_options(option_def options[]);

#endif
