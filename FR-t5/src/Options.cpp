/************************************************************
 *
 * options -- an interface for getopt_long
 *
 * Copyright (C) 1999 Sebastian Will
 *
 * wills@informatik.uni-muenchen.de
 *
 ************************************************************/
 
//============================================================================
//
// Changed by Wu Aiping (wuaiping@moon.ibp.ac.cn), Aug, 2007
//
//============================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		OPTIONS.cpp
				ABSTRACT:	THREADING OPTIONS
				VERSION:	0.01
				AUTHOR:	Miao Zhichao
				DATE:		2012.03.19 Mon
				CONTACT:	chichaumiau@gmail.com

				NOTICE: This is free software and the source code is freely
				available. You are free to redistribute or modify under the
				conditions that (1) this notice is not removed or modified
				in any way and (2) any modified versions of the program are
				also available for free.
						** Absolutely no Warranty **
				***********************************************************/
#include "Options.hpp"


using namespace std;



/* prototype */
/* decode a string according to arg_type, for internal use */ 
bool decode_argument(void *argument, int arg_type, char *optarg);
int count_opts(option_def *options);
char *sprint_option_name(char *buf,option_def *options,int i);
char *sprint_option_name_opt(char *buf,option_def *options,int i);

char buf[512]; /* buffer used for constructing strings */


/************************************************************
 *
 * process_options()
 * interface to process options 
 *
 * returns TRUE   if all required options are read
 *         FALSE  if a required option is missing,
 *                or unknown options occured 
 *
 * last entry in definitions array should contain only 0's
 *   (in fact longname==0 and shortname==0 is tested)
 *
 ************************************************************/
bool
process_options(int argc, char *argv[], option_def *options) {
  int i,j,k;        /* counter */
  int num_opts; /* number of options in options[] */
  
  char c;
  
  char *short_opts; 
  struct option *long_opts;
  bool *is_set;

  int index;

  num_opts=count_opts(options);

  /* alloc sufficient memory */
  short_opts = (char *)malloc(num_opts*2 + 1);
  long_opts  = (struct option *)malloc((num_opts+1)*sizeof(struct option));
  is_set     = (bool *)malloc(num_opts*sizeof(bool)); 
			
  /* generate short options string and long options struct */
  for (i=0,j=0,k=0; i<num_opts; ++i) {
    /* short options */
    if (options[i].shortname) {
      short_opts[j] = options[i].shortname;
      ++j;
      if (options[i].argument!=0) { /* with argument */
	short_opts[j]=':';
	++j;
	if (options[i].deflt!=0) { /* with optional argument */
	  short_opts[j]=':';
	  ++j;
	}
      }
    }
    /* long options */
    if (options[i].longname) {
      long_opts[k].name = options[i].longname;
      if (options[i].argument==0) 
	long_opts[k].has_arg = no_argument;
      else {
	if (options[i].deflt==0)
	  long_opts[k].has_arg = required_argument;
	else
	  long_opts[k].has_arg = optional_argument;
      }
      long_opts[k].flag = &index;
      long_opts[k].val  = i;
      ++k;
    }
  }
  
  /* clear option flags */
  for (i=0,j=0,k=0; i<num_opts; ++i) {
    if (options[i].flag) *(options[i].flag) = FALSE;
    is_set[i] = FALSE;
  }

  /* set default values */
  for (i=0,j=0,k=0; i<num_opts; ++i) {
    is_set[i] = decode_argument(options[i].argument,
				options[i].arg_type,
				options[i].deflt);
  }

  /* main loop to process options */
  while ((c=getopt_long(argc,argv,short_opts,long_opts,0))!=EOF) {
    switch (c) {      
    case '?':
      return FALSE;
    case ':':
      return FALSE;
    default:
      if (c!=0) { /* short option */
	/* find option index */
	for (i=0; i<num_opts && options[i].shortname != c; ++i)
	  ;
	assert(i < num_opts);
	index = i;
      } /* else: long option */
      
      if (options[index].flag) *(options[index].flag)=TRUE;
      
      is_set[index] |= decode_argument(options[index].argument,
				       options[index].arg_type,
				       optarg);
    }
  }

  /* read no option arguments */
  for (i=0; i<num_opts; ++i)
    if (options[i].longname==0 && options[i].shortname==0)
      if (optind<argc) {
	if (options[i].flag) *options[i].flag=TRUE;
	is_set[i] |= decode_argument(options[i].argument,
				     options[i].arg_type,
				     argv[optind]);
	optind++;
      }


  /* are mandatory option arguments set */
  for (i=0; i<num_opts; ++i)
    if (options[i].argument != 0
	&& is_set[i]        == FALSE
	&& options[i].flag  == 0) 
      {
	fprintf(stderr,"Mandatory option and/or argument missing: ");
	if (options[i].longname)
	  fprintf(stderr,"--%s\n",options[i].longname);
	else if (options[i].shortname)
	  fprintf(stderr,"-%c\n",options[i].shortname);
	else
	  fprintf(stderr,"<%s>\n",
		  options[i].argname?options[i].argname:"param");
	return FALSE;
      }

  /* free allocated memory */
  free(short_opts);
  free(long_opts);
  free(is_set);
  
  return TRUE;
}


void
print_options(option_def options[]) {
  int i;        /* counter */
  int num_opts; /* number of options in options[] */

  num_opts=count_opts(options);

  for (i=0; i < num_opts; ++i) {
    if (!(options[i].flag!=0 
	  && *(options[i].flag)==FALSE 
	  && options[i].argument==0
	  )
	) {
      printf("  %-30s = ", sprint_option_name(buf,options,i));
      
      if (options[i].flag!=0
	  && options[i].argument==0 
	  && *(options[i].flag) ) printf("ON");
      else {
	if (options[i].flag==0 || *options[i].flag) {
	  if (options[i].argument)
	    switch (options[i].arg_type) {
	    case O_ARG_STRING:
	      printf ("%s",*((char **)options[i].argument));
	      break;
	    case O_ARG_INT:
	      printf ("%d",*((int*)(options[i].argument)));
	      break;
	    case O_ARG_FLOAT:
	      printf ("%f",*((float*)(options[i].argument)));
	      break;
	    case O_ARG_DOUBLE:
	      printf ("%f",*((double*)(options[i].argument)));
	      break;
	    default:
	      printf ("has unknown type");
	  }
	  else printf ("ON");
	}
      }
      printf("\n");
    }
  }
}

/************************************************************
 * void
 * print_usage()
 *
 * prints a standard usage string suited for short help output 
 *
 ************************************************************/
void
print_usage(char *progname, option_def options[]) {
  int i;        /* counter */
  int num_opts; /* number of options in options[] */

  num_opts=count_opts(options);
  
  printf("%s ", progname);
   
  for (i=0; i < num_opts; ++i) {
    /* options and no options*/
    printf("%s",sprint_option_name_opt(buf,options,i));
  }
}

void
print_help(char *progname, option_def options[]) {
  int i;
  int num_opts = count_opts(options);

  printf("USAGE: "); print_usage(progname, options); printf("\n");

  for (i=0; i<num_opts; i++) {
    printf("    %-30s",sprint_option_name(buf,options,i)); 
    
    if (options[i].description)
      printf("%s",options[i].description);

    printf("\n");
  }
}

/************************************************************/


char *sprint_option_name(char *buf,option_def *options,int i) {
  char *start=buf;
  if (options[i].shortname) buf += sprintf(buf,"-%c",options[i].shortname);
  if (options[i].shortname && options[i].longname) buf += sprintf(buf,", ");
  if (options[i].longname) buf += sprintf(buf,"--%s ",options[i].longname);
  
  if (options[i].argument) 
    buf += sprintf(buf,"<%s>", options[i].argname?options[i].argname:"param");
  return start;
}

char *sprint_option_name_opt(char *buf,option_def *options,int i) {
  char *start=buf;
  bool mandatory = options[i].flag==0 && options[i].deflt==0 ;

  buf += sprintf(buf," ");

  if (!mandatory) buf += sprintf(buf,"[");
  
  if (options[i].shortname) buf += sprintf(buf,"-%c",options[i].shortname);
  if (options[i].shortname && options[i].longname) buf += sprintf(buf,",");
  if (options[i].longname) buf += sprintf(buf,"--%s",options[i].longname);
      
  if (options[i].argument) {
    
    buf += sprintf(buf,"<%s", options[i].argname?options[i].argname:"param");
    if (options[i].deflt) buf += sprintf(buf,"(%s)",options[i].deflt);
    buf += sprintf(buf,">");
  }
  
  if (!mandatory) buf += sprintf(buf,"]");
  return start;
}

bool
decode_argument(void *argument, int arg_type, char *optarg) {
  if (optarg == NULL) return FALSE;

  if (argument == 0) {
    fprintf(stderr,"process_options: no argument variable\n");
    exit(-1);
  }
  
  switch(arg_type) {
  case O_ARG_STRING:
    *((char **)argument) = optarg;
    break;
  case O_ARG_INT:
    sscanf(optarg, "%d", (int *)argument);
    break;
  case O_ARG_FLOAT:
    sscanf(optarg, "%f", (float *)argument);
    break;
  case O_ARG_DOUBLE:
    sscanf(optarg, "%lf", (double *)argument);
    break;
  default:
    fprintf(stderr,"process_options: unknown argument type\n");
    exit(-1);
    break;
  }
  return TRUE;;
}

/* count entries in options */
int count_opts(option_def *options) {
  int i;
  for (i=0; !(options[i].argument == NULL && options[i].flag == NULL)
	 ; ++i)
    ;
  return i;
}
