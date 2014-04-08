/*
  File autogenerated by gengetopt version 2.22.2
  generated with the following command:
  gengetopt --unamed-opts 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include "getopt.h"

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "Fit models to network data";

const char *gengetopt_args_info_usage = "Usage: Stochastic Block Models [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "      --git-version             detailed version description  (default=off)",
  "      --assume_N_nodes=INT      Pre-create N nodes (0 to N-1), which may be \n                                  left with zero degree  (default=`0')",
  "      --stringIDs               string IDs in the input  (default=off)",
  "      --seed=INT                seed to drand48() and gsl_rng_set  \n                                  (default=`0')",
  "      --GT=STRING               The ground truth, one line per community.",
  "  -K, --K=INT                   Number of clusters, K  (default=`-1')",
  "  -i, --iterations=INT          How many iterations  (default=`10000')",
  "      --metroK.algo=INT         Use the simple Metropolis move on K  \n                                  (default=`1')",
  "      --metro1Comm1Edge.algo=INT\n                                Use the simple Metropolis move on K  \n                                  (default=`1')",
  "      --NearbyGibbs.algo=INT    Gibbs updated on Nearby comms  (default=`1')",
  "      --Simplest1Node.algo=INT    (default=`0')",
  "      --AnySM.algo=INT            (default=`1')",
  "      --SharedSM.algo=INT         (default=`1')",
  "      --M3.algo=INT               (default=`1')",
  "      --m.iidBernoulli=FLOAT    A simpler model for the edges. Default is off \n                                  (-1)  (default=`-1')",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_FLOAT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->git_version_given = 0 ;
  args_info->assume_N_nodes_given = 0 ;
  args_info->stringIDs_given = 0 ;
  args_info->seed_given = 0 ;
  args_info->GT_given = 0 ;
  args_info->K_given = 0 ;
  args_info->iterations_given = 0 ;
  args_info->metroK_algo_given = 0 ;
  args_info->metro1Comm1Edge_algo_given = 0 ;
  args_info->NearbyGibbs_algo_given = 0 ;
  args_info->Simplest1Node_algo_given = 0 ;
  args_info->AnySM_algo_given = 0 ;
  args_info->SharedSM_algo_given = 0 ;
  args_info->M3_algo_given = 0 ;
  args_info->m_iidBernoulli_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->git_version_flag = 0;
  args_info->assume_N_nodes_arg = 0;
  args_info->assume_N_nodes_orig = NULL;
  args_info->stringIDs_flag = 0;
  args_info->seed_arg = 0;
  args_info->seed_orig = NULL;
  args_info->GT_arg = NULL;
  args_info->GT_orig = NULL;
  args_info->K_arg = -1;
  args_info->K_orig = NULL;
  args_info->iterations_arg = 10000;
  args_info->iterations_orig = NULL;
  args_info->metroK_algo_arg = 1;
  args_info->metroK_algo_orig = NULL;
  args_info->metro1Comm1Edge_algo_arg = 1;
  args_info->metro1Comm1Edge_algo_orig = NULL;
  args_info->NearbyGibbs_algo_arg = 1;
  args_info->NearbyGibbs_algo_orig = NULL;
  args_info->Simplest1Node_algo_arg = 0;
  args_info->Simplest1Node_algo_orig = NULL;
  args_info->AnySM_algo_arg = 1;
  args_info->AnySM_algo_orig = NULL;
  args_info->SharedSM_algo_arg = 1;
  args_info->SharedSM_algo_orig = NULL;
  args_info->M3_algo_arg = 1;
  args_info->M3_algo_orig = NULL;
  args_info->m_iidBernoulli_arg = -1;
  args_info->m_iidBernoulli_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->git_version_help = gengetopt_args_info_help[2] ;
  args_info->assume_N_nodes_help = gengetopt_args_info_help[3] ;
  args_info->stringIDs_help = gengetopt_args_info_help[4] ;
  args_info->seed_help = gengetopt_args_info_help[5] ;
  args_info->GT_help = gengetopt_args_info_help[6] ;
  args_info->K_help = gengetopt_args_info_help[7] ;
  args_info->iterations_help = gengetopt_args_info_help[8] ;
  args_info->metroK_algo_help = gengetopt_args_info_help[9] ;
  args_info->metro1Comm1Edge_algo_help = gengetopt_args_info_help[10] ;
  args_info->NearbyGibbs_algo_help = gengetopt_args_info_help[11] ;
  args_info->Simplest1Node_algo_help = gengetopt_args_info_help[12] ;
  args_info->AnySM_algo_help = gengetopt_args_info_help[13] ;
  args_info->SharedSM_algo_help = gengetopt_args_info_help[14] ;
  args_info->M3_algo_help = gengetopt_args_info_help[15] ;
  args_info->m_iidBernoulli_help = gengetopt_args_info_help[16] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->assume_N_nodes_orig));
  free_string_field (&(args_info->seed_orig));
  free_string_field (&(args_info->GT_arg));
  free_string_field (&(args_info->GT_orig));
  free_string_field (&(args_info->K_orig));
  free_string_field (&(args_info->iterations_orig));
  free_string_field (&(args_info->metroK_algo_orig));
  free_string_field (&(args_info->metro1Comm1Edge_algo_orig));
  free_string_field (&(args_info->NearbyGibbs_algo_orig));
  free_string_field (&(args_info->Simplest1Node_algo_orig));
  free_string_field (&(args_info->AnySM_algo_orig));
  free_string_field (&(args_info->SharedSM_algo_orig));
  free_string_field (&(args_info->M3_algo_orig));
  free_string_field (&(args_info->m_iidBernoulli_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->git_version_given)
    write_into_file(outfile, "git-version", 0, 0 );
  if (args_info->assume_N_nodes_given)
    write_into_file(outfile, "assume_N_nodes", args_info->assume_N_nodes_orig, 0);
  if (args_info->stringIDs_given)
    write_into_file(outfile, "stringIDs", 0, 0 );
  if (args_info->seed_given)
    write_into_file(outfile, "seed", args_info->seed_orig, 0);
  if (args_info->GT_given)
    write_into_file(outfile, "GT", args_info->GT_orig, 0);
  if (args_info->K_given)
    write_into_file(outfile, "K", args_info->K_orig, 0);
  if (args_info->iterations_given)
    write_into_file(outfile, "iterations", args_info->iterations_orig, 0);
  if (args_info->metroK_algo_given)
    write_into_file(outfile, "metroK.algo", args_info->metroK_algo_orig, 0);
  if (args_info->metro1Comm1Edge_algo_given)
    write_into_file(outfile, "metro1Comm1Edge.algo", args_info->metro1Comm1Edge_algo_orig, 0);
  if (args_info->NearbyGibbs_algo_given)
    write_into_file(outfile, "NearbyGibbs.algo", args_info->NearbyGibbs_algo_orig, 0);
  if (args_info->Simplest1Node_algo_given)
    write_into_file(outfile, "Simplest1Node.algo", args_info->Simplest1Node_algo_orig, 0);
  if (args_info->AnySM_algo_given)
    write_into_file(outfile, "AnySM.algo", args_info->AnySM_algo_orig, 0);
  if (args_info->SharedSM_algo_given)
    write_into_file(outfile, "SharedSM.algo", args_info->SharedSM_algo_orig, 0);
  if (args_info->M3_algo_given)
    write_into_file(outfile, "M3.algo", args_info->M3_algo_orig, 0);
  if (args_info->m_iidBernoulli_given)
    write_into_file(outfile, "m.iidBernoulli", args_info->m_iidBernoulli_orig, 0);
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  FIX_UNUSED (field);
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_FLOAT:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (
  int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "git-version",	0, NULL, 0 },
        { "assume_N_nodes",	1, NULL, 0 },
        { "stringIDs",	0, NULL, 0 },
        { "seed",	1, NULL, 0 },
        { "GT",	1, NULL, 0 },
        { "K",	1, NULL, 'K' },
        { "iterations",	1, NULL, 'i' },
        { "metroK.algo",	1, NULL, 0 },
        { "metro1Comm1Edge.algo",	1, NULL, 0 },
        { "NearbyGibbs.algo",	1, NULL, 0 },
        { "Simplest1Node.algo",	1, NULL, 0 },
        { "AnySM.algo",	1, NULL, 0 },
        { "SharedSM.algo",	1, NULL, 0 },
        { "M3.algo",	1, NULL, 0 },
        { "m.iidBernoulli",	1, NULL, 0 },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVK:i:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          cmdline_parser_print_version ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'K':	/* Number of clusters, K.  */
        
        
          if (update_arg( (void *)&(args_info->K_arg), 
               &(args_info->K_orig), &(args_info->K_given),
              &(local_args_info.K_given), optarg, 0, "-1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "K", 'K',
              additional_error))
            goto failure;
        
          break;
        case 'i':	/* How many iterations.  */
        
        
          if (update_arg( (void *)&(args_info->iterations_arg), 
               &(args_info->iterations_orig), &(args_info->iterations_given),
              &(local_args_info.iterations_given), optarg, 0, "10000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "iterations", 'i',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          /* detailed version description.  */
          if (strcmp (long_options[option_index].name, "git-version") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->git_version_flag), 0, &(args_info->git_version_given),
                &(local_args_info.git_version_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "git-version", '-',
                additional_error))
              goto failure;
          
          }
          /* Pre-create N nodes (0 to N-1), which may be left with zero degree.  */
          else if (strcmp (long_options[option_index].name, "assume_N_nodes") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->assume_N_nodes_arg), 
                 &(args_info->assume_N_nodes_orig), &(args_info->assume_N_nodes_given),
                &(local_args_info.assume_N_nodes_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "assume_N_nodes", '-',
                additional_error))
              goto failure;
          
          }
          /* string IDs in the input.  */
          else if (strcmp (long_options[option_index].name, "stringIDs") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->stringIDs_flag), 0, &(args_info->stringIDs_given),
                &(local_args_info.stringIDs_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "stringIDs", '-',
                additional_error))
              goto failure;
          
          }
          /* seed to drand48() and gsl_rng_set.  */
          else if (strcmp (long_options[option_index].name, "seed") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->seed_arg), 
                 &(args_info->seed_orig), &(args_info->seed_given),
                &(local_args_info.seed_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "seed", '-',
                additional_error))
              goto failure;
          
          }
          /* The ground truth, one line per community..  */
          else if (strcmp (long_options[option_index].name, "GT") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->GT_arg), 
                 &(args_info->GT_orig), &(args_info->GT_given),
                &(local_args_info.GT_given), optarg, 0, 0, ARG_STRING,
                check_ambiguity, override, 0, 0,
                "GT", '-',
                additional_error))
              goto failure;
          
          }
          /* Use the simple Metropolis move on K.  */
          else if (strcmp (long_options[option_index].name, "metroK.algo") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->metroK_algo_arg), 
                 &(args_info->metroK_algo_orig), &(args_info->metroK_algo_given),
                &(local_args_info.metroK_algo_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "metroK.algo", '-',
                additional_error))
              goto failure;
          
          }
          /* Use the simple Metropolis move on K.  */
          else if (strcmp (long_options[option_index].name, "metro1Comm1Edge.algo") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->metro1Comm1Edge_algo_arg), 
                 &(args_info->metro1Comm1Edge_algo_orig), &(args_info->metro1Comm1Edge_algo_given),
                &(local_args_info.metro1Comm1Edge_algo_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "metro1Comm1Edge.algo", '-',
                additional_error))
              goto failure;
          
          }
          /* Gibbs updated on Nearby comms.  */
          else if (strcmp (long_options[option_index].name, "NearbyGibbs.algo") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->NearbyGibbs_algo_arg), 
                 &(args_info->NearbyGibbs_algo_orig), &(args_info->NearbyGibbs_algo_given),
                &(local_args_info.NearbyGibbs_algo_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "NearbyGibbs.algo", '-',
                additional_error))
              goto failure;
          
          }
          /* .  */
          else if (strcmp (long_options[option_index].name, "Simplest1Node.algo") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->Simplest1Node_algo_arg), 
                 &(args_info->Simplest1Node_algo_orig), &(args_info->Simplest1Node_algo_given),
                &(local_args_info.Simplest1Node_algo_given), optarg, 0, "0", ARG_INT,
                check_ambiguity, override, 0, 0,
                "Simplest1Node.algo", '-',
                additional_error))
              goto failure;
          
          }
          /* .  */
          else if (strcmp (long_options[option_index].name, "AnySM.algo") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->AnySM_algo_arg), 
                 &(args_info->AnySM_algo_orig), &(args_info->AnySM_algo_given),
                &(local_args_info.AnySM_algo_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "AnySM.algo", '-',
                additional_error))
              goto failure;
          
          }
          /* .  */
          else if (strcmp (long_options[option_index].name, "SharedSM.algo") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->SharedSM_algo_arg), 
                 &(args_info->SharedSM_algo_orig), &(args_info->SharedSM_algo_given),
                &(local_args_info.SharedSM_algo_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "SharedSM.algo", '-',
                additional_error))
              goto failure;
          
          }
          /* .  */
          else if (strcmp (long_options[option_index].name, "M3.algo") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->M3_algo_arg), 
                 &(args_info->M3_algo_orig), &(args_info->M3_algo_given),
                &(local_args_info.M3_algo_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "M3.algo", '-',
                additional_error))
              goto failure;
          
          }
          /* A simpler model for the edges. Default is off (-1).  */
          else if (strcmp (long_options[option_index].name, "m.iidBernoulli") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->m_iidBernoulli_arg), 
                 &(args_info->m_iidBernoulli_orig), &(args_info->m_iidBernoulli_given),
                &(local_args_info.m_iidBernoulli_given), optarg, 0, "-1", ARG_FLOAT,
                check_ambiguity, override, 0, 0,
                "m.iidBernoulli", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */




  cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
