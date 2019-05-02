/* ===================================================================== */
/* to convert reflection files in other formats to CNS format

   written by: Paul Adams 4-98 

   copyright Yale University */
/* ===================================================================== */

/* include files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

/* macros */

#define MAXLINE 256
#define REFBLOCK 1000
#define MISSING -9999
#define DTREK_NOMEAS -1
#define HKL_NOMEAS -9998

/* structures */

/* an index - hkl */
struct index {
  int h;
  int k;
  int l; 
};
/* a measurement - iobs, sigi, fobs, sigf */
struct measurement {
  double i;
  double si;
  double f;
  double sf;
};
/* a reflection record - index and measurements (mean, observation) */
struct reflection {
  struct index index;
  struct measurement mean;
  struct measurement obs;
};

/* function definitions */
void usage(int);
int make_space();
int read_header(FILE *, char *);
int read_reflections(FILE *, char *, int);
int get_hkl_refs(FILE *, int *, int *, int *, 
		 double *, double *, double *, double *);
int calc_fobs(struct reflection *);
int write_cns_file(FILE *, int, char *, char *, char *);

/* global variables */
int space=0;
struct reflection *records;
struct reflection *reflection;
char *data_type;

/* main */

main(int argc, char *argv[]) {

  FILE *input, *output;
  int nref, i, write_mean;
  char  file_type[MAXLINE];
  char    in_name[MAXLINE];
  char   out_name[MAXLINE];

  /* defaults */
  write_mean = 1;
  strcpy(file_type,"hkl");
  strcpy(in_name,"");
  strcpy(out_name,"");

  setbuf(stdout,NULL);

  /* parse command line arguments */

  if ( (argc - 1) >= 1 ) {
    if ( strcmp(argv[1],"-h") == 0 ||
	 strcmp(argv[1],"-help") == 0 ) {
      usage(0);
    }
  }

  if ( (argc - 1) < 2 ) {
    usage(1);
  }
  
  for (i=1; i<argc; i++) {
    if ( strcmp(argv[i],"-f") == 0 ||
	 strcmp(argv[i],"-format") == 0 ) {
      if ( sscanf(argv[++i],"%s",file_type) < 0 ) {
	usage(2);
      }
      else {
	if ( strcmp(file_type,"hkl") != 0 &&
	     strcmp(file_type,"dtrek") != 0 ) {
	  usage(2);
	}
	printf("Type of reflection file: %s\n",file_type);
      }
    }
    else if ( strcmp(argv[i],"-nomean") == 0 ||
	      strcmp(argv[i],"-nm") == 0 ) {
      write_mean = 0;
      printf("Mean intensities will not be written\n");
    }
    else if ( strcmp(in_name,"") == 0 ) {
      if ( sscanf(argv[i],"%s",in_name) < 0 ) {
	usage(3);
      }
      else {
	printf("Input reflection file: %s\n",in_name);
      }
    }
    else if ( strcmp(out_name,"") == 0 ) {
      if ( sscanf(argv[i],"%s",out_name) < 0 ) {
	usage(3);
      }
      else {
	printf("Output reflection file: %s\n",out_name);
      }
    }
  }

  /* open input file */
  if ( ( input = fopen(in_name,"r") ) == NULL ) {
    printf("to_cns: cannot open file %s\n",in_name);
    exit(1); 
  }

  /* open output file */
  if ( ( output = fopen(out_name,"w") ) == NULL ) {
    printf("to_cns: cannot open file %s\n",out_name);
    exit(1); 
  }

  /* read header of input file */

  read_header(input,file_type);

  /* get records from file */

  nref = read_reflections(input,file_type,write_mean);
    
  /* write reflections to CNS format file */

  write_cns_file(output,nref,file_type,in_name,out_name);

  /* free up record space */

  free(records);

}

void usage(int error) {
  if ( error == 1 ) {
    fprintf(stderr,"to_cns: too few arguments\n");
  }
  if ( error == 2 ) {
    fprintf(stderr,"to_cns: mangled file type specifier\n");
    fprintf(stderr,"to_cns: supported formats -> hkl dtrek\n");
  }
  if ( error == 3 ) {
    fprintf(stderr,"to_cns: mangled file name \n");
  }
  fprintf(stderr,
	  "usage: to_cns [-f hkl|dtrek] [-nomean] in-file out-file\n");
  fprintf(stderr,
	  "defaults: format is hkl\n");
  fprintf(stderr,
	  "          mean I written for anomalous d*trek files\n");
  if ( error > 0 ) {
    exit (1);
  }
  else {
    exit (0);
  }
}

int make_space() {

  /* allocate an initial block of space */
  if ( space == 0 ) {
    if ( ! ( records = (struct reflection *) 
	     malloc(REFBLOCK * sizeof(struct reflection)) ) ) {
      printf("to_cns: could not allocate memory for records\n");
      printf("        current allocation= 0 bytes\n");
      return 0;
    }
    reflection = records;
    ++space;
  }
  else {
    /* if the reflection space is exhausted allocate some more */
    if ( reflection-records == (space * REFBLOCK)) {
      if ( ! ( records = (struct reflection *) 
	       realloc(records, (space + 1) * REFBLOCK 
		       * sizeof(struct reflection)) ) ) {
	printf("to_cns: could not allocate memory for records\n");
     /* added "int" ATB  6/18/08 */
	printf("        current allocation= %d bytes\n",
	       (int) ((space * REFBLOCK) * sizeof(struct reflection)));
	return 0;
      }
      /* make sure point points to the first element of the new space */
      reflection = records + (space * REFBLOCK);
      ++space;
    }
  }
  return 1;
}

int read_header(FILE *stream, char *file_type) {

  int i, nindex, ndata, nextra, ninfo;
  char  line[MAXLINE];
  char index[MAXLINE];
  char  data[MAXLINE];
  char extra[MAXLINE];
  char  info[MAXLINE];

  data_type = "unknown";

  if ( strcmp(file_type,"hkl") == 0 ) {
    fscanf(stream,"%d",&nindex);
    fscanf(stream,"%d",&ndata);
    fgets(line,MAXLINE,stream);
    fgets(line,MAXLINE,stream);
    if ( ndata > 0 ) {
      printf("to_cns: input file format error\n");
      printf("to_cns: file does not appear to be in hkl format\n");
      exit(1);
    }
  }
  else if ( strcmp(file_type,"dtrek") == 0 ) {
    fscanf(stream,"%d %d %d %d\n",&nindex,&ndata,&nextra,&ninfo);
    if ( nindex != 3 || ndata < 0 ) {
      printf("to_cns: input file format error\n");
      printf("to_cns: file does not appear to be in dtrek format\n");
      exit(1);
    }
    for (i=0;i<ninfo;i++) {
      fgets(info,MAXLINE,stream);
    }
    for (i=0;i<nindex;i++) {
      fscanf(stream,"%s",index);
    }
    data_type = "nonanomalous";
    for (i=0;i<ndata;i++) {
      fscanf(stream,"%s",data);
      if ( strcmp(data,"fIntensity+") == 0 ||
	   strcmp(data,"fIntensity-") == 0 ||
	   strcmp(data,"fSigmaI+") == 0 ||
	   strcmp(data,"fSigmaI-") == 0 ) {
	data_type = "anomalous";
      }
    }
    for (i=0;i<nextra;i++) {
      fscanf(stream,"%s",extra);
    }
  }
  return 1;
}

int read_reflections(FILE *stream, char *file_type, int write_mean) {

  int nref;
  int h, k, l;
  double mean_i, mean_s, plus_i, plus_s, minus_i, minus_s;

  if ( strcmp(file_type,"dtrek") == 0 ) {

    while ( ! feof(stream) ) {

      if ( strcmp(data_type,"anomalous") == 0 ) {
	if ( ! fscanf(stream,"%d %d %d %lf %lf %lf %lf %lf %lf",
		      &h, &k, &l,
		      &mean_i, &mean_s,
		      &plus_i, &plus_s,
		      &minus_i, &minus_s) == 9 ) {
	  printf("to_cns: error reading reflection record\n");
	  exit(1);
	}
	if ( ! feof(stream) ) {
	  if ( plus_i != DTREK_NOMEAS || plus_s != DTREK_NOMEAS ) {
	    
	    if ( ! make_space() ) {exit(1);}
	    
	    (reflection->index).h = h;
	    (reflection->index).k = k;
	    (reflection->index).l = l;
	    if ( write_mean ) {
	      (reflection->mean).i = mean_i;
	      (reflection->mean).si = mean_s;
	    }
	    else {
	      (reflection->mean).i = MISSING;
	      (reflection->mean).si = MISSING;
	    }
	    (reflection->obs).i = plus_i;
	    (reflection->obs).si = plus_s;
	    calc_fobs(reflection);
	    ++reflection;
	  }
	  if ( minus_i != DTREK_NOMEAS || minus_s != DTREK_NOMEAS ) {
	    
	    if ( ! make_space() ) {exit(1);}
	    
	    (reflection->index).h = -h;
	    (reflection->index).k = -k;
	    (reflection->index).l = -l;
	    if ( write_mean ) {
	      (reflection->mean).i = mean_i;
	      (reflection->mean).si = mean_s;
	    }
	    else {
	      (reflection->mean).i = MISSING;
	      (reflection->mean).si = MISSING;
	    }
	    (reflection->obs).i = minus_i;
	    (reflection->obs).si = minus_s;
	    calc_fobs(reflection);
	    ++reflection;
	  }
	}
      }

      if ( strcmp(data_type,"nonanomalous") == 0 ) {
	if ( ! fscanf(stream,"%d %d %d %lf %lf",
		      &h, &k, &l,
		      &mean_i, &mean_s) == 5 ) {
	  printf("to_cns: error reading reflection record\n");
	  exit(1);
	}
	
	if ( ! feof(stream) ) {
	  if ( ! make_space() ) {exit(1);}
	
	  (reflection->index).h = h;
	  (reflection->index).k = k;
	  (reflection->index).l = l;
	  (reflection->mean).i = MISSING;
	  (reflection->mean).si = MISSING;
	  (reflection->obs).i = mean_i;
	  (reflection->obs).si = mean_s;
	  calc_fobs(reflection);
	  ++reflection;
	}
      }
    }
  }
  else if ( strcmp(file_type,"hkl") == 0 ) {

    while ( ! feof(stream) ) {

      if ( ! feof(stream) ) {
	if ( ! (get_hkl_refs(stream, &h, &k, &l,
			     &plus_i, &plus_s, 
			     &minus_i, &minus_s)) ) {
	  printf("to_cns: error reading reflection record\n");
	  exit(1);
	}
	
	if ( plus_i != HKL_NOMEAS && minus_i != HKL_NOMEAS ) {
	  data_type = "anomalous";
	}
	
	if ( plus_i != HKL_NOMEAS || plus_s != HKL_NOMEAS ) {
	  
	  if ( ! make_space() ) {exit(1);}
	  
	  (reflection->index).h = h;
	  (reflection->index).k = k;
	  (reflection->index).l = l;
	  (reflection->mean).i = MISSING;
	  (reflection->mean).si = MISSING;
	  (reflection->obs).i = plus_i;
	  (reflection->obs).si = plus_s;
	  calc_fobs(reflection);
	  ++reflection;
	}
	if ( minus_i != HKL_NOMEAS || minus_s != HKL_NOMEAS ) {
	  
	  if ( ! make_space() ) {exit(1);}
	  
	  (reflection->index).h = -h;
	  (reflection->index).k = -k;
	  (reflection->index).l = -l;
	  (reflection->mean).i = MISSING;
	  (reflection->mean).si = MISSING;
	  (reflection->obs).i = minus_i;
	  (reflection->obs).si = minus_s;
	  calc_fobs(reflection);
	  ++reflection;
	}
      }
    
      if ( strcmp(data_type,"unknown") == 0 ) {
	data_type = "nonanomalous";
      }
    }
  }

  nref = reflection - records;

  printf("Diffraction data is of type: %s\n",data_type);
  printf("Number of reflections read: %d\n",nref);

  return nref;

}

int get_hkl_refs(FILE *stream, int *h, int *k, int *l,
		 double *plus_i, double *plus_s,
		 double *minus_i, double *minus_s) {
  int i, nvars;
  char line[MAXLINE];
  char word[MAXLINE];
/* changes to initialization ATB 6/18/08 */
  for (i=0;i<MAXLINE;i++) {
    line[i] = 0;
    word[i] = '\0';
  }

  fgets(line,MAXLINE,stream);

  strncpy(word,&line[0],4);
  *h = atoi(word);

  strncpy(word,&line[4],4);
  *k = atoi(word);

  strncpy(word,&line[8],4);
  *l = atoi(word);

  strncpy(word,&line[12],8);
  *plus_i = atof(word);

  strncpy(word,&line[20],8);
  *plus_s = atof(word);

  strncpy(word,&line[28],8);
  *minus_i = atof(word);

  strncpy(word,&line[36],8);
  *minus_s = atof(word);

  nvars = 3;

  if ( ( *plus_i == 0.0 && *plus_s ==  0.0 ) ||
       ( *plus_i == 0.0 && *plus_s == -1.0 ) ) {
    *plus_i = HKL_NOMEAS;
    *plus_s = HKL_NOMEAS;
  }
  else {
    nvars += 2;
  }
  if ( *minus_i == 0.0 && *minus_s == 0.0 ) {
    *minus_i = HKL_NOMEAS;
    *minus_s = HKL_NOMEAS;
  }
  else {
    nvars += 2;
  }

  return nvars;

}

int calc_fobs(struct reflection *reflection) {

  if ( (reflection->obs).i >= 0 ) {
    (reflection->obs).f = sqrt((reflection->obs).i);
    if ( (reflection->obs).si < (reflection->obs).i ) {
      (reflection->obs).sf = (reflection->obs).f - 
	               sqrt(((reflection->obs).i) - 
                             (reflection->obs).si);
    }
    else {
      (reflection->obs).sf = (reflection->obs).f;
    }
  }
  else {
    (reflection->obs).f = 0;
    (reflection->obs).sf = 0;
  }

  if ( (reflection->mean).i >= 0 ) {
    (reflection->mean).f = sqrt((reflection->mean).i);
    if ( (reflection->mean).si < (reflection->mean).i ) {
      (reflection->mean).sf = (reflection->mean).f - 
	               sqrt(((reflection->mean).i) - 
                             (reflection->mean).si);
    }
    else {
      (reflection->mean).sf = (reflection->mean).f;
    }
  }
  else {
    (reflection->mean).f = 0;
    (reflection->mean).sf = 0;
  }
  return 1;
}

int write_cns_file(FILE *stream, int nref, 
		   char *file_type, char *in_name, char *out_name) {

  int i_mean=0, nprint=0;

  fprintf(stream,"{* CNS file %s converted from %s file %s \n",
	  out_name,file_type,in_name);
  fprintf(stream,"   Converted by program to_cns *}\n");

  reflection = records;
  while ( (reflection-records) < nref ) {
    if ( (reflection->mean).i != MISSING || 
	 (reflection->mean).si != MISSING ) {
      i_mean = 1;
    }
    if ( (reflection->obs).i != MISSING ||
	 (reflection->mean).i != MISSING ) {
      nprint++;
    }
    ++reflection;
  }

  printf("Number of reflections to be written: %d\n",nprint);
  
  fprintf(stream,"NREFlections= %d\n",nprint);

  if ( strcmp(data_type,"anomalous") == 0 ) {
    fprintf(stream,"ANOMalous= TRUE\n");
  }
  else {
    fprintf(stream,"ANOMalous= FALSe\n");
  }

  fprintf(stream,"DECLare NAME=IOBS  DOMAin=RECIprocal TYPE=REAL END\n");
  fprintf(stream,"DECLare NAME=SIGI  DOMAin=RECIprocal TYPE=REAL END\n");
  fprintf(stream,"DECLare NAME=FOBS  DOMAin=RECIprocal TYPE=REAL END\n");
  fprintf(stream,"DECLare NAME=SIGMA DOMAin=RECIprocal TYPE=REAL END\n");

  if ( i_mean ) {
    fprintf(stream,"DECLare NAME=IOBS_MEAN  DOMAin=RECIprocal TYPE=REAL END\n");
    fprintf(stream,"DECLare NAME=SIGI_MEAN  DOMAin=RECIprocal TYPE=REAL END\n");
    fprintf(stream,"DECLare NAME=FOBS_MEAN  DOMAin=RECIprocal TYPE=REAL END\n");
    fprintf(stream,"DECLare NAME=SIGMA_MEAN DOMAin=RECIprocal TYPE=REAL END\n");
  }

  printf("Writing reflections ");

  reflection = records;
  while ( (reflection-records) < nref ) {
    if ( (reflection->obs).i != MISSING || 
	 (reflection->mean).i != MISSING ) {
      fprintf(stream,"INDEx= %4d %4d %4d ",
	      (reflection->index).h,
	      (reflection->index).k,
	      (reflection->index).l);
      if ( (reflection->obs).i != MISSING ) {
	fprintf(stream,"IOBS= %12.2f SIGI= %12.2f\n",
		(reflection->obs).i,
		(reflection->obs).si);
      }
      else {
	fprintf(stream,"                      ");
	fprintf(stream,"IOBS=          0.00 SIGI=          0.00\n");
      }
      fprintf(stream,"                      ");
      fprintf(stream,"FOBS= %10.2f SIGMA= %10.2f\n",
	      (reflection->obs).f,
	      (reflection->obs).sf);

      if ( i_mean ) {
	if ( (reflection->mean).i != MISSING ) {
	  fprintf(stream,"                      ");
	  fprintf(stream,"IOBS_MEAN= %12.2f SIGI_MEAN= %12.2f\n",
		  (reflection->mean).i,
		  (reflection->mean).si);
	}
	else {
	  fprintf(stream,"                      ");
	  fprintf(stream,"IOBS_MEAN=         0.00 SIGI_MEAN=         0.00\n");
	}
	fprintf(stream,"                      ");
	fprintf(stream,"FOBS_MEAN= %10.2f SIGMA_MEAN= %10.2f\n",
		(reflection->mean).f,
		(reflection->mean).sf);
      }
    }
    if ( (reflection-records) % 1000 == 0 ) {
      if ( (reflection-records) / 1000 == 50 ) {
	printf(".\n");
      }
      else {
	printf(".");
      }
    }
    ++reflection;
  }
  printf("\nClosing CNS reflection file\n");
  return 1;
}
