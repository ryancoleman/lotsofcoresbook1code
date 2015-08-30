
#include <mpi.h> //#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


#ifndef FALSE
#define FALSE (1 == 0)
#endif

#ifndef TRUE
#define TRUE (1 == 1)
#endif

#ifndef ASSERT
#define ASSERT(lcond,msg,ival) \
{ if (!(lcond)) { fprintf(stderr,"%s %d \n",msg,ival); exit(0); }  }
#endif

#ifndef doloop
#define doloop(i,istart,iend) for(i=(istart); i <= (iend); i++)
#endif



#ifndef STREQUAL
#define STREQUAL(s1,s2) (strcmp(s1,s2) == 0)
#endif

#ifndef STRCPY
#define STRCPY(idest,isrc) { strcpy(idest,isrc); }
#endif

/*
#define CLOCK0  get_current_time

#ifndef CLOCK0
#define CLOCK0(dummy)  (((double)clock())/( (double) CLOCKS_PER_SEC))
#endif				
*/

#ifndef STRING
#define STRING char*
#endif


#ifndef MAXROUTINES
#define MAXROUTINES 256
#endif

#ifndef MAXLEVELS
#define MAXLEVELS 8
#endif

#ifndef real8
#define real8 double
#endif

#ifndef integer
#define integer int
#endif

#define dictstart(i) dictstart_[(i)-1]
#define dicttotal(i) dicttotal_[(i)-1]
#define dictcount(i) dictcount_[(i)-1]
#define dictname(i) dictname_[(i)-1]
#define lastroutine(i) lastroutine_[(i)-1]

#ifndef logical
#define logical int
#endif

#ifndef STRING_LENGTH
#define STRING_LENGTH 128
#endif



#include "profinit.h"



static real8 dictstart(MAXROUTINES+1);
static real8 dicttotal(MAXROUTINES+1);
static integer dictcount(MAXROUTINES+1);
static STRING dictname(MAXROUTINES+1); 
static STRING lastroutine(MAXLEVELS+1);

static integer nroutine = 0; 
static integer nlevels = 0;

/*
#include <sys/time.h>
static double get_current_time()
{
   static struct timeval time_val;

   gettimeofday( &time_val, NULL );

   return(time_val.tv_sec + time_val.tv_usec /(1000.0 * 1000.0));
}
*/  
void profinit()
{

    integer i;

    nroutine = 0;
    doloop(i, 1, MAXROUTINES) {
        dictname(i) = (char *) malloc( MAXLEVELS*(STRING_LENGTH+1) );
        STRCPY( dictname(i), " ");
	dictstart(i) = 0.0;
	dictcount(i) = 0;
	dicttotal(i) = 0.0;
    };

    nlevels = 0;
    doloop(i, 1, MAXLEVELS) {
        lastroutine(i) = (char *) malloc( i*(STRING_LENGTH+1) );
        STRCPY( lastroutine(i)," ");
    };

    return;
}


void profstart(STRING rname)
{


    integer j, i, ipos;
    logical found;
    STRING name;

    /* ======= start execution =  */
//  printf("profstart %s\n", rname);

    name = rname;

    nlevels = nlevels + 1;
    ASSERT((1 <= nlevels) && (nlevels <= MAXLEVELS),
	   " ** profstart: invalid nlevels ", nlevels);

    if(nlevels > 1){
        STRCPY( lastroutine(nlevels), lastroutine(nlevels-1));
        strcat( lastroutine(nlevels), " " );
        strcat( lastroutine(nlevels), name );
    }else{
        STRCPY( lastroutine(nlevels), name );
    }
    found = FALSE;
    doloop(j, 1, nroutine) {
	i = nroutine - j + 1;	/* count down loop  */


	found = STREQUAL(dictname(i), lastroutine(nlevels));
	if (found) {
	    ipos = i;
	    break;
	};
    };

    if (!found) {
	nroutine = nroutine + 1;
	ASSERT(nroutine <= MAXROUTINES,
	       " ** profstart: nroutine > MAXROUTINES", nroutine);

	ipos = nroutine;
	STRCPY(dictname(ipos), lastroutine(nlevels));
	dictcount(ipos) = 0;
	dicttotal(ipos) = 0.0;

    };

    dictstart(ipos) = MPI_Wtime(); //=CLOCK0();
    dictcount(ipos) = dictcount(ipos) + 1;



    return;
}


void profend(STRING rname)
{

    integer j, i, ipos;
    logical found;
    STRING name;

    real8 tend;


    /* ======= start execution =  */
//  printf("profend %s\n", rname);

    tend = MPI_Wtime();   //=CLOCK0();

  
    if(nlevels > 1){
        name = (char *) malloc( nlevels*(STRING_LENGTH+1) );
        STRCPY( name, lastroutine(nlevels-1));
        strcat( name, " ");
        strcat( name, rname);
    }else{
        name = (char *) malloc( nlevels*(STRING_LENGTH+1) );
        STRCPY( name, rname);
    }



    ASSERT((1 <= nlevels) && (nlevels <= MAXLEVELS),
	   " ** profend: invalid nlevels ", nlevels);


    if (!STREQUAL(name, lastroutine(nlevels))) {
        fprintf(stderr," ** profend: name != lastroutine \n");
        fprintf(stderr," ** profend: name = %s \n",name );
        for (j=nlevels; j >= 1; j--) {
	   fprintf(stderr,"lastroutine(%d) = %s ", j,lastroutine(j));
	   };
	exit(0);
    };


    found = FALSE;
    doloop(j, 1, nroutine) {
	i = nroutine - j + 1;	/* count down loop  */


	found = STREQUAL(dictname(i), name);
	if (found) {
	    ipos = i;
	    break;
	};
    };



    if (!found) {
        fprintf(stderr,"** profend: routine %s not found\n",name);
	exit(0);
    };

    dicttotal(ipos) = dicttotal(ipos) + (tend - dictstart(ipos));

    nlevels = nlevels - 1;



    return;
}


void profstat()
{



    integer i, j;


//  const char *fname = "profstat.dat";
    const int MAX_FOUT_LENGTH = 30;
    char fname[MAX_FOUT_LENGTH]; // profstat%d.dat, myrank
    const char *mode = "w+";


    FILE *outdev;

    logical compute_max = FALSE;
    logical separate_output = TRUE;

    int myrank, size, rank;
    real8 maxtotal[MAXROUTINES];
    int maxtrank[MAXROUTINES];
    int maxtcount[MAXROUTINES];
    STRING allname[MAXROUTINES];
    int allroutine;
    int found;

//  MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    snprintf(fname, MAX_FOUT_LENGTH, "profstat%d.dat", myrank);

//  if(compute_max || !separate_output){
//      if(myrank == 0){
//          // initialization
//          MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//          outdev = fopen(fname, mode );
//          printf("rank %%d / total %%f secs / %%d times / name %%s\n");
//          fprintf(outdev, "rank %%d / total %%f secs / %%d times / name %%s\n");
//
//          if(compute_max){
//              allroutine = nroutine;
//              doloop(i, 1, allroutine) {
//                  allname[i-1] = (char *) malloc( MAXLEVELS*(STRING_LENGTH+1) );
//                  STRCPY(allname[i-1], dictname(i));
//                  maxtotal[i-1] = dicttotal(i);
//                  maxtrank[i-1] = 0; //myrank
//                  maxtcount[i-1] = dictcount(i);
//              };
//          }
//          if(!separate_output){
//              doloop(i, 1, nroutine){
//                  printf("%d / %lf / %d / %s\n",
//          	 	myrank,dicttotal(i),dictcount(i),dictname(i));
//                  fprintf(outdev,"%d / %lf / %d / %s\n",
//          	 	myrank,dicttotal(i),dictcount(i),dictname(i));
//              };
//          }
//          // recieve information from each process
//          //     update maximum of each routine
//          //     print out each of them
//          for(rank = 1; rank < size; rank++){
//              MPI_Recv(&nroutine, 1, MPI_INT, rank, MAXROUTINES+1, MPI_COMM_WORLD, &status);
//              MPI_Recv(&dicttotal(1), nroutine, MPI_DOUBLE, rank, MAXROUTINES+2, MPI_COMM_WORLD, &status);
//              MPI_Recv(&dictcount(1), nroutine, MPI_INT, rank, MAXROUTINES+3, MPI_COMM_WORLD, &status);
//              doloop(i, 1, nroutine) {
//                  MPI_Recv(dictname(i), MAXLEVELS*(STRING_LENGTH+1), MPI_CHAR, rank, i, MPI_COMM_WORLD, &status);
//                  if(compute_max){
//                      doloop(j, 1, allroutine) {
//                          if(found = STREQUAL(allname[j-1], dictname(i))){
//                              if(maxtotal[j-1] < dicttotal(i)){
//                                  maxtotal[j-1] = dicttotal(i);
//                                  maxtrank[j-1] = rank;
//                                  maxtcount[j-1] = dictcount(i);
//                              }
//                              break;
//                          }
//                      };
//                      if(!found){
//                          allroutine = allroutine + 1;
//                          ASSERT(allroutine <= MAXROUTINES,
//                             " ** profstart: allroutine > MAXROUTINES", allroutine);
//                          allname[allroutine-1] = (char *) malloc( MAXLEVELS*(STRING_LENGTH+1) );
//                          STRCPY(allname[allroutine-1], dictname(i));
//                          maxtotal[allroutine-1] = dicttotal(i);
//                          maxtrank[allroutine-1] = rank;
//                          maxtcount[allroutine-1] = dictcount(i);
//                      }
//                  }
//                  if(!separate_output){
//                      doloop(i, 1, nroutine){
//                          printf("%d / %lf / %d / %s\n",
//            	      	myrank,dicttotal(i),dictcount(i),dictname(i));
//                          fprintf(outdev,"%d / %lf / %d / %s\n",
//            	      	myrank,dicttotal(i),dictcount(i),dictname(i));
//                      };
//                  }
//              };
//          }
//          if(compute_max){
//              doloop(i, 1, allroutine) {
//                  printf("%d / %lf / %d / %s\n",
//              	 	maxtrank[i-1],maxtotal[i-1],maxtcount[i-1],allname[i-1]);
//                  fprintf(outdev, "%d / %lf / %d / %s\n",
//              	 	maxtrank[i-1],maxtotal[i-1],maxtcount[i-1],allname[i-1]);
//              };
//          }
//          fclose(outdev);
//      }else{ 
//          // other processes send information to root
//          MPI_Send(&nroutine, 1, MPI_INT, 0, MAXROUTINES+1, MPI_COMM_WORLD);
//          MPI_Send(&dicttotal(1), nroutine, MPI_DOUBLE, 0, MAXROUTINES+2, MPI_COMM_WORLD);
//          MPI_Send(&dictcount(1), nroutine, MPI_INT, 0, MAXROUTINES+3, MPI_COMM_WORLD);
//          doloop(i, 1, nroutine) {
//              MPI_Send(dictname(i), MAXLEVELS*(STRING_LENGTH+1), MPI_CHAR, 0, i, MPI_COMM_WORLD);
//          };
//          if(separate_output){
//              outdev = fopen(fname, mode );
//              printf("rank %%d / total %%f secs / %%d times / name %%s\n");
//              fprintf(outdev, "rank %%d / total %%f secs / %%d times / name %%s\n");
//       
//              doloop(i, 1, nroutine){
//                  printf("%d / %lf / %d / %s\n",
//              	myrank,dicttotal(i),dictcount(i),dictname(i));
//                  fprintf(outdev,"%d / %lf / %d / %s\n",
//              	myrank,dicttotal(i),dictcount(i),dictname(i));
//              };
//              fclose(outdev);
//          }
//      }
//  }else{ // output on its own
        outdev = fopen(fname, mode );
        printf("rank %%d / total %%f secs / %%d times / name %%s\n");
        fprintf(outdev, "rank %%d / total %%f secs / %%d times / name %%s\n");
  
        doloop(i, 1, nroutine){
            printf("%d / %lf / %d / %s\n",
        	myrank,dicttotal(i),dictcount(i),dictname(i));
            fprintf(outdev,"%d / %lf / %d / %s\n",
         	myrank,dicttotal(i),dictcount(i),dictname(i));
        };
        fclose(outdev);
//  }

    doloop( i,1,MAXROUTINES) {
      if (dictname(i) != NULL) {
	 free( dictname(i) );
         };
    };

    doloop(i, 1, MAXLEVELS) {
	if (lastroutine(i) != NULL) {
          free( lastroutine(i) );
	};
    };

    return;
}

void profinit_()
{
  profinit();
}

void profstart_( STRING str)
{
  profstart(str);
}

void profend_( STRING str)
{
  profend(str);
}

void profstat_()
{
  profstat();
}
