/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   main program to drive LAMMPS
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  int myid;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  clock_t t1 = clock();
  struct timeval newtime,gtime;
  gettimeofday(&gtime, NULL);

  LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
  lammps->input->file();
  delete lammps;
  clock_t t2 = clock();
      double totaltime = ((double)(t2 - t1))/CLOCKS_PER_SEC;
           printf("\n---- Total time taken by LAMMPS: %e \n", totaltime);
   if(myid == 0)
   {
           gettimeofday(&newtime, NULL);
            // Get difference in time
        double time1 = static_cast<double>(gtime.tv_sec) +
                        static_cast<double>(gtime.tv_usec)/(1.0e6);
        double time2 = static_cast<double>(newtime.tv_sec) +
                        static_cast<double>(newtime.tv_usec)/(1.0e6);
        double difference = time2 - time1;
           double total_time = ((double)(t2 - t1))/CLOCKS_PER_SEC;
           printf("\n---- Total time taken by LAMMPS: %e \n", total_time);
           printf("\n---- Total time taken by LAMMPS (timer): %e \n", difference);
            struct rusage data;
          getrusage(RUSAGE_SELF, &data);
           printf("*** MEMORY USED: %ld ***\n", data.ru_maxrss);

   }
  MPI_Finalize();
}
