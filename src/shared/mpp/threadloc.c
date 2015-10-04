/*********************************************************************/
/*                                                                   */
/*                   GNU General Public License                      */
/*                                                                   */
/* This file is part of the Flexible Modeling System (FMS).          */
/*                                                                   */
/* FMS is free software; you can redistribute it and/or modify it    */
/* under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version.                               */
/*                                                                   */
/* FMS is distributed in the hope that it will be useful,            */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/* GNU General Public License for more details.                      */
/*                                                                   */
/* You should have received a copy of the GNU General Public License */
/* along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  */
/*                                                                   */
/*********************************************************************/

/* Fortran-callable routine for returning the MLD ("brick"?) where
   this thread/process is located. */
#include <stdio.h>
#include <unistd.h>
#ifdef use_libMPI
#include <mpi.h>
#endif
int pe, npes;

#ifdef __sgi
#include <sys/pmo.h>
#include <sys/types.h>
#include <sys/stat.h>

extern pmo_handle_t *mpi_sgi_mld;
extern int mpi_sgi_dsm_ppm;

int find_nodenum(int mynodedev);

int mld_id_() {
/*   pmo_handle_t mymld; */
/*   int mynodedev; */
/*   int mymemorynode; */
#define SIZE 1000000
  int array[SIZE];
  pm_pginfo_t pginfo_buf;
  int thisdev, thisnode;

  bzero( array, sizeof(array) ); /* zero to force allocation */

  __pm_get_page_info( array, 1, &pginfo_buf, 1 );
  thisdev = pginfo_buf.node_dev;
  thisnode = find_nodenum(thisdev);
  return thisnode;
}

int find_nodenum(int mynodedev) {
  int i;
  struct stat sbuf;
  char buff[80];
  for (i=0; ;i++)	{
    sprintf(buff,"/hw/nodenum/%d",i);
    stat(buff, &sbuf);
    if (sbuf.st_ino	== mynodedev)
      return(i);
  }

}
#else
int mld_id_() { /* dummy routine for portability */
  return 0;
}
#endif /* sgi */

#ifdef test_threadloc
void main(int argc, char **argv) {
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &pe );
  MPI_Comm_size( MPI_COMM_WORLD, &npes );
#ifdef _OPENMP
#pragma omp parallel
 {
   int thrnum = omp_get_thread_num();
   printf( "pe=%d thrnum=%d mld=%d\n", pe, thrnum, mld_id_() );
 }
#endif
  printf( "pe=%d mld=%d\n", pe, mld_id_() );
  MPI_Finalize();
}
#endif
