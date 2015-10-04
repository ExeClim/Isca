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

#if defined(__sgi) || defined(__aix) || defined(__SX)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define RUSAGE_SELF      0         /* calling process */
#define RUSAGE_CHILDREN  -1        /* terminated child processes */

#ifdef __aix
int memuse(void)
#else
int memuse_(void)
#endif
{
 struct rusage my_rusage;
 int iret;

 my_rusage.ru_maxrss = 0;
 iret = getrusage(RUSAGE_SELF,&my_rusage);
 iret = (int) my_rusage.ru_maxrss;
 return iret;
 /* printf("Max memuse in KB is %d \n",iret); */
}
#endif
