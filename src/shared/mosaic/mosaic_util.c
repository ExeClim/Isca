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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef use_libMPI 
#include <mpi.h>
#endif
#include "mosaic_util.h"
#include "constant.h"

#define HPI (0.5*M_PI)
#define TPI (2.0*M_PI)
#define TOLORENCE (1.e-6)
#define EPSLN (1.e-10)
#define EPSLN10 (1.e-10)


/***********************************************************
    void error_handler(char *str)
    error handler: will print out error message and then abort
***********************************************************/

void error_handler(const char *msg)
{
  fprintf(stderr, "FATAL Error: %s\n", msg );
#ifdef use_libMPI      
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(1);
#endif  
}; /* error_handler */

/*********************************************************************

   int nearest_index(double value, const double *array, int ia)

   return index of nearest data point within "array" corresponding to "value".
   if "value" is outside the domain of "array" then nearest_index = 0
   or = size(array)-1 depending on whether array(0) or array(ia-1) is
   closest to "value"

   Arguments:
     value:  arbitrary data...same units as elements in "array"
     array:  array of data points  (must be monotonically increasing)
     ia   :  size of array.

 ********************************************************************/
int nearest_index(double value, const double *array, int ia)
{
  int index, i;
  int keep_going;

  for(i=1; i<ia; i++){
    if (array[i] < array[i-1]) 
      error_handler("nearest_index: array must be monotonically increasing"); 
  }
  if (value < array[0] )
    index = 0;
  else if ( value > array[ia-1]) 
    index = ia-1;
  else
    {
      i=0;
      keep_going = 1;
      while (i < ia && keep_going) {
	i = i+1;
	if (value <= array[i]) {
	  index = i;
	  if (array[i]-value > value-array[i-1]) index = i-1;
	  keep_going = 0;
	}
      }
    }
  return index;

};

/******************************************************************/

void tokenize(const char * const string, const char *tokens, unsigned int varlen,
	      unsigned int maxvar, char * pstring, unsigned int * const nstr)
{
  size_t i, j, nvar, len, ntoken;
  int found, n;
  
  nvar = 0; j = 0;
  len = strlen(string);
  ntoken = strlen(tokens);
  /* here we use the fact that C array [][] is contiguous in memory */
  if(string[0] == 0)error_handler("Error from tokenize: to-be-parsed string is empty");
  
  for(i = 0; i < len; i ++){
    if(string[i] != ' ' && string[i] != '\t'){
      found = 0;
      for(n=0; n<ntoken; n++) {
	if(string[i] == tokens[n] ) {
	  found = 1;
	  break;
	}
      }
      if(found) {
	if( j != 0) { /* remove :: */
	  *(pstring + (nvar++)*varlen + j) = 0;
	  j = 0;
	  if(nvar >= maxvar) error_handler("Error from tokenize: number of variables exceeds limit");
	}
      }
      else {
        *(pstring + nvar*varlen + j++) = string[i];
        if(j >= varlen ) error_handler("error from tokenize: variable name length exceeds limit during tokenization");
      }
    }
  }
  *(pstring + nvar*varlen + j) = 0;
  
  *nstr = ++nvar;

}

/*******************************************************************************
  double maxval_double(int size, double *data)
  get the maximum value of double array
*******************************************************************************/
double maxval_double(int size, const double *data)
{
  int n;
  double maxval;

  maxval = data[0];
  for(n=1; n<size; n++){
    if( data[n] > maxval ) maxval = data[n];
  }

  return maxval;
  
}; /* maxval_double */


/*******************************************************************************
  double minval_double(int size, double *data)
  get the minimum value of double array
*******************************************************************************/
double minval_double(int size, const double *data)
{
  int n;
  double minval;

  minval = data[0];
  for(n=1; n<size; n++){
    if( data[n] < minval ) minval = data[n];
  }

  return minval;
  
}; /* minval_double */

/*******************************************************************************
  double avgval_double(int size, double *data)
  get the average value of double array
*******************************************************************************/
double avgval_double(int size, const double *data)
{
  int n;
  double avgval;

  avgval = 0;
  for(n=0; n<size; n++) avgval += data[n];
  avgval /= size;
  
  return avgval;
  
}; /* avgval_double */


/*******************************************************************************
  void latlon2xyz
  Routine to map (lon, lat) to (x,y,z)
******************************************************************************/
void latlon2xyz(int size, const double *lon, const double *lat, double *x, double *y, double *z)
{
  int n;
  
  for(n=0; n<size; n++) {
    x[n] = cos(lat[n])*cos(lon[n]);
    y[n] = cos(lat[n])*sin(lon[n]);
    z[n] = sin(lat[n]);
  }

} /* latlon2xyz */

/*------------------------------------------------------------
       void xyz2laton(np, p, xs, ys)
   Transfer cartesian coordinates to spherical coordinates
   ----------------------------------------------------------*/
void xyz2latlon( int np, const double *x, const double *y, const double *z, double *lon, double *lat)
{

  double xx, yy, zz;
  double dist, sinp;
  int i;

  for(i=0; i<np; i++) {
    xx = x[i];
    yy = y[i];
    zz = z[i];
    dist = sqrt(xx*xx+yy*yy+zz*zz);
    xx /= dist;
    yy /= dist;
    zz /= dist;
    
    if ( fabs(xx)+fabs(yy)  < EPSLN ) 
       lon[i] = 0;
     else
       lon[i] = atan2(yy, xx);
     lat[i] = asin(zz);
    
     if ( lon[i] < 0.) lon[i] = 2.*M_PI + lon[i];
  }

} /* xyz2latlon */

/*------------------------------------------------------------------------------
  double box_area(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
  return the area of a lat-lon grid box. grid is in radians.
  ----------------------------------------------------------------------------*/
double box_area(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
{
  double dx = ur_lon-ll_lon;
  double area;
  
  if(dx > M_PI)  dx = dx - 2.0*M_PI;
  if(dx < -M_PI) dx = dx + 2.0*M_PI;

  return (dx*(sin(ur_lat)-sin(ll_lat))*RADIUS*RADIUS ) ;
  
}; /* box_area */


/*------------------------------------------------------------------------------
  double poly_area(const x[], const y[], int n)
  obtains area of input polygon by line integrating -sin(lat)d(lon)
  Vertex coordinates must be in degrees.
  Vertices must be listed counter-clockwise around polygon.
  grid is in radians.
  ----------------------------------------------------------------------------*/
double poly_area_dimensionless(const double x[], const double y[], int n)
{
  double area = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1, lat2;
    double dy, dat;
    
    lat1 = y[ip];
    lat2 = y[i];
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;
    if (dx==0.0) continue;
    
    if ( fabs(lat1-lat2) < SMALL_VALUE) /* cheap area calculation along latitude */
      area -= dx*sin(0.5*(lat1+lat2));
    else {
#ifdef fix_truncate      
      dy = 0.5*(lat1-lat2);
      dat = sin(dy)/dy;
      area -= dx*sin(0.5*(lat1+lat2))*dat;
#else
      area += dx*(cos(lat1)-cos(lat2))/(lat1-lat2);
#endif
    }
  }
  if(area < 0)
    return (-area/(4*M_PI));
  else
    return (area/(4*M_PI));

}; /* poly_area */

double poly_area(const double x[], const double y[], int n)
{
  double area = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1, lat2;
    double dy, dat;
    
    lat1 = y[ip];
    lat2 = y[i];
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;
    if (dx==0.0) continue;
    
    if ( fabs(lat1-lat2) < SMALL_VALUE) /* cheap area calculation along latitude */
      area -= dx*sin(0.5*(lat1+lat2));
    else {
#ifdef fix_truncate      
      dy = 0.5*(lat1-lat2);
      dat = sin(dy)/dy;
      area -= dx*sin(0.5*(lat1+lat2))*dat;
#else
      area += dx*(cos(lat1)-cos(lat2))/(lat1-lat2);
#endif
    }
  }
  if(area < 0)
     return -area*RADIUS*RADIUS;
  else  
     return area*RADIUS*RADIUS;

}; /* poly_area */

double poly_area_no_adjust(const double x[], const double y[], int n)
{
  double area = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1, lat2;

    lat1 = y[ip];
    lat2 = y[i];
    if (dx==0.0) continue;
    
    if ( fabs(lat1-lat2) < SMALL_VALUE) /* cheap area calculation along latitude */
      area -= dx*sin(0.5*(lat1+lat2));
    else
      area += dx*(cos(lat1)-cos(lat2))/(lat1-lat2);
  }
  if(area < 0)
     return area*RADIUS*RADIUS;
  else
     return area*RADIUS*RADIUS;
}; /* poly_area_no_adjust */

int delete_vtx(double x[], double y[], int n, int n_del)
{
  for (;n_del<n-1;n_del++) {
    x[n_del] = x[n_del+1];
    y[n_del] = y[n_del+1];
  }
  
  return (n-1);
} /* delete_vtx */

int insert_vtx(double x[], double y[], int n, int n_ins, double lon_in, double lat_in)
{
  int i;

  for (i=n-1;i>=n_ins;i--) {
    x[i+1] = x[i];
    y[i+1] = y[i];
  }
  
  x[n_ins] = lon_in;
  y[n_ins] = lat_in;
  return (n+1);
} /* insert_vtx */

void v_print(double x[], double y[], int n)
{
  int i;

  for (i=0;i<n;i++) printf(" %20g   %20g\n", x[i], y[i]);
} /* v_print */

int fix_lon(double x[], double y[], int n, double tlon)
{
  double x_sum, dx;
  int i, nn = n, pole = 0;

  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLORENCE) pole = 1;
  if (0&&pole) {
    printf("fixing pole cell\n");
    v_print(x, y, nn);
    printf("---------");
  }

  /* all pole points must be paired */
  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLORENCE) {
    int im=(i+nn-1)%nn, ip=(i+1)%nn;

    if (y[im]==y[i] && y[ip]==y[i]) {
      nn = delete_vtx(x, y, nn, i);
      i--;
    } else if (y[im]!=y[i] && y[ip]!=y[i]) {
      nn = insert_vtx(x, y, nn, i, x[i], y[i]);
      i++;
    }
  }
  /* first of pole pair has longitude of previous vertex */
  /* second of pole pair has longitude of subsequent vertex */
  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLORENCE) {
    int im=(i+nn-1)%nn, ip=(i+1)%nn;

    if (y[im]!=y[i]) x[i] = x[im];
    if (y[ip]!=y[i]) x[i] = x[ip];
  }

  if (nn) x_sum = x[0]; else return(0);
  for (i=1;i<nn;i++) {
    double dx = x[i]-x[i-1];

    if      (dx < -M_PI) dx = dx + TPI;
    else if (dx >  M_PI) dx = dx - TPI;
    x_sum += (x[i] = x[i-1] + dx);
  }

  dx = (x_sum/nn)-tlon;
  if      (dx < -M_PI) for (i=0;i<nn;i++) x[i] += TPI;
  else if (dx >  M_PI) for (i=0;i<nn;i++) x[i] -= TPI;

  if (0&&pole) {
    printf("area=%g\n", poly_area(x, y,nn));
    v_print(x, y, nn);
    printf("---------");
  }

  return (nn);
} /* fix_lon */


/*------------------------------------------------------------------------------
  double great_circle_distance()
  computes distance between two points along a great circle
  (the shortest distance between 2 points on a sphere)
  returned in units of meter
  ----------------------------------------------------------------------------*/
double great_circle_distance(double *p1, double *p2)
{
  double dist, beta;
  
  /* This algorithm is not accurate for small distance 
  dist = RADIUS*ACOS(SIN(p1[1])*SIN(p2[1]) + COS(p1[1])*COS(p2[1])*COS(p1[0]-p2[0]));
  */
  beta = 2.*asin( sqrt( sin((p1[1]-p2[1])/2.)*sin((p1[1]-p2[1])/2.) + 
                               cos(p1[1])*cos(p2[1])*(sin((p1[0]-p2[0])/2.)*sin((p1[0]-p2[0])/2.)) ) );
  dist = RADIUS*beta;
  return dist;

}; /* great_circle_distance */


/* Compute the great circle area of a polygon on a sphere */
double great_circle_area(int n, const double *x, const double *y, const double *z) {
  int i;
  double pnt0[3], pnt1[3], pnt2[3];
  double sum, area;
  
  /* sum angles around polygon */
  sum=0.0;
  for ( i=0; i<n; i++) {
    /* points that make up a side of polygon */
    pnt0[0] = x[i];
    pnt0[1] = y[i];
    pnt0[2] = z[i];
    pnt1[0] = x[(i+1)%n];
    pnt1[1] = y[(i+1)%n];
    pnt1[2] = z[(i+1)%n];
    pnt2[0] = x[(i+2)%n];
    pnt2[1] = y[(i+2)%n];
    pnt2[2] = z[(i+2)%n];    

    /* compute angle for pnt1 */
    sum += spherical_angle(pnt1, pnt2, pnt0);

  }
  area = (sum - (n-2.)*M_PI) * RADIUS* RADIUS;
  return area;
}

/*------------------------------------------------------------------------------
  double spherical_angle(const double *p1, const double *p2, const double *p3)
           p3
         /
        /
       p1 ---> angle
         \
          \
           p2
 -----------------------------------------------------------------------------*/
double spherical_angle(const double *v1, const double *v2, const double *v3)
{
  double angle;
  double px, py, pz, qx, qy, qz, ddd;

  /* vector product between v1 and v2 */
  px = v1[1]*v2[2] - v1[2]*v2[1];
  py = v1[2]*v2[0] - v1[0]*v2[2];
  pz = v1[0]*v2[1] - v1[1]*v2[0];
  /* vector product between v1 and v3 */
  qx = v1[1]*v3[2] - v1[2]*v3[1];
  qy = v1[2]*v3[0] - v1[0]*v3[2];
  qz = v1[0]*v3[1] - v1[1]*v3[0];
    
  ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
  if ( ddd <= 0.0 )
    angle = 0. ;
  else {
    ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd);
    if ( abs(ddd)>1.) {
      angle = 2.*atan(1.0);    /* 0.5*pi */
      /*FIX (lmh) to correctly handle co-linear points (angle near pi or 0) */
      if (ddd < 0.)
	angle = 4.*atan(1.0); /*should be pi*/
      else  
	angle = 0.;
    }
    else 
      angle = acos( ddd );
  }
  
  return angle;
}; /* spherical_angle */

/*------------------------------------------------------------------------------
  double spherical_excess_area(p_lL, p_uL, p_lR, p_uR) 
  get the surface area of a cell defined as a quadrilateral 
  on the sphere. Area is computed as the spherical excess
  [area units are m^2]
  ----------------------------------------------------------------------------*/
double spherical_excess_area(const double* p_ll, const double* p_ul,
			     const double* p_lr, const double* p_ur, double radius)
{
  double area, ang1, ang2, ang3, ang4;
  double v1[3], v2[3], v3[3];

  /*   S-W: 1   */  
  latlon2xyz(1, p_ll, p_ll+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_lr, p_lr+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_ul, p_ul+1, v3, v3+1, v3+2);
  ang1 = spherical_angle(v1, v2, v3);

  /*   S-E: 2   */  
  latlon2xyz(1, p_lr, p_lr+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_ur, p_ur+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_ll, p_ll+1, v3, v3+1, v3+2);
  ang2 = spherical_angle(v1, v2, v3);

  /*   N-E: 3   */  
  latlon2xyz(1, p_ur, p_ur+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_ul, p_ul+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_lr, p_lr+1, v3, v3+1, v3+2);
  ang3 = spherical_angle(v1, v2, v3);
  
  /*   N-W: 4   */  
  latlon2xyz(1, p_ul, p_ul+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_ur, p_ur+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_ll, p_ll+1, v3, v3+1, v3+2);
  ang4 = spherical_angle(v1, v2, v3);

  area = (ang1 + ang2 + ang3 + ang4 - 2.*M_PI) * radius* radius;

  return area;
  
}; /* spherical_excess_area */


/*----------------------------------------------------------------------
    void vect_cross(e, p1, p2)
    Perform cross products of 3D vectors: e = P1 X P2
    -------------------------------------------------------------------*/
    
void vect_cross(const double *p1, const double *p2, double *e )
{
  
  e[0] = p1[1]*p2[2] - p1[2]*p2[1];
  e[1] = p1[2]*p2[0] - p1[0]*p2[2];
  e[2] = p1[0]*p2[1] - p1[1]*p2[0];

}; /* vect_cross */


/*----------------------------------------------------------------------
    double* vect_cross(p1, p2)
    return cross products of 3D vectors: = P1 X P2
    -------------------------------------------------------------------*/
    
double *cross(const double *p1, const double *p2)
{
  double *e=NULL;

  e = (double *)malloc(3*sizeof(double));
  e[0] = p1[1]*p2[2] - p1[2]*p2[1];
  e[1] = p1[2]*p2[0] - p1[0]*p2[2];
  e[2] = p1[0]*p2[1] - p1[1]*p2[0];

  return e;
  
}; /* cross */

double dot(const double *p1, const double *p2)
{

  return( p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2] );

}


double metric(const double *p) {
  return (sqrt(p[0]*p[0] + p[1]*p[1]+p[2]*p[2]) );
}


/* ----------------------------------------------------------------
   make a unit vector
   --------------------------------------------------------------*/
void normalize_vect(double *e)
{
  double pdot;
  int k;

  pdot = e[0]*e[0] + e[1] * e[1] + e[2] * e[2];
  pdot = sqrt( pdot ); 

  for(k=0; k<3; k++) e[k] /= pdot;
};


/*------------------------------------------------------------------
  void unit_vect_latlon(int size, lon, lat, vlon, vlat)

  calculate unit vector for latlon in cartesian coordinates

  ---------------------------------------------------------------------*/
void unit_vect_latlon(int size, const double *lon, const double *lat, double *vlon, double *vlat)
{
  double sin_lon, cos_lon, sin_lat, cos_lat;
  int n;
  
  for(n=0; n<size; n++) {
    sin_lon = sin(lon[n]);
    cos_lon = cos(lon[n]);
    sin_lat = sin(lat[n]);
    cos_lat = cos(lat[n]);
    
    vlon[3*n] = -sin_lon;
    vlon[3*n+1] =  cos_lon;
    vlon[3*n+2] =  0.;
    
    vlat[3*n]   = -sin_lat*cos_lon;
    vlat[3*n+1] = -sin_lat*sin_lon;
    vlat[3*n+2] =  cos_lat;
  }
}; /* unit_vect_latlon */


/* Intersect a line and a plane
   Intersects between the plane ( three points ) (entries in counterclockwise order)
   and the line determined by the endpoints l1 and l2 (t=0.0 at l1 and t=1.0 at l2)
   returns true if the two intersect and the output variables are valid
   outputs p containing the coordinates in the tri and t the coordinate in the line
   of the intersection.
   NOTE: the intersection doesn't have to be inside the tri or line for this to return true
*/
int intersect_tri_with_line(const double *plane, const double *l1, const double *l2, double *p,
			     double *t) {

  double M[3*3], inv_M[3*3];
  double V[3];
  double X[3];
  int is_invert=0;
  
  const double *pnt0=plane;
  const double *pnt1=plane+3;
  const double *pnt2=plane+6;

  /* To do intersection just solve the set of linear equations for both
     Setup M
  */
  M[0]=l1[0]-l2[0]; M[1]=pnt1[0]-pnt0[0]; M[2]=pnt2[0]-pnt0[0];
  M[3]=l1[1]-l2[1]; M[4]=pnt1[1]-pnt0[1]; M[5]=pnt2[1]-pnt0[1];
  M[6]=l1[2]-l2[2]; M[7]=pnt1[2]-pnt0[2]; M[8]=pnt2[2]-pnt0[2];


  /* Invert M */
  is_invert = invert_matrix_3x3(M,inv_M);
  if (!is_invert) return 0;

  /* Set variable holding vector */
  V[0]=l1[0]-pnt0[0];
  V[1]=l1[1]-pnt0[1];
  V[2]=l1[2]-pnt0[2];

  /* Calculate solution */
  mult(inv_M, V, X);

  /* Get answer out */
  *t=X[0];
  p[0]=X[1];
  p[1]=X[2];

  return 1;
}


void mult(double m[], double v[], double out_v[]) {

  out_v[0]=m[0]*v[0]+m[1]*v[1]+m[2]*v[2];
  out_v[1]=m[3]*v[0]+m[4]*v[1]+m[5]*v[2];
  out_v[2]=m[6]*v[0]+m[7]*v[1]+m[8]*v[2];

}


/* returns 1 if matrix is inverted, 0 otherwise */
int invert_matrix_3x3(double m[], double m_inv[]) {

  const double det =  m[0] * (m[4]*m[8] - m[5]*m[7])
                     -m[1] * (m[3]*m[8] - m[5]*m[6])
                     +m[2] * (m[3]*m[7] - m[4]*m[6]);
#ifdef test_invert_matrix_3x3
  printf("det = %g\n", det);
#endif  
  if (fabs(det) < EPSLN10 ) return 0;

  const double deti = 1.0/det;

  m_inv[0] = (m[4]*m[8] - m[5]*m[7]) * deti;
  m_inv[1] = (m[2]*m[7] - m[1]*m[8]) * deti;
  m_inv[2] = (m[1]*m[5] - m[2]*m[4]) * deti;

  m_inv[3] = (m[5]*m[6] - m[3]*m[8]) * deti;
  m_inv[4] = (m[0]*m[8] - m[2]*m[6]) * deti;
  m_inv[5] = (m[2]*m[3] - m[0]*m[5]) * deti;

  m_inv[6] = (m[3]*m[7] - m[4]*m[6]) * deti;
  m_inv[7] = (m[1]*m[6] - m[0]*m[7]) * deti;
  m_inv[8] = (m[0]*m[4] - m[1]*m[3]) * deti;

  return 1;
}

void deleteNode(struct Node *list)
{
  struct Node *temp1=NULL, *temp2=NULL;
  
  temp1 = list->Next;
  if( temp1 ) {
    temp2 = temp1->Next;
    free(temp1);
    temp1 = temp2;
  }
}

void initNode(struct Node *node)
{
    node->x = 0;
    node->y = 0;
    node->z = 0;
    node->u = 0;  
    node->intersect = 0;
    node->inbound = 0;
    node->isInside = 0;
    node->Next = NULL;
    node->initialized=0;
    
}
  
void addBeg(struct Node *list, double x, double y, double z, int intersect, double u, int inbound)
{
  struct Node *temp=NULL;

  if(list == NULL) error_handler("addBeg: list is NULL");
  
  if(list->initialized) {
    temp = (struct Node *)malloc(sizeof(struct Node));
    temp->Next = list;
  }
  else{
    temp=list;
    temp->Next=NULL;
  }
  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->u = u;
  temp->intersect = intersect;
  temp->inbound = inbound;
  temp->initialized=1;
  temp->isInside = 0;
}

void addEnd(struct Node *list, double x, double y, double z, int intersect, double u, int inbound)
{

  struct Node *temp1=NULL, *temp2=NULL;

  if(list == NULL) error_handler("addBeg: list is NULL");
  
  if(list->initialized) {
    temp1 = (struct Node *)malloc(sizeof(struct Node));
    /* (x,y,z) might already in the list */ 
    temp2 = list;
    while (temp2) {
      if(fabs(temp2->x-x) <EPSLN10 && fabs(temp2->y-y) <EPSLN10 && fabs(temp2->z-z) <EPSLN10)
	return;
      temp2=temp2->Next;
    }
    temp2 = list;
    while(temp2->Next != NULL)  
      temp2=temp2->Next;  
  
    /* Append at the end of the list.  */
    temp2->Next=temp1;
  }
  else {
    temp1 = list;
  }

  temp1->x = x;
  temp1->y = y;
  temp1->z = z;
  temp1->u = u;
  temp1->intersect = intersect;
  temp1->inbound = inbound;
  temp1->Next = NULL;
  temp1->initialized=1;
  temp1->isInside = 0;
}

int length(struct Node *list)  
{  
   struct Node *cur_ptr=NULL;  
   int count=0;  
  
   cur_ptr=list;  
  
   while(cur_ptr)  
   {
     if(!cur_ptr->initialized) break;
      cur_ptr=cur_ptr->Next;
      count++;  
   }  
   return(count);  
}  

int sameNode(struct Node node1, struct Node node2)
{
  if( fabs(node1.x-node2.x)<EPSLN10 && fabs(node1.y-node2.y)<EPSLN10 && fabs(node1.z-node2.z)<EPSLN10 )
    return 1;
  else
    return 0;
}
      

void addNode(struct Node *list, struct Node inNode)
{

  addEnd(list, inNode.x, inNode.y, inNode.z, inNode.intersect, inNode.u, inNode.inbound);
  
}

struct Node *getNode(struct Node *list, struct Node inNode)
{
  struct Node *thisNode=NULL;
  struct Node *temp=NULL;

  temp = list;
  while( temp ) {
    if( sameNode( *temp, inNode ) ) {
      thisNode = temp;
      temp = NULL;
      break;
    }
    temp = temp->Next;
  }

  return thisNode;
}

struct Node *getNextNode(struct Node *list)
{
  return list->Next;
}

void copyNode(struct Node *node_out, struct Node node_in)
{

  node_out->x = node_in.x;
  node_out->y = node_in.y;
  node_out->z = node_in.z;
  node_out->u = node_in.u;
  node_out->intersect = node_in.intersect;
  node_out->inbound   = node_in.inbound;
  node_out->Next = NULL;
  node_out->initialized = node_in.initialized;
  node_out->isInside = node_in.isInside;
}

void printNode(struct Node *list, char *str)
{
  struct Node *temp;

  if(list == NULL) error_handler("printNode: list is NULL");
  if(str) printf("  %s \n", str);
  temp = list;
  while(temp) {
    if(temp->initialized) {
      printf(" (x, y, z, interset, inbound, isInside) = (%g,%g,%g,%d,%d,%d)\n",
	     temp->x, temp->y, temp->z, temp->intersect, temp->inbound, temp->isInside);
    }
    temp = temp->Next;
  }
  printf("\n");
}

int isHeadNode(struct Node *list, struct Node nodeIn)
{
  if(sameNode(*list, nodeIn))
    return 1;
  else
    return 0;
}
  
void removeNode(struct Node *list, struct Node nodeIn)
{

  struct Node *temp1=NULL;
  struct Node *temp2=NULL;
  int found=0;

  if( isHeadNode(list, nodeIn) )
    error_handler("removeNode: can not remove the head node of the list");
  
    /* search the list to find the point has the same value as nodeIn */
  temp1 = list;
  temp2 = list->Next;
  found = 0;
  while (temp2) {
    if(sameNode(*temp2, nodeIn)) {
      found = 1;
      temp1->Next = temp2->Next;
      free(temp2);
      break;
    }
    temp1 = temp2;
    temp2 = temp2->Next;
  }
  if(!found) error_handler("removeNode: nodeIn is not in the list");    
}


/* The following insert a intersection after non-intersect point (x2,y2,z2), if the point
   after (x2,y2,z2) is an intersection, if u is greater than the u value of the intersection,
   insert after, otherwise insert before
*/
void insertAfter(struct Node *list, double x, double y, double z, int intersect, double u, int inbound,
		 double x2, double y2, double z2)
{
  struct Node *temp1=NULL, *temp2=NULL;
  struct Node *temp;
  int found=0;

  /* (x,y,z) might in the list. Just set intersect to 2 to means both vertice and intersection */ 
  temp2 = list;
  while (temp2) {
    if(fabs(temp2->x-x) <EPSLN10 && fabs(temp2->y-y) <EPSLN10 && fabs(temp2->z-z) <EPSLN10) {
      if( temp2->intersect == 0) temp2->intersect = 2;
      temp2->u = u;
      return;
    }	
    temp2=temp2->Next;
  }

  
  temp1 = list;
  found = 0;
  while ( temp1 ) {
    if( temp1->x == x2 && temp1->y == y2 && temp1->z == z2 ) {
      found = 1;
      break;
    }
    temp1=temp1->Next;
  }
  if (!found) error_handler("inserAfter: non-intersect point (x,y,z) is not found in the list");
  temp2 = temp1->Next;
  while ( temp2 ) {
    if( temp2->intersect == 1 ) {
      if( temp2->u > u ) {
	break;
      }
    }
    else
      break;
    temp1 = temp2;
    temp2 = temp2->Next;
  }

  /* insert after temp1 */
  temp = (struct Node *)malloc(sizeof(struct Node));
  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->u = u;
  temp->intersect = intersect;
  temp->inbound = inbound;
  temp->isInside = 0;
  temp->initialized = 1;
  temp1->Next = temp;
  temp->Next = temp2;
  
}

double gridArea(struct Node *grid) {
  double x[20], y[20], z[20];
  struct Node *temp=NULL;
  double area;
  int n;
  
  temp = grid;
  n = 0;
  while( temp ) {
    x[n] = temp->x;
    y[n] = temp->y;
    z[n] = temp->z;
    n++;
    temp = temp->Next;
  }

  area = great_circle_area(n, x, y, z);

  return area;
  
}

int isIntersect(struct Node node) {

  return node.intersect;

}

void resetIntersectValue(struct Node *list ) {
  struct Node *temp=NULL;
  
  temp = list;
  while(temp) {
    if(temp->intersect == 2) temp->intersect = 0;
    temp=temp->Next;
  }
  
}

int getInbound( struct Node node )
{
  return node.inbound;
}

struct Node *getLast(struct Node *list)
{
  struct Node *temp1;

  temp1 = list;
  if( temp1 ) {
    while( temp1->Next ) {
      temp1 = temp1->Next;
    }
  }

  return temp1;
}


int getFirstInbound( struct Node *list, struct Node *nodeOut)
{
  struct Node *temp=NULL;

  temp=list;

  while(temp) {
    if( temp->inbound == 2 ) {
      copyNode(nodeOut, *temp);
      return 1;
    }
    temp=temp->Next;
  }

  return 0;
}

void getCoordinate(struct Node node, double *x, double *y, double *z)
{


  *x = node.x;
  *y = node.y;
  *z = node.z;

}

void setInbound(struct Node *interList, struct Node *list)
{

  struct Node *temp1=NULL, *temp2=NULL;
  int prev_is_inside, next_is_inside;
  
  /* for each point in interList, search through list to decide the inbound value the interList point */
  /* For each inbound point, the prev node should be outside and the next is inside. */
  if(length(interList) == 0) return;
  
  temp1 = interList;
  
  while(temp1) {
    temp2 = list;
    prev_is_inside = -1;
    while(temp2) {
      if(sameNode(*temp1, *temp2)) {
	if( prev_is_inside == -1 ) { /* -1 means temp2 is the first element in list */
	  temp2 = getLast(list);
	  prev_is_inside = isInside(temp2);
	}
	temp2 = temp2->Next;
	if(temp2==NULL) temp2 = list;
	next_is_inside = isInside(temp2);
	if( !prev_is_inside && next_is_inside )
	  temp1->inbound = 2;
	else
	  temp1->inbound = 1;
	break;
      }
      prev_is_inside = isInside(temp2);
      temp2 = temp2->Next;
    }
    if(prev_is_inside == -1) error_handler("mosaic_util.c: temp1 is not inside the list");
    temp1 = temp1->Next;
  }
  
}

int isInside(struct Node *node) {
  if( node->isInside || node->intersect)
    return 1;
  else
    return 0;
}

  

int insidePolygon( double x1, double y1, double z1, int npts, double *x2, double *y2, double *z2)
{
  int i, ip;
  double pnt0[3], pnt1[3], pnt2[3];
  double anglesum;
  
  anglesum = 0;

  pnt0[0] = x1;
  pnt0[1] = y1;
  pnt0[2] = z1;
  
  for (i=0;i<npts;i++) {
    ip = (i+1)%npts;
    pnt1[0] = x2[i];
    pnt1[1] = y2[i];
    pnt1[2] = z2[i];
    pnt2[0] = x2[ip];
    pnt2[1] = y2[ip];
    pnt2[2] = z2[ip];
    anglesum += spherical_angle(pnt0, pnt2, pnt1);
  
  }
#ifdef debug_test_create_xgrid 
  printf("anglesum is %g\n", anglesum);
#endif
  if( fabs(anglesum - 2*M_PI ) < EPSLN10 )
    return 1;
  else if( anglesum > 0 ) {/* The point maybe a vertice */
    for (i=0;i<npts;i++) {
      if( fabs(x2[i]-x1) < EPSLN10 && fabs(y2[i]-y1) < EPSLN10 && fabs(z2[i]-z1) < EPSLN10 )
	return 1;
    }
    return 0;
  }
  else
    return 0;

}
