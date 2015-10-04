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

/***********************************************************************
                      mosaic_util.h
    This header file provide some utilities routine that will be used in many tools.
    
    contact: GFDL.Climate.Model.Info@noaa.gov
***********************************************************************/
#ifndef MOSAIC_UTIL_H_
#define MOSAIC_UTIL_H_

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)
#define SMALL_VALUE ( 1.e-10 )
struct Node{
  double x, y, z, u;
  int intersect; /* indicate if this point is an intersection, 0 = no, 1= yes, 2=both intersect and vertices */ 
  int inbound;      /* -1 uninitialized, 0 coincident, 1 outbound, 2 inbound */
  int initialized; /* = 0 means empty list */
  int isInside;   /* = 1 means one point is inside the other polygon */
  struct Node *Next;
};


void error_handler(const char *msg);
int nearest_index(double value, const double *array, int ia);
int lon_fix(double *x, double *y, int n_in, double tlon);
double minval_double(int size, const double *data);
double maxval_double(int size, const double *data);
double avgval_double(int size, const double *data);
void latlon2xyz(int size, const double *lon, const double *lat, double *x, double *y, double *z); 
void xyz2latlon(int size, const double *x, const double *y, const double *z, double *lon, double *lat);
double box_area(double ll_lon, double ll_lat, double ur_lon, double ur_lat);
double poly_area(const double lon[], const double lat[], int n);
double poly_area_dimensionless(const double lon[], const double lat[], int n);
double poly_area_no_adjust(const double x[], const double y[], int n);
int fix_lon(double lon[], double lat[], int n, double tlon);
void tokenize(const char * const string, const char *tokens, unsigned int varlen,
	      unsigned int maxvar, char * pstring, unsigned int * const nstr);
double great_circle_distance(double *p1, double *p2);
double spherical_excess_area(const double* p_ll, const double* p_ul,
			     const double* p_lr, const double* p_ur, double radius);
void vect_cross(const double *p1, const double *p2, double *e );
double spherical_angle(const double *v1, const double *v2, const double *v3);
void normalize_vect(double *e);
void unit_vect_latlon(int size, const double *lon, const double *lat, double *vlon, double *vlat);
double great_circle_area(int n, const double *x, const double *y, const double *z);
double * cross(const double *p1, const double *p2);
double dot(const double *p1, const double *p2);
int intersect_tri_with_line(const double *plane, const double *l1, const double *l2, double *p,
			     double *t);
int invert_matrix_3x3(double m[], double m_inv[]);
void mult(double m[], double v[], double out_v[]);

double gridArea(struct Node *grid);

void addBeg(struct Node *list, double x, double y, double z, int intersect, double u, int inbound);

void addEnd(struct Node *list, double x, double y, double z, int intersect, double u, int inbound);

int getInbound( struct Node node );
int length(struct Node *list);
struct Node *getNextNode(struct Node *list);
struct Node *getNode(struct Node *list, struct Node inNode);
int sameNode(struct Node node1, struct Node node2);
void deleteNode(struct Node *list);
void initNode(struct Node *node);
void copyNode(struct Node *node_out, struct Node node_in);
void assignNode(struct Node *node, double x, double y, double z, int intersect, double u, int inbound);
void removeNode(struct Node *list, struct Node nodeIn);
void resetIntersectValue(struct Node *list );
void addNode(struct Node *list, struct Node nodeIn);
void printNode(struct Node *list, char *str);
void insertAfter(struct Node *list, double x, double y, double z, int intersect, double u, int inbound,
		 double x2, double y2, double z2);
int isIntersect(struct Node node);
int getFirstInbound( struct Node *list, struct Node *nodeOut);
void getCoordinate(struct Node node, double *x, double *y, double *z);
double metric(const double *p);
int insidePolygon( double x1, double y1, double z1, int npts, double *x2, double *y2, double *z2);
struct Node *getLast(struct Node *list);
void setInbound(struct Node *interList, struct Node *list);
int isInside(struct Node *node);
#endif
