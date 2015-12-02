#include "stdafx.h"
#include "cgrid.h"

int CGrid::SPTR_BUF = 1024;
int CGrid::STRU_BUF = 4096;

////////////////////////////////////////////////
//...auxilliary procedure -- zero of grid nodes;
void CGrid::zero_grid()
{
   delete_struct(hit);
   delete_struct(geom);
   delete_struct(geom_ptr);
   delete_struct(cond);
   delete_struct(cond_ptr);
   delete_struct(X);
   delete_struct(Y);
   delete_struct(Z);
   delete_struct(nX);
   delete_struct(nY);
   delete_struct(nZ);
   N         = 0; 
   N1        = 0; 
   N2        = 0; 
   buf_count = 0;
   buf_X     = 0;
   buf_Y     = 0;
   buf_Z     = 0;
   for (; pp && N_par > 0; N_par--) delete_struct(pp[N_par-1]);
   delete_struct(pp);
}

/////////////////////////////////////
//...bufferization of the grid nodes;
int  CGrid::bufferization (int buf_size, Num_State mm)
{
   int m = 1, i;
   int * new_hit = (int *)new_struct((N+buf_size)*sizeof(int));
   if (mm == OK_STATE && hit) memcpy(new_hit, hit, N*sizeof(int));
   if (! new_hit) m = 0; else {
       delete_struct(hit); hit = new_hit;
   }
   double * new_X = (double *)new_struct((N+buf_size)*sizeof(double));
   if (mm == OK_STATE && X) memcpy(new_X, X, N*sizeof(double));
   if (!    new_X) m = 0; else {
       delete_struct(X); X = new_X;
   }
   double * new_Y = (double *)new_struct((N+buf_size)*sizeof(double));
   if (mm == OK_STATE && Y) memcpy(new_Y, Y, N*sizeof(double));
   if (!    new_Y) m = 0; else {
       delete_struct(Y); Y = new_Y;
   }
   double * new_Z = (double *)new_struct((N+buf_size)*sizeof(double));
   if (mm == OK_STATE && Z) memcpy(new_Z, Z, N*sizeof(double));
   if (!    new_Z) m = 0; else {
       delete_struct(Z); Z = new_Z;
   }
   double * new_nX = (double *)new_struct((N+buf_size)*sizeof(double));
   if (mm == OK_STATE && nX) memcpy(new_nX, nX, N*sizeof(double));
   if (!    new_nX) m = 0; else {
       delete_struct(nX); nX = new_nX;
   }
   double * new_nY = (double *)new_struct((N+buf_size)*sizeof(double));
   if (mm == OK_STATE && nY) memcpy(new_nY, nY, N*sizeof(double));
   if (!    new_nY) m = 0; else {
       delete_struct(nY); nY = new_nY;
   }
   double * new_nZ = (double *)new_struct((N+buf_size)*sizeof(double));
   if (mm == OK_STATE && nZ) memcpy(new_nZ, nZ, N*sizeof(double));
   if (!    new_nZ) m = 0; else {
       delete_struct(nZ); nZ = new_nZ;
   }
   for (i = 0; pp && i < N_par; i++) {
       double * new_pp = (double *)new_struct((N+buf_size)*sizeof(double));
       if (mm == OK_STATE && pp[i]) memcpy(new_pp, pp[i], N*sizeof(double));
       if (!    new_pp) m = 0; else {
           delete_struct(pp[i]); pp[i] = new_pp;
       }
   }
   if (m) buf_count = buf_size;
   return(m);
}

//////////////////////////////////////////////////////
//...bufferization of the X-array for multilevel grid;
int  CGrid::bufferizat_X (int buf_size, Num_State mm)
{
   int m = 1;
   double * new_X = (double *)new_struct((N+buf_size)*sizeof(double));
   if (mm == OK_STATE && X) memcpy(new_X, X, N*sizeof(double));
   if (!    new_X) m = 0; else {
       delete_struct(X); X = new_X;
   }
   if (m) buf_X = buf_size;
   return(m);
}

//////////////////////////////////////////////////////
//...bufferization of the Y-array for multilevel grid;
int  CGrid::bufferizat_Y (int buf_size, Num_State mm)
{
   int m = 1;
   double * new_Y = (double *)new_struct((N1+buf_size)*sizeof(double));
   if (mm == OK_STATE && Y) memcpy(new_Y, Y, N1*sizeof(double));
   if (!    new_Y) m = 0; else {
       delete_struct(Y); Y = new_Y;
   }
   if (m) buf_Y = buf_size;
   return(m);
}

//////////////////////////////////////////////////////
//...bufferization of the Z-array for multilevel grid;
int  CGrid::bufferizat_Z (int buf_size, Num_State mm)
{
   int m = 1;
   double * new_Z = (double *)new_struct((N2+buf_size)*sizeof(double));
   if (mm == OK_STATE && Z) memcpy(new_Z, Z, N2*sizeof(double));
   if (!    new_Z) m = 0; else {
       delete_struct(Z); Z = new_Z;
   }
   if (m) buf_Z = buf_size;
   return(m);
}

//////////////////////////////////////////////////
//...adding with bufferization of the grid values;
int CGrid::add_new_point(double X, double Y, double Z, double nX, double nY, double nZ, double * pp)
{
   if (buf_count || bufferization(N_buf)) {
       CGrid::X[N]     = X;
       CGrid::Y[N]     = Y;
       CGrid::Z[N]     = Z;
       CGrid::nX[N]    = nX;
       CGrid::nY[N]    = nY;
       CGrid::nZ[N]    = nZ;
       CGrid::hit[N++] = -1; set_param(pp); buf_count--;
       return(1);
   }
   return(0);
}

///////////////////////////////////////////////////////////////////
//...adding with bufferization of the X-values for multilevel grid;
int CGrid::add_new_point_X(double val, double eps)
{
   if (! buf_X && ! bufferizat_X(N_buf)) return(0);
   if (eps) {
      int    m = 0, j;
      while (m < N && X[m] < val) m++;
      if ((m == N || fabs(X[m] - val) > eps*(1.+fabs(X[m]))) &&
          (m == 0 || fabs(X[m-1]-val) > eps*(1.+fabs(X[m-1])))) {
          for (j = N; j > m; j--) X[j] = X[j-1];
          X[m] = val;
          N++; buf_X--;
      }
   }
   else {
      X[N++] = val; buf_X--;
   }
   return(1);
}

///////////////////////////////////////////////////////////////////
//...adding with bufferization of the Y-values for multilevel grid;
int CGrid::add_new_point_Y(double val, double eps)
{
   if (! buf_Y && ! bufferizat_Y(N_buf)) return(0);
   if (eps) {
      int    m = 0, j;
      while (m < N1 && Y[m] < val) m++;
      if ((m == N1 || fabs(Y[m] - val) > eps*(1.+fabs(Y[m]))) &&
          (m == 0  || fabs(Y[m-1]-val) > eps*(1.+fabs(Y[m-1])))) {
          for (j = N1; j > m; j--) Y[j] = Y[j-1];
          Y[m] = val;
          N1++; buf_Y--;
      }
   }
   else {
      Y[N1++] = val; buf_Y--;
   }
   return(1);
}

///////////////////////////////////////////////////////////////////
//...adding with bufferization of the Z-values for multilevel grid;
int CGrid::add_new_point_Z(double val, double eps)
{
   if (! buf_Z && ! bufferizat_Z(N_buf)) return(0);
   if (eps) {
      int    m = 0, j;
      while (m < N2 && Z[m] < val) m++;
      if ((m == N2 || fabs(Z[m] - val) > eps*(1.+fabs(Z[m]))) &&
          (m == 0  || fabs(Z[m-1]-val) > eps*(1.+fabs(Z[m-1])))) {
          for (j = N2; j > m; j--) Z[j] = Z[j-1];
          Z[m] = val;
          N2++; buf_Z--;
      }
   }
   else {
      Z[N2++] = val; buf_Z--;
   }
   return(1);
}

////////////////////////////////////////////////////////////////
//...adding with bufferization of the X-values for grid lattice;
int CGrid::add_new_lattice_X(double val, double eps)
{
   if (! buf_X && ! bufferizat_X(N_buf)) return(0);
   int    m = 0, j;
   while (m < N && X[m] < val) m++;
   if (m != N && fabs(X[m] - val) <= eps) X[m] = (X[m]+val)/2.; else   
	if (m != 0 && fabs(X[m-1]-val) <= eps) X[m-1] = (X[m-1]+val)/2.; 
	else {
		 for (j = N; j > m; j--) X[j] = X[j-1];
		 X[m] = val;
		 N++; buf_X--;
   }
   return(1);
}

/////////////////////////////////////////////////////////////////
//...adding with bufferization of the Y-values for grid lattice;
int CGrid::add_new_lattice_Y(double val, double eps)
{
   if (! buf_Y && ! bufferizat_Y(N_buf)) return(0);
   int    m = 0, j;
   while (m < N1 && Y[m] < val) m++;
   if (m != N1 && fabs(Y[m] - val) <= eps) Y[m] = (Y[m]+val)/2.; else   
	if (m != 0  && fabs(Y[m-1]-val) <= eps) Y[m-1] = (Y[m-1]+val)/2.; 
	else {
		for (j = N1; j > m; j--) Y[j] = Y[j-1];
		Y[m] = val;
		N1++; buf_Y--;
   }
   return(1);
}

////////////////////////////////////////////////////////////////
//...adding with bufferization of the Z-values for grid lattice;
int CGrid::add_new_lattice_Z(double val, double eps)
{
   if (! buf_Z && ! bufferizat_Z(N_buf)) return(0);
   int    m = 0, j;
   while (m < N2 && Z[m] < val) m++;
   if (m != N2 && fabs(Z[m] - val) <= eps) Z[m] = (Z[m]+val)/2.; else   
	if (m != 0  && fabs(Z[m-1]-val) <= eps) Z[m-1] = (Z[m-1]+val)/2.; 
	else {
		for (j = N2; j > m; j--) Z[j] = Z[j-1];
		Z[m] = val;
		N2++; buf_Z--;
   }
   return(1);
}

///////////////////////////////////////////////////
//...adding new table of parameters for grid nodes;
void CGrid::add_params(int m, int id_bufferizat)
{
   if (m > 0) {
       double ** new_pp = (double **)new_struct((m+N_par)*sizeof(double *));
       if (pp) memcpy(new_pp, pp, N_par*sizeof(double *));
       if (new_pp) {
           delete_struct(pp); pp = new_pp; N_par += m;
       }
       if (id_bufferizat) bufferization(N_buf);
   }
}

/////////////////////////////////
//...setting value of parameters;
void CGrid::set_param(double * pp)
{
   if  (CGrid::pp && pp) 
   for (int m = N_par-1; m >= 0; m--) 
        set_param(m, N-1, pp[m]);
}

/////////////////////////////////////////////////////////
//           OPERATIONS WITH NODES TOPOLOGY            //
/////////////////////////////////////////////////////////
//////////////////////////////////////////////
//...��������� ������������� ��������� ������;
int CGrid::stru_install(int k, int mm_facet[6][5]) 
{
	int mm[8], m = 0, i, j;
	if  (geom && geom_ptr)
	for (m = geom[(i = geom_ptr[k])+1]-3, j = 0; j < m; j++) mm[j] = geom[i+5+j];
	switch (geom[i]) {
		case GL_BOXS: {
			m = 6;
			mm_facet[0][0] = 4; mm_facet[0][1] = mm[0]; mm_facet[0][2] = mm[4]; mm_facet[0][3] = mm[7]; mm_facet[0][4] = mm[3];
			mm_facet[1][0] = 4; mm_facet[1][1] = mm[1]; mm_facet[1][2] = mm[2]; mm_facet[1][3] = mm[6]; mm_facet[1][4] = mm[5];
			mm_facet[2][0] = 4; mm_facet[2][1] = mm[0]; mm_facet[2][2] = mm[1]; mm_facet[2][3] = mm[5]; mm_facet[2][4] = mm[4];
			mm_facet[3][0] = 4; mm_facet[3][1] = mm[3]; mm_facet[3][2] = mm[7]; mm_facet[3][3] = mm[6]; mm_facet[3][4] = mm[2];
			mm_facet[4][0] = 4; mm_facet[4][1] = mm[0]; mm_facet[4][2] = mm[3]; mm_facet[4][3] = mm[2]; mm_facet[4][4] = mm[1];
			mm_facet[5][0] = 4; mm_facet[5][1] = mm[4]; mm_facet[5][2] = mm[5]; mm_facet[5][3] = mm[6]; mm_facet[5][4] = mm[7];
		}  break;
		case GL_PENTA: {
			m = 5;
			mm_facet[0][0] = 4; mm_facet[0][1] = mm[0]; mm_facet[0][2] = mm[1]; mm_facet[0][3] = mm[4]; mm_facet[0][4] = mm[3];
			mm_facet[1][0] = 4; mm_facet[1][1] = mm[1]; mm_facet[1][2] = mm[2]; mm_facet[1][3] = mm[5]; mm_facet[1][4] = mm[4];
			mm_facet[2][0] = 4; mm_facet[2][1] = mm[2]; mm_facet[2][2] = mm[0]; mm_facet[2][3] = mm[3]; mm_facet[2][4] = mm[5];
			mm_facet[3][0] = 3; mm_facet[3][1] = mm[0]; mm_facet[3][2] = mm[2]; mm_facet[3][3] = mm[1];
			mm_facet[4][0] = 3; mm_facet[4][1] = mm[3]; mm_facet[4][2] = mm[4]; mm_facet[4][3] = mm[5];
		}	break;
		case GL_TETRA: {
			m = 4;
			mm_facet[0][0] = 3; mm_facet[0][1] = mm[0]; mm_facet[0][2] = mm[2]; mm_facet[0][3] = mm[1];
			mm_facet[1][0] = 3; mm_facet[1][1] = mm[1]; mm_facet[1][2] = mm[2]; mm_facet[1][3] = mm[3];
			mm_facet[2][0] = 3; mm_facet[2][1] = mm[2]; mm_facet[2][2] = mm[0]; mm_facet[2][3] = mm[3];
			mm_facet[3][0] = 3; mm_facet[3][1] = mm[0]; mm_facet[3][2] = mm[1]; mm_facet[3][3] = mm[3];
		}  break;
		case GL_PYRAMID: {
			m = 5;
			mm_facet[0][0] = 3; mm_facet[0][1] = mm[0]; mm_facet[0][2] = mm[1]; mm_facet[0][3] = mm[4];
			mm_facet[1][0] = 3; mm_facet[1][1] = mm[1]; mm_facet[1][2] = mm[2]; mm_facet[1][3] = mm[4];
			mm_facet[2][0] = 3; mm_facet[2][1] = mm[2]; mm_facet[2][2] = mm[3]; mm_facet[2][3] = mm[4];
			mm_facet[3][0] = 3; mm_facet[3][1] = mm[3]; mm_facet[3][2] = mm[0]; mm_facet[3][3] = mm[4];
			mm_facet[4][0] = 4; mm_facet[4][1] = mm[0]; mm_facet[4][2] = mm[3]; mm_facet[4][3] = mm[2]; mm_facet[4][4] = mm[1];
		}  break;
		default: m = 0;
	}
	return(m);
}

///////////////////////////////////////////////////////////////////////////////////
//...��������������� �������, ������������ ����� ������ ������ ������ ������ �����;
int CGrid::link_number(Topo * link, int id_arc)
{
  int k = 0, m = id_arc;
  if (link) {
      while (0 <= id_arc && ++id_arc <= link[0] && 0 <= link[id_arc]) k++;

      if (id_arc == link[0]+1) {
          id_arc = 0;
          while (++id_arc <= m && 0 <= link[id_arc]) k++;
      }
  }
  return(k);
}

/////////////////////////////////////////////////////////
//...��������� �� ��������� �������� � ������� ���������;
int CGrid::geom_element(int k)
{
  int i = 0, l;
  if (geom && geom[0] >= k && k >= 0) {
      if (geom_ptr)  i = geom_ptr[k]; else
      for (i = 1, l = 0; l < k; l++) i += geom[++i]+1;
  }
  return(i);
}

/////////////////////////////////////////
//...����������� �������������� ��������;
void CGrid::quick_sort(int first, int last, int * index)
{
   int i = first, j = last, median = index[(first+last)/2];
   do {
      while (index[i] < median) i++;
      while (index[j] > median) j--;
      if (i <= j) {
			swap(index[i], index[j]); i++; j--;
      }
   } 
	while (i <= j);

   if (first < j) quick_sort(first, j, index);
   if (i  < last) quick_sort(i,  last, index);
}

void CGrid::quick_sort_inverse(int first, int last, int * index)
{
   int i = first, j = last, median = index[(first+last)/2];
   do {
      while (index[i] > median) i++;
      while (index[j] < median) j--;
      if (i <= j) {
			swap(index[i], index[j]); i++; j--;
      }
   } 
	while (i <= j);

   if (first < j) quick_sort_inverse(first, j, index);
   if (i  < last) quick_sort_inverse(i,  last, index);
}

int CGrid::binary_search(int first, int last, int elem, int * index)
{
   int median;
   do if (elem <  index[median = (first+last)/2]) last = median-1; else	first = median+1;
	while (elem != index[median] && first <= last);
   return(elem == index[median] ? median : -1);
}

int CGrid::binary_search_inverse(int first, int last, int elem, int * index)
{
   int median;
   do if (elem >  index[median = (first+last)/2]) last = median-1; else	first = median+1;
	while (elem != index[median] && first <= last);
   return(elem == index[median] ? median : -1);
}

///////////////////////////////////////////////////
//...������������ ����� ����� �� ������� ���������;
void CGrid::grid_add(CGrid * nd)
{
	if (bufferization(nd->N)) {
		memcpy(hit+N, nd->hit,  nd->N*sizeof(int));
		memcpy( X+N,  nd->X, nd->N*sizeof(double));
		memcpy( Y+N,  nd->Y, nd->N*sizeof(double));
		memcpy( Z+N,  nd->Z, nd->N*sizeof(double));
		memcpy(nX+N, nd->nX, nd->N*sizeof(double));
		memcpy(nY+N, nd->nY, nd->N*sizeof(double));
		memcpy(nZ+N, nd->nZ, nd->N*sizeof(double));
		if (  pp && nd->pp && N_par == nd->N_par)
		for (int i = 0; pp && i < N_par; i++) {
			memcpy(pp[i]+N, nd->pp[i], N*sizeof(double));
		}
		buf_count -= nd->N; N += nd->N; nd->zero_grid();
	}
}

//////////////////////////////////////////////////////////////////
//...auxilliary function for adding information to the grid nodes;
int CGrid::grid_add(double *& X0,  double *& Y0,double *& Z0, double *& nX0, double *& nY0, 
                    double *& nZ0, int N0, int *& gm, int id_block)
{
	if (N0 <= 0) return(1);

///////////////////////////////////////////
//...������������ ������� ����� ��p�������;
	double * new_X  = X0 ? (double *)new_struct((N+N0)*sizeof(double)) : NULL,
			* new_Y  = Y0 ? (double *)new_struct((N+N0)*sizeof(double)) : NULL,
			* new_Z  = Z0 ? (double *)new_struct((N+N0)*sizeof(double)) : NULL,
			* new_nX = nX0 ? (double *)new_struct((N+N0)*sizeof(double)) : NULL,
			* new_nY = nY0 ? (double *)new_struct((N+N0)*sizeof(double)) : NULL,
			* new_nZ = nZ0 ? (double *)new_struct((N+N0)*sizeof(double)) : NULL;
	int   * new_hit = hit ? (int *)new_struct((N+N0)*sizeof(int)) : NULL, 
		   * new_gm, k, j, l, shift = 2;

///////////////////////////////////////////////////////////////////////////////
//...������������ ����� �������, ������������ ���������, � ������������ ������;
	if (! geom || ! gm) new_gm = NULL; 
	else {
		for (k = j = 0; k < geom[0]; k++) j += geom[(j += shift)];
		for (k = l = 0; k <   gm[0]; k++) l +=   gm[(l += shift)];
		new_gm  = (int *)new_struct((j+l+1)*sizeof(int));
	}

//////////////////////////////////////////////////////
//...��������� �������� ����� ����� � ����� ���������;
	if (X && new_X) memcpy(new_X, X, N*sizeof(double));
	if (Y && new_Y) memcpy(new_Y, Y, N*sizeof(double));
	if (Z && new_Z) memcpy(new_Z, Z, N*sizeof(double));
	if (nX && new_nX) memcpy(new_nX, nX, N*sizeof(double));
	if (nY && new_nY) memcpy(new_nY, nY, N*sizeof(double));
	if (nZ && new_nZ) memcpy(new_nZ, nZ, N*sizeof(double));
	if (hit && new_hit) memcpy(new_hit, hit, N*sizeof(int));

	if (X0 && new_X) memcpy(new_X+N, X0, N0*sizeof(double));
	if (Y0 && new_Y) memcpy(new_Y+N, Y0, N0*sizeof(double));
	if (Z0 && new_Z) memcpy(new_Z+N, Z0, N0*sizeof(double));
	if (nX0 && new_nX) memcpy(new_nX+N, nX0, N0*sizeof(double));
	if (nY0 && new_nY) memcpy(new_nY+N, nY0, N0*sizeof(double));
	if (nZ0 && new_nZ) memcpy(new_nZ+N, nZ0, N0*sizeof(double));
	if (new_hit) for (k = 0; k < N0; k++) new_hit[k+N] = id_block;

/////////////////////////////////////////////////
//...��������� �������� ��������� � ����� ������;
	if (new_gm) {
		memcpy(new_gm, geom, (j+1)*sizeof(int));
		memcpy(new_gm+j+1, gm+1, l*sizeof(int));
		new_gm[0] += gm[0];
	}

//////////////////////
//...��������� ������;
	if (! new_X && X && X0 || 
		! new_Y && Y && Y0 || 
		! new_Z && Z && Z0 ||
		! new_nX && nX && nX0 ||
		! new_nY && nY && nY0 ||
		! new_nZ && nZ && nZ0 ||
		! new_gm && geom && gm || 
		! new_hit && hit) {
		delete_struct(X);  delete_struct(new_X);
		delete_struct(Y);  delete_struct(new_Y);
		delete_struct(Z);  delete_struct(new_Z);
		delete_struct(nX); delete_struct(new_nX);
		delete_struct(nY); delete_struct(new_nY);
		delete_struct(nZ); delete_struct(new_nZ); delete_struct(new_hit);
		delete_struct(gm); delete_struct(new_gm);
		return(0);
	}
	else {
		delete_struct(X0);  delete_struct(X);
		delete_struct(Y0);  delete_struct(Y);
		delete_struct(Z0);  delete_struct(Z);
		delete_struct(nX0); delete_struct(nX);
		delete_struct(nY0); delete_struct(nY);
		delete_struct(nZ0); delete_struct(nZ); delete_struct(hit);
		delete_struct(gm);  delete_struct(geom);
	}
//////////////////////////////////////////////////////////////
//...������� ����� ������� � ��������� � ������� �� ���������;
	N   += N0;
	X    = new_X;
	Y    = new_Y;
	Z    = new_Z;
	nX   = new_nX;
	nY   = new_nY;
	nZ   = new_nZ;
	hit  = new_hit;
	geom = new_gm;

	return(1);
}

//////////////////////////////////////////////
//...isometrical transformation of grid nodes;
void CGrid::grid_iso(int k0, int k1, double * P0,
                     double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX)
{
  if (0 <= k0 && k1 <= N) {
      double P[3];
      for (int k = k0; k < k1; k++) {
           P[0] = X[k];
           P[1] = Y[k];
           P[2] = Z[k];  point_iso(P, P0, CZ, SZ, CY, SY, CX, SX);
           X[k] = P[0];
           Y[k] = P[1];
           Z[k] = P[2];
           P[0] = nX[k];
           P[1] = nY[k];
           P[2] = nZ[k]; point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
           nX[k] = P[0];
           nY[k] = P[1];
           nZ[k] = P[2];
      }
  }
}

/////////////////////////////////////////////////
//...auxilliary function for triangle grid nodes;
void CGrid::grid_tria_3D_refine(int N_elem, double * P, double * pp, int id_block)
{
  double d, px1, px2, py1, py2, pz1, pz2;
  int    i, j, k, l, m, shift = N, * new_geom;
//////////////////
//...size of geom;
  for (l  = 1, i = 0; geom && i < geom[0]; i++) 
       l += geom[++l]+1;

  if (N_elem > 0 && P && bufferization(j = (N_elem+1)*(N_elem+2)/2) &&
     (new_geom = (int *)new_struct((m = --l+N_elem*(N_elem+8)-3)*sizeof(int))) != NULL) {
      if (geom) 
          memcpy(new_geom,  geom, (l+1)*sizeof(int)); 
      delete_struct (geom); geom = new_geom; N += j;

/////////////////////////////////////////////////////////
//...shifting vector homothetic triangles in subdivision;
      px1 = (P[3]-P[0])*(d = 1./N_elem); px2 = (P[6]-P[3])*d;
      py1 = (P[4]-P[1])* d;              py2 = (P[7]-P[4])*d;
      pz1 = (P[5]-P[2])* d;              pz2 = (P[8]-P[5])*d;

//////////////////////////////////////////////////////
//...filling points -- apices of homothetic triangles;
      for ( X [j = shift] = P[0],  Y[j] = P[1],   Z[j] = P[2],
           nX [j]         = P[9], nY[j] = P[10], nZ[j] = P[11],
           hit[j]         = -1, geom[0] += 2*N_elem-1, i = 1; i <= N_elem; i++) {
           geom[++l] = GL_TRIANGLES;
           geom[++l] = 3;
           geom[++l] = j;
           geom[++l] = j+1;
           geom[++l] = (m = j+N_elem-i+1)+1;
           X[j+1]    = X[j]+px1;
           Y[j+1]    = Y[j]+py1;
           Z[j+1]    = Z[j]+pz1;
           X[m+1]    = X[j+1]+px2;
           Y[m+1]    = Y[j+1]+py2;
           Z[m+1]    = Z[j+1]+pz2;
           nX[j+1]   = nX[m+1] = P[9];
           nY[j+1]   = nY[m+1] = P[10];
           nZ[j+1]   = nZ[m+1] = P[11];
           hit[j+1]  = id_block;
           hit[m+1]  = id_block;
           if (i < N_elem) {
               geom[++l] = GL_QUAD_STRIP;
               geom[++l] = (N_elem-i+1)*2;
               geom[++l] = ++j;
               geom[++l] = ++m;
               for (k = i; k < N_elem; k++) {
                    X[j+1]    = X[j]+px1;
                    Y[j+1]    = Y[j]+py1;
                    Z[j+1]    = Z[j]+pz1;
                    X[m+1]    = X[j+1]+px2;
                    Y[m+1]    = Y[j+1]+py2;
                    Z[m+1]    = Z[j+1]+pz2;
                    nX[j+1]   = nX[m+1] = P[9];
                    nY[j+1]   = nY[m+1] = P[10];
                    nZ[j+1]   = nZ[m+1] = P[11];
                    hit[j+1]  = id_block;
                    hit[m+1]  = id_block;
                    geom[++l] = ++j;
                    geom[++l] = ++m;
               }
               ++j; 
           }
      }

////////////////////////
//...parameters setting;
		if (pp)
      for (i = shift;   i <  N; i++) 
      for (m = N_par-1; m >= 0; m--) set_param(m, i, pp[m]);

      buf_count = 0; return;
  }
//zero_grid();
}

////////////////////////////////////////
//...auxilliary function for grid nodes;
void CGrid::grid_quad_3D_refine(int N_elem1, int N_elem2, double * P, double * pp, int id_block)
{//...������� ������� ���������� ����� �������
  double d, px1, px2, px3, py1, py2, py3, pz1, pz2, pz3, pnt[6];
  int    i, j, k, l, m, shift = N, * new_geom;
//////////////////
//...size of geom;
  for (l  = 1, i = 0; geom && i < geom[0]; i++) 
       l += geom[++l]+1;

  if (N_elem1 > 0 && N_elem2 > 0 && P && bufferization(j = (N_elem1+1)*(N_elem2+1)) &&
     (new_geom = (int *)new_struct((m = --l+N_elem1*(N_elem2+2)*2+1)*sizeof(int))) != NULL) {
      if (geom) 
          memcpy(new_geom,  geom, (l+1)*sizeof(int)); 
      delete_struct (geom); geom = new_geom; N += j;

///////////////////////////////////////////////////////////
//...shifting vector homothetic quadrangles in subdivision;
      px1 = (P[9] -P[0])*(d = 1./N_elem2); px2 = (P[6]-P[3])*d;
      py1 = (P[10]-P[1])* d;               py2 = (P[7]-P[4])*d;
      pz1 = (P[11]-P[2])* d;               pz2 = (P[8]-P[5])*d;

////////////////////////////////////////////////////////
//...filling points -- apices of homothetic quadrangles;
      pnt[0] = P[0];
      pnt[1] = P[1];
      pnt[2] = P[2];
      pnt[3] = P[3];
      pnt[4] = P[4];
      pnt[5] = P[5];
      for (j = shift, i = 0; i <= N_elem2; i++) {
          px3 = (pnt[3]-(X[j] = pnt[0]))*(d = 1./N_elem1); 
          py3 = (pnt[4]-(Y[j] = pnt[1]))* d;               
          pz3 = (pnt[5]-(Z[j] = pnt[2]))* d;               
          nX [j] = P[12]; 
          nY [j] = P[13]; 
          nZ [j] = P[14];
          hit[j] = id_block; 
          for (k = 0; k < N_elem1; k++) { 
               X[j+1]   = X[j]+px3;
               Y[j+1]   = Y[j]+py3;
               Z[j+1]   = Z[j]+pz3;
               nX[j+1]  = P[12];
               nY[j+1]  = P[13];
               nZ[j+1]  = P[14];
               hit[j+1] = id_block;
               j++;
          }
          pnt[0] += px1; pnt[3] += px2;
          pnt[1] += py1; pnt[4] += py2;
          pnt[2] += pz1; pnt[5] += pz2;
          j++;
      }
      for (geom[0]  += N_elem2, j = shift, i = 1; i <= N_elem2; i++) {
           geom[++l] = GL_QUAD_STRIP;
           geom[++l] = (N_elem1+1)*2;
           geom[++l] = j;
           geom[++l] = m = j+N_elem1+1;
           for (k = 0; k < N_elem1; k++) {
                geom[++l] = ++j;
                geom[++l] = ++m;
           }
           ++j;
      }

////////////////////////
//...parameters setting;
		if (pp)
      for (i = shift;   i <  N; i++) 
      for (m = N_par-1; m >= 0; m--) set_param(m, i, pp[m]);

      buf_count = 0;
  }
}

//////////////////////////////////////////////////////////
//...auxilliary function for grid nodes in infitity block;
void CGrid::grid_quad_3D_refine_spec(int N_elem1, int N_elem2, double * P, double * pp, int id_block)
{
  double d, px1, px2, px3, py1, py2, py3, pz1, pz2, pz3, pnt[6];
  int    i, j, k, l, m, shift = N, * new_geom;
//////////////////
//...size of geom;
  for (l  = 1, i = 0; geom && i < geom[0]; i++) 
       l += geom[++l]+1;

  if (N_elem1 > 0 && N_elem2 > 0 && P && bufferization(j = (N_elem1+1)*(N_elem2+1)) &&
     (new_geom = (int *)new_struct((m = --l+N_elem1*(N_elem2+2)*2+1)*sizeof(int))) != NULL) {
      if (geom) 
          memcpy(new_geom,  geom, (l+1)*sizeof(int)); 
      delete_struct (geom); geom = new_geom; N += j;

///////////////////////////////////////////////////////////
//...shifting vector homothetic quadrangles in subdivision;
      px1 = (P[3]-P[6])*(d = 1./N_elem2); px2 = (P[0]-P[9 ])*d;
      py1 = (P[4]-P[7])* d;               py2 = (P[1]-P[10])*d;
      pz1 = (P[5]-P[8])* d;               pz2 = (P[2]-P[11])*d;

      double homotetic = P[15], //...coefficient of homotetic expansion;
             L_inf     = P[16]; //...initial step for homotetic expansion;

		if (L_inf == 0.) N_elem1 = N_elem2 = 1;

////////////////////////////////////////////////////////
//...filling points -- apices of homothetic quadrangles;
      pnt[0] = P[6];
      pnt[1] = P[7];
      pnt[2] = P[8];
      pnt[3] = P[9];
      pnt[4] = P[10];
      pnt[5] = P[11];
      for (j = shift, i = 0; i <= N_elem2; i++) {
          px3 = (pnt[3]-(X[j] = pnt[0]))*d; 
          py3 = (pnt[4]-(Y[j] = pnt[1]))*d;               
          pz3 = (pnt[5]-(Z[j] = pnt[2]))*d; 
          if (L_inf > 0.) {
              double f = L_inf/sqrt(px3*px3+py3*py3+pz3*pz3);
              px3 *= f; 
              py3 *= f; 
              pz3 *= f; 
          }
          nX [j] = P[12]; 
          nY [j] = P[13]; 
          nZ [j] = P[14];
          hit[j] = id_block; 
          for (k = 0; k < N_elem1; k++) { 
               X[j+1]   = X[j]+px3; px3 += homotetic*px3;
               Y[j+1]   = Y[j]+py3; py3 += homotetic*py3;
               Z[j+1]   = Z[j]+pz3; pz3 += homotetic*pz3;
               nX[j+1]  = P[12];
               nY[j+1]  = P[13];
               nZ[j+1]  = P[14];
               hit[j+1] = id_block;
               j++;
          }
          pnt[0] += px1; pnt[3] += px2;
          pnt[1] += py1; pnt[4] += py2;
          pnt[2] += pz1; pnt[5] += pz2;
          j++;
      }
      for (geom[0]  += N_elem2, j = shift, i = 1; i <= N_elem2; i++) {
           geom[++l] = GL_QUAD_STRIP;
           geom[++l] = (N_elem1+1)*2;
           geom[++l] = j;
           geom[++l] = m = j+N_elem1+1;
           for (k = 0; k < N_elem1; k++) {
                geom[++l] = ++j;
                geom[++l] = ++m;
           }
           ++j;
      }

////////////////////////
//...parameters setting;
		if (pp)
      for (i = shift;   i <  N; i++) 
      for (m = N_par-1; m >= 0; m--) set_param(m, i, pp[m]);

      buf_count = 0;
  }
}

/////////////////////////////////////////////////////////////////////
//...construction grid lattice for arbitrary structure of inclusions;
void CGrid::grid_lattice(CGrid * nd, double * par, double eps)
{
	if (nd && nd->N) {
		int k, j, mX, mY, mZ;
		for (par[6] = par[7] = 0., k = 0; k < nd->N; k++) {
			add_new_lattice_X(nd->X[k], eps);
			add_new_lattice_Y(nd->Y[k], eps);
			add_new_lattice_Z(nd->Z[k], eps);
			par[6] += 4./3.*M_PI*nd->get_param(0, k)*sqr(nd->get_param(0, k));
			par[7] += 4./3.*M_PI*nd->get_param(1, k)*sqr(nd->get_param(1, k));
		}	par[6] /= (par[1]-par[0])*(par[3]-par[2])*(par[5]-par[4]);
			par[7] /= (par[1]-par[0])*(par[3]-par[2])*(par[5]-par[4]); par[7] -= par[6];

///////////////////////////////////////////
//...����������� ������� ����� �����������;
		for (j = 1; j < N;  j++) X[j-1] = (X[j-1]+X[j])/2.;
		for (j = 1; j < N1; j++) Y[j-1] = (Y[j-1]+Y[j])/2.;
		for (j = 1; j < N2; j++) Z[j-1] = (Z[j-1]+Z[j])/2.;

/////////////////////////////
//...��������� ����� �������;
		add_new_lattice_X(par[0]); X[N-1]  = par[1];
		add_new_lattice_Y(par[2]); Y[N1-1] = par[3];
		add_new_lattice_Z(par[4]); Z[N2-1] = par[5];

////////////////////////////////////////////////////
//...���������� ������ ������, ���������� ���������;
		if (nd->hit)
		for (k = 0; k < nd->N;  k++) {
			mX = mY = mZ = 0;
			while (mX < N  && X[mX] < nd->X[k]) mX++;
			while (mY < N1 && Y[mY] < nd->Y[k]) mY++;
			while (mZ < N2 && Z[mZ] < nd->Z[k]) mZ++;
			nd->hit[k] = (mX-1)+(N-1)*((mY-1)+(N1-1)*(mZ-1));
		}
	}
}

////////////////////////////////////////////////
//...converts block structure from NASTRAN file;
int CGrid::converts_nas(char * id_NODES, unsigned long & count, unsigned long upper_limit, int id_long)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long count_beg = count, ppos_cur = count_beg;
		const int STR_SIZE = 200;
		int  m1, m2, m3, m4, m5, m6, m7, m8, j, m, l, l0, N_geom, part = id_long ? 16 : 8;
		char one_line[STR_SIZE+1], * pchar, temp = '\0';

////////////////////////////////////////
//...�������������� ��������� ���������;
		zero_grid(); N_geom = 0;
		while (ppos_cur < upper_limit) {
			if ((pchar = strstr(id_NODES+ppos_cur, "\xA")) != NULL)
			count = (unsigned long)(pchar-id_NODES)+1; else
			count = upper_limit;

			memcpy(one_line, id_NODES+ppos_cur, min((unsigned long)STR_SIZE, count-ppos_cur));
			one_line[min((unsigned long)STR_SIZE, count-ppos_cur)] = '\x0';
			ppos_cur = count;

///////////////////////////////////////
//...,������������ ��������� ���������;
			if (! ::strncmp(one_line, "CBEAM",  5) || 
				 ! ::strncmp(one_line, "CROD",   4)) l0 = 2; else
			if (! ::strncmp(one_line, "CTRIA3", 6)) l0 = 3;	else
			if (! ::strncmp(one_line, "CQUAD4", 6)) {
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); m3 = atoi(one_line+part*5); swap(temp, one_line[part*6]); 
				swap(temp, one_line[part*7]); m4 = atoi(one_line+part*6); swap(temp, one_line[part*7]); 
				if (m1 == m4) {                             m4 = -1;} else
				if (m2 == m1) { m1 = m2; m2 = m3; m3 = m4;  m4 = -1;} else
				if (m3 == m2) { swap(m1, m3); swap(m2, m4); m4 = -1;} else
				if (m4 == m3) { m3 = m2; m2 = m1; m1 = m4;  m4 = -1;}
				l0 = 2;
				if (m3 >= 0) l0++;
				if (m4 >= 0) l0++;
			}
			else
			if (! ::strncmp(one_line, "CTETRA",   6)) l0 = 4; else
			if (! ::strncmp(one_line, "CHEXA",    5)) l0 = 8; else
			if (! ::strncmp(one_line, "CPENTA",   6)) l0 = 6; else
			if (! ::strncmp(one_line, "CPYRAMID", 8)) l0 = 5; else l0 = -5;
			l0 += 5;

///////////////////////////////////////////////
//...��������� ��������� �� ��������� ��������;
			geom_ptr_add(l0-2, N_geom);
		}

/////////////////////////////////////////////////////////////////////
//...��������� ������ ��������� ���������, �������� �� ����� �������;
		ppos_cur = count_beg;
		if ((geom = (int *)new_struct((geom_ptr[N_geom]+1)*sizeof(int))) != NULL)
		while (ppos_cur < upper_limit) {

			if ((pchar = strstr(id_NODES+ppos_cur, "\xA")) != NULL)
			count = (unsigned long)(pchar-id_NODES)+1; else
			count = upper_limit;

			memcpy(one_line, id_NODES+ppos_cur, min((unsigned long)STR_SIZE, count-ppos_cur));
			one_line[min((unsigned long)STR_SIZE, count-ppos_cur)] = '\x0';
			ppos_cur = count;

/////////////////////////////////////
//...,����������� ��������� ��������;
			if (! ::strncmp(one_line, "CBEAM", 5)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = GL_LINE_STRIP;
				geom[l+1] = 5;
				geom[l+2] = -1;
				geom[l+3] = -m;
				geom[l+4] = -j;
				geom[l+5] = m1;
				geom[l+6] = m2;
			}
			else
			if (! ::strncmp(one_line, "CROD", 4)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = GL_LINE_STRIP;
				geom[l+1] = 5;
				geom[l+2] = -1;
				geom[l+3] = -m;
				geom[l+4] = -j;
				geom[l+5] = m1;
				geom[l+6] = m2;
			}
			else
			if (! ::strncmp(one_line, "CTRIA3", 6)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); m3 = atoi(one_line+part*5); swap(temp, one_line[part*6]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = GL_TRIANGLES;
				geom[l+1] = 6;
				geom[l+2] = -1;
				geom[l+3] = -m;
				geom[l+4] = -j;
				geom[l+5] = m1;
				geom[l+6] = m2;
				geom[l+7] = m3;
		  }
		  else
		  if (! ::strncmp(one_line, "CQUAD4", 6)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); m3 = atoi(one_line+part*5); swap(temp, one_line[part*6]); 
				swap(temp, one_line[part*7]); m4 = atoi(one_line+part*6); swap(temp, one_line[part*7]); 
				if (m1 == m4) {                             m4 = -1;} else
				if (m2 == m1) { m1 = m2; m2 = m3; m3 = m4;  m4 = -1;} else
				if (m3 == m2) { swap(m1, m3); swap(m2, m4); m4 = -1;} else
				if (m4 == m3) { m3 = m2; m2 = m1; m1 = m4;  m4 = -1;}
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = m4 >= 0 ? GL_QUADS : (m3 > 0 ? GL_TRIANGLES : GL_LINE_STRIP);
				geom[l+1] = m4 >= 0 ? 7 : (m3 > 0 ? 6 : 5);
				geom[l+2] = -1;
				geom[l+3] = -m;
				geom[l+4] = -j;
				geom[l+5] = m1;
				geom[l+6] = m2;
				l0 = 7;
				if (m3 >= 0) geom[l+l0++] = m3;
				if (m4 >= 0) geom[l+l0++] = m4;
		  }
		  else
		  if (! ::strncmp(one_line, "CTETRA", 6)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); m3 = atoi(one_line+part*5); swap(temp, one_line[part*6]); 
				swap(temp, one_line[part*7]); m4 = atoi(one_line+part*6); swap(temp, one_line[part*7]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = GL_TETRA;
				geom[l+1] = 7;
				geom[l+2] = -1;
				geom[l+3] = -m;
				geom[l+4] = -j;
				geom[l+5] = m1;
				geom[l+6] = m2;
				geom[l+7] = m3;
				geom[l+8] = m4;
			}
			else
			if (! ::strncmp(one_line, "CHEXA", 5)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); m3 = atoi(one_line+part*5); swap(temp, one_line[part*6]); 
				swap(temp, one_line[part*7]); m4 = atoi(one_line+part*6); swap(temp, one_line[part*7]); 
				swap(temp, one_line[part*8]); m5 = atoi(one_line+part*7); swap(temp, one_line[part*8]); 
				swap(temp, one_line[part*9]); m6 = atoi(one_line+part*8); swap(temp, one_line[part*9]); 
				if (ppos_cur < upper_limit) {//...���������� ��������� ������;
					if ((pchar = strstr(id_NODES+ppos_cur, "\xA")) != NULL)
					count = (unsigned long)(pchar-id_NODES)+1; else
					count = upper_limit;

					memcpy(one_line, id_NODES+ppos_cur, min((unsigned long)STR_SIZE, count-ppos_cur));
					one_line[min((unsigned long)STR_SIZE, count-ppos_cur)] = '\x0';
					ppos_cur = count;
						
					swap(temp, one_line[part*2]); m7 = atoi(one_line+part);	 swap(temp, one_line[part*2]); 
					swap(temp, one_line[part*3]); m8 = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				}
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = GL_BOXS;
				geom[l+1] = 11;
				geom[l+2]  = -1;
				geom[l+3]  = -m;
				geom[l+4]  = -j;
				geom[l+5]  = m1;
				geom[l+6]  = m2;
				geom[l+7]  = m3;
				geom[l+8]  = m4;
				geom[l+9]  = m5;
				geom[l+10] = m6;
				geom[l+11] = m7;
				geom[l+12] = m8;
			}
			else
			if (! ::strncmp(one_line, "CPENTA", 6)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); m3 = atoi(one_line+part*5); swap(temp, one_line[part*6]); 
				swap(temp, one_line[part*7]); m4 = atoi(one_line+part*6); swap(temp, one_line[part*7]); 
				swap(temp, one_line[part*8]); m5 = atoi(one_line+part*7); swap(temp, one_line[part*8]); 
				swap(temp, one_line[part*9]); m6 = atoi(one_line+part*8); swap(temp, one_line[part*9]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = GL_PENTA;
				geom[l+1] = 9;
				geom[l+2] = -1;
				geom[l+3] = -m;
				geom[l+4] = -j;
				geom[l+5] = m1;
				geom[l+6] = m2;
				geom[l+7] = m3;
				geom[l+8] = m4;
				geom[l+9] = m5;
				geom[l+10] = m6;
			}
			else
			if (! ::strncmp(one_line, "CPYRAMID", 8)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*4]); m1 = atoi(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); m2 = atoi(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); m3 = atoi(one_line+part*5); swap(temp, one_line[part*6]); 
				swap(temp, one_line[part*7]); m4 = atoi(one_line+part*6); swap(temp, one_line[part*7]); 
				swap(temp, one_line[part*8]); m5 = atoi(one_line+part*7); swap(temp, one_line[part*8]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = geom_ptr[geom[0]++];
				geom[l]   = GL_PYRAMID;
				geom[l+1] = 8;
				geom[l+2] = -1;
				geom[l+3] = -m;
				geom[l+4] = -j;
				geom[l+5] = m1;
				geom[l+6] = m2;
				geom[l+7] = m3;
				geom[l+8] = m4;
				geom[l+9] = m5;
			}
		}
	}
	return(1);
}

void CGrid::converts_nas(char * ch_NODES, int id_long)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_nas  (id_NODES, count, upper_limit, id_long);
  delete_struct(id_NODES);
}

void CGrid::converts_nas(const char * ch_NODES, int id_long)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_nas  (id_NODES, count, upper_limit, id_long);
  delete_struct(id_NODES);
}

///////////////////////////////////////////////////
//...reading boundary conditions from NASTRAN file;
int CGrid::condit_nas(char * id_NODES, unsigned long & count, unsigned long upper_limit, int id_long)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long count_beg = count, ppos_cur = count_beg;
		const int STR_SIZE = 200;
		int  m1, m2, j, m, l, l0, N_cond, part = id_long ? 16 : 8;
		char one_line[STR_SIZE+1], * pchar, temp = '\0';

////////////////////////////////////////
//...�������������� ��������� ���������;
		N_cond = 0;
		while (ppos_cur < upper_limit) {
			if ((pchar = strstr(id_NODES+ppos_cur, "\xA")) != NULL)
			count = (unsigned long)(pchar-id_NODES)+1; else
			count = upper_limit;

			memcpy(one_line, id_NODES+ppos_cur, min((unsigned long)STR_SIZE, count-ppos_cur));
			one_line[min((unsigned long)STR_SIZE, count-ppos_cur)] = '\x0';
			ppos_cur = count;

////////////////////////////////////
//...,������������ �������� �������;
			if (! ::strncmp(one_line, "SPCD",   4)) l0 = 4; else
			if (! ::strncmp(one_line, "PLOAD4", 6)) l0 = 6; else l0 = 0;

///////////////////////////////////////////////
//...��������� ��������� �� ��������� ��������;
			cond_ptr_add(l0-2, N_cond);
		}

///////////////////////////////////////////////////////////////////
//...��������� ������ ��������� �������, �������� �� ����� �������;
		ppos_cur = count_beg;
		if ((cond = (int *)new_struct((cond_ptr[N_cond]+1)*sizeof(int))) != NULL)
		while (ppos_cur < upper_limit) {

			if ((pchar = strstr(id_NODES+ppos_cur, "\xA")) != NULL)
			count = (unsigned long)(pchar-id_NODES)+1; else
			count = upper_limit;

			memcpy(one_line, id_NODES+ppos_cur, min((unsigned long)STR_SIZE, count-ppos_cur));
			one_line[min((unsigned long)STR_SIZE, count-ppos_cur)] = '\x0';
			ppos_cur = count;

////////////////////////////////////
//...,����������� ��������� �������;
		  if (! ::strncmp(one_line, "SPCD", 4)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); m1 = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = cond_ptr[cond[0]++];
				cond[l]   = GL_SPDC;
				cond[l+1] = 2;
				cond[l+2] = m;
				cond[l+3] = m1;
			}
			else
			if (! ::strncmp(one_line, "PLOAD4", 6)) {
				swap(temp, one_line[part*2]); m  = atoi(one_line+part);   swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*3]); j  = atoi(one_line+part*2); swap(temp, one_line[part*3]); 
				swap(temp, one_line[part*8]); m1 = atoi(one_line+part*7); swap(temp, one_line[part*8]); 
				swap(temp, one_line[part*9]); m2 = atoi(one_line+part*8); swap(temp, one_line[part*9]); 
/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
				l = cond_ptr[cond[0]++];
				cond[l]   = GL_PLOAD4;
				cond[l+1] = 4;
				cond[l+2]  = m;
				cond[l+3]  = j;
				cond[l+4]  = m1;
				cond[l+5]  = m2;
			}
		}
	}
	return(1);
}

void CGrid::condit_nas(char * ch_NODES, int id_long)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  condit_nas		(id_NODES, count, upper_limit, id_long);
  delete_struct(id_NODES);
}

void CGrid::condit_nas(const char * ch_NODES, int id_long)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  condit_nas		(id_NODES, count, upper_limit, id_long);
  delete_struct(id_NODES);
}

/////////////////////////////////////
//...reading nodes from NASTRAN file;
int CGrid::nodes_nas(char * id_NODES, unsigned long & count, unsigned long upper_limit, int id_long, int id_status)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long count_beg = count, ppos_cur = count_beg;
		const int STR_SIZE = 200;
		int  k, part = id_long ? 16 : 8;
		char one_line[STR_SIZE+1], * pchar, temp = '\0';
		double X, Y, Z;

///////////////////////////////////////////
//...���������� ���������� ����� ���������;
		while (ppos_cur < upper_limit) {

///////////////////////////
//...,�������� ���� ������;
			if ((pchar = strstr(id_NODES+ppos_cur, "\xA")) != NULL)
			count = (unsigned long)(pchar-id_NODES)+1; else
			count = upper_limit;

			memcpy(one_line, id_NODES+ppos_cur, min((unsigned long)STR_SIZE, count-ppos_cur));
			one_line[min((unsigned long)STR_SIZE, count-ppos_cur)] = '\x0';
			ppos_cur = count;

////////////////////////////////////////
//...,���������� ���������� ������ ����;
			if (! ::strncmp(one_line, "GRID", 4)) {
				swap(temp, one_line[part*2]); k = atoi(one_line+part);			 swap(temp, one_line[part*2]); 
				swap(temp, one_line[part*4]); X = user_strtod(one_line+part*3); swap(temp, one_line[part*4]); 
				swap(temp, one_line[part*5]); Y = user_strtod(one_line+part*4); swap(temp, one_line[part*5]); 
				swap(temp, one_line[part*6]); Z = user_strtod(one_line+part*5); swap(temp, one_line[part*6]); 

///////////////////////////////////////////////////////////////////////////
//...������� ���������� ���� � ��� ����� �� ���������� ������������� �����;
				if (add_new_point (X, Y, Z)) hit[N-1] = k;
				else return(0);
			}
		}

///////////////////////////////////////////////////////////////////////////
//...��������� ������� ����� (� ��������� ��������) �� ��������� ���������;
		if (id_status == OK_STATE && hit && N) {
			int * h_mask = NULL, N_h_mask, N_i_mask, l, j;

			for (N_h_mask = N_i_mask = hit[k = 0]; k < N; k++) 
			if  (N_h_mask < hit[k]) N_h_mask = hit[k]; else
			if  (N_i_mask > hit[k]) N_i_mask = hit[k];
			if ((h_mask = (int *)new_struct((N_h_mask-N_i_mask+1)*sizeof(int))) != NULL)
			for (k = 0; k < N; k++) h_mask[hit[k]-N_i_mask] = k;

			if (geom && geom_ptr && h_mask)
			for (k = 0; k < geom[0]; k++)
			for (l = geom_ptr[k], j = 5; j <= geom[l+1]+1; j++)
				geom[l+j] = h_mask[geom[l+j]-N_i_mask];

			if (cond && cond_ptr && h_mask)
			for (k = 0; k < cond[0]; k++)
				if (cond[l = cond_ptr[k]] == GL_PLOAD4) {
					cond[l+4] = h_mask[cond[l+4]-N_i_mask];
					cond[l+5] = h_mask[cond[l+5]-N_i_mask];
				}
				else
				if (cond[l] == GL_SPDC) cond[l+3] = h_mask[cond[l+3]-N_i_mask];

			delete_struct(h_mask);
		}

////////////////////////////////////////////////////////////////////////////////
//...��������� ������� ��������� (� ��������� ��������) �� ���������� ���������;
		if (id_status == OK_STATE && cond && cond_ptr) {
			int * b_mask = NULL, N_b_mask, N_l_mask, l, j;

			for (l = cond_ptr[k = 0];	k < cond[0]; l = cond_ptr[++k]) 
			if (cond[l] == GL_PLOAD4) {
				if (! k)			N_b_mask = N_l_mask = cond[l+3];	else
				if (N_b_mask < cond[l+3]) N_b_mask = cond[l+3]; else
				if (N_l_mask > cond[l+3]) N_l_mask = cond[l+3];
			}
			if (geom && geom_ptr && (b_mask = (int *)new_struct((N_b_mask-N_l_mask+1)*sizeof(int))) != NULL)
			for (k = 0; k < geom[0]; k++)  
				if ((j = -geom[geom_ptr[k]+2]) <= N_b_mask && N_l_mask <= j) b_mask[j-N_l_mask] = k;

			if (b_mask)
			for (k = 0; k < cond[0]; k++)
				if (cond[l = cond_ptr[k]] == GL_PLOAD4) cond[l+3] = b_mask[cond[l+3]-N_l_mask];

			delete_struct(b_mask);
		}
	}
	return(1);
}

void CGrid::nodes_nas(char * ch_NODES, int id_long)
{
  unsigned long count, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_nas (id_NODES, count = 0, upper_limit, id_long);
  condit_nas   (id_NODES, count = 0, upper_limit, id_long);
  nodes_nas    (id_NODES, count = 0, upper_limit, id_long);
  delete_struct(id_NODES);
}

void CGrid::nodes_nas(const char * ch_NODES, int id_long)
{
  unsigned long count, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_nas (id_NODES, count = 0, upper_limit, id_long);
  condit_nas   (id_NODES, count = 0, upper_limit, id_long);
  nodes_nas    (id_NODES, count = 0, upper_limit, id_long);
  delete_struct(id_NODES);
}

///////////////////////////////////////////////////
//...converts block structure from ABAQUS_INP file;
int CGrid::converts_inp(char * id_NODES, unsigned long & count, unsigned long upper_limit)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long count_beg, ppos_cur = count, upper, upper_element;
		const int STR_SIZE = 250;
		int  m[21], j, j_max = 19, l, l0, N_geom, elem, N_part, default_phase = -1;;
		char one_line[STR_SIZE+1], part_name[STR_SIZE+1], * pchar, * ppos;

////////////////////////////////////////
//...�������������� ��������� ���������;
		zero_grid(); N_geom = 0;

/////////////////////////////////////
//...������������� ��� ������� *Part;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name="); count_beg = ppos_cur;
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");

////////////////////////////////////////////////////////////////////
//...���������� ��� ������� *Part � ���������� ��� ������� ��������;
			ONE_LINE(id_NODES, pchar, ppos_cur, upper, part_name, STR_SIZE);
			if ((pchar = strstr(part_name, ",")) != NULL) pchar[0] = 0;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");

//////////////////////////////////////////////////////////
//...���������� ��������� ���� ��������� �� ������� *Part;
			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);

///////////////////////////////////////
//...,������������ ��������� ���������;
				l0 = element_type(elem, id_NODES+ppos_cur);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

/////////////////////////////////////////////////
//...������������� ��� �������� ������ ���������;
				if (l0 >= 3)
				while (ppos_cur < upper_element) {
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

///////////////////////////////////////////////
//...��������� ��������� �� ��������� ��������;
					geom_ptr_add(l0-2, N_geom);
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}

///////////////////////////////////////////////////////////////////////////////////////////////
//...��������� ������ ��������� ���������, �������� �� ����� �������, ������� ������ ���������;
		N_part = 0;
		if ((ppos_cur = count_beg) < upper_limit && geom_ptr && (geom = (int *)new_struct(geom_ptr[N_geom]*sizeof(int))) != NULL)
			while ((upper = ppos_cur)  < upper_limit) {                                   
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); count = upper; upper -= sizeof("*End Part");
			N_part++;

			ONE_LINE(id_NODES, pchar, ppos_cur, upper, part_name, STR_SIZE);
			if ((pchar = strstr(part_name, ",")) != NULL) pchar[0] = 0;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");

			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				l0 = element_type(elem, id_NODES+ppos_cur);				
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

/////////////////////////////////////////////////
//...������������� ��� �������� ������ ���������;
				if (l0 >= 3)
				while (ppos_cur < upper_element) {
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

////////////////////////////////////
//...����������� ��������� ��������;
					m[0] = atoi(one_line);
					for (ppos = strstr(one_line, ","), j = 5; ppos && j < l0; ppos = j == j_max ? one_line : strstr(ppos+1, ","),  j++) {
						m[j-4] = atoi(ppos+1);
						if (j == j_max) ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
					}

/////////////////////////////////////////////////////
//...������� ������� � ��������� (������� ���������);
					l = geom_ptr[geom[0]++];
					geom[l]   = elem;
					geom[l+1] = l0-2;
					geom[l+2] = -N_part;
					geom[l+3] = -m[0];
					geom[l+4] = default_phase;
					for (j = 5; j < l0; j++)
					geom[l+j] = m[j-4];
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}
	}
	return(1);
}

void CGrid::converts_inp(char * ch_NODES)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_inp (id_NODES, count, upper_limit);
  delete_struct(id_NODES);
}

void CGrid::converts_inp(const char * ch_NODES)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_inp (id_NODES, count, upper_limit);
  delete_struct(id_NODES);
}

//////////////////////////////////////////////////////
//...reading boundary conditions from ABAQUS_INP file;
int CGrid::condit_inp(char * id_NODES, unsigned long & count, unsigned long upper_limit, int max_phase, int id_status)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long count_beg, count_part, count_length, ppos_cur = count, upper, upper_element;
		const int STR_SIZE = 250;
		int  elem, j_node = -1, j_elem = -1, j, k, l, l0, N_cond, k0, N_part;
		char one_line[STR_SIZE+1], part_name[STR_SIZE+1], * pchar, * ppos, temp = '\0';

////////////////////////////////////////
//...�������������� ��������� ���������;
		unsigned long N_buf = 50, buf_incr = 20, current_list, 
						* elements_list = (unsigned long *)new_struct((N_buf*2+1)*sizeof(unsigned long));
		N_cond = 0;

/////////////////////////////////////
//...������������� ��� ������� *Part;
		ppos_cur = count;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name="); count_beg = ppos_cur;
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");

//////////////////////////////////////////////////////////////////////
//...���������� ��� ������� *Part � ���������� ������ ��������� �����;
			ONE_LINE(id_NODES, pchar, ppos_cur, upper, part_name, STR_SIZE);
			if ((pchar = strstr(part_name, ",")) != NULL) pchar[0] = 0; count_part = ppos_cur;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");

/////////////////////////////////////
//...������������ ������������ �����;
			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////
//...,���������� ��� ����;
				if (strstr(one_line, "internal")) elem = ERR_STATE; 
				if (strstr(one_line, "generate")) elem = GL_NODES_GENERATE; else elem = GL_NODES;
				l0 = 0;

////////////////////////////////////////////////////
//...������������� ��� ��������� ����� ������� ����;
				if (elem == GL_NODES || elem == GL_NODES_GENERATE) {
					while (ppos_cur < upper_element) {
						ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
////////////////////////
//...,������������ ����;
						ppos = one_line; l0++; 
						while ((ppos = strstr(ppos+1, ",")) != NULL && (k0 = atoi(ppos+1)) > 0) l0++;
					}
///////////////////////////////////////////////
//...��������� ��������� �� ��������� ��������;
					cond_ptr_add(l0+1, N_cond);
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
			}

/////////////////////////////////////////
//...������������ ������������ ���������;
			ppos_cur = count_part;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");

			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////////
//...,���������� ��� ��������;
				if (strstr(one_line, "internal")) elem = ERR_STATE; 
				if (strstr(one_line, "generate")) elem = GL_ELEMENTS_GENERATE; else elem = GL_ELEMENTS; 
				l0 = 0;

///////////////////////////////////////////////////////
//...������������� ��� ��������� ��������� ������� ����;
				if (elem == GL_ELEMENTS || elem == GL_ELEMENTS_GENERATE) {
					while (ppos_cur < upper_element) {
						ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
////////////////////////////////////
//...������������ ��������� �������;
						ppos = one_line; l0++;
						while ((ppos = strstr(ppos+1, ",")) != NULL && (k0 = atoi(ppos+1)) > 0) l0++;
					}
///////////////////////////////////////////////
//...��������� ��������� �� ��������� ��������;
					cond_ptr_add(l0+1, N_cond);
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}

/////////////////////////////////////////////////////
//...��������� ������� �� �������� ����� � ���������;
		ppos_cur = count;
		if (cond_ptr && (cond = (int *)new_struct(cond_ptr[N_cond]*sizeof(int))) != NULL && geom && geom_ptr && hit) {
                                                                
///////////////////////////////////////////////
//...������� ������ ���������� � �������������;
			FILE * TST = NULL;
			if (id_status == OK_STATE) {
				TST = fopen("LayersCondit.sta", "w");
				fprintf(TST, "Layers condition, N_loading_sets = %d: \n", N_cond);
			}

/////////////////////////////////////
//...������������� ��� ������� *Part;
			ppos_cur = count_beg; N_part = 0;
			while ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
				N_part++;
////////////////////////////////////
//...���� ������ ������������ �����;
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, part_name, STR_SIZE);
				if ((pchar = strstr(part_name, ",")) != NULL) pchar[0] = 0; count_part = ppos_cur;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");

/////////////////////////////////////////
//...������� ������� ������������� �����;
				int * node_pernum = NULL, node_max = 0, node_min = 0, k_beg, k_end;
				for (k_end = k_beg = 0; k_beg < N; k_beg++) {
					if (! k_beg || hit[k_beg] < hit[k_beg-1]) k_end++; 
					if (! k_beg || hit[k_beg] < hit[k_beg-1] && N_part == k_end) break;
				}
				if (k_beg < N) {
					node_max = hit[k_beg]; node_min = node_max; 
					for (k_end = k_beg+1; k_end < N; k_end++) 
						if (hit[k_end] < hit[k_end-1]) break;

					for (k = k_beg+1; k < k_end; k++) {
						if (hit[k] < node_min) node_min = hit[k];
						if (hit[k] > node_max) node_max = hit[k];
					}
					node_pernum = (int *)new_struct((node_max-node_min+1)*sizeof(int)); 
					for (k = k_beg; k < k_end; k++) node_pernum[hit[k]-node_min] = k;
				}

//////////////////////////////////
//...��������� ������������ �����;
				if (node_pernum)
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////
//...,���������� ��� ����;
					if (strstr(one_line, "internal")) elem = ERR_STATE; 
					if (strstr(one_line, "generate")) elem = GL_NODES_GENERATE; else elem = GL_NODES;

/////////////////////////////////////////////////////////////////////
//...������������� ��� ���� ������� ���� � ��������� ���������� �����;
					if (elem == GL_NODES || elem == GL_NODES_GENERATE) {
						if (TST) {
							if ((pchar = strstr(one_line, ",")) != NULL) pchar[0] = '\x0';
							else one_line[strlen(one_line)-1] = '\x0';
							fprintf(TST, "%s, j_node = %i\n", one_line, j_node);
						}
						l = cond_ptr[cond[0]++];
						cond[l++] = elem; l0 = l;
						cond[l++] = 1;
						cond[l++] = j_node--;
						while (ppos_cur < upper_element) {
							ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

////////////////////////////////////////////////////////////////////////////////////////////
//...����������� ������� ����� �� ������ (���������� ��������� ����� ������� �������������);
							cond[l++] = node_pernum[atoi(ppos = one_line)-node_min];
							cond[l0]++;

							while ((ppos = strstr(ppos+1, ",")) != NULL && (k0 = atoi(ppos+1)) > 0) {
								cond[l++] = node_pernum[k0-node_min];
								cond[l0]++;
							}
							if (elem == GL_NODES_GENERATE) cond[l-1] = k0;
						}
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
				}
				delete_struct(node_pernum);

/////////////////////////////////////////////
//...������� ������� ������������� ���������;
				int * elem_pernum = NULL, elem_max = 0, elem_min = 0;
				for (k_beg = 0; k_beg < geom[0]; k_beg++) 
					if (N_part+geom[geom_ptr[k_beg]+2] == 0) break;
				if (k_beg < geom[0]) {
					elem_max = -geom[geom_ptr[k_beg]+3]; elem_min = elem_max; 
					for (k = k_beg; k < geom[0]; k++) if (N_part+geom[geom_ptr[k]+2] == 0) {
						if (-geom[geom_ptr[k]+3] < elem_min) elem_min = -geom[geom_ptr[k]+3];
						if (-geom[geom_ptr[k]+3] > elem_max) elem_max = -geom[geom_ptr[k]+3];
					}
					elem_pernum = (int *)new_struct((elem_max-elem_min+1)*sizeof(int));
					for (k = k_beg; k < geom[0]; k++) 
						if (N_part+geom[geom_ptr[k]+2] == 0) elem_pernum[-geom[geom_ptr[k]+3]-elem_min] = k;
				}

//////////////////////////////////////
//...��������� ������������ ���������;
				ppos_cur = count_part; current_list = elements_list[0]; 
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");

				if (elem_pernum)
				while ((upper_element = ppos_cur) < upper) {
					if (elements_list[0] == N_buf) {
						unsigned long * new_elements_list = elements_list; elements_list = (unsigned long *)new_struct(((N_buf += buf_incr)*2+1)*sizeof(unsigned long));
						memcpy(elements_list, new_elements_list, (new_elements_list[0]*2+1)*sizeof(unsigned long)); delete_struct(new_elements_list);
					}
					elements_list[elements_list[0]*2+1] = ppos_cur;

					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

/////////////////////////////
//...���������� ��� ��������;
					if (strstr(one_line, "internal")) elem = ERR_STATE; 
					if (strstr(one_line, "generate")) elem = GL_ELEMENTS_GENERATE; else elem = GL_ELEMENTS; 

/////////////////////////////////////////////////////////////////////////
//...������������� ��� �������� ������� ���� � ��������� ��������� �����;
					if (elem == GL_ELEMENTS || elem == GL_ELEMENTS_GENERATE) {
						if ((pchar = strstr(one_line, ", generate")) != NULL || (pchar = strstr(one_line, "\xD")) != NULL || (pchar = strstr(one_line, "\xA")) != NULL) {
							swap(pchar[0], temp);
							elements_list[elements_list[0]*2+2] = strlen(one_line);
							elements_list[0]++;

							if (TST) fprintf(TST, "%s, j_elem = %i\n", one_line, j_elem);
							swap(pchar[0], temp);
						}
						l = cond_ptr[cond[0]++];
						cond[l++] = elem; l0 = l;
						cond[l++] = 1;
						cond[l++] = j_elem--;
						while (ppos_cur < upper_element) {
							ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

////////////////////////////////////////////////////////////////////////////////////////////////
//...����������� ������� ��������� �� ������ (���������� ��������� ����� ������� �������������);
							cond[l++] = elem_pernum[atoi(ppos = one_line)-elem_min];
							cond[l0]++;

							while ((ppos = strstr(ppos+1, ",")) != NULL && (k0 = atoi(ppos+1)) > 0) {
								cond[l++] = elem_pernum[k0-elem_min];
								cond[l0]++;
							}
							if (elem == GL_ELEMENTS_GENERATE) cond[l-1] = k0;
						}
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
				}
				delete_struct(elem_pernum);

//////////////////////////////////////////////////////////
//...������� ������ �� �������� *Solid Section � ��������;
				ppos_cur = count_part;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");

				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////////////////////////////////////////
//...���� ��� ������� � ������� ������ �� �������� � ��������;
					if ( (pchar = strstr(one_line, ", material")) != NULL) {
						pchar[0] = temp; count_length = strlen(one_line);
						
						for (count = current_list; count < elements_list[0]; count++)  
							if (elements_list[count*2+2] == count_length && ! strncmp(one_line, id_NODES+elements_list[count*2+1], count_length)) break;

						if (geom && geom_ptr && cond && cond_ptr) {
							for (l = cond_ptr[k = 0];	k < cond[0]; l = cond_ptr[++k]) 
							if (cond[l] == (int)GL_ELEMENTS) {
								if (cond[l+2]+(int)count+1 == 0)
									for (j = cond[l+1]; j > 1; j--) if (cond[l+1+j] < N)
										geom[geom_ptr[cond[l+1+j]]+4] = cond[l+2];
							}
							else
							if (cond[l] == (int)GL_ELEMENTS_GENERATE) {
								if (cond[l+2]+(int)count+1 == 0)
									for (j = cond[l+3]; j <= cond[l+4]; j += cond[l+5]) if (j < N)
										geom[geom_ptr[j]+4] = cond[l+2];
							}
						}
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");
				}

//////////////////////////////////////////////////////////
//...������� ������ �� �������� *Shell Section � ��������;
				ppos_cur = count_part;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Shell Section, elset=");

				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
					
//////////////////////////////////////////////////////////////
//...���� ��� ������� � ������� ������ �� �������� � ��������;
					if ( (pchar = strstr(one_line, ", material")) != NULL) {
						pchar[0] = temp; count_length = strlen(one_line);
						
						for (count = current_list; count < elements_list[0]; count++)  
							if (elements_list[count*2+2] == count_length && ! strncmp(one_line, id_NODES+elements_list[count*2+1], count_length)) break;

						if (geom && geom_ptr && cond && cond_ptr) {
							for (l = cond_ptr[k = 0];	k < cond[0]; l = cond_ptr[++k]) 
							if (cond[l] == (int)GL_ELEMENTS) {
								if (cond[l+2]+(int)count+1 == 0)
									for (j = cond[l+1]; j > 1; j--) if (cond[l+1+j] < N)
										geom[geom_ptr[cond[l+1+j]]+4] = cond[l+2];
							}
							else
							if (cond[l] == (int)GL_ELEMENTS_GENERATE) {
								if (cond[l+2]+(int)count+1 == 0)
									for (j = cond[l+3]; j <= cond[l+4]; j += cond[l+5]) if (j < N)
										geom[geom_ptr[j]+4] = cond[l+2];
							}
						}
						swap(pchar[0], temp);
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Shell Section, elset=");
				}

/////////////////////////////////////////////////////////
//...������� ������ �� �������� *Beam Section � ��������;
				ppos_cur = count_part;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Beam Section, elset=");

				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
					
//////////////////////////////////////////////////////////////
//...���� ��� ������� � ������� ������ �� �������� � ��������;
					if ( (pchar = strstr(one_line, ", material")) != NULL) {
						pchar[0] = temp; count_length = strlen(one_line);
						
						for (count = current_list; count < elements_list[0]; count++)  
							if (elements_list[count*2+2] == count_length && ! strncmp(one_line, id_NODES+elements_list[count*2+1], count_length)) break;

						if (geom && geom_ptr && cond && cond_ptr) {
							for (l = cond_ptr[k = 0];	k < cond[0]; l = cond_ptr[++k]) 
							if (cond[l] == (int)GL_ELEMENTS) {
								if (cond[l+2]+(int)count+1 == 0)
									for (j = cond[l+1]; j > 1; j--) if (cond[l+1+j] < N)
										geom[geom_ptr[cond[l+1+j]]+4] =  cond[l+2];
							}
							else
							if (cond[l] == (int)GL_ELEMENTS_GENERATE) {
								if (cond[l+2]+(int)count+1 == 0)
									for (j = cond[l+3]; j <= cond[l+4]; j += cond[l+5]) if (j < N)
										geom[geom_ptr[j]+4] = cond[l+2];
							}
						}
						swap(pchar[0], temp);
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Beam Section, elset=");
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
			}
			if (TST) fclose(TST);
		}
		delete_struct(elements_list);
	}
	return(1);
}

void CGrid::condit_inp(char * ch_NODES)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  condit_inp	(id_NODES, count, upper_limit);
  delete_struct(id_NODES);
}

void CGrid::condit_inp(const char * ch_NODES)
{
  unsigned long count    = 0, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  condit_inp	(id_NODES, count, upper_limit);
  delete_struct(id_NODES);
}

////////////////////////////////////////
//...reading nodes from ABAQUS_INP file;
int CGrid::nodes_inp(char * id_NODES, unsigned long & count, unsigned long upper_limit)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long ppos_cur = count, upper, upper_element;
		const int STR_SIZE = 250;
		char one_line[STR_SIZE+1], * pchar, * ppos;
		double X, Y, Z;
		int  k, l, j, l0, N_part = 0, nodes_part;
		
/////////////////////////////////////
//...������������� ��� ������� *Part;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Node");
			N_part++; nodes_part = N; 

///////////////////////////////////////////
//...���������� ���������� ����� ���������;
			if ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

				while (ppos_cur < upper_element) {
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
/////////////////////////////
//...,���������� ������ ����;
					k = atoi(one_line); X = Y = Z = 0.;
					if ((ppos = strstr(one_line, ",")) != NULL) {
						X = user_strtod(ppos+1);
						if ((ppos = strstr(ppos+1, ",")) != NULL) {
							Y = user_strtod(ppos+1);
							if ((ppos = strstr(ppos+1, ",")) != NULL) {
								Z = user_strtod(ppos+1);
							}
						}
					}
///////////////////////////////////////////////////
//...��������� ���������� ���� � ��� ������� �����;
					if (add_new_point(X, Y, Z)) hit[N-1] = k;
					else return(0);
				}
			}

/////////////////////////////////////////
//...������� ������� ������������� �����;
			int * node_pernum = NULL, node_max = 0, node_min = 0;
			node_max = hit[nodes_part]; node_min = node_max; 
			for (k = nodes_part; k < N; k++) {
				if (hit[k] < node_min) node_min = hit[k];
				if (hit[k] > node_max) node_max = hit[k];
			}
			node_pernum = (int *)new_struct((node_max-node_min+1)*sizeof(int)); 
			for (k = nodes_part; k < N; k++) node_pernum[hit[k]-node_min] = k;

////////////////////////////////////////////////////////
//...��������� ���� � ��������� �� ���������� ���������;
			if (geom && geom_ptr) {
				for (k = 0; k < geom[0]; k++) 
					if (geom[(l = geom_ptr[k])+2]+N_part == 0) {
						l0 = geom[l+1]+2;
						for (j = 5; j < l0; j++) if (geom[l+j] <= node_max) 
							geom[l+j] = node_pernum[geom[l+j]-node_min]; 
					}
			}
			delete_struct(node_pernum);

			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}
	}
	return(1);
}

void CGrid::nodes_inp(char * ch_NODES, int max_phase)
{
  unsigned long count, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_inp (id_NODES, count = 0, upper_limit);
  nodes_inp    (id_NODES, count = 0, upper_limit);
  condit_inp   (id_NODES, count = 0, upper_limit, max_phase);
  delete_struct(id_NODES);
}

void CGrid::nodes_inp(const char * ch_NODES, int max_phase)
{
  unsigned long count, upper_limit;
  char        * id_NODES = read_struct_ascii(ch_NODES);
  if         (! id_NODES) return;
  user_Count   (id_NODES, 0, upper_limit, '\x0');
  converts_inp (id_NODES, count = 0, upper_limit);
  nodes_inp    (id_NODES, count = 0, upper_limit);
  condit_inp   (id_NODES, count = 0, upper_limit, max_phase);
  delete_struct(id_NODES);
}

/////////////////////////////////////////////////////////
//...reading structure from text file (X, Y, Z, rad, ll);
int CGrid::stru_in(char * id_STRU, unsigned long & count, unsigned long upper_limit, double * par)
{
	if (id_STRU && par && count < upper_limit) {
		unsigned long ppos_cur = count, upper;
		const int STR_SIZE = 250;
		char one_line[STR_SIZE+1], * pchar, * ppos;
		double X, Y, Z, pp[2];
		int  k = 0, M;

////////////////////////////////
//...���������� ����� ���������;
		ONE_LINE(id_STRU, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
		M = atoi(one_line);

////////////////////////////////
//...���������� ������� �������;
		ONE_LINE(id_STRU, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
		par[0] = MAX_HIT; par[1] = MIN_HIT;
		par[0] = user_strtod(one_line);
		if ((ppos = strstr(one_line, ",")) != NULL) {
			par[1] = user_strtod(ppos+1);
		}
		ONE_LINE(id_STRU, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
		par[2] = MAX_HIT; par[3] = MIN_HIT;
		par[2] = user_strtod(one_line);
		if ((ppos = strstr(one_line, ",")) != NULL) {
			par[3] = user_strtod(ppos+1);
		}
		ONE_LINE(id_STRU, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
		par[4] = MAX_HIT; par[5] = MIN_HIT;
		par[4] = user_strtod(one_line);
		if ((ppos = strstr(one_line, ",")) != NULL) {
			par[5] = user_strtod(ppos+1);
		}

//////////////////////////////////////////////////////////
//...���������� ���������� ������� ��������� � �� �������;
		zero_grid(); add_params(2);
		while ((upper = ppos_cur) < upper_limit && k++ < M) {
				ONE_LINE(id_STRU, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
				X = Y = Z = pp[0] = pp[1] = 0.;
				X = user_strtod(one_line);
				if ((ppos = strstr(one_line, ",")) != NULL) {
					Y = user_strtod(ppos+1);
					if ((ppos = strstr(ppos+1, ",")) != NULL) {
						Z = user_strtod(ppos+1);
						if ((ppos = strstr(ppos+1, ",")) != NULL) {
							pp[0] = user_strtod(ppos+1);
							if ((ppos = strstr(ppos+1, ",")) != NULL) {
								pp[1] = user_strtod(ppos+1);
							}
						}
					}
				}
//////////////////////////////////////////////////////////////
//...��������� ���������� ������� � ������� � ��������� �����;
				if (add_new_point(X, Y, Z, 0., 0., 1., pp)) hit[N-1] = k;
				else return(0);
		}
	}
	return(1);
}

void CGrid::stru_in(char * ch_STRU, double * par)
{
  unsigned long count, upper_limit;
  char        * id_STRU = read_struct_ascii(ch_STRU);
  if         (! id_STRU) return;
  user_Count   (id_STRU, 0, upper_limit, '\x0');
  stru_in      (id_STRU, count = 0, upper_limit, par);
  delete_struct(id_STRU);
}

void CGrid::stru_in(const char * ch_STRU, double * par)
{
  unsigned long count, upper_limit;
  char        * id_STRU = read_struct_ascii(ch_STRU);
  if         (! id_STRU) return;
  user_Count   (id_STRU, 0, upper_limit, '\x0');
  stru_in      (id_STRU, count = 0, upper_limit, par);
  delete_struct(id_STRU);
}

///////////////////////////////////////////////////////////////////////
//...reproducing grid in Surfer format with isometrical transformation;
void CGrid::grid_out(FILE * GR, double fi, double theta, double fX, int m_axis, int shift, int m_contour, int element)
{
	if (! GR) return;
	double P[3], X, Y = M_PI/180., CZ = cos(fi    *= Y), SZ = sin(fi),
											 CY = cos(theta *= Y), SY = sin(theta),
											 CX = cos(fX    *= Y), SX = sin(fX);
	int    k, l, j, k1, k2;
	switch (m_axis) {
         case AXIS_X: k1 = 1; k2 = 2; break;
         case AXIS_Y: k1 = 2; k2 = 0; break;
         default    : k1 = 0; k2 = 1; break;
	}
	for (l = 0, k = min(geom[0], element); k > 0; k--) l += geom[l += (m_contour ? 1 : 2)];
	for (k = element; k < geom[0]; k++) if (! m_contour) {
      switch (geom[(++l)++]) {
         case GL_LINE_STRIP: {
              fprintf(GR, "  %d  0\n", geom[l]);
              for (j = shift; j < geom[l]; j++) {
                  set_point(P, geom[l+j+1], CZ, SZ, CY, SY, CX, SX);
                  fprintf (GR, " % e  % e\n", P[k1], P[k2]);
              }
              l += geom[l];
         }    break;
         case GL_TRIANGLE_STRIP: {
              fprintf(GR, "  %i  0\n", j = geom[l]);
              for (j = shift; j < geom[l]; j++) {
                  set_point(P, geom[l+j+1], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }
              fprintf(GR, "  %i  0\n", j = (geom[l]+1-shift)/2);
              for (j = shift; j < geom[l]; j += 2) {
                  set_point(P, geom[l+j+1], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }
              if ((j = (geom[l]-shift)/2) > 1) {
                  fprintf(GR, "  %i  0\n", j);
                  for (j = 1+shift; j < geom[l]; j += 2) {
                      set_point(P, geom[l+j+1], CZ, SZ, CY, SY, CX, SX);
                      fprintf(GR, " % e  % e\n", P[k1], P[k2]);
                  }
              }
              l += geom[l];
         }    break;
         case GL_QUAD_STRIP: {
              for (j = shift; j < geom[l]-3; j += 2) {
                  set_point(P, geom[l+j+1], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, "  5  0\n % e  % e\n", X = P[k1], Y = P[k2]);

                  set_point(P, geom[l+j+2], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);

                  set_point(P, geom[l+j+4], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);

                  set_point(P, geom[l+j+3], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n % e  % e\n", P[k1], P[k2], X, Y);
              }
              l += geom[l];
         }    break;
         case GL_TRIANGLE_FAN: {
              for (j = shift; j < geom[l]-2; j++) {
                  set_point(P, geom[l+1+shift], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, "  4  0\n % e  % e\n", X = P[k1], Y = P[k2]);

                  set_point(P, geom[l+j+2], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);

                  set_point(P, geom[l+j+3], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n % e  % e\n", P[k1], P[k2], X, Y);
              }
              l += geom[l];
         }    break;
         case GL_QUADS:
         case GL_TRIANGLES: {
              fprintf(GR, "  %i  0\n", geom[l]-shift+1);

              set_point(P, geom[l+1+shift], CZ, SZ, CY, SY, CX, SX);
              fprintf(GR, " % e  % e\n", X = P[k1], Y = P[k2]);

              for (j = 2+shift; j <= geom[l]; j++) {
                  set_point(P, geom[l+j], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }
              fprintf(GR, " % e  % e\n", X, Y);
              l += geom[l];
         }    break;
         case GL_BOXS: {
              fprintf(GR, "  5  0\n");

              set_point(P, geom[l+1+shift], CZ, SZ, CY, SY, CX, SX);
              fprintf(GR, " % e  % e\n", X = P[k1], Y = P[k2]);

              for (j = 2+shift; j <= 4+shift; j++) {
                  set_point(P, geom[l+j], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }
              fprintf(GR, " % e  % e\n", X, Y);

              fprintf(GR, "  4  0\n");
              for (j = 1+shift; j <= 2+shift; j++) {
                  set_point(P, geom[l+4*j], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              set_point(P, geom[l+5+shift], CZ, SZ, CY, SY, CX, SX);
              fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              fprintf(GR, " % e  % e\n", X, Y);

              fprintf(GR, "  4  0\n");
              for (j = 1+shift; j <= 2+shift; j++) {
                  set_point(P, geom[l+4*j-1], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }
              for (j = 2+shift; j >= 1+shift; j--) {
                  set_point(P, geom[l+4*j-2], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              fprintf(GR, "  2  0\n");
              for (j = 5+shift; j <= 6+shift; j++) {
                  set_point(P, geom[l+j], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              fprintf(GR, "  2  0\n");
              for (j = 7+shift; j <= 8+shift; j++) {
                  set_point(P, geom[l+j], CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              l += geom[l];
         }    break;
         case FV_BOXS: {
              fprintf(GR, "  5  0\n");

              set_point(P, geom[l+1+shift], 1+shift, CZ, SZ, CY, SY, CX, SX);
              fprintf(GR, " % e  % e\n", X = P[k1], Y = P[k2]);

              for (j = 2+shift; j <= 4+shift; j++) {
                  set_point(P, geom[l+1], j, CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }
              fprintf(GR, " % e  % e\n", X, Y);

              fprintf(GR, "  4  0\n");
              for (j = 1+shift; j <= 2+shift; j++) {
                  set_point(P, geom[l+1], 4*j, CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              set_point(P, geom[l+1], 5+shift, CZ, SZ, CY, SY, CX, SX);
              fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              fprintf(GR, " % e  % e\n", X, Y);

              fprintf(GR, "  4  0\n");
              for (j = 1+shift; j <= 2+shift; j++) {
                  set_point(P, geom[l+1], 4*j-1, CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }
              for (j = 2+shift; j >= 1+shift; j--) {
                  set_point(P, geom[l+1], 4*j-2, CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              fprintf(GR, "  2  0\n");
              for (j = 5+shift; j <= 6+shift; j++) {
                  set_point(P, geom[l+1], j, CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              fprintf(GR, "  2  0\n");
              for (j = 7+shift; j <= 8+shift; j++) {
                  set_point(P, geom[l+1], j, CZ, SZ, CY, SY, CX, SX);
                  fprintf(GR, " % e  % e\n", P[k1], P[k2]);
              }

              l += geom[l];
         }    break;
      }
	} 
	else {
      fprintf(GR, "%3i  0\n", geom[(++l)]);
      for (j = 0; j < geom[l]; j++) {
           set_point(P, geom[l+j+1], CZ, SZ, CY, SY, CX, SX);
           fprintf (GR, " % e  % e\n", P[k1], P[k2]);
      }
      l += geom[l];
	}
}

////////////////////////////////
//...��������� �������� �������;
void CGrid::grid_out(char * GR_FILE, double fi, double theta, double fX, int m_axis, int shift, int m_contour, int element)
{
	if (! GR_FILE) return;
	FILE * GR = fopen(GR_FILE, "w");
	grid_out(GR, fi, theta, fX, m_axis, shift, m_contour, element);
	fclose(GR);
}
void CGrid::grid_out(const char * GR_FILE, double fi, double theta, double fX, int m_axis, int shift, int m_contour, int element)
{
	if (! GR_FILE) return;
	FILE * GR = fopen(GR_FILE, "w");
	grid_out(GR, fi, theta, fX, m_axis, shift, m_contour, element);
	fclose(GR);
}

///////////////////////////////////////////////////////////////////////////
//...function for reproducing centers of multilevel grids in Surfer format;
void CGrid::TestGrid(FILE * GR, double mtb, double fi, double theta, double fX, int m_axis, int id_normal)
{
	if (! GR) return;
	double P[3], PN[3], Y = M_PI/180., CZ = cos(fi    *= Y), SZ = sin(fi),
												 CY = cos(theta *= Y), SY = sin(theta),
												 CX = cos(fX    *= Y), SX = sin(fX);
	int    k1, k2, i, j, l, m, n;
	switch (m_axis) {
         case AXIS_X: k1 = 1; k2 = 2; break;
         case AXIS_Y: k1 = 2; k2 = 0; break;
         default    : k1 = 0; k2 = 1; break;
	}
	for (i = 0; i < N;  i++)
	for (j = 0; j < N1 || j == 0 && N1 == 0; j++)
	for (l = 0; l < N2 || l == 0 && N2 == 0; l++) if (! hit || hit[i+(j+l*N1)*N]) {
       set_point (P, i, m = N1 == 0 ? i : j, n = N2 == 0 ? i : l, CZ, SZ, CY, SY, CX, SX);
       fprintf(GR, "  5  0\n % e  % e\n % e  % e\n % e  % e\n % e  % e\n % e  % e\n",
               P[k1]-mtb, P[k2]-mtb,
               P[k1]+mtb, P[k2]-mtb,
               P[k1]+mtb, P[k2]+mtb,
               P[k1]-mtb, P[k2]+mtb,
               P[k1]-mtb, P[k2]-mtb);
       if (id_normal && nX && nY && nZ) {
           PN[0] = nX[i];
           PN[1] = nY[m];
           PN[2] = nZ[n];
           point_iso<double>(PN,   NULL, CZ, SZ, CY, SY, CX, SX);
           fprintf(GR, "  2  0\n % e  % e\n % e  % e\n",
                   P[k1],                             P[k2],
                   P[k1]+12.*mtb*PN[k1],              P[k2]+12.*mtb*PN[k2]);
       }
	}
}

////////////////////////////////
//...��������� �������� �������;
void CGrid::TestGrid(char * GR_FILE, double mtb, double fi, double theta, double fX, int m_axis, int id_normal)
{
	if (! GR_FILE) return;
	FILE * GR = fopen(GR_FILE, "w");
	TestGrid(GR, mtb, fi, theta, fX, m_axis, id_normal);
	fclose(GR);
}
void CGrid::TestGrid(const char * GR_FILE, double mtb, double fi, double theta, double fX, int m_axis, int id_normal)
{
	if (! GR_FILE) return;
	FILE * GR = fopen(GR_FILE, "w");
	TestGrid(GR, mtb, fi, theta, fX, m_axis, id_normal);
	fclose(GR);
}
