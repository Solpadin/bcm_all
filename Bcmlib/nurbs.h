/*============================================*/
/*                   NURBS                    */
/*============================================*/
#ifndef ___NURBS___
#define ___NURBS___

#include "ccells.h"

///////////////////////////////////////////////////////////////
//...�������, ����������� ���������� ������ � ����������� � 2D;
void nurbs2D(double t, CMap * mp, double * P, int id_map = 0);
void nurbs2D(double u, double v, CMap * mp, double * P, int id_map = 0);

///////////////////////////////////////////////////////////////
//...�������, ����������� ���������� ������ � ����������� � 3D;
void nurbs3D(double t, CMap * mp, double * P, int id_map = 0);
void nurbs3D(double u, double v, CMap * mp, double * P, int id_map = 0);

//////////////////////////////////////////////////////////////////////////////
//...�������, ����������� ���������� ������ � ����������� � 3D, � ����� �����;
void nurbs3D(double t, CMap * mp, double * P, double * ort, int id_map = 0);
void nurbs3D(double u, double v, CMap * mp, double * P, double * ort, int id_map = 0);

////////////////////////////////////////////////////////////////////////////////////////
//...�������, ����������� ���������� ������ � ����������� � 3D, �����, � ����� ��������;
void nurbs3D(double t, CMap * mp, double * P, double * ort, double * R_t, int id_map = 0);
void nurbs3D(double u, double v, CMap * mp, double * P, double * ort, double * R_t, int id_map = 0);

///////////////////////////////////////////////////////
//...�������, ����������� �������� ����������� �������;
void Nurbs3D(double u, double v, CMap * mp, double * P, int id_map = 0);
void Nurbs3D(double u, double v, CMap * mp, double * P, double * ort, int id_map = 0);

////////////////////////////////////////////////////
//...���������� ���������������� �������� �� ������;
CMap * get_pm_2D(Cells * ce, int Max_N);

/////////////////////////////////////////////////////////////////////////
//...�������������� ���������������� �������� � ������� ���������� �����;
void pm_sph_convert(CMap * pm, CMap * mp, CMap * sph);

///////////////////////////////////////////////////////////////////////////////////
//...������� �������� �� ������ � �� ������� ������ �� ���������� ����������� � 3D;
void Trassa3D(double & u, double & v, CMap * mp, double * P, double * ort, int MAX_NEWTON = 0);
void Trassa3D_skel(double & u, CMap * mp, double * P, double * ort, int MAX_NEWTON = 0);

///////////////////////
//...���������� ������;
double nurbs3D_u_par(double v, CMap * mp, double * P, double eps = 1e-1, int M_iter = 40, int N_ini = 10, double fnorm = 2.);
double nurbs3D_v_par(double u, CMap * mp, double * P, double eps = 1e-1, int M_iter = 40, int N_ini = 10, double fnorm = 2.);

/////////////////////////////////////////////////////
//...����������� ����������� ������� �� ������� IGES;
CMap * Nurbs3D_IGES	(char * ch_IGES, int IGES_cod, unsigned long id_num = 1);
void Trassa_write		(char * ch_TRASSA, CMap ** trassa, int max_trassa);
int  Trassa_read		(char * ch_TRASSA, CMap ** trassa, int MAX_TRASSA, int * max_trassa = NULL);
int  Trassa_convert	(char * ch_TRASSA, char * ch_TRASSA_OUT, CMap ** trassa, int MAX_TRASSA, int cod = 128);
void Trassa_link		(char * ch_TRASSA, char * ch_TRASSA_OUT, CMap ** trassa, int MAX_TRASSA, double f = 1.);
void Trassa_link_skel(char * ch_TRASSA, char * ch_TRASSA_OUT, CMap ** trassa, int MAX_TRASSA, double f = 1.);
CCells* Trassa_read	(char * ch_TRASSA);
CCells* Trassa_read_skel(char * ch_TRASSA);
void track_igs			(char * ch_TRACK,  CGrid * nd_track);
void marker_igs		(char * ch_MARKER, CGrid * nd_marker);
#endif
