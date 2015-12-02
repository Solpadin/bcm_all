/*============================================*/
/*                  PRO_FITTING               */
/*============================================*/
#ifndef ___PRO_FITTING___
#define ___PRO_FITTING___
#define ___PRO_FITTING_cpp___
#ifndef ___ABRIDGE_PROFILE_MODE___
#include "cgrid.h"
#else
#include "utils.h"
#endif

///////////////////////////////////////////////////
//...����� ����� ��������� ������������ � �� �����;
int    GetFITTINGSampleCount(void);
char * GetFITTINGSampleName (int N_sm);
char * GetFITTINGGOSTName   (int N_sm);

///////////////////////////////////////////////////////////////////////////////////
//...������� ����������� ��������� ��� ���� ������ ��������� ������������ ��������;
  void * CreateFITTINGContext(int N_sm);

////////////////////////////////////////////////////////////////////////
//...�p����p�������� ������������� ������������ ��p���� ��� ������������;
#ifndef ___ABRIDGE_PROFILE_MODE___
  int fitting_init(void * context, CGrid * block_nd = NULL);
#endif

/////////////////////////////////////////////
//...p������� ������� ��p�� ��� ������������;
int fitting_solv(void * context, int & sec, int & hund, int id_solv);

/////////////////////////////////////////////////////////////////
//...������������ ������� ��� �������������� � �������� ��������;

#endif
