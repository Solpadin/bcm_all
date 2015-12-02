/*============================================*/
/*                  PRO_SKIN                  */
/*============================================*/
#ifndef ___PRO_SKIN___
#define ___PRO_SKIN___
#define ___PRO_SKIN2D_cpp___
#define ___nPRO_SKIN3D_cpp___
#include "cgrid.h"

/////////////////////////////////////////////////////////////
//...����� 2D ����������� �������� ��� ����-����� � �� �����;
  int    GetSK2DSampleCount(void);
  char * GetSK2DSampleName (int N_sm);

///////////////////////////////////////////////////////////////////////
//...������� ����������� ��������� ��� ����������� �������� ����������;
  void * CreateSK2DContext(int N_sm);
////////////////////////////////////////////////////////////////////////////
//...��������� ���������� ���������� ����-����� ��� ����� ������������ ����;
  int SetFasaHmg(void * context, double nju, double G1, double G2, double C0);

/////////////////////////////////////////////////////////////////
//...��������� ���������� ����-����� ��� ����� ������������ ����;
  int SetSkin(void * context, double skin);

///////////////////////////////////////////////////////////////////////
//...�p����p�������� ������������� ������������ ��p���� ��� ����-�����;
  int skin3D_init(void * context, CGrid * block_nd = NULL);
  int skin2D_init(void * context, CGrid * block_nd = NULL);

///////////////////////////////////////////
//...p������� ������� ��p�� ��� ����-�����;
int skin2D_solv(void * context, int & sec, int & hund, int id_solv);

//////////////////////////////////////////////////////////
//...���������� ��� 2D ������� (������ ������������ ����);
char * GetFuncSKName(int num);
enum Func_SK {
    _FK1 = 0,  
    _FK2,  
    _FK3,  
    _FK4,  
    _FK5,  
    _FK6,
    _FK7,  
    _FK8,
    _FK9,  
    _FK10,
    _NUM_SK_FUNCTIONS
};
#define NUM_SK_FUNCTIONS _NUM_SK_FUNCTIONS
inline int GetFuncSKCount  (void) { return(NUM_SK_FUNCTIONS);}
inline int GetFuncSKDefault(void) { return(_FK1);}

////////////////////////////////////////
//...������� ����������� ��� ����-�����;
void solver_SK(void * context, CGrid * nd, double * F, double * par = NULL, int id_F = GetFuncSKDefault());

#endif
