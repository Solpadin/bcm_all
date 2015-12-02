/*============================================*/
/*                  PRO_WIRE                  */
/*============================================*/
#ifndef ___PRO_WIRE___
#define ___PRO_WIRE___
#define ___PRO_WIRE2D_cpp___
#ifndef ___ABRIDGE_PROFILE_MODE___
#include "cgrid.h"
#else
#include "utils.h"
#endif

/////////////////////////////////////////////////////////
//...����� 2D ����������� �������� ��� Bars25 � �� �����;
  int    GetWIRESampleCount(void);
  char * GetWIRESampleName (int N_sm);
  char * GetWIREGOSTName   (int N_sm);

///////////////////////////////////////////////////
//...����� ����� ��������� ������������ � �� �����;
int    GetREDUCEDWIRESampleCount(void);
char * GetREDUCEDWIRESampleName (int N_sm);
char * GetREDUCEDWIREGOSTName   (int N_sm);

///////////////////////////////////////////////////////////////////////
//...������� ����������� ��������� ��� ����������� �������� ����������;
  void * CreateWIREContext(int N_sm);
  void * CreateREDUCEDWIREContext(int N_sm);

////////////////////////////////////////////////////////////
//...���������� ������ ��������������� ������� �� ���������;
int GetREDUCEDWIRENumber(void * context);

/////////////////////////////////////////////////////////////////////
//...�p����p�������� ������������� ������������ ��p���� ��� ��������;
#ifndef ___ABRIDGE_PROFILE_MODE___
int wire2D_init (void * context, CGrid * block_nd = NULL);
#endif

///////////////////////////////////////////
//...p������� ������� ��p�� ��� ����-�����;
int wire2D_solv(void * context, int & sec, int & hund, int id_solv);

/////////////////////////////////////////////////////
//...���������� ��� 2D ������� (������ ��� ��������);
char * GetFuncWIREName(int num);
enum Func_Wire {
    _FWIRE1 = 0,  
    _FWIRE2,  
    _FWIRE3,  
    _FWIRE4,  
    _FWIRE5,  
    _FWIRE6,
    _NUM_WIRE_FUNCTIONS
};
#define NUM_WIRE_FUNCTIONS _NUM_WIRE_FUNCTIONS
inline int GetFuncWIRECount  (void) { return(NUM_WIRE_FUNCTIONS);}
inline int GetFuncWIREDefault(void) { return(_FWIRE1);}

//////////////////////////////////////
//...������� ����������� ��� ��������;
#ifndef ___ABRIDGE_PROFILE_MODE___
void solver_wire2D(void * context, CGrid * nd, double * F, double * par = NULL, int id_F = GetFuncWIREDefault());
#endif

/////////////////////////////////////////////////////////////////
//...������������ ������� ��� �������������� � �������� ��������;
double get_full_square(void * context);
double get_square(void * context);
double get_inertia_Y(void * context);
double get_inertia_Z(void * context);
double get_inertia_YZ(void * context);
double get_static_Y(void * context);
double get_static_Z(void * context);
double get_torsion_J(void * context);
double get_effect_max_E(void * context);
double get_effect_min_E(void * context);
double get_effect_max_At(void * context);
double get_E_core(void * context);
double get_E_wire(void * context);
double get_nju_core(void * context);
double get_nju_wire(void * context);
double get_diameter(void * context);
double get_last_diameter(void * context);
double get_diameter(void * context, int layer);
int    get_last_num_layer(void * context);
int    get_num_layer(void * context, int layer);
double get_last_twist_angle(void * context);
double get_twist_angle(void * context, int layer);
double get_effect_GJ(void * context);
double get_effect_EI(void * context);

#endif
