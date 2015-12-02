#include "stdafx.h"

#include "shapes.h"
#include "cacou3d.h"

#define Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CAcou3D::NUM_KAPPA = 8;
int CAcou3D::NUM_SHEAR = 13;
int CAcou3D::NUM_HESS  = 18;

//////////////////////////////////
//...initialization of the blocks;
int CAcou3D::block_shape_init(Block<complex> & B, Num_State id_free)
{
	int	k, m, set_cmpl = 1; //...set_cmpl == 2 -- real case for extern problem;
   if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
	if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<complex>;
/*		if ((B.type & ERR_CODE) == ZOOM_BLOCK) {
			B.shape->add_shape(CreateShape<double>(AU3D_ZOOM_SHAPE));
			if (B.mp[7] < 0.) B.shape->set_cmpl(set_cmpl);
		}
		else
		if ((B.type & ERR_CODE) == ELLI_BLOCK)	B.shape->add_shape(CreateShape<double>(AU3D_ELLI_SHAPE));
		else*/												B.shape->add_shape(CreateShape<complex>(AU3D_POLY_SHAPE));

////////////////////////
//...setting parameters;
		/*if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(param[1]*fabs(B.mp[7]), 0., 0., 0.);
		else*/												B.shape->set_shape(param[1]*fabs(B.mp[7]));

//////////////////////////////////
//...multipoles number adaptation;
		if (B.link[NUM_PHASE] == -2) //...another phase of media!!!
			B.shape->init1(UnPackInts(get_param(0), 1), solver.id_norm*2, 3); else
		/*if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[7] < 0.) B.shape->init1(UnPackInts(get_param(0), 1), solver.id_norm*2, 1); else
		if ((B.type & ERR_MASK) == SUB_UGOLOK_BLOCK)			    B.shape->init1(UnPackInts(get_param(0), 1), solver.id_norm*2, 1);
		else*/																	 B.shape->init1(UnPackInts(get_param(0)),		solver.id_norm*2, 1); 

//////////////////////////////////////////////////////////////////////////////
//...setting acselerator, local system of coordinate and init parametrization;
		B.shape->set_local(B.mp+1);
		B.shape->release  ();
   }

//////////////////////////////////////
//...setting frequency and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
         B.shape->set_local(B.mp+1);
			if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(param[1]*fabs(B.mp[7]), 0., 0., 0.);
			else												B.shape->set_shape(param[1]*fabs(B.mp[7]));
		}
		else
		if (id_free == OK_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->set_potential(solver.hh[k = (int)(&B-B.B)][0][m], m*2);
		else
		if (id_free == NO_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->get_potential(solver.hh[k = (int)(&B-B.B)][0][solver.id_norm+m], m*2);
		else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
	}
	return(B.shape != NULL);
}

//////////////////////////////////////////////
//...realization of complex acoustic pressure;
void  CAcou3D::jump1(double * P, int i, int m)
{
	m += solver.id_norm;
	B[i].shape->parametrization(P, 1);
	B[i].shape->cpy(0, B[i].shape->FULL(solver.hh[i][0][m]));
}

////////////////////////////////////////////////////////////////////
//...realization of normal derivatives of complex acoustic pressure;
void  CAcou3D::jump2(double * P, int i, int m)
{
	m += solver.id_norm;
	B[i].shape->parametrization_grad(P);
	B[i].shape->deriv_N();
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
}

////////////////////////////////////////
//...realization of acoustic attmidance;
void CAcou3D::jump3(double * P, int i, int m)
{
	m += solver.id_norm;
	B[i].shape->parametrization_grad(P);
	B[i].shape->deriv_N();
	B[i].shape->admittance(0, comp(P[6], P[7]));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
}

/////////////////////////////
//...вспомогательные функции;
void CAcou3D::jump_cmpl_cpy(double * p_cpy, double * p_cmp, complex * hh, int NN, int id_memset)
{
	if (id_memset) memset (hh, 0, sizeof(complex)*NN);
	int j;
	for (j = 0; p_cpy && j < NN; j++) hh[j] += comp(p_cpy[j]);
	for (j = 0; p_cmp && j < NN; j++) hh[j] += comp(0., p_cmp[j]);
}

void CAcou3D::jump_conj_cpy(int i, int m, int k, complex adm)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m += solver.id_norm;
	k += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) solver.hh[i][0][m][j] = conj(adm*solver.hh[i][0][k][j]);
}

//////////////////////////////////////////
//...realization vibro displacements (Ux);
void CAcou3D::jump1_common_x(double * P, int i, int m)
{
	//double G0 = 1./get_param(NUM_SHEAR),
	//		 C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), (complex *)NULL, 0., 0.);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), (complex *)NULL, 0., 0.);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), (complex *)NULL, 0., 0.);

	//B[i].shape->adm_xx(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), C0);
	//B[i].shape->adm_xy(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), C0);
	//B[i].shape->adm_xz(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), C0);

	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->p_cpy, 1., G0);
}

////////////////////////////////////////////////////
//...additional vibro compression displacement (Ux);
void CAcou3D::jump1_compress_x(double * P, int i, int m)
{
	//double C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//B[i].shape->adm_xx(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), -C0);
	//B[i].shape->adm_xy(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), -C0);
	//B[i].shape->adm_xz(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), -C0);
}

//////////////////////////////////////////
//...realization vibro displacements (Uy);
void CAcou3D::jump1_common_y(double * P, int i, int m)
{
	//double G0 = 1./get_param(NUM_SHEAR),
	//		 C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), (complex *)NULL, 0., 0.);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), (complex *)NULL, 0., 0.);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), (complex *)NULL, 0., 0.);

	//B[i].shape->adm_xy(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), C0);
	//B[i].shape->adm_yy(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), C0);
	//B[i].shape->adm_yz(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), C0);

	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->p_cpy, 1., G0);
}

////////////////////////////////////////////////////
//...additional vibro compression displacement (Uy);
void CAcou3D::jump1_compress_y(double * P, int i, int m)
{
	//double C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//B[i].shape->adm_xy(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), -C0);
	//B[i].shape->adm_yy(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), -C0);
	//B[i].shape->adm_yz(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), -C0);
}

//////////////////////////////////////////
//...realization vibro displacements (Uz);
void CAcou3D::jump1_common_z(double * P, int i, int m)
{
	//double G0 = 1./get_param(NUM_SHEAR),
	//		 C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), (complex *)NULL, 0., 0.);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), (complex *)NULL, 0., 0.);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), (complex *)NULL, 0., 0.);

	//B[i].shape->adm_xz(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), C0);
	//B[i].shape->adm_yz(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), C0);
	//B[i].shape->adm_zz(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), C0);

	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->p_cpy, 1., G0);
}

////////////////////////////////////////////////////
//...additional vibro compression displacement (Uz);
void CAcou3D::jump1_compress_z(double * P, int i, int m)
{
	//double C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//B[i].shape->adm_xz(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), C0);
	//B[i].shape->adm_yz(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), C0);
	//B[i].shape->adm_zz(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), C0);
}

///////////////////////////////////////////////////////////////
//...realization normal derivative of vibro displacements (Ux);
void CAcou3D::jump2_common_x(double * P, int i, int m)
{
	//double G0 = 1./get_param(NUM_SHEAR),
	//		 C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 0), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[3]*G0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[4]*G0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[5]*G0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 0., C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 0., C0);
}

//////////////////////////////////////////////////////////////////////////
//...additional normal derivative of vibro compression displacements (Ux);
void CAcou3D::jump2_compress_x(double * P, int i, int m)
{
	//double C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 0), 1., -C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 1., -C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 1., -C0);
}

///////////////////////////////////////////////////////////////
//...realization normal derivative of vibro displacements (Uy);
void CAcou3D::jump2_common_y(double * P, int i, int m)
{
	//double G0 = 1./get_param(NUM_SHEAR),
	//		 C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 0., C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[3]*G0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[4]*G0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[5]*G0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 0., C0);
}

//////////////////////////////////////////////////////////////////////////
//...additional normal derivative of vibro compression displacements (Uy);
void CAcou3D::jump2_compress_y(double * P, int i, int m)
{
	//double C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 1., -C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2), 1., -C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 1., -C0);
}

///////////////////////////////////////////////////////////////
//...realization normal derivative of vibro displacements (Uz);
void CAcou3D::jump2_common_z(double * P, int i, int m)
{
	//double G0 = 1./get_param(NUM_SHEAR),
	//		 C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 0., C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 0., C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[3]*G0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[4]*G0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[5]*G0);
}

//////////////////////////////////////////////////////////////////////////
//...additional normal derivative of vibro compression displacements (Uz);
void CAcou3D::jump2_compress_z(double * P, int i, int m)
{
	//double C0 = 1./get_param(NUM_SHEAR-1); m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 1., -C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 1., -C0);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2), 1., -C0);
}

////////////////////////////////////////////////////////////////
//...realization of surface forces for vibro displacements (Px);
void CAcou3D::jump4_common_x(double * P, int i, int m)
{
	//double C0 = get_param(NUM_SHEAR)/get_param(NUM_SHEAR-1)*2.; m += solver.id_norm;
	//int	 num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 0), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[3]*2.);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[4]);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[5]);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[4]);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[5]);
}

//////////////////////////////////////////////////////
//...additional vibro compression surface forces (Px);
void CAcou3D::jump4_compress_x(double * P, int i, int m)
{
	//double C0 = -get_param(NUM_SHEAR  )/get_param(NUM_SHEAR-1)*2., 
	//	 alpha = -get_param(NUM_SHEAR+1)/get_param(NUM_SHEAR+2)*4.; m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 0), 1., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[3]*alpha);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 1., C0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[3]*alpha);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 1., C0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[3]*alpha);
}

////////////////////////////////////////////////////////////////
//...realization of surface forces for vibro displacements (Pyx);
void CAcou3D::jump4_common_y(double * P, int i, int m)
{
	//double C0 = get_param(NUM_SHEAR)/get_param(NUM_SHEAR-1)*2.; m += solver.id_norm;
	//int	 num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 0., C0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[3]);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[3]);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[4]*2.);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[5]);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 0., C0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[5]);
}

//////////////////////////////////////////////////////
//...additional vibro compression surface forces (Py);
void CAcou3D::jump4_compress_y(double * P, int i, int m)
{
	//double C0 = -get_param(NUM_SHEAR  )/get_param(NUM_SHEAR-1)*2., 
	//	 alpha = -get_param(NUM_SHEAR+1)/get_param(NUM_SHEAR+2)*4.; m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1), 1., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[4]*alpha);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2), 1., C0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[4]*alpha);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 1., C0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[4]*alpha);
}

////////////////////////////////////////////////////////////////
//...realization of surface forces for vibro displacements (Pz);
void CAcou3D::jump4_common_z(double * P, int i, int m)
{
	//double C0 = get_param(NUM_SHEAR)/get_param(NUM_SHEAR-1)*2.; m += solver.id_norm;
	//int	 num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 0., C0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[3]);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 0., C0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[4]);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2), 0., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[3]);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[4]);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[5]*2.);
}

//////////////////////////////////////////////////////
//...additional vibro compression surface forces (Pz);
void CAcou3D::jump4_compress_z(double * P, int i, int m)
{
	//double C0 = -get_param(NUM_SHEAR  )/get_param(NUM_SHEAR-1)*2., 
	//	 alpha = -get_param(NUM_SHEAR+1)/get_param(NUM_SHEAR+2)*4.; m += solver.id_norm;
	//int num_hess = NUM_HESS+solver.id_norm;
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 0), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0), 1., C0);
	//B[i].shape->adm_x     (B[i].shape->FULL(solver.hh[i][0][m], 0, 0), P[5]*alpha);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1), 1., C0);
	//B[i].shape->adm_y     (B[i].shape->FULL(solver.hh[i][0][m], 0, 1), P[5]*alpha);
	//B[i].shape->admittance(B[i].shape->FULL(solver.hh[i][0][m], 0, 2), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2), 1., C0);
	//B[i].shape->adm_z     (B[i].shape->FULL(solver.hh[i][0][m], 0, 2), P[5]*alpha);
}

///////////////////////////////////////////////
//...realization normal derivatives of hessian;
void CAcou3D::hessian_deriv_N(double * P, int i)
{
	int num_hess = NUM_HESS+solver.id_norm;

	B[i].shape->cpy_xx();
	B[i].shape->deriv_N();
	B[i].shape->cpy_xx();
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 0));

	B[i].shape->cpy_xy();
	B[i].shape->deriv_N();
	B[i].shape->cpy_xy();
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 1));

	B[i].shape->cpy_yy();
	B[i].shape->deriv_N();
	B[i].shape->cpy_yy();
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], 0, 2));

	B[i].shape->cpy_xz();
	B[i].shape->deriv_N();
	B[i].shape->cpy_xz();
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 0));

	B[i].shape->cpy_yz();
	B[i].shape->deriv_N();
	B[i].shape->cpy_yz();
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 1));

	B[i].shape->cpy_zz();
	B[i].shape->deriv_N();
	B[i].shape->cpy_zz();
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], 0, 2));
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CAcou3D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	//for (int j = 0; j < solver.dim[i]; j++) {
	//	complex P[] = {solver.hh[i][0][m  ][j], 
	//						solver.hh[i][0][m+1][j], 
	//						solver.hh[i][0][m+2][j] };
	//	B[i].shape->norm_local(P);
	//	solver.hh[i][0][m  ][j] = P[0];
	//	solver.hh[i][0][m+1][j] = P[1];
	//	solver.hh[i][0][m+2][j] = P[2];
	//}
}

void CAcou3D::jump_make_common(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	//for (int j = 0; j < solver.dim[i]; j++) {
	//	complex P[] = { solver.hh[i][0][m  ][j], 
	//						 solver.hh[i][0][m+1][j], 
	//						 solver.hh[i][0][m+2][j] };
	//	B[i].shape->norm_common(P);
	//	solver.hh[i][0][m  ][j] = P[0];
	//	solver.hh[i][0][m+1][j] = P[1];
	//	solver.hh[i][0][m+2][j] = P[2];
	//}
}

///////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для акустической фазы пространства;
Num_State CAcou3D::gram1_phase1(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		complex h, hh, attm;
		double  f, P[8];
		int m = solver.id_norm, id_zero;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);
	
////////////////////////////////////////////////////////////////////////////////////////////
//...realization of Zommerfield or common boundary condition for external or internal block;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			id_zero =   nd->get_param(1, l) == MIN_HIT;
			attm = comp(nd->get_param(0, l),  nd->get_param(1, l))*comp(0., get_param(NUM_KAPPA));
			hh   = comp(nd->get_param(2, l),  nd->get_param(3, l));
			f = (id_zero ? 1.:1./(1.+norm(attm)))*nd->get_param(4, l);

			if (id_zero) attm = comp(0., MIN_HIT);
			if (hh.y == MAX_HIT) hh.y = 0.; //...for SONEX absorber identification;
			if (hh.y == MIN_HIT) {          //...for VN corrector identification;
				 hh.y = hh.x; 
				 hh.x = nd->get_param(1, l); 
				 attm = comp(0.);
			}
			if (! id_zero) hh *= comp(0., -get_param(NUM_KAPPA+1)); //..frequency dependence condition!!!
			P[6] = real(attm); //...передаем аттмиданс поглощения через точку;
			P[7] = imag(attm); 

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));

///////////////////////////////////////////////////////////
//...вычисляем необходимые моменты коллокационного вектора;
			B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
			if (id_zero) {
				P[3] = P[4] = 0.; P[5] = 1.;
				jump1(P, i, 0);
				jump_conj_cpy(i, 2, 0);

				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m], f);
				if (abs(h = hh) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+2], B[i].shape->get_NN(), h*f);
			}
			else {
				jump3(P, i, 1);
				jump_conj_cpy(i, 3, 1);

				solver.to_equationDD(i, solver.hh[i][0][m+3], solver.hh[i][0][m+1], f);
				if (abs(h = hh) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+3], B[i].shape->get_NN(), h*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для механической фазы пространства;
Num_State CAcou3D::gram1_phase2(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), nu1 = get_param(NUM_SHEAR+1), hx, hy, hz, p4, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			hy = nd->get_param(1, l); 
			hz = nd->get_param(2, l);
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == MIN_HIT || p4 == 2.) {
			  hy = nd->nY[l]*hx;
			  hz = nd->nZ[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));

/////////////////////////////////////////////
//...jump of all neaded displacement moments;
			B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
			B[i].shape->parametrization_hess(P, 1);
			if (p4 == MIN_HIT || p4 == MAX_HIT) {
				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 
				jump1_common_z(P, i, 2); 
			}
			else
			if (p4 == (double)(SKEWS_BND-SPECIAL_BND) || p4 == (double)(UPTAKE_BND-SPECIAL_BND)) {
				hessian_deriv_N(P, i);
				jump4_common_x(P, i, 0); 
				jump4_common_y(P, i, 1); 
				jump4_common_z(P, i, 2); 

				jump1_common_x(P, i, 6); 
				jump1_common_y(P, i, 7); 
				jump1_common_z(P, i, 8); 
			}
			else
			if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				hessian_deriv_N(P, i);
				jump4_common_x(P, i, 0); 
				jump4_common_y(P, i, 1); 
				jump4_common_z(P, i, 2); 
			}
			B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
			B[i].shape->parametrization_hess(P, 1);
			if (p4 == MIN_HIT || p4 == MAX_HIT) {
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);
				jump1_compress_z(P, i, 2); solver.admittance(i, 2, G1);
			}
			else
			if (p4 == (double)(SKEWS_BND-SPECIAL_BND) || p4 == (double)(UPTAKE_BND-SPECIAL_BND)) {
				hessian_deriv_N(P, i);
				jump4_compress_x(P, i, 0); solver.admittance(i, 5, 0., 0, P[3]);
				jump4_compress_y(P, i, 1); solver.admittance(i, 5, 1., 1, P[4]); 
				jump4_compress_z(P, i, 2); solver.admittance(i, 5, 1., 2, P[5]);

				solver.admittance(i, 0, 1., 5, -P[3]); //...вычисляем касательную силу и нормальное перемещение;
				solver.admittance(i, 1, 1., 5, -P[4]);
				solver.admittance(i, 2, 1., 5, -P[5]);

				jump1_compress_x(P, i, 6); solver.admittance(i, 4, 0., 6, P[3]*G1);
				jump1_compress_y(P, i, 7); solver.admittance(i, 4, 1., 7, P[4]*G1);
				jump1_compress_z(P, i, 8); solver.admittance(i, 4, 1., 8, P[5]*G1);

				if (p4 == (double)(UPTAKE_BND-SPECIAL_BND)) { //...вычисляем поглощение;
					solver.admittance(i, 5, 1., 4, comp(0., get_param(NUM_SHEAR+3)*(1.-nu1)/(.5-nu1)));
					jump_conj_cpy	  (i, 6, 5);
				}
			}
			else
				if (0.<= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				hessian_deriv_N(P, i);
				jump4_compress_x(P, i, 0);
				jump4_compress_y(P, i, 1); 
				jump4_compress_z(P, i, 2);
			}
			jump_make_common(i, 0);

////////////////////////////
//...composition functional;
			if (p4 == (double)(SKEWS_BND-SPECIAL_BND))
			solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
			if (p4 == (double)(UPTAKE_BND-SPECIAL_BND))
			solver.to_equationDD(i, solver.hh[i][0][m+6], solver.hh[i][0][m+5], f);
			solver.to_equationDD(i, solver.hh[i][0][m],   solver.hh[i][0][m],   f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
			solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);

			if (p4 == MIN_HIT || p4 == (double)(SKEWS_BND-SPECIAL_BND) || p4 == MAX_HIT) f *= G1;
			if (fabs(hx) > EE && p4 == (double)(SKEWS_BND-SPECIAL_BND))
				solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*f);
			if (fabs(hx) > EE && p4 != (double)(SKEWS_BND-SPECIAL_BND))
				solver.to_equationHH(i, 0, solver.hh[i][0][m], hx*f);
			if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+1], hy*f);
			if (fabs(hz) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+2], hz*f);

		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе акустическая - акустическая фаза;
Num_State CAcou3D::transfer1_phase1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      complex f0 = comp(1.), g1 = comp(.5), g2 = comp(-.5/(get_param(NUM_KAPPA))), g0 = comp(-1.);
      double  f, P[6];
      int m = solver.id_norm;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(complex));
				}

///////////////////////////////////////////////////////
//...jump of acoustic pressure and velocity (P*g1+V*g2);
				B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
				jump1(P, i, 0); 
				jump2(P, i, 1); 
				solver.admittance(i, 1, g2, 0, g1); jump_conj_cpy(i, 3, 1);
				solver.admittance(i, 0, f0, 1, g0); jump_conj_cpy(i, 2, 0);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_KAPPA)); 
				jump1(P, k, 0); 
				jump2(P, k, 1); 
				solver.admittance(k, 1, g2, 0, g1); jump_conj_cpy(k, 2, 0);
				solver.admittance(k, 0, f0, 1, g0); jump_conj_cpy(k, 3, 1);

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[i][0][m], f);
				solver.to_transferTD(i, j, solver.hh[k][0][m+3], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+1], f);
			}
			break;
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе механическая - механическая фаза;
Num_State CAcou3D::transfer1_phase2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], f0 = 1., g1 = G1*.5, g2 = G1*.5, g0 = -G1;
      int m = solver.id_norm;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(complex));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
				B[i].shape->parametrization_hess(P, 1);
				hessian_deriv_N (P, i);

				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 
				jump1_common_z(P, i, 2); 

				jump2_common_x(P, i, 3); 
				jump2_common_y(P, i, 4); 
				jump2_common_z(P, i, 5); 

				B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
				B[i].shape->parametrization_hess(P, 1);
				hessian_deriv_N (P, i);

				jump1_compress_x(P, i, 0);
				jump1_compress_y(P, i, 1); 
				jump1_compress_z(P, i, 2); 
				
				jump2_compress_x(P, i, 3); solver.admittance(i, 0, g1, 3, g2); solver.admittance(i, 3, g0, 0, f0); 
				jump2_compress_y(P, i, 4); solver.admittance(i, 1, g1, 4, g2); solver.admittance(i, 4, g0, 1, f0); 
				jump2_compress_z(P, i, 5); solver.admittance(i, 2, g1, 5, g2); solver.admittance(i, 5, g0, 2, f0); 
				jump_make_common(i, 0);
				jump_make_common(i, 3);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+4)); 
				B[k].shape->parametrization_hess(P, 1);
				hessian_deriv_N (P, k);

				jump1_common_x(P, k, 0); 
				jump1_common_y(P, k, 1); 
				jump1_common_z(P, k, 2); 

				jump2_common_x(P, k, 3); 
				jump2_common_y(P, k, 4); 
				jump2_common_z(P, k, 5); 

				B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+3)); 
				B[k].shape->parametrization_hess(P, 1);
				hessian_deriv_N (P, k);

				jump1_compress_x(P, k, 0); 
				jump1_compress_y(P, k, 1); 
				jump1_compress_z(P, k, 2); 

				jump2_compress_x(P, k, 3); solver.admittance(k, 0, g1, 3, g2); solver.admittance(k, 3, g0, 0, f0); 
				jump2_compress_y(P, k, 4); solver.admittance(k, 1, g1, 4, g2); solver.admittance(k, 4, g0, 1, f0); 
				jump2_compress_z(P, k, 5); solver.admittance(k, 2, g1, 5, g2); solver.admittance(k, 5, g0, 2, f0); 
				jump_make_common(k, 0);
				jump_make_common(k, 3);

/////////////////////////////////////////////////////////////////////////////////////////////
//...сшивка функций и нормальных производных методом наименьших квадратов (в шкале давлений);
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferTT(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
			}
			break;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе акустическая - механическая фаза;
Num_State CAcou3D::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double G1 = get_param(NUM_KAPPA+2)/get_param(NUM_KAPPA), f, P[6], 
			    g2 = 1./get_param(NUM_KAPPA);
      int m = solver.id_norm;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

				f = nd->get_param(0, l);

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(complex));
				}

///////////////////////////////////////////////
//...jump of all neaded moments for this block;
				if (B[i].link[NUM_PHASE] == -1) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
					jump1(P, i, 0); 
					jump2(P, i, 1); solver.admittance(i, 1, g2); 
					jump_conj_cpy(i, 2, 0);
					jump_conj_cpy(i, 3, 1); 
				}
				if (B[i].link[NUM_PHASE] == -2) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, i);
					
					jump4_common_x(P, i, 6); 
					jump4_common_y(P, i, 7); 
					jump4_common_z(P, i, 8); 

					jump1_common_x(P, i, 9); 
					jump1_common_y(P, i, 10); 
					jump1_common_z(P, i, 11); 

					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, i);

					jump4_compress_x(P, i, 6); solver.admittance(i, 5, 0., 6, P[3]); 
					jump4_compress_y(P, i, 7); solver.admittance(i, 5, 1., 7, P[4]); 
					jump4_compress_z(P, i, 8); solver.admittance(i, 5, 1., 8, P[5]); 

					jump1_compress_x(P, i,  9); solver.admittance(i, 4, 0.,  9, P[3]*G1);
					jump1_compress_y(P, i, 10); solver.admittance(i, 4, 1., 10, P[4]*G1);
					jump1_compress_z(P, i, 11); solver.admittance(i, 4, 1., 11, P[5]*G1);
				}

///////////////////////////////////////////////////////
//...jump of all neaded moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
	
				if (B[k].link[NUM_PHASE] == -1) {
					B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_KAPPA)); 
					jump1(P, k, 0); 
					jump2(P, k, 1); solver.admittance(k, 1, g2); 
					jump_conj_cpy(k, 2, 0);
					jump_conj_cpy(k, 3, 1); 
				}
				if (B[k].link[NUM_PHASE] == -2) {
					B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[k].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, k);

					jump4_common_x(P, k, 6); 
					jump4_common_y(P, k, 7); 
					jump4_common_z(P, k, 8); 

					jump1_common_x(P, k, 9); 
					jump1_common_y(P, k, 10); 
					jump1_common_z(P, k, 11); 

					B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[k].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, k);

					jump4_compress_x(P, k, 6); solver.admittance(k, 5, 0., 6, P[3]); 
					jump4_compress_y(P, k, 7); solver.admittance(k, 5, 1., 7, P[4]); 
					jump4_compress_z(P, k, 8); solver.admittance(k, 5, 1., 8, P[5]); 

					jump1_compress_x(P, k,  9); solver.admittance(k, 4, 0.,  9, P[3]*G1);
					jump1_compress_y(P, k, 10); solver.admittance(k, 4, 1., 10, P[4]*G1);
					jump1_compress_z(P, k, 11); solver.admittance(k, 4, 1., 11, P[5]*G1);
				}

///////////////////////////////////////////////////////////////////////////////////
//...сшивка перемещений и давлений методом наименьших квадратов (в шкале давлений);
				if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -2) {
					jump_make_common (k, 6);
					solver.admittance(k, 6, 1., 5, -nd->nX[l]); //...вычисляем касательную силу;
					solver.admittance(k, 7, 1., 5, -nd->nY[l]);
					solver.admittance(k, 8, 1., 5, -nd->nZ[l]);
					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+5], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[i][0][m],   f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+8], solver.hh[k][0][m+8], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+7], solver.hh[k][0][m+7], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+6], solver.hh[k][0][m+6], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+4], solver.hh[k][0][m+4], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+4], f);
				}
				if (B[k].link[NUM_PHASE] == -1 && B[i].link[NUM_PHASE] == -2) {
					jump_make_common (i, 6);
					solver.admittance(i, 6, 1., 5, -nd->nX[l]); //...вычисляем касательную силу;
					solver.admittance(i, 7, 1., 5, -nd->nY[l]);
					solver.admittance(i, 8, 1., 5, -nd->nZ[l]);
					solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+1], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+6], solver.hh[i][0][m+6], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+7], solver.hh[i][0][m+7], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+8], solver.hh[i][0][m+8], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+2], solver.hh[k][0][m], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m], f);
				}
			}
			break;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе акустическая - механическая фаза энергетическим методом;
Num_State CAcou3D::transfer3(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double G1 = get_param(NUM_KAPPA+2)/get_param(NUM_KAPPA), f, P[6],
 			    g2 = 1./get_param(NUM_KAPPA);
     int	 m  = solver.id_norm;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

				f = nd->get_param(0, l);

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(complex));
				}

///////////////////////////////////////////////
//...jump of all neaded moments for this block;
				if (B[i].link[NUM_PHASE] == -1) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
					jump1(P, i, 0); 
					jump2(P, i, 1); solver.admittance(i, 1, g2); 
					jump_conj_cpy(i, 2, 0);
					jump_conj_cpy(i, 3, 1); 
				}
				if (B[i].link[NUM_PHASE] == -2) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, i);
					
					jump4_common_x(P, i, 6); 
					jump4_common_y(P, i, 7); 
					jump4_common_z(P, i, 8); 

					jump1_common_x(P, i, 9); 
					jump1_common_y(P, i, 10); 
					jump1_common_z(P, i, 11); 

					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, i);

					jump4_compress_x(P, i, 6); solver.admittance(i, 5, 0., 6, P[3]); 
					jump4_compress_y(P, i, 7); solver.admittance(i, 5, 1., 7, P[4]); 
					jump4_compress_z(P, i, 8); solver.admittance(i, 5, 1., 8, P[5]); 
					jump_make_common(i, 6);

					jump1_compress_x(P, i,  9); solver.admittance(i,  9, G1); solver.admittance(i, 4, 0.,  9, P[3]);
					jump1_compress_y(P, i, 10); solver.admittance(i, 10, G1); solver.admittance(i, 4, 1., 10, P[4]); 
					jump1_compress_z(P, i, 11); solver.admittance(i, 11, G1); solver.admittance(i, 4, 1., 11, P[5]); 
					jump_make_common(i, 9);
				}

///////////////////////////////////////////////////////
//...jump of all neaded moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
	
				if (B[k].link[NUM_PHASE] == -1) {
					jump1(P, k, 0); 
					jump2(P, k, 1); solver.admittance(k, 1, g2); 
					jump_conj_cpy(k, 2, 0);
					jump_conj_cpy(k, 3, 1); 
				}
				if (B[k].link[NUM_PHASE] == -2) {
					B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[k].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, k);

					jump4_common_x(P, k, 6); 
					jump4_common_y(P, k, 7); 
					jump4_common_z(P, k, 8); 

					jump1_common_x(P, k, 9); 
					jump1_common_y(P, k, 10); 
					jump1_common_z(P, k, 11); 

					B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[k].shape->parametrization_hess(P, 1);
					hessian_deriv_N (P, k);

					jump4_compress_x(P, k, 6); solver.admittance(k, 5, 0., 6, P[3]);
					jump4_compress_y(P, k, 7); solver.admittance(k, 5, 1., 7, P[4]);
					jump4_compress_z(P, k, 8); solver.admittance(k, 5, 1., 8, P[5]);
					jump_make_common(k, 6);

					jump1_compress_x(P, k,  9); solver.admittance(k,  9, G1); solver.admittance(k, 4, 0.,  9, P[3]);
					jump1_compress_y(P, k, 10); solver.admittance(k, 10, G1); solver.admittance(k, 4, 1., 10, P[4]); 
					jump1_compress_z(P, k, 11); solver.admittance(k, 11, G1); solver.admittance(k, 4, 1., 11, P[5]); 
					jump_make_common(k, 9);
				}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -2) {
					solver.to_equationER(i, solver.hh[i][0][m+2],  solver.hh[i][0][m+1],  f);
					solver.to_equationEL(k, solver.hh[k][0][m+11], solver.hh[k][0][m+8], -f);
					solver.to_equationEL(k, solver.hh[k][0][m+10], solver.hh[k][0][m+7], -f);
					solver.to_equationEL(k, solver.hh[k][0][m+9],  solver.hh[k][0][m+6], -f);
				}
				if (B[k].link[NUM_PHASE] == -1 && B[i].link[NUM_PHASE] == -2) {
					solver.to_equationER(i, solver.hh[i][0][m+9],  solver.hh[i][0][m+6],  f);
					solver.to_equationER(i, solver.hh[i][0][m+10], solver.hh[i][0][m+7],  f);
					solver.to_equationER(i, solver.hh[i][0][m+11], solver.hh[i][0][m+8],  f);
					solver.to_equationEL(k, solver.hh[k][0][m+2],  solver.hh[k][0][m+1], -f);
				}

/////////////////////////////////////////////////////
//...сшивка перемещений методом наименьших квадратов;
				if (B[i].link[NUM_PHASE] == -1 && B[k].link[NUM_PHASE] == -2) {
					solver.admittance(k, 6, 1., 5, -nd->nX[l]); //...вычисляем касательную силу;
					solver.admittance(k, 7, 1., 5, -nd->nY[l]);
					solver.admittance(k, 8, 1., 5, -nd->nZ[l]);
					solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+4], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+3], solver.hh[i][0][m+1], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+8], solver.hh[k][0][m+8], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+7], solver.hh[k][0][m+7], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+6], solver.hh[k][0][m+6], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+4], solver.hh[k][0][m+4], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+4], f);
				}
				if (B[k].link[NUM_PHASE] == -1 && B[i].link[NUM_PHASE] == -2) {
					solver.admittance(i, 6, 1., 5, -nd->nX[l]); //...вычисляем касательную силу;
					solver.admittance(i, 7, 1., 5, -nd->nY[l]);
					solver.admittance(i, 8, 1., 5, -nd->nZ[l]);
					solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+1], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+6], solver.hh[i][0][m+6], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+7], solver.hh[i][0][m+7], f);
					solver.to_transferTT(i, j, solver.hh[i][0][m+8], solver.hh[i][0][m+8], f);
					solver.to_transferTD(i, j, solver.hh[k][0][m+3], solver.hh[k][0][m+1], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+1], f);
				}
			}
			break;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии для акустической фазы;
Num_State CAcou3D::gram4_phase1(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double f, P[8];
		int m = solver.id_norm, id_zero;
		complex attm, hh, h;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////////////////////////////
//...realization of Zommerfield or common boundary condition for external or internal block;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			id_zero =   nd->get_param(1, l) == MIN_HIT;
			attm = comp(nd->get_param(0, l),  nd->get_param(1, l))*comp(0., get_param(NUM_KAPPA));
			hh   = comp(nd->get_param(2, l),  nd->get_param(3, l));
			f = (id_zero ? 1.:1./(1.+norm(attm)))*nd->get_param(4, l);

			if (id_zero) attm = comp(0., MIN_HIT);
			if (hh.y == MAX_HIT) hh.y = 0.; //...for SONEX absorber identification;
			if (hh.y == MIN_HIT) {          //...for VN corrector identification;
				 hh.y = hh.x; 
				 hh.x = nd->get_param(1, l); 
				 attm = comp(0.);
			}
			if (! id_zero) hh *= comp(0., -get_param(NUM_KAPPA+1)); //..frequency dependence condition!!!
			P[6] = real(attm); //...передаем аттмиданс поглощения через точку;
			P[7] = imag(attm); 

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));

///////////////////////////////////////////////////////////
//...вычисляем необходимые моменты коллокационного вектора;
			B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
			jump1(P, i, 0);
			jump_conj_cpy(i, 2, 0);

			if (id_zero) { //...граничное условие -- давление;
				P[3] = P[4] = 0.; P[5] = 1.;
				jump2(P, i, 1);
			}
			else //...граничное условие -- скорость или поглощение;
				jump3(P, i, 1);

////////////////////////////////////////////////////
//...граничное условие методом наименьших квадратов;
			if (id_zero) {
				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m], f);
				if (abs(h = hh) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+2], B[i].shape->get_NN(), h*f);
			}
				
//////////////////////////////
//...энергетические слагаемые;
			solver.to_equationEE(i, solver.hh[i][0][m+2], solver.hh[i][0][m+1], f);
			if (! id_zero && abs(h = hh) > EE)
				solver.to_equationEH(i, 0, solver.hh[i][0][m+2], h*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии для механической фазы;
Num_State CAcou3D::gram4_phase2(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_KAPPA+2), nu1 = get_param(NUM_SHEAR+1), hx, hy, hz, p4, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);

			hx = nd->get_param(0, l); 
			hy = nd->get_param(1, l); 
			hz = nd->get_param(2, l);
			p4 = nd->get_param(3, l);
			f  = nd->get_param(4, l);

			if (p4 == MIN_HIT || p4 == 2.) {
			  hy = nd->nY[l]*hx;
			  hz = nd->nZ[l]*hx;
			  hx *= nd->nX[l];
			}
			else
			if (p4 == 0.) hx = hy = hz = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));

/////////////////////////////////////////////
//...jump of all neaded displacement moments;
			B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
			B[i].shape->parametrization_hess(P, 1);
			hessian_deriv_N(P, i);

			jump1_common_x(P, i, 0); 
			jump1_common_y(P, i, 1); 
			jump1_common_z(P, i, 2); 

			jump4_common_x(P, i, 3);
			jump4_common_y(P, i, 4); 
			jump4_common_z(P, i, 5);

			B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
			B[i].shape->parametrization_hess(P, 1);
			hessian_deriv_N (P, i);

			jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
			jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);
			jump1_compress_z(P, i, 2); solver.admittance(i, 2, G1);

			jump4_compress_x(P, i, 3);
			jump4_compress_y(P, i, 4); 
			jump4_compress_z(P, i, 5);

			if (p4 == (double)(SKEWS_BND-SPECIAL_BND) || p4 == (double)(UPTAKE_BND-SPECIAL_BND)) {
				solver.admittance(i, 5, 0., 0, P[3]);
				solver.admittance(i, 5, 1., 1, P[4]); 
				solver.admittance(i, 5, 1., 2, P[5]);

				solver.admittance(i, 4, 0., 6, P[3]*G1);
				solver.admittance(i, 4, 1., 7, P[4]*G1);
				solver.admittance(i, 4, 1., 8, P[5]*G1);

				if (p4 == (double)(UPTAKE_BND-SPECIAL_BND)) { //...вычисляем поглощение;
					solver.admittance(i, 5, 1., 4, comp(0., get_param(NUM_SHEAR+3)*(1.-nu1)/(.5-nu1)));
					jump_conj_cpy	  (i, 6, 5);
				}
			}
			jump_make_common(i, 0);
			jump_make_common(i, 3);

////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов;
			if (p4 == (double)(SKEWS_BND-SPECIAL_BND))
			solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
			if (p4 == (double)(UPTAKE_BND-SPECIAL_BND))
			solver.to_equationDD(i, solver.hh[i][0][m+6], solver.hh[i][0][m+5], f);
			if (p4 == MIN_HIT || p4 == MAX_HIT) {
				solver.to_equationDD(i, solver.hh[i][0][m],	solver.hh[i][0][m],	f);
				solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
				solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);

				if (fabs(hx) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m], hx*G1*f);
				if (fabs(hy) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+1], hy*G1*f);
				if (fabs(hz) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+2], hz*G1*f);
			}
				
/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+3], (f *= G1));
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4], f);
			solver.to_equationEE(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5], f);

			if (1. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				if (fabs(hx) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m], -hx*f);
				if (fabs(hy) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m+1], -hy*f);
				if (fabs(hz) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m+2], -hz*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе акустическая - акустическая фаза энергетическим методом;
Num_State CAcou3D::transfer4_phase1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N && ! src) {
		int j, l, m = solver.id_norm;
		double f, P[6];
      complex D_p[4], h = comp(0.);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				f = nd->get_param(0, l);
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(complex));
				}

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
				B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
				jump1(P, i, 0); 
				jump2(P, i, 1);
				jump_conj_cpy(i, 2, 0);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_KAPPA)); 
				jump1(P, k, 0); 
				jump2(P, k, 1); 
				jump_conj_cpy(k, 2, 0);

//////////////////////////////////////////////////
//...сшивка давлений методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m], f);
				solver.to_transferTT(i, j, solver.hh[i][0][m+2], solver.hh[i][0][m], f);
				solver.to_transferTD(i, j, solver.hh[k][0][m+2], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m], f);

//////////////////////////////
//...энергетические слагаемые;
				solver.to_equationER(i, solver.hh[i][0][m+2], solver.hh[i][0][m+1],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+2], solver.hh[k][0][m+1], -f);
			}
			break;
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода на границе механическая - механическая фаза энергетическим методом;
Num_State CAcou3D::transfer4_phase2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      int    j, l, m = solver.id_norm;
      double f, P[6], G1 = get_param(NUM_KAPPA+2);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				f = nd->get_param(0, l);
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(complex));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
				B[i].shape->parametrization_hess(P, 1);
				hessian_deriv_N(P, i);

				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 
				jump1_common_z(P, i, 2); 

				jump4_common_x(P, i, 3);
				jump4_common_y(P, i, 4); 
				jump4_common_z(P, i, 5);

				B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
				B[i].shape->parametrization_hess(P, 1);
				hessian_deriv_N (P, i);
				
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);
				jump1_compress_z(P, i, 2); solver.admittance(i, 2, G1);

				jump4_compress_x(P, i, 3);
				jump4_compress_y(P, i, 4); 
				jump4_compress_z(P, i, 5);
				jump_make_common(i, 0);
				jump_make_common(i, 3);
				
////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
	
				B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+4)); 
				B[k].shape->parametrization_hess(P, 1);
				hessian_deriv_N(P, k);

				jump1_common_x(P, k, 0); 
				jump1_common_y(P, k, 1); 
				jump1_common_z(P, k, 2); 

				jump4_common_x(P, k, 3);
				jump4_common_y(P, k, 4); 
				jump4_common_z(P, k, 5);

				B[k].shape->set_shape(B[k].shape->get_R(), get_param(NUM_SHEAR+3)); 
				B[k].shape->parametrization_hess(P, 1);
				hessian_deriv_N (P, k);

				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);
				jump1_compress_z(P, k, 2); solver.admittance(k, 2, G1);

				jump4_compress_x(P, k, 3);
				jump4_compress_y(P, k, 4); 
				jump4_compress_z(P, k, 5);
				jump_make_common(k, 0);
				jump_make_common(k, 3);

#ifndef ___PRESSURE_SEWING___
//////////////////////////////////////////////////
//...сшивка давлений методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
#else
////////////////////////////////////////////////////////////////////////
//...сшивка перемещений методом наименьших квадратов (в шкале давлений);
				double g = f/sqr(get_param(NUM_KAPPA));
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], g);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], g);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], g);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], g);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], g);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], g);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], g);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], g);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], g);
#endif
/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+3],  f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+3], -f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+4], -f);

				solver.to_equationER(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+2], solver.hh[k][0][m+5], -f);
			}
			break;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////
//...istalling parameters of multi-phase media;
void CAcou3D::set_fasa_hmg(double Hz, double Ro0, double C0, double Ro1, double nju1, double G1)
{
  if (size_of_param() > NUM_KAPPA+1 && size_of_param() > NUM_SHEAR+4) {
      param[5] = Hz;
		param[NUM_KAPPA-2] = Ro0;						//...кг/м^3;
		param[NUM_KAPPA-1] = C0;						//...м/сек;
		param[NUM_KAPPA] = 2.*M_PI/C0*Hz;			//...1/м;
		param[NUM_KAPPA+1] = Ro0*2.*M_PI*Hz;		//...кГ*сек/м^4;
		param[NUM_KAPPA+2] = Ro0*sqr(2.*M_PI*Hz);	//...кГ/м^4;
      param[NUM_SHEAR-2] = Ro1;						//...кг/м^3;
      param[NUM_SHEAR-1] = Ro1*sqr(2.*M_PI*Hz);	//...кГ/м^4;
      param[NUM_SHEAR]	 = G1;						//...кГ/м^2;
      param[NUM_SHEAR+1] = nju1;
      param[NUM_SHEAR+2] = .25/(nju1-1.);
      param[NUM_SHEAR+3] = sqrt((.5-nju1)/(1.-nju1)*Ro1/G1)*2.*M_PI*Hz; //...1/м;
      param[NUM_SHEAR+4] = sqrt(Ro1/G1)*2.*M_PI*Hz;								   //...1/м;

//		double L = 2.1, ggg = param[NUM_KAPPA+1]; //...оценка звукового давления (в кГ/м^2) для 1 мм/сек;
//		ggg /= param[NUM_KAPPA];
//		ggg /= sin(L*param[NUM_KAPPA])*1000.;
//		ggg *= cos(L*param[NUM_KAPPA]);
  }
}
void CAcou3D::set_fasa_hmg(double Hz, double Ro0, double C0)
{
	if (size_of_param() > NUM_KAPPA+1) { //...акустическая среда;
      param[5] = Hz;
		param[NUM_KAPPA-2] = Ro0;						//...кг/м^3;
		param[NUM_KAPPA-1] = C0;						//...м/сек;
		param[NUM_KAPPA] = 2.*M_PI/C0*Hz;			//...1/м;
		param[NUM_KAPPA+1] = Ro0*2.*M_PI*Hz;		//...кГ*сек/м^4;
		param[NUM_KAPPA+2] = Ro0*sqr(2.*M_PI*Hz);	//...кГ/м^4;
	}
}
void CAcou3D::set_fasa_hmg(double Hz, double Ro1, double nju1, double G1)
{
  if (size_of_param() > NUM_SHEAR+4) { //...упругое тело;
      param[5] = Hz;
      param[NUM_SHEAR-2] = Ro1;						//...кг/м^3;
      param[NUM_SHEAR-1] = Ro1*sqr(2.*M_PI*Hz);	//...кГ/м^4;
      param[NUM_SHEAR]	 = G1;						//...кГ/м^2;
      param[NUM_SHEAR+1] = nju1;
      param[NUM_SHEAR+2] = .25/(nju1-1.);
      param[NUM_SHEAR+3] = sqrt((.5-nju1)/(1.-nju1)*Ro1/G1)*2.*M_PI*Hz; //...1/м;
      param[NUM_SHEAR+4] = sqrt(Ro1/G1)*2.*M_PI*Hz;								   //...1/м;
  }
}

//////////////////////////////////////////////////
//...counting header for solving acoustic problem;
Num_State CAcou3D::counting_header(Num_State Num_counting)
{
	int N_elem = UnPackInts(get_param(3)), i, j, k, n_rhs = max(2, src ? (int)src[0] : 1);
	char msg[201];
	
	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");
		sprintf(msg, "CAcou3D sample: N_sm = %d, N_mpl = %d, N_elem = %d, Hz = %g", N, 
				 UnPackInts(get_param(0)), N_elem, get_param(5));
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num_counting){
			case BASIC_COMPUT: Message("FEM Blocks...");	break;
		}
		Message(" ");
	}

//////////////////////////////////////////////
//...установливаем связи на границе включений;
	if (NUM_PHASE > 0 && ! solver.mode(NO_PHASE))
	for (k = 0; k < N; k++) if (B[k].link[0] >= NUM_PHASE) {
		for (j = B[k].link[0]; j > 0; j--) 
		if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE]) {
			B[k].link[j] = -i+SRF_STATE;
			link_add(B[k], i);
			link_add(B[k], SRF_STATE, NULL_STATE);
		}
		if (B[k].link[NUM_PHASE] <= -3) B[k].link[NUM_PHASE] = -1; //...коррекция парметров включений;
		if (B[k].mp && 1) //...разворачиваем мультиполи по оси X (make_common() описывает положение локальной системы координат);
			 B[k].mp[5] = M_PI_2;
	}

////////////////////////////////////
//...устанавливаем параметры задачи;
	set_fasa_hmg(param[5], param[NUM_KAPPA-2], param[NUM_KAPPA-1]);
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials !!!
	solver.n += 2+NUM_HESS;		  //<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);

/////////////////////////////////////////////////////////
//...делаем перенумерацию структуры и задаем размерность;
	if (! solver.struct_permutat(solver.id_change == EXTERN_STATE ? NULL_STATE : OK_STATE) || 
		 ! solver.inverse_index()) {
		return ERR_STATE;
	}
	for (k = 0; k < solver.N; k++)
		solver.set_dimension(k, freedom_block(k));

	solver.struct_init();

	if (solver.mode(FULLY_MODE)) { 
		solver.test_struct("1", 0);
		solver.test_struct("2", 1);
	}
   return OK_STATE;
}

//////////////////////////////////////
//...вычмсление невязки междк блоками;
void CAcou3D::block_descrap(char * OUT_FILE)
{
   int N_elem = (int)get_param(3);
	char msg[201];

////////////////////////////////////////
//...auxilliary grids for discrete norm;
   CGrid * bnd = CreateNodes();

	CGrid * block_bnd = CreateNodes();
			  block_bnd->add_params(3);

	CGrid * gauss_bnd = CreateNodes(GRID_QG_NODES);
			  gauss_bnd->add_params(1);

///////////////////////////////////////////
//...вычисление среднеквадраимчной невязки;
	FILE * OUT = OUT_FILE ? fopen(OUT_FILE, "w") : NULL;
	sprintf(msg, "Block descrapency...");
	Message(msg); fprintf(OUT, "%s\n", msg);

   double pp[6], P0[6], Po[12];
	int  k, l, i, j, j_surf, m;

	for (k = 0; k < N; k++) if (B[k].bar && B[k].link) {
		sprintf(msg, "block %4i: ", k);
		Message(msg); fprintf(OUT, "%s\n", msg);

		for (i = 0; i < B[k].link[0]; i++) if ((j = B[k].link[i+1]) >= 0) {
			bnd->zero_grid(); 
			m = block_comput(bnd, k, j, sqr(get_param(4)), j_surf, 1);
       
/////////////////////////////////
//...накапливаем граничные точки;
			if (bnd->geom)
			for (l = 0; l <  bnd->geom[0]; l++) {
				int num  = bnd->geom_element(l), num_n = bnd->geom[num+1],
					num_f = num_n+num, cnt = 0;

				if (bnd->geom[num] == GL_TRIANGLES) {
					P0[3] = bnd->nX[bnd->geom[num+2]];
					P0[4] = bnd->nY[bnd->geom[num+2]];
					P0[5] = bnd->nZ[bnd->geom[num+2]];
					for (; num < num_f; num++) {
						Po[cnt++] = bnd->X[bnd->geom[num+2]];
						Po[cnt++] = bnd->Y[bnd->geom[num+2]];
						Po[cnt++] = bnd->Z[bnd->geom[num+2]];
					}
					gauss_bnd->facet_QG(Po, N_elem, NULL_STATE, NULL_STATE);
					if (NUM_PHASE >= i+1) { //...коррекция квадратур;
						if (m && bar && bar->graph && -j_surf+SRF_STATE < bar->graph[0]) gauss_bnd->QG_tria_surface(bar->ce[-j_surf+SRF_STATE]->mp, Po);
						for (int lp = 0; lp < gauss_bnd->N; lp++) {
							pp[0] = gauss_bnd->get_param(0, lp);
							block_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp);
						}
					}
					gauss_bnd->add_buffer(gauss_bnd->N);
				}
				else 
				if (bnd->geom[num] == GL_QUAD_STRIP) {
					P0[3] = bnd->nX[bnd->geom[num+2]];
					P0[4] = bnd->nY[bnd->geom[num+2]];
					P0[5] = bnd->nZ[bnd->geom[num+2]];
					for (; num < num_f-3; num += 2) {
						Po[0] = bnd->X[bnd->geom[num+2]];
						Po[1] = bnd->Y[bnd->geom[num+2]];
						Po[2] = bnd->Z[bnd->geom[num+2]];
						Po[3] = bnd->X[bnd->geom[num+3]];
						Po[4] = bnd->Y[bnd->geom[num+3]];
						Po[5] = bnd->Z[bnd->geom[num+3]];
						Po[6] = bnd->X[bnd->geom[num+4]];
						Po[7] = bnd->Y[bnd->geom[num+4]];
						Po[8] = bnd->Z[bnd->geom[num+4]];
						Po[9] = bnd->X[bnd->geom[num+5]];
						Po[10] = bnd->Y[bnd->geom[num+5]];
						Po[11] = bnd->Z[bnd->geom[num+5]];

						gauss_bnd->facet_QG(Po, N_elem, OK_STATE, NULL_STATE);
						if (NUM_PHASE >= i+1) { //...коррекция квадратур;
							if (m && bar && bar->graph && -j_surf+SRF_STATE < bar->graph[0]) gauss_bnd->QG_quad_surface(bar->ce[-j_surf+SRF_STATE]->mp, NULL, Po);
							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								pp[0] = gauss_bnd->get_param(0, lp);
								block_bnd->add_new_point(gauss_bnd->X[lp], gauss_bnd->Y[lp], gauss_bnd->Z[lp], P0[3], P0[4], P0[5], pp); 
							}
						}
						gauss_bnd->add_buffer(gauss_bnd->N);
					}
				}
			}

/////////////////////////////////////
//...выисление невязки между блоками;
			if (NUM_PHASE >= i) {
				double F[6], sum = 0., norm = 0.;
				for (int lp = 0; lp < block_bnd->N; lp++) if (block_bnd->hit[lp]) {

					memset(F, 0, 6*sizeof(double));
					GetFuncAllValues(block_bnd->X[lp], block_bnd->Y[lp], block_bnd->Z[lp], F,   k, SPL_VALUE);
					GetFuncAllValues(block_bnd->X[lp], block_bnd->Y[lp], block_bnd->Z[lp], F+3, j, SPL_VALUE);

					sum  += (sqr(F[0]-F[3])+sqr(F[1]-F[4])+sqr(F[2]-F[5]))*block_bnd->get_param(0, lp);
					norm += (sqr(F[0])+sqr(F[1])+sqr(F[2]))*block_bnd->get_param(0, lp);
				}
				sum  = sqrt(sum);
				norm = sqrt(norm);

				sprintf(msg, "            block_bnd->N = %i  sum = %g", block_bnd->N, sum/(1.+norm));
				Message(msg); fprintf(OUT, "%s\n", msg);
			}
			block_bnd->add_buffer(block_bnd->N);
		}
	}
	sprintf(msg, "");
	Message(msg); fprintf(OUT, "%s\n", msg);
	if (OUT) fclose(OUT);

   delete block_bnd;
	delete gauss_bnd;
   delete bnd;
}

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
void CAcou3D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 0., 0., 1.}, L1 = 0.8, L2 = 1.1, L = 3., Vn = .001, f;
	complex zz;

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp) {
		int m = solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(complex));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
		case SPL_VALUE: { //...SPL in dB;
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				if (B[i].link[NUM_PHASE] == -1) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
					jump1(P, i, 0);
					zz = B[i].shape->potential(solver.hh[i][0][m], id_variant*2);
				}
				else if (B[i].link[NUM_PHASE] == -2) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N(P, i);

					jump4_common_x(P, i, 0); 
					jump4_common_y(P, i, 1); 
					jump4_common_z(P, i, 2); 

					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N(P, i);

					jump4_compress_z(P, i, 0);
					jump4_compress_z(P, i, 1);
					jump4_compress_z(P, i, 2);
					jump_make_common(i, 0);
					zz = B[i].shape->potential(solver.hh[i][0][m], id_variant*2);
				}
				F[0] = 20.*max(log10(abs(zz)*M_SQRT2*.25e5), 0.); 
		}		break;
      case PRESSURE_VALUE: { //...complex pressure (tzz);
				if (B[i].link[NUM_PHASE] == -1) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
					jump1(P, i, 0);
					zz = B[i].shape->potential(solver.hh[i][0][m], id_variant*2);
				}
				else
				if (B[i].link[NUM_PHASE] == -2) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N(P, i);

					jump4_common_z(P, i, 0); 
					jump4_common_z(P, i, 1); 
					jump4_common_z(P, i, 2); 

					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[i].shape->parametrization_hess(P, 1);
					hessian_deriv_N(P, i);

					jump4_compress_z(P, i, 0);
					jump4_compress_z(P, i, 1);
					jump4_compress_z(P, i, 2);
					jump_make_common(i, 0);
					zz = B[i].shape->potential(solver.hh[i][0][m+2], id_variant*2);
				}
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case FLUX_X_VALUE: { //...complex velocity (Vx and Ux);
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				if (B[i].link[NUM_PHASE] == -1) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
					jump2(P, i, 0);
					zz = B[i].shape->potential(solver.hh[i][0][m], id_variant*2)/comp(0., -get_param(NUM_KAPPA+1));
				}
				else
				if (B[i].link[NUM_PHASE] == -2) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[i].shape->parametrization_hess(P, 1);

					jump1_common_x(P, i, 0); 
					jump1_common_y(P, i, 1); 
					jump1_common_z(P, i, 2); 

					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[i].shape->parametrization_hess(P, 1);

					jump1_compress_x(P, i, 0);
					jump1_compress_y(P, i, 1);
					jump1_compress_z(P, i, 2);
					jump_make_common(i, 0);

					zz = B[i].shape->potential(solver.hh[i][0][m], id_variant*2);
//					zz = zz*comp(0., 2.*M_PI*get_param(5));
				}
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case FLUX_Y_VALUE: { //...complex velocity (Vy and Uy);
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				if (B[i].link[NUM_PHASE] == -1) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
					jump2(P, i, 0);
					zz = B[i].shape->potential(solver.hh[i][0][m], id_variant*2)/comp(0., -get_param(NUM_KAPPA+1));
				}
				else
				if (B[i].link[NUM_PHASE] == -2) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[i].shape->parametrization_hess(P, 1);

					jump1_common_x(P, i, 0); 
					jump1_common_y(P, i, 1); 
					jump1_common_z(P, i, 2); 

					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[i].shape->parametrization_hess(P, 1);

					jump1_compress_x(P, i, 0);
					jump1_compress_y(P, i, 1);
					jump1_compress_z(P, i, 2);
					jump_make_common(i, 0);

					zz = B[i].shape->potential(solver.hh[i][0][m+1], id_variant*2);
					zz = zz*comp(0., 2.*M_PI*get_param(5));
				}
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case FLUX_Z_VALUE: { //...complex velocity (Vz and Uz);
				P[5] = 1.; P[3] = P[4] = 0.;
				B[i].shape->norm_local(P+3);
				if (B[i].link[NUM_PHASE] == -1) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_KAPPA)); 
					jump2(P, i, 0);
					zz = B[i].shape->potential(solver.hh[i][0][m], id_variant*2)/comp(0., -get_param(NUM_KAPPA+1));
				}
				else
				if (B[i].link[NUM_PHASE] == -2) {
					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+4)); 
					B[i].shape->parametrization_hess(P, 1);

					jump1_common_x(P, i, 0); 
					jump1_common_y(P, i, 1); 
					jump1_common_z(P, i, 2); 

					B[i].shape->set_shape(B[i].shape->get_R(), get_param(NUM_SHEAR+3)); 
					B[i].shape->parametrization_hess(P, 1);

					jump1_compress_x(P, i, 0);
					jump1_compress_y(P, i, 1);
					jump1_compress_z(P, i, 2);
					jump_make_common(i, 0);

					zz = B[i].shape->potential(solver.hh[i][0][m+2], id_variant*2);
					zz = zz*comp(0., 2.*M_PI*get_param(5));
				}
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case SPL_ANALYT_RIGID_VALUE: { //...SPL in dB, analytical solution, rigid wall;
				Z = X;
				double  K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1)),
						kkk = get_param(NUM_KAPPA+2),
						kp0 = get_param(NUM_KAPPA), kp = get_param(NUM_SHEAR+3);
				double C0 = cos(kp0*(L-L2))/(K*kp), 
						 B0 = kp0*sin(kp0*(L-L2))/kkk;
				double A0 = K*kp*(B0*sin(kp*(L2-L1))+C0*cos(kp*(L2-L1))),
						 BB = kkk/kp0*(B0*cos(kp*(L2-L1))-C0*sin(kp*(L2-L1)));
				complex D = comp(0., -get_param(NUM_KAPPA+1))*Vn/(kp0*(A0*sin(kp0*L1)+BB*cos(kp0*L1)));

				if (Z < L1) zz = D*(A0*cos(kp0*(Z-L1))+BB*sin(kp0*(Z-L1))); else 
				if (Z < L2) zz = D*K*kp*(-B0*sin(kp*(Z-L2))+C0*cos(kp*(Z-L2))); else 
								zz = D*cos(kp0*(Z-L));
				F[0] = 20.*max(log10(abs(zz)*M_SQRT2*.25e5), 0.);
		}		break;
      case PRESS_ANALYT_RIGID_VALUE: { //...complex pressure (tzz), analytical solution, rigid wall;
				Z = X;
				double  K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1)),
						kkk = get_param(NUM_KAPPA+2),
						kp0 = get_param(NUM_KAPPA), kp = get_param(NUM_SHEAR+3);
				double C0 = cos(kp0*(L-L2))/(K*kp), 
						 B0 = kp0*sin(kp0*(L-L2))/kkk;
				double A0 = K*kp*(B0*sin(kp*(L2-L1))+C0*cos(kp*(L2-L1))),
						 BB = kkk/kp0*(B0*cos(kp*(L2-L1))-C0*sin(kp*(L2-L1)));
				complex D = comp(0., -get_param(NUM_KAPPA+1))*Vn/(kp0*(A0*sin(kp0*L1)+BB*cos(kp0*L1)));

				if (Z < L1) zz = D*(A0*cos(kp0*(Z-L1))+BB*sin(kp0*(Z-L1))); else 
				if (Z < L2) zz = D*K*kp*(-B0*sin(kp*(Z-L2))+C0*cos(kp*(Z-L2))); else 
								zz = D*cos(kp0*(Z-L));
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case FLUX_Z_ANALYT_RIGID_VALUE: { //...complex velocity (Vz), analytical solution, rigid wall;
				double  K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1)),
						kkk = get_param(NUM_KAPPA+2),
						kp0 = get_param(NUM_KAPPA), kp = get_param(NUM_SHEAR+3);
				double C0 = cos(kp0*(L-L2))/(K*kp), 
						 B0 = kp0*sin(kp0*(L-L2))/kkk;
				double A0 = K*kp*(B0*sin(kp*(L2-L1))+C0*cos(kp*(L2-L1))),
						 BB = kkk/kp0*(B0*cos(kp*(L2-L1))-C0*sin(kp*(L2-L1)));
				complex D = Vn/(kp0*(A0*sin(kp0*L1)+BB*cos(kp0*L1)));

				if (Z < L1) zz = D*kp0*(-A0*sin(kp0*(Z-L1))+BB*cos(kp0*(Z-L1))); else 
				if (Z < L2) zz = D*kkk*(B0*cos(kp*(Z-L2))+C0*sin(kp*(Z-L2))); else 
								zz = D*kp0*(-sin(kp0*(Z-L)));
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case FLUX_X_ANALYT_RIGID_VALUE: { //...complex velocity (Vz), analytical solution, rigid wall;
				double  K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1)),
						kkk = get_param(NUM_KAPPA+2),
						kp0 = get_param(NUM_KAPPA), kp = get_param(NUM_SHEAR+3);
				double C0 = cos(kp0*(L-L2))/(K*kp), 
						 B0 = kp0*sin(kp0*(L-L2))/kkk;
				double A0 = K*kp*(B0*sin(kp*(L2-L1))+C0*cos(kp*(L2-L1))),
						 BB = kkk/kp0*(B0*cos(kp*(L2-L1))-C0*sin(kp*(L2-L1)));
				complex D = Vn/(kp0*(A0*sin(kp0*L1)+BB*cos(kp0*L1)));

				if (X < L1) zz = D*kp0*(-A0*sin(kp0*(X-L1))+BB*cos(kp0*(X-L1))); else 
				if (X < L2) zz = D*kkk*(B0*cos(kp*(X-L2))+C0*sin(kp*(X-L2))); else 
								zz = D*kp0*(-sin(kp0*(X-L)));
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case SPL_ANALYT_ABSORB_VALUE: { //...SPL in dB, analytical solution, full absorption;
				double   K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1)),
						 kkk = get_param(NUM_KAPPA+2),
						 kp0 = get_param(NUM_KAPPA), kp = get_param(NUM_SHEAR+3);
				complex C0 = comp(cos(kp0*L2), -sin(kp0*L2))/(K*kp), 
						  B0 = kp0*comp(-sin(kp0*L2), -cos(kp0*L2))/kkk;
				complex A0 = K*kp*(B0*sin(kp*(L2-L1))+C0*cos(kp*(L2-L1))),
						  BB = kkk/kp0*(B0*cos(kp*(L2-L1))-C0*sin(kp*(L2-L1)));
				complex  D = comp(0., -get_param(NUM_KAPPA+1))*Vn/(kp0*(A0*sin(kp0*L1)+BB*cos(kp0*L1))), zz;

				if (Z < L1) zz = D*(A0*cos(kp0*(Z-L1))+BB*sin(kp0*(Z-L1))); else 
				if (Z < L2) zz = D*K*kp*(-B0*sin(kp*(Z-L2))+C0*cos(kp*(Z-L2))); else 
								zz = D*comp(cos(kp0*Z), -sin(kp0*Z));
				F[0] = 20.*max(log10(abs(zz)*M_SQRT2*.25e5), 0.);
		}		break;
      case PRESS_ANALYT_ABSORB_VALUE: { //...complex pressure (tzz), analytical solution, full absorption;
				double   K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1)),
						 kkk = get_param(NUM_KAPPA+2),
						 kp0 = get_param(NUM_KAPPA), kp = get_param(NUM_SHEAR+3);
				complex C0 = comp(cos(kp0*L2), -sin(kp0*L2))/(K*kp), 
						  B0 = kp0*comp(-sin(kp0*L2), -cos(kp0*L2))/kkk;
				complex A0 = K*kp*(B0*sin(kp*(L2-L1))+C0*cos(kp*(L2-L1))),
						  BB = kkk/kp0*(B0*cos(kp*(L2-L1))-C0*sin(kp*(L2-L1)));
				complex  D = comp(0., -get_param(NUM_KAPPA+1))*Vn/(kp0*(A0*sin(kp0*L1)+BB*cos(kp0*L1))), zz;

				if (Z < L1) zz = D*(A0*cos(kp0*(Z-L1))+BB*sin(kp0*(Z-L1))); else 
				if (Z < L2) zz = D*K*kp*(-B0*sin(kp*(Z-L2))+C0*cos(kp*(Z-L2))); else 
								zz = D*comp(cos(kp0*Z), -sin(kp0*Z));
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case FLUX_Z_ANALYT_ABSORB_VALUE: { //...complex velocity (Vz), analytical solution, full absorption;
				double   K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1)),
						 kkk = get_param(NUM_KAPPA+2),
						 kp0 = get_param(NUM_KAPPA), kp = get_param(NUM_SHEAR+3);
				complex C0 = comp(cos(kp0*L2), -sin(kp0*L2))/(K*kp), 
						  B0 = kp0*comp(-sin(kp0*L2), -cos(kp0*L2))/kkk;
				complex A0 = K*kp*(B0*sin(kp*(L2-L1))+C0*cos(kp*(L2-L1))),
						  BB = kkk/kp0*(B0*cos(kp*(L2-L1))-C0*sin(kp*(L2-L1)));
				complex  D = Vn/(kp0*(A0*sin(kp0*L1)+BB*cos(kp0*L1))), zz;

				if (Z < L1) zz = D*kp0*(-A0*sin(kp0*(Z-L1))+BB*cos(kp0*(Z-L1))); else 
				if (Z < L2) zz = D*kkk*(B0*cos(kp*(Z-L2))+C0*sin(kp*(Z-L2))); else 
								zz = D*kp0*comp(-sin(kp0*Z), -cos(kp0*Z));
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
		case ANALYT_VALUE: { //...аналитическое решение ONO BOX с полным поглощением на концевом торце (давление);
				double A_inv = 1./.6, B_inv = 1./.5, tx = .25*A_inv, ty = .25*B_inv, length = 2.5;
				int N = id_variant;
				//zz = Box_Pressure_Stair(tx, ty, A_inv, B_inv, length, N, N, X-.6, Y-.5, Z, 
				//								get_param(NUM_KAPPA), get_param(NUM_KAPPA-2)*get_param(NUM_KAPPA-1));
				//f  = abs(zz)*.25e5*M_SQRT2;
				//F[0] = 20.*log10(f >  1./*EE_ker*/ ? f : 1.);
		}		break;
		case ANALYT2VALUE: { //...аналитическое решение ONO BOX с полным поглощением на концевом торце (комплексная скорость Vz);
				double A_inv = 1./.6, B_inv = 1./.5, tx = .25*A_inv, ty = .25*B_inv, length = 2.5;
				int N = id_variant;
				//zz = Box_VelocityZ_Stair(tx, ty, A_inv, B_inv, length, N, N, X-.6, Y-.5, Z, 
				//								get_param(NUM_KAPPA), get_param(NUM_KAPPA-2)*get_param(NUM_KAPPA-1));
				//F[0] = real(zz);
				//F[1] = imag(zz);
		}		break;
		case ANALYT3VALUE: { //...аналитическое решение ONO BOX (давление);
				double A_inv = 1./.6, B_inv = 1./.5, tx = .25*A_inv, ty = .25*B_inv, length = 2.5;
				int N = id_variant;
				//zz = Box_Pressure_Stair_abs(tx, ty, A_inv, B_inv, length, N, N, X-.6, Y-.5, Z, 
				//								get_param(NUM_KAPPA), get_param(NUM_KAPPA-2)*get_param(NUM_KAPPA-1));
				//f  = abs(zz)*.25e5*M_SQRT2;
				//F[0] = 20.*log10(f >  1./*EE_ker*/ ? f : 1.);
		}		break;
		case ANALYT4VALUE: { //...аналитическое решение ONO BOX (комплексная скорость Vz);
				double A_inv = 1./.6, B_inv = 1./.5, tx = .25*A_inv, ty = .25*B_inv, length = 2.5;
				//int N = id_variant;
				//zz = Box_VelocityZ_Stair_abs(tx, ty, A_inv, B_inv, length, N, N, X-.6, Y-.5, Z, 
				//								get_param(NUM_KAPPA), get_param(NUM_KAPPA-2)*get_param(NUM_KAPPA-1));
				//F[0] = real(zz);
				//F[1] = imag(zz);
		}		break;
		case ANALYT5VALUE: { //...аналитическое решение в трубе с излучателем (давление);
				double t0 = .25, RR_inv = 1./.55, length = 2.5;
				int N = id_variant;
				//zz = Cyl_Pressure_Stair(t0, RR_inv, length, F+3, N, sqrt(sqr(X)+sqr(Y)), Z, 
				//								get_param(12)*get_param(6), get_param(10)*get_param(11));
				//f  = abs(zz)*.25e5*M_SQRT2;
				//F[0] = 20.*log10(f >  1./*EE_ker*/ ? f : 1.);
		}		break;
		case ANALYT6VALUE: { //...аналитическое решение в трубе с излучателем (комплексная скорость Vz);
				double t0 = .25, RR_inv = 1./.55, length = 2.5;
				int N = id_variant;
				//zz = Cyl_VelocityZ_Stair(t0, RR_inv, length, F+3, N, sqrt(sqr(X)+sqr(Y)), Z, 
				//								get_param(NUM_KAPPA), get_param(NUM_KAPPA-2)*get_param(NUM_KAPPA-1));
				//F[0] = real(zz);
				//F[1] = imag(zz);
		}		break;
////////////////////////////////////////////////
//...для тестирования акустического виброкласса;
      case 10000: { //...SPL in dB, analytical solution, rigid wall, simple acoustic;
				zz = comp(0., -get_param(NUM_KAPPA+1)/get_param(NUM_KAPPA))*cos(get_param(NUM_KAPPA)*(Z-L))/sin(get_param(NUM_KAPPA)*L)*Vn;
				F[0] = 20.*max(log10(abs(zz)*M_SQRT2*.25e5), 0.);
		}		break;
      case 10003: { //...complex pressure (tzz), analytical solution, rigid wall, simple acoustic;
				Z = X;
				zz = comp(0., -get_param(NUM_KAPPA+1)/get_param(NUM_KAPPA))*cos(get_param(NUM_KAPPA)*(Z-L))/sin(get_param(NUM_KAPPA)*L)*Vn;
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case 10006: { //...complex velocity (Vz), analytical solution, rigid wall, simple acoustic;
				zz = -sin(get_param(NUM_KAPPA)*(Z-L))/sin(get_param(NUM_KAPPA)*L)*Vn;
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case 20000: { //...SPL in dB, analytical solution, full absorption, simple acoustic;
				zz = get_param(NUM_KAPPA-2)*get_param(NUM_KAPPA-1)*comp(cos(get_param(NUM_KAPPA)*Z), sin(get_param(NUM_KAPPA)*Z))*Vn;
				F[0] = 20.*max(log10(abs(zz)*M_SQRT2*.25e5), 0.);
		}		break;
      case 20003: { //...complex pressure (tzz), analytical solution, full absorption, simple acoustic;
				zz = get_param(NUM_KAPPA-2)*get_param(NUM_KAPPA-1)*comp(cos(get_param(NUM_KAPPA)*Z), sin(get_param(NUM_KAPPA)*Z))*Vn;
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case 20006: { //...complex velocity (Vz), analytical solution, full absorption, simple acoustic;
				zz = comp(cos(get_param(NUM_KAPPA)*Z), sin(get_param(NUM_KAPPA)*Z))*Vn;
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case 30000: { //...SPL in dB, analytical solution, rigid wall, simple vibro-velocity;
				double K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1));
				zz = comp(0., sqrt(get_param(NUM_SHEAR-2)*K))*cos(get_param(NUM_SHEAR+3)*(Z-L))/sin(get_param(NUM_SHEAR+3)*L)*Vn;
				F[0] = 20.*max(log10(abs(zz)*M_SQRT2*.25e5), 0.);
		}		break;
      case 30003: { //...complex pressure (tzz), analytical solution, rigid wall, simple vibro-velocity;
				double K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1));
				zz = comp(0., sqrt(get_param(NUM_SHEAR-2)*K))*cos(get_param(NUM_SHEAR+3)*(Z-L))/sin(get_param(NUM_SHEAR+3)*L)*Vn;
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case 30006: { //...complex velocity (Vz), analytical solution, rigid wall, simple vibro-velocity;
				Z = X;
				zz = -sin(get_param(NUM_SHEAR+3)*(Z-L))/sin(get_param(NUM_SHEAR+3)*L)*Vn;
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case 40000: { //...SPL in dB, analytical solution, rigid wall, simple vibro-pressure;
				zz = -sin(get_param(NUM_SHEAR+3)*(Z-L))/sin(get_param(NUM_SHEAR+3)*L);
				F[0] = 20.*max(log10(abs(zz)*M_SQRT2*.25e5), 0.);
		}		break;
      case 40003: { //...complex pressure (tzz), analytical solution, rigid wall, simple vibro-pressure;
				zz = -sin(get_param(NUM_SHEAR+3)*(Z-L))/sin(get_param(NUM_SHEAR+3)*L);
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case 40006: { //...complex velocity (Vz), analytical solution, rigid wall, simple vibro-pressure;
				double K = 2.*(1.-get_param(NUM_SHEAR+1))*get_param(NUM_SHEAR)/(1.-2.*get_param(NUM_SHEAR+1));
				zz = comp(0., 1./sqrt(get_param(NUM_SHEAR-2)*K))*cos(get_param(NUM_SHEAR+3)*(Z-L))/sin(get_param(NUM_SHEAR+3)*L);
				F[0] = real(zz);
				F[1] = imag(zz);
		}		break;
      case PROCESSOR_VALUE: {
				solver.hh[i][0][m][0] = comp(1.);
//				F[0] = (int)real(B[i].shape->potential(solver.hh[i][0][m], id_variant*2));
				F[0] = (int)real(solver.hh[i][0][1][0]);
		}		break;
		default: F[0] = i; F[1] = 0.;
		}
		B[i].shape->make_common(P);
		B[i].shape->norm_common(P+3);
	}
}

////////////////////////////////////////////////////////////////////////////////////
//...подготовка данных для визуализации в формате Surfer (для акустического класса);
void CAcou3D::GetSurferFormat(FILE * SURF, FILE * SURF1, FILE * SURF2, CGrid * nd, Num_Value _FMF, int id_variant, int id_axis, int iparam)
{
  size_t res;
  if (SURF && SURF1 && nd && nd->N > 0 && nd->N1 > 0) {
		short int i0 = (short int)nd->N,
					 j0 = (short int)nd->N1;
		double * out_F = (double *)new_struct(2004*sizeof(double)), X, Y, Z,
					min1F = 0., max1F = 1.,
					min2F = 0., max2F = 1.;
		int hit;

		if (iparam == SPECIAL_STATE) { //...симметричная (удвоенная относительно правого конца) область; 
			if (id_axis == AXIS_X || id_axis == AXIS_Z) nd->X[nd->N-1]  = nd->X[(nd->N-1)/2]*2.-nd->X[0]; 
			if (id_axis == AXIS_Y || id_axis == AXIS_Z) nd->Y[nd->N1-1] = nd->Y[(nd->N1-1)/2]*2.-nd->Y[0]; 
		}
		else
		if (iparam == SPECIAL3STATE) { //...симметричная (удвоенная относительно правого конца X) область; 
			if (id_axis == AXIS_X || id_axis == AXIS_Z) nd->X[nd->N-1]  = nd->X[(nd->N-1)/2]*2.-nd->X[0]; 
			if (id_axis == AXIS_Y)							  nd->Y[nd->N1-1] = nd->Y[(nd->N1-1)/2]*2.-nd->Y[0]; 
		}
		else
		if (iparam == SPECIAL2STATE) { //...симметричная (удвоенная относительно левого конца) область; 
			if (id_axis == AXIS_X || id_axis == AXIS_Z) nd->X[0] = nd->X[(nd->N-1)/2]*2.-nd->X[nd->N-1]; 
			if (id_axis == AXIS_Y || id_axis == AXIS_Z) nd->Y[0] = nd->Y[(nd->N1-1)/2]*2.-nd->Y[nd->N1-1]; 
		}
		else
		if (iparam == SPECIAL4STATE) { //...симметричная (удвоенная относительно левого конца X) область; 
			if (id_axis == AXIS_X || id_axis == AXIS_Z) nd->X[0] = nd->X[(nd->N-1)/2]*2.-nd->X[nd->N-1]; 
			if (id_axis == AXIS_Y)							  nd->Y[0] = nd->Y[(nd->N1-1)/2]*2.-nd->Y[nd->N1-1]; 
		}
		res = fwrite("DSBB", sizeof(char)*4,  1, SURF);
		res = fwrite(& i0,   sizeof(short int), 1, SURF); res = fwrite(& j0,         sizeof(short int), 1, SURF);
		res = fwrite(nd->X,  sizeof(double),  1, SURF); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF);
		res = fwrite(nd->Y,  sizeof(double),  1, SURF); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF);
		res = fwrite(& min1F, sizeof(double), 1, SURF); res = fwrite(& max1F,        sizeof(double),  1, SURF);

		res = fwrite("DSBB", sizeof(char)*4,  1, SURF1);
		res = fwrite(& i0,   sizeof(short int), 1, SURF1); res = fwrite(& j0,         sizeof(short int), 1, SURF1);
		res = fwrite(nd->X,  sizeof(double),  1, SURF1); res = fwrite(nd->X+nd->N-1,  sizeof(double),  1, SURF1);
		res = fwrite(nd->Y,  sizeof(double),  1, SURF1); res = fwrite(nd->Y+nd->N1-1, sizeof(double),  1, SURF1);
		res = fwrite(& min2F, sizeof(double), 1, SURF1); res = fwrite(& max2F,        sizeof(double),  1, SURF1);
		if (iparam == SPECIAL_STATE) { 
			if (id_axis == AXIS_X || id_axis == AXIS_Z) nd->X[nd->N-1]  = nd->X[0]; 
			if (id_axis == AXIS_Y || id_axis == AXIS_Z) nd->Y[nd->N1-1] = nd->Y[0]; 
		}
		else
		if (iparam == SPECIAL2STATE) {
			if (id_axis == AXIS_X || id_axis == AXIS_Z) nd->X[0] = nd->X[nd->N-1]; 
			if (id_axis == AXIS_Y || id_axis == AXIS_Z) nd->Y[0] = nd->Y[nd->N1-1]; 
		}
		min1F = min2F = MAX_HIT;
		max1F = max2F = MIN_HIT;
		//if (_FMF == ANALYT5VALUE || 
		//	 _FMF == ANALYT6VALUE) Fourier_Bessel_utils(out_F+3, 2000);
		for (int j = 0; j < nd->N1; j++)
		for (int i = 0; i < nd->N;  i++) {
			  if (id_axis == AXIS_X) {
					Y = nd->X[i];
					Z = nd->Y[j];
					X = nd->N2 == 1 ? nd->Z[0] : 0.;
			  }
			  else
			  if (id_axis == AXIS_Y) {
					Z = nd->X[i];
					X = nd->Y[j];
					Y = nd->N2 == 1 ? nd->Z[0] : 0.;
			  }
			  else
			  if (id_axis == AXIS_Z) {
					X = nd->X[i];
					Y = nd->Y[j];
					Z = nd->N2 == 1 ? nd->Z[0] : 0.;
			  }
			  else
			  if (id_axis == AXIS_SPH) {
					Z = nd->N2 == 1 ? nd->Z[0] : 1.;
					X = cos(nd->X[i])*sin(nd->Y[j])*Z*1.;
					Y = sin(nd->X[i])*sin(nd->Y[j])*Z*2.;
					Z = cos(nd->Y[j])*Z*3.;
			  }

			  float ff = NOT_HIT;
			  if (nd->hit) hit = nd->hit[i+j*nd->N];
			  if (hit != -1)  {
					GetFuncAllValues(X, Y, Z, out_F, hit, _FMF, id_variant);

					ff = (float)out_F[0];
					if (_FMF < 0) ff = (float)hit;
					res = fwrite(& ff, sizeof(float), 1, SURF);
					if (min1F > out_F[0]) min1F = out_F[0];
					if (max1F < out_F[0]) max1F = out_F[0];

					ff = (float)out_F[1];
					if (_FMF < 0) ff = (float)hit;
					res = fwrite(& ff, sizeof(float), 1, SURF1);
					if (min2F > out_F[1]) min2F = out_F[1];
					if (max2F < out_F[1]) max2F = out_F[1];
			  }
			  else {
					res = fwrite(& ff, sizeof(float), 1, SURF);
					res = fwrite(& ff, sizeof(float), 1, SURF1);
			  }
		}
		delete_struct(out_F);

//////////////////////////////////////////////////////////////
//...перезапись максимального и минимального значения функции;
		res = fseek(SURF, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
		res = fwrite(& min1F, sizeof(double), 1, SURF);
		res = fwrite(& max1F, sizeof(double), 1, SURF);
		res = fseek(SURF, 0L, SEEK_END);

		res = fseek(SURF1, sizeof(char)*4+sizeof(short int)*2+sizeof(double)*4, SEEK_SET);
      res = fwrite(& min2F, sizeof(double), 1, SURF1);
      res = fwrite(& max2F, sizeof(double), 1, SURF1);
      res = fseek(SURF1, 0L, SEEK_END);
  }
}
#undef  Message
