#include "stdafx.h"

#include "shapes.h"
#include "ccohes2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CCohes2D::NUM_ADHES = 3;
int CCohes2D::NUM_SHEAR = 6;
int CCohes2D::NUM_SHIFT = 5;
int CCohes2D::MAX_PHASE = 2;
int CCohes2D::NUM_HESS  = 8;
int CCohes2D::regul = 1;

//////////////////////////////////
//...initialization of the blocks;
int CCohes2D::block_shape_init(Block<double> & B, Num_State id_free)
{
	int k, m;
   if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));

		extern int gradient_model;
		if (gradient_model)	{//...using gradient displacements;
			if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->add_shape(CreateShape<double>(SK2D_ELLI_SHAPE));
			else											   B.shape->add_shape(CreateShape<double>(SK2D_POLY_SHAPE/*SK2D_BEAMZ_SHAPE*/));
		}

////////////////////////
//...setting parameters;
		B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
		if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
		else												B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
      
		B.shape->init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 2);
		if (B.link[NUM_PHASE] == -2) //...another degree of multipoles for nclusion!!!
		B.shape->init1(UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, 2);

//////////////////////////////////////////////////////////////////////////////
//...setting acselerator, local system of coordinate and init parametrization;
         B.shape->set_local(B.mp+1);
         B.shape->release  ();
   }

///////////////////////////////////////////////
//...setting cohesion parameter and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...������������� ������� � ������ �����������;
//			if (B.link[NUM_PHASE] == -1) B.mp[4] = M_PI/3.;

			B.shape->set_local(B.mp+1);
			B.shape->set_shape(NUM_MPLS, get_param(1)*fabs( B.mp[7]));
			if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			else												B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
		}
		else
		if (id_free == OK_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->set_potential(solver.hh[k = (int)(&B-B.B)][0][m], m);
		else
		if (id_free == NO_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->get_potential(solver.hh[k = (int)(&B-B.B)][0][solver.id_norm+m], m);
		else
		if (id_free == NULL_STATE) //...������������� ����������� (� ������ �������� �������, ��������);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

//////////////////////////////////////////////
//...realization of common displacements (Rx);
void CCohes2D::jump1_common_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	B[i].shape->cpy_x     (0, B[i].shape->deriv);
	B[i].shape->cpy       (0, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G0, P[0]*alpha*G0);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*alpha*G0, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	jump_admittance	(1, i, m-solver.id_norm, 0.);
	B[i].shape->adm_xy(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), C0);
	B[i].shape->adm_xx(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), C0);
	B[i].shape->adm	(1, ptr, -G0);
}

////////////////////////////////////////////////////////////////////
//...additional inclusion of cohesion compression displacement (ux);
void CCohes2D::jump1_compress_x(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;
	B[i].shape->adm_xy(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -C0);
	B[i].shape->adm_xx(1, B[i].shape->FULL(solver.hh[i][0][m], 1), -C0);
}

//////////////////////////////////////////////
//...realization of common displacements (Ry);
void CCohes2D::jump1_common_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	B[i].shape->cpy_y     (0, B[i].shape->deriv);
	B[i].shape->cpy       (0, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G0, P[1]*alpha*G0);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*alpha*G0, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	jump_admittance	(1, i, m-solver.id_norm, 0.);
	B[i].shape->adm_xy(1, B[i].shape->FULL(solver.hh[i][0][m], 1), C0);
	B[i].shape->adm_yy(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), C0);
	B[i].shape->adm	(1, ptr, -G0);
}

////////////////////////////////////////////////////////////////////
//...additional inclusion of cohesion compression displacement (uy);
void CCohes2D::jump1_compress_y(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;
	B[i].shape->adm_xy(1, B[i].shape->FULL(solver.hh[i][0][m], 1), -C0);
	B[i].shape->adm_yy(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -C0);
}

////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (Rx);
void CCohes2D::jump2_common_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[4]);

	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*alpha*G0, 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]*(2.*alpha+1.)*G0);
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]*(1.+alpha)*G0);
	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1]*alpha*G0, 0.);
	B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[4]*alpha*G0);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1), 0., C0);
	B[i].shape->adm_x     (1, ptr, -P[3]*G0);
	B[i].shape->adm_y     (1, ptr, -P[4]*G0);
	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 0., C0);
}

//////////////////////////////////////////////////////////////////////////////////////////
//...additional inclusion of normal derivatives of cohesion compression displacement (ux);
void CCohes2D::jump2_compress_x(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;

	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1),		B[i].shape->FULL(solver.hh[i][0][num_hess], 1),	  1., -C0);
	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
}

////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (Ry);
void CCohes2D::jump2_common_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	B[i].shape->cpy_yy	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[3]);

	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*alpha*G0, 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]*(1.+alpha)*G0);
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]*(2.*alpha+1.)*G0);
	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[0]*alpha*G0, 0.);
	B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[3]*alpha*G0);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 0., C0);
	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 1), 0., C0);
	B[i].shape->adm_x     (1, ptr, -P[3]*G0);
	B[i].shape->adm_y     (1, ptr, -P[4]*G0);
}

//////////////////////////////////////////////////////////////////////////////////////////
//...additional inclusion of normal derivatives of cohesion compression displacement (uy);
void CCohes2D::jump2_compress_y(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;

	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 1), 1., -C0);
}

/////////////////////////////////////////////////////////////////
//...realization of surface forces for common displacements (Px);
void CCohes2D::jump4_common_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;
	B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[4]);

	B[i].shape->admittance(0, B[i].shape->deriv, NULL, alpha, 0.);
	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0], 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]);
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]*(1.+alpha));

	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1], 0.);
	B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[3]*(-2.*alpha-1.));
	B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[4]*( 1.+alpha));

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1), 0., C0);
	B[i].shape->adm_x     (1, ptr, -P[3]*2.);
	B[i].shape->adm_y     (1, ptr, -P[4]);
	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 0., C0);
	B[i].shape->adm_x     (1, ptr, -P[4]);
}

//////////////////////////////////////////////////////////////////////////////
//...additional inclusion of surface forces for compression displacement (ux);
void CCohes2D::jump4_compress_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = get_param(NUM_SHEAR+1+shift)/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;

	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1), 1., -C0);
	B[i].shape->adm_x     (1, ptr, P[3]*alpha);
	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
	B[i].shape->adm_y     (1, ptr, P[3]*alpha);
}

/////////////////////////////////////////////////////////////////
//...realization of surface forces for common displacements (Py);
void CCohes2D::jump4_common_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;
	B[i].shape->cpy_yy	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[3]);

	B[i].shape->admittance(0, B[i].shape->deriv, NULL, alpha, 0.);
	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1], 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]*(1.+alpha));
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]);
	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[0], 0.);
	B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[4]*(-2.*alpha-1.));
	B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[3]*( 1.+alpha));

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 0., C0);
	B[i].shape->adm_y     (1, ptr, -P[3]);
	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 1), 0., C0);
	B[i].shape->adm_x     (1, ptr, -P[3]);
	B[i].shape->adm_y     (1, ptr, -P[4]*2.);
}

//////////////////////////////////////////////////////////////////////////////
//...additional inclusion of surface forces for compression displacement (uy);
void CCohes2D::jump4_compress_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = get_param(NUM_SHEAR+1+shift)/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;

	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
	B[i].shape->adm_x     (1, ptr, P[4]*alpha);
	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 1), 1., -C0);
	B[i].shape->adm_y     (1, ptr, P[4]*alpha);
}

///////////////////////////////////////////////
//...realization normal derivatives of hessian;
void CCohes2D::hessian_deriv_N(int k, double * P, int i)
{
	int num_hess = NUM_HESS+solver.id_norm;
	B[i].shape->set_norm_cs(P);

	B[i].shape->cpy_xx(k);
	B[i].shape->deriv_N(k);
	B[i].shape->cpy_xx(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k));

	B[i].shape->cpy_xy(k);
	B[i].shape->deriv_N(k);
	B[i].shape->cpy_xy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k, 1));

	B[i].shape->cpy_yy(k);
	B[i].shape->deriv_N(k);
	B[i].shape->cpy_yy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], k));
}

///////////////////////////////////////
//...composition of collocation vector;
void CCohes2D::jump_admittance(int l, int i, int m, double adm_re, int k, double adm_im)
{
	m += solver.id_norm;
	if (m >= 0 && adm_re != 1.) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),    NULL, adm_re, 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, adm_re, 0.);
	}
	if (k >= 0) {
		k += solver.id_norm;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),	 B[i].shape->FULL(solver.hh[i][0][k], l),	   1., adm_im);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][k], l, 1), 1., adm_im);
	}
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CCohes2D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0. };
		B[i].shape->norm_local(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

void CCohes2D::jump_make_common(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0.};
		B[i].shape->norm_common(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...������������ ������� ����� � ������ ����������� ������� (������� �������������� ������);
Num_State CCohes2D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, f, P[6], TX, TY, hx;
      int	 m  = solver.id_norm, id_dir, k, j, shift;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			f  = nd->get_param(3, l);

			for (k  = (int)nd->get_param(4, l), j = 0; j < solver.JR[i][0]; j++) 
			if ( k == solver.JR[i][j+solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = 0.;       P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...��������� ��������� ������� �������������� ������ � �������� ��������������� �����;
				TX = TY = hx = 0.;
				switch (abs(id_dir = (int)nd->get_param(2, l))) {
					case 1: TX =  AX; hx = -G1*AX; break;
					case 2: TX = -AX; hx =  G1*AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

////////////////////////////////////////////////////////////////////////////////////
//...������� ������ ��� ������������ ������������ ���� ������� ���������� ���������;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(1, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(1, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

				if (fabs(hx) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m+4],  hx*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m+5],  hx*f);

				  solver.to_equationHL(k, 0, solver.hh[k][0][m+4], -hx*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+5], -hx*f);
				}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...��������� ������� ��� ������������ ������������ ���� ����������� (������������� � ������ ������������� ���������);
				if (regul) {
					if (id_dir == 1 || id_dir == 2) {
						solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
						solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*f*.5);

						solver.to_equationDD(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_equationHH(k, 0, solver.hh[k][0][m+4], -hx*f*.5);
					}
					if (id_dir == 3 || id_dir == 4) {
						solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
						solver.to_equationDD(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], f);
					}
					//����������� ������������� ������ ������???
				}

//////////////////////////////////////////////////////////
//...������ ����������� ���� ������� ���������� ���������;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

/////////////////////////////////////////////////
//...�������������� ��������� (������������� G1);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
	
				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////
//...������������ ������ �������� � ������ ����������� ������� �� ������� ���;
Num_State CCohes2D::transfer3(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], Ai, Ak, Bi, Bk;
      int m = solver.id_norm, shift;

////////////////////////////////////////
//...������������ ������������� �������;
		Ai = Ak = get_param(NUM_ADHES)*.5-(Bi = Bk = get_param(NUM_ADHES+1)*.5);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = 0.;        P[5] = 0.;
				f = nd->get_param(0, l);
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_common_x (P, i, 4);
					jump2_common_y (P, i, 5);
				}
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_compress_x(P, i, 4); solver.admittance(i, 4, G1);
					jump2_compress_y(P, i, 5); solver.admittance(i, 5, G1); 
				}
				jump_make_common(i, 0);
				jump_make_common(i, 2);
				jump_make_common(i, 4);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_common_x (P, k, 4);
					jump2_common_y (P, k, 5);
				}
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_compress_x(P, k, 4); solver.admittance(k, 4, G1); 
					jump2_compress_y(P, k, 5); solver.admittance(k, 5, G1); 
				}
				jump_make_common(k, 0);
				jump_make_common(k, 2);
				jump_make_common(k, 4);

//////////////////////////////////////////////////////////////////////////
//...������ ������� � ���������� ����������� ������� ���������� ���������;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

////////////////////////////////////////////////////////////////
//...��������� ������������� ������� ������� (������������� G1);
				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], -f*(Ai*sqr(P[3])+Bi*sqr(P[4])));
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], -f*(Ak*sqr(P[3])+Bk*sqr(P[4])));

					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+5], -f*Ai*Bi*P[3]*P[4]*2.);
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+5], -f*Ak*Bk*P[3]*P[4]*2.);

					solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], -f*(Ai*sqr(P[4])+Bi*sqr(P[3])));
					solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], -f*(Ak*sqr(P[4])+Bk*sqr(P[3])));
				}

/////////////////////////////////////////////////
//...�������������� ��������� (������������� G1);
				if (B[i].link[NUM_PHASE] != B[k].link[NUM_PHASE] && B[i].link[NUM_PHASE] < -1) f = -f;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CCohes2D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
//		nd->TestGrid("nodes.bln", 0.0005, 0., 0., 0., AXIS_Z, 1);
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			
			hx = nd->get_param(0, l);
			hy = nd->get_param(1, l);
			p3 = nd->get_param(2, l);
			f  = nd->get_param(3, l);

			if (p3 == MIN_HIT || p3 == 2.) {
			  hy  = nd->nY[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p3 == NUMS_BND) { //...����������� ������ -- ��������� ����������;
				double nju = get_param(NUM_SHEAR+1), G1 = .5/get_param(NUM_SHEAR), 
						 AAA = -nju*G1, BBB = (1.-nju)*G1;
				hx = nd->X[l]*BBB;
				hy = nd->Y[l]*AAA; 
			}
			else 
			if (p3 == 0.) hx = hy = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////////////////////////
//...jump of all neaded displacement moments;
			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
			B[i].shape->parametrization_hess(P, 1);

			hessian_deriv_N(1, P, i);
			jump1_common_x (P, i, 0);
			jump1_common_y (P, i, 1);

			jump4_common_x (P, i, 2);
			jump4_common_y (P, i, 3);

			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
			B[i].shape->parametrization_hess(1, P, 1);

			hessian_deriv_N (1, P, i);
			jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
			jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

			jump4_compress_x(P, i, 2); 
			jump4_compress_y(P, i, 3); 

			jump_make_common(i, 0);
			jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////////////////
//...��������� ������� ������� ���������� ��������� (������������ ������������);
			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.); 

				if (fabs(hx) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);

				if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+5], hy*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
			}

///////////////////////////////////////////////////
//...��������� ������������ ���� (�� ���� �������);
			solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
			solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.); 
			solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
			solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);

/////////////////////////////////////////////////
//...�������������� ��������� (������������� G1);
			solver.to_equationEE(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
			solver.to_equationEE(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...������������ ������ �������� � ������ ����������� �������;
Num_State CCohes2D::transfer4(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6];
      int m = solver.id_norm, shift;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = 0.;        P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

/////////////////////////////////////////////////////////////////////////////
//...������ ������� ������� ���������� ��������� (������������ ������������);
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(1, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(1, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

////////////////////////////////////////////////////////////////////
//...������ ������� ������� ���������� ��������� (����������� ����);
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

/////////////////////////////////////////////////
//...�������������� ��������� (������������� G1);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////
//...�������������� ��� �� �������� ������ ����� ��� ������������� ������;
Num_State CCohes2D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
      int l, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;
      double alpha = get_param(NUM_SHEAR+shift+1)/(.5-get_param(NUM_SHEAR+shift+1)), 
					 G0 = get_param(NUM_SHEAR+shift), f, P[6], UX, UY, RX, RY;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		nd->TestGrid("nodes.bln", 0.0007, 0., 0., 0., AXIS_Z, 1);
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[0][0][num], 0, solver.dim[0]*sizeof(double));

///////////////////////////////////////////////////////////////////////////////////////////////
//...������������ ��������� �� ������ ���������� � ������������ ���������� (��� ��������� ���);
			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
 			B[i].shape->parametrization_grad(0, P, 1);
			B[i].shape->parametrization_grad(1, P, 2);

			jump1_common_x(P, i, 0); 
			jump1_common_y(P, i, 1); 

			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
			B[i].shape->parametrization_grad(1, P, 2);

			jump1_compress_x(P, i, 0);
			jump1_compress_y(P, i, 1);
			jump_make_common(i, 0); //...��������� � ����� ������� ���������;

//////////////////////////////////
//...��������� ������ �����������;
			RX = B[i].shape->potential(solver.hh[i][0][m],   0);
			RY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			
			K[3] -= RX*nd->nX[l]*f;
			K[4] -=(RX*nd->nY[l]+RY*nd->nX[l])*f*.5;
			K[5] -= RY*nd->nY[l]*f;

			RX = B[i].shape->potential(solver.hh[i][0][m],   1);
			RY = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			
			K[9]  -= RX*nd->nX[l]*f;
			K[10] -=(RX*nd->nY[l]+RY*nd->nX[l])*f*.5;
			K[11] -= RY*nd->nY[l]*f;

////////////////////////////////////////
//...��������� ������������ �����������;
			jump_admittance(1, i, 0, 0.);
			jump_admittance(1, i, 1, 0.);

			UX = B[i].shape->potential(solver.hh[i][0][m],   0);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			
			K[0] -= G0*(UX*nd->nX[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
			K[1] -= G0*(UX*nd->nY[l]+UY*nd->nX[l])*f;
			K[2] -= G0*(UY*nd->nY[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;

			UX = B[i].shape->potential(solver.hh[i][0][m],   1);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			
			K[6]  -= G0*(UX*nd->nX[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
			K[7]  -= G0*(UX*nd->nY[l]+UY*nd->nX[l])*f;
			K[8]  -= G0*(UY*nd->nY[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CCohes2D::set_fasa_hmg(double nju1, double nju2, double G1, double G2, double C1, double C2)
{
	if (size_of_param() > NUM_SHEAR+3+NUM_SHIFT) {
		param[NUM_SHEAR-1] = C1;
		param[NUM_SHEAR-1+NUM_SHIFT] = C2;
		param[NUM_SHEAR] = G1;
		param[NUM_SHEAR+NUM_SHIFT] = G2;
		param[NUM_SHEAR+1] = nju1;
		param[NUM_SHEAR+1+NUM_SHIFT] = nju2;
		param[NUM_SHEAR+2] = sqrt((.5-nju1)/(1.-nju1)*C1/G1);
		param[NUM_SHEAR+2+NUM_SHIFT] = sqrt((.5-nju2)/(1.-nju2)*C2/G2);
		param[NUM_SHEAR+3] = sqrt(C1/G1);
		param[NUM_SHEAR+3+NUM_SHIFT] = sqrt(C2/G2);
	}
}

//////////////////////////////////////////////////////
//...counting kernel for solving double plane problem;
Num_State CCohes2D::computing_header(Num_Comput Num)
{
	int i, j, k, elem, id_dir, n_rhs = 2;
	char msg[201];

	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");
		sprintf(msg, "CCohes2D sample: N_sm = %d, N_mpl = %d, C1 = %g, C2 = %g", N,
				  UnPackInts(get_param(NUM_MPLS)), get_param(NUM_SHEAR-1), get_param(NUM_SHEAR-1+NUM_SHIFT));
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num){
			case   BASIC_COMPUT: Message("Analytical Blocks..."); break;
			case MAPPING_COMPUT: Message("FEM Blocks...");			break;
		}
		Message(" ");
	}

//////////////////////////////////////////////
//...������������� ���������� ���������� ����;
	if (NUM_PHASE > 0)
	for (k = 0; k < N; k++) if (B[k].link[0] >= NUM_PHASE)
	for (j = B[k].link[0]; j > 0; j--) 
	if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE])
	if (/*B[k].link[NUM_PHASE] == -1 &&*/ B[k].mp) { 
		int m, m1, m2;
		m  = stru.geom[(i = stru.geom_ptr[k])+1]-2;
		m1 = stru.geom[i+4+(j-1)%m];
		m2 = stru.geom[i+4+j%m];
		double X0 = stru.X[m1], Y0 = stru.Y[m1], 
				 X1 = stru.X[m2], Y1 = stru.Y[m2];
		B[k].type = ELLI_BLOCK;
		B[k].mp[4] = arg0(comp(X1-X0, Y1-Y0))-M_PI_2;
		B[k].mp[5] = 0.;
		B[k].mp[6] = 0.;
	}

//////////////////////////////////
//...���������� ������� ���������;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials!!!
	solver.n += 10;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);

////////////////////////////////////////////////////////////////////
//...��������� ����� �� ������������� ������ � �� ������� ���������;
	if (solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; SetGeomBounding(par); d_parm[0] = par[1]-par[0]; d_parm[1] = par[3]-par[2];
		for (k = 0; k < N; k++) SkeletonBounding(B[k], par);
		for (k = 0; k < N; k++) if (B[k].link && B[k].link[NUM_PHASE] == -1) {
			i = 0; while ((elem =  geom_plink_2D(B[k], i, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
			i = j = 0; while ((elem = block_plink_2D(B[k], i, j, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
		}
	}
	LinkPhase2D(MAX_PHASE);

/////////////////////////////////////////////////////////
//...������ ������������� ��������� � ������ �����������;
	if (! solver.struct_permutat(solver.id_change == EXTERN_STATE ? NULL_STATE : /*NULL*/OK_STATE) || 
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

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
void CCohes2D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp && B[i].link[0] >= NUM_PHASE) {
		int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
			case DISPL_VALUE: {
//////////////////////////
//...common displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_compress_x(P, i, 0);
				jump1_compress_y(P, i, 1);
//				jump_make_common(i, 0); //...������ B[i].shape->norm_common(F);

/////////////////////////////////////
//...calculation common displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				//if ((solv%ENERGY_SOLVING) && 1) F[(solv%ENERGY_SOLVING)-1] -=X;
				F[1] = F[0]-X;
			}  break;
			case DISPL_CLASSIC_VALUE: {
//////////////////////////////////
//...displacement (classic field);
				B[i].shape->parametrization_grad(0, P, 1);

				jump1_common_x(P, i, 0); jump_admittance(1, i, 0, 0.);
				jump1_common_y(P, i, 1); jump_admittance(1, i, 1, 0.);

////////////////////////////////////////
//...calculation classical displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				if ((solv%ENERGY_SOLVING) && 0) F[(solv%ENERGY_SOLVING)-1] -=X;
			}  break;
			case DISPL_COHESION_VALUE: {
///////////////////////////////////
//...displacement (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_common_x (P, i, 0); jump_admittance(0, i, 0, 0.);
				jump1_common_y (P, i, 1); jump_admittance(0, i, 1, 0.);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_compress_x(P, i, 0);
				jump1_compress_y(P, i, 1);

///////////////////////////////////////
//...calculation cohesion displacement;
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_X_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0);
				jump4_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_Y_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0);
				jump4_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_X_CLASSIC_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
//////////////////////////////////////
//...surface stresses (classic field);
				B[i].shape->parametrization_hess(0, P, 1);

				jump4_common_x(P, i, 0); jump_admittance(1, i, 0, 0.);
				jump4_common_y(P, i, 1); jump_admittance(1, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_Y_CLASSIC_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
//////////////////////////////////////
//...surface stresses (classic field);
				B[i].shape->parametrization_hess(0, P, 1);

				jump4_common_x(P, i, 0); 
				jump4_common_y(P, i, 1); 

/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_X_COHESION_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////////////////////////
//...surface stresses (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0); jump_admittance(0, i, 0, 0.);
				jump4_compress_y(P, i, 1); jump_admittance(0, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_Y_COHESION_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////////////////////////
//...surface stresses (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0); jump_admittance(0, i, 0, 0.);
				jump4_compress_y(P, i, 1); jump_admittance(0, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case NORMAL_X_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case NORMAL_Y_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case DILAT_VALUE: {//...divergency of displacements;
				double G = (.5-get_param(NUM_SHEAR+1+shift))/((1.-get_param(NUM_SHEAR+1+shift))*get_param(NUM_SHEAR+shift));
///////////////////////////////////////
//...classical and pressure potentials;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
 				B[i].shape->parametrization_grad(0, P, 1);
				B[i].shape->parametrization_grad(1, P, 2);

				B[i].shape->adm_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 0), 1.); 
				B[i].shape->adm_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), 1.); 
				B[i].shape->adm_x(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 0), -1.); 
				B[i].shape->adm_y(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -1.); 

/////////////////////////////////
//...divergency of displacements;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant)*G;
			}	break;
			case ROTATION_VALUE: {//...rotation of displacements;
				double G = .5/get_param(NUM_SHEAR+shift);
////////////////////////////////////
//...classical and shear potentials;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
 				B[i].shape->parametrization_grad(0, P, 1);
				B[i].shape->parametrization_grad(1, P, 2);

				B[i].shape->adm_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 0), -1.); 
				B[i].shape->adm_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), 1.); 
				B[i].shape->adm_y(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 0), 1.); 
				B[i].shape->adm_x(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -1.); 

///////////////////////////////
//...rotation of displacements;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant)*G;
			}	break;
			case ENERGY_VALUE: {//...local energy in cell;
				double e11, e12, e22, theta, u1, u2, G1 = get_param(NUM_SHEAR+shift), 
						 G2 = get_param(NUM_SHEAR+1+shift)*G1/(.5-get_param(NUM_SHEAR+1+shift)), C1 = get_param(NUM_SHEAR-1+shift);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				jump1_common_x (P, i, 4); jump_admittance(0, i, 4, 0.);
				jump1_common_y (P, i, 5); jump_admittance(0, i, 5, 0.);

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//...��������� ���������, ���������� ������� ���������� � ����������� ���� (� ������ �������� ���������);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				e11 =	F[0];
				e12 =	F[1];

				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common(F);

				e22 = F[1];
				e12 =(F[0]+e12)*.5;
				theta = e11+e22;

				u1 = B[i].shape->potential(solver.hh[i][0][m+4], id_variant);
				u2 = B[i].shape->potential(solver.hh[i][0][m+5], id_variant);

//////////////////
//...local energy;
				F[0] = G1*(e11*e11+2.*e12*e12+e22*e22)+(G2*theta*theta+C1*(u1*u1+u2*u2))*.5;
			}	break;
			case STRAIN_X_VALUE: {
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////////
//...calculation common dilatation tensor (exx and exy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common(F); double exy = F[0];

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				F[1] = (F[1]+exy)*.5;
			}  break;
			case STRAIN_Y_VALUE: {
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////////
//...calculation common dilatation tensor (exx and exy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F); double exy = F[1];

				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common(F);

				F[0] = (F[0]+exy)*.5;
			}  break;
			default : F[0] = i; F[1] = 0.;
     }
  }
}

///////////////////////////////////////////
//...calculaion energy by block functional;
void CCohes2D::GetEnergyValue(int k, double * energy)
{
   if (solver.dim && solver.dim[k] && solver.JR && solver.JR[k] && solver.hh && solver.hh[k][0].GetMatrix() && 
							solver.TL && solver.TL[k]) {
		double ** TL = solver.TL[k][solver.JR[k][solver.JR_DIAG]-solver.JR_SHIFT].GetMatrix(), * h = solver.hh[k][0][solver.id_norm];
		int i, j, NN = solver.dim[k];

///////////////////////
//...��������� �������;
		block_shape_init(B[k], NO_STATE);
		if (h && TL) { 
			double measure = .5/get_param(NUM_SHEAR);
			for (i = 0; i < NN; i++)
			for (j = 0; j < NN; j++)
				energy[-B[k].link[NUM_PHASE]-1] += TL[i][j]*h[i]*h[j]*measure;
		}
	}
}

/////////////////////////////////////////////////////////////////////
//...���������� ������ ���������� ���� (������������ �������� �����);
void CCohes2D::TakeLayerModel(double L, double H, double l, double nju1, double nju2, double G1, double G2, double l1, double l2, double A)
{
	double P1[3] = {-L, -H*.5, 0.},
			 P2[3] = {-l, -H*.5, 0.},
			 P3[3] = {-l,  H*.5, 0.},
			 P4[3] = {-L,  H*.5, 0.},
			 P5[3] = {-l, -H*.5, 0.},
			 P6[3] = { l, -H*.5, 0.},
			 P7[3] = { l,  H*.5, 0.},
			 P8[3] = {-l,  H*.5, 0.};

/////////////////////////////////////
//...�������� ������� �� ���� ������;
	CCells * ce; 
   init_blocks(3);

	B[0].type = ELLI_BLOCK;
	B[0].bar	 = new CCells;
	ce = new CCells;
	ce->get_quad_facet(P1, P2, P3, P4, ZERO_STATE);
	B[0].bar->bar_add(ce);
	set_block3D(B[0], SPHERE_GENUS, L-l, 0, -L);
	set_link (B[0], NUM_PHASE+1);
	B[0].link[2] = -1+SRF_STATE;
	B[0].link[NUM_PHASE+1] = 1;
	B[0].mp[7] = 0.; //...��������� ����� �� ���� ������;

	B[1].type = ELLI_BLOCK;
	B[1].bar	 = new CCells;
	ce = new CCells;
	ce->get_quad_facet(P5, P6, P7, P8, ZERO_STATE);
	B[1].bar->bar_add(ce);
	set_block3D(B[1], SPHERE_GENUS, l);
	set_link (B[1], NUM_PHASE+2);
	B[1].link[2] = -2+SRF_STATE;
	B[1].link[4] = SRF_STATE;
	B[1].link[NUM_PHASE]	 = -2;
	B[1].link[NUM_PHASE+1] = 2;
	B[1].link[NUM_PHASE+2] = 0;
	B[1].mp[7] = 0.; //...��������� ����� � ������ ������;

	B[2].type = ELLI_BLOCK;
	B[2].bar	 = new CCells;
	ce = new CCells; P1[0] = -P1[0]; P2[0] = -P2[0]; P3[0] = -P3[0]; P4[0] = -P4[0];
	ce->get_quad_facet(P2, P1, P4, P3, ZERO_STATE);
	B[2].bar->bar_add(ce);
	set_block3D(B[2], SPHERE_GENUS, L-l, 0, L);
	set_link (B[2], NUM_PHASE+1);
	B[2].link[4] = -1+SRF_STATE;
	B[2].link[NUM_PHASE+1] = 1;
	B[2].mp[7] = 0.; //...��������� ����� �� ���� ������;

////////////////////////////////////
//...������������� ��������� ������;
	double C1, C2, k1 = (1.-nju1)/(.5-nju1)*G1, k2 = (1.-nju2)/(.5-nju2)*G2;
	set_mpls(PackInts(1, 1)); //...multipoles degree;
	set_quad(PackInts(4, 2)); //...quadrature degree;
	set_normaliz(1.);
	set_param(NUM_ADHES, A);
	set_lagrange(1e5);		  //...Lagrange corfficient for LSM;
	set_fasa_hmg(nju1, nju2, G1, G2, C1 = l1 ? G1/sqr(l1) : 1e33, C2 = l2 ? G2/sqr(l2) : 1e33);

/////////////////////////////////////////////////
//...�������������� ���������� � ������� �������;
	solver.set_blocks(N, 1); //<==== number of saved potentials !!!
	solver.n += 10;//<==== number of additional auxilliary arrays!!!
   shapes_init(INITIAL_STATE);
	B[0].shape->init1(1, 0, solver.id_norm, 2);//...� ���� ������ ������� ��������;
	B[1].shape->init1(1, 0, solver.id_norm, 2);
	B[2].shape->init1(1, 0, solver.id_norm, 2);
	shapes_init(NULL_STATE);
	int k;
	for (k = 0; k < solver.N;  k++)
	solver.set_dimension(k, freedom_block(k));
   solver.struct_init();

   double par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	for (k = 0; k < N; k++)
		SkeletonBounding(B[k], par);
	d_parm[0] = par[1] - par[0];
	d_parm[1] = par[3] - par[2];

////////////////////////////////////////////////
//...������ ������������ �������������� �������;
	double kk1 = sqrt(C1/k1), t1 = exp(-kk1*(L-l)), tt1 = sqr(t1), al1 = kk1*(1.+tt1)/(1.-tt1),
			 kk2 = sqrt(C2/k2), t2 = exp(-kk2*l),		tt2 = sqr(t2), al2 = kk2*(1.+tt2)/(1.-tt2),
			 AA = k2*al1+k1*al2-A*al1*al2,
			 BB = k2-k1-A*al2,
			 CC = k2-k1+A*al1,
			 QQ = AA*(k2*(L-l)+k1*l+A)-BB*CC,
			 HM = k2*L*AA/QQ,
			 HD = k1*L*AA/QQ,
			 DM = k2*L*BB/QQ*2.*t1/(1.-tt1),
			 DD = k1*L*CC/QQ*2.*t2/(tt2-1.);

////////////////////////////////////////////////////////////
//...������� ������������ � ������������� ���������-�������;
	B[0].shape->set_R(1.);
	B[0].shape->A[0][0] = -L*4.*G1*(1.-nju1)/(3.-4.*nju1);
	B[0].shape->A[0][2] = HM*2.*G1*(1.-nju1)/(1.-2.*nju1);
	B[0].shape->A[0][7] = DM*C1/kk1;
	B[1].shape->set_R(1.);
	B[1].shape->A[0][2] = HD*2.*G2*(1.-nju2)/(1.-2.*nju2);
	B[1].shape->A[0][7] = DD*C2/kk2;
	B[2].shape->set_R(1.);
	B[2].shape->A[0][0] =  L*4.*G1*(1.-nju1)/(3.-4.*nju1);
	B[2].shape->A[0][2] = HM*2.*G1*(1.-nju1)/(1.-2.*nju1);
	B[2].shape->A[0][7] = DM*C1/kk1;
}

//////////////////////////////////////////////////////////////////////////
//...����������� ������ ���� �� ������ ������ ������������ �������� �����;
double CCohes2D::TakeLayer_E1(double ff)
{
	double nj1 = get_param(NUM_SHEAR+1), 
			 nj2 = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT),
			 k1 = (1.-nj1)/(.5-nj1)*G1,
			 k2 = (1.-nj2)/(.5-nj2)*G2,
			 C1 = get_param(NUM_SHEAR-1),
			 C2 = get_param(NUM_SHEAR-1+NUM_SHIFT),
			 A  = get_param(NUM_ADHES),
			 B  = get_param(NUM_ADHES+1);
////////////////////////////////////////////////////////////
//...����� ������������ ������ ���������� ����������/������;
	double kk1 = sqrt(C1/k1), t1 = ff ? exp(-kk1*(1./ff-1.)) : 0., tt1 = sqr(t1), al1 = kk1*(1.+tt1)/(1.-tt1),
			 kk2 = sqrt(C2/k2), t2 = exp(-kk2),	tt2 = sqr(t2), al2 = kk2*(1.+tt2)/(1.-tt2),
			 AA = k2*al1+k1*al2-A*al1*al2,
			 BB = k2-k1-A*al2,
			 CC = k2-k1+A*al1,
			 QQ = AA*(k2*(1.-ff)+k1*ff+A*ff)-BB*CC*ff,
			 kk_eff = k1*k2*AA/QQ, GG_eff;
////////////////////////////////////////////////////
//...��������� ����� ��� ������������ ������ ������;
			 kk1 = sqrt(C1/G1); t1 = ff ? exp(-kk1*(1./ff-1.)) : 0.; tt1 = sqr(t1); al1 = kk1*(1.+tt1)/(1.-tt1);
			 kk2 = sqrt(C2/G2); t2 = exp(-kk2);	tt2 = sqr(t2); al2 = kk2*(1.+tt2)/(1.-tt2);
			 AA = G2*al1+G1*al2-B*al1*al2;
			 BB = G2-G1-B*al2;
			 CC = G2-G1+B*al1;
			 QQ = AA*(G2*(1.-ff)+G1*ff+B*ff)-BB*CC*ff;
			 GG_eff = G1*G2*AA/QQ;
	return(3.*kk_eff-4.*GG_eff)*GG_eff/(kk_eff-GG_eff);
}


////////////////////////////////////////////////////////////////////////////
//...����������� ������ ������ �� ������ ������ ������������ �������� �����;
double CCohes2D::TakeLayer_G1(double ff)
{
	double nj1 = get_param(NUM_SHEAR+1), 
			 nj2 = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT),
			 C1 = get_param(NUM_SHEAR-1),
			 C2 = get_param(NUM_SHEAR-1+NUM_SHIFT),
			 A  = get_param(NUM_ADHES),
			 B  = get_param(NUM_ADHES+1);
////////////////////////////////////////////////////
//...��������� ����� ��� ������������ ������ ������;
	double kk1 = sqrt(C1/G1), t1 = ff ? exp(-kk1*(1./ff-1.)) : 0., tt1 = sqr(t1), al1 = kk1*(1.+tt1)/(1.-tt1),
			 kk2 = sqrt(C2/G2), t2 = exp(-kk2),	tt2 = sqr(t2), al2 = kk2*(1.+tt2)/(1.-tt2),
			 AA = G2*al1+G1*al2-B*al1*al2,
			 BB = G2-G1-B*al2,
			 CC = G2-G1+B*al1,
			 QQ = AA*(G2*(1.-ff)+G1*ff+B*ff)-BB*CC*ff,
			 GG_eff = G1*G2*AA/QQ;
	return(GG_eff);
}

////////////////////////////////////////////////////////////////////
//...���������� ������ ��� ������������ ��������� (������ ��������);
double CCohes2D::TakeEshelby_volm_two(double ff)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
			 C_I = get_param(NUM_SHEAR-1+NUM_SHIFT), 
			 C_M = get_param(NUM_SHEAR-1), AA = get_param(NUM_ADHES),
			 K_I = mu_I/(1.-2.*nu_I), K_M = mu_M/(1.-2.*nu_M);
	if (C_I == 0. || C_M == 0.) { //...������������ ���������� ������;
		return(K_M+ff*(K_I-K_M)/(1.+(1.-ff)*(K_I-K_M)/(K_M+mu_M)));
	}
	double ku_I = (1.-nu_I)/(.5-nu_I)*mu_I,
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, RR2 = 1./pow(ff, 1./3.), 
			 kk_I = sqrt(C_I/ku_I), tt_I = exp(-2.*kk_I),
			 kk_M = sqrt(C_M/ku_M), tt_D = exp(kk_M*(RR2-1.)),
			 HH1  = ((1.+tt_I)*kk_I-(1.-tt_I))*.5, //...���������� ������� ���� ������������ !!!
			 HH2  = ((1.-tt_I)*(sqr(kk_I)+3.)-3.*(1.+tt_I)*kk_I)*.5,
			 JJP1 =  (kk_M-1.), JJP2 = ((sqr(kk_M)+3.)-3.*kk_M), JHP1 =  (kk_M*RR2-1.)*ff*tt_D, 
			 JJM1 = -(kk_M+1.), JJM2 = ((sqr(kk_M)+3.)+3.*kk_M), JHM1 = -(kk_M*RR2+1.)*ff/tt_D,
			 matr[7][8] = {
					{ 1., -HH1, -1., -1., JJP1, JJM1, 0., 0.},
					{ 0., -HH2,  0.,  2., JJP2, JJM2, 0., 0.},
					{ AA, -ku_I*HH1-AA*(HH2+HH1), 0.,0., ku_M*JJP1, ku_M*JJM1, 0., 0.},
					{ K_I, (mu_I+3.*ku_I)/4.*HH1, -K_M, mu_M, -(mu_M+3.*ku_M)/4.*JJP1, -(mu_M+3.*ku_M)/4.*JJM1, 0., 0.},
					{ 0., 0., 0., 0., JHP1, JHM1,  0., 0.},
					{ 0., 0., K_M, -mu_M*ff, 0., 0., -1., 0.}, //...����������� ������ � ������������� ����;
					{ 0., 0., 1., ff, 0., 0., 0., 1.},			 //...����������� ���� ������� � ������ �����;
	};

//////////////////////////////////////////////////////////////////
//...������ ������� �������� ��������� A0, C0, A1, B1, C1, D1, KH;
	int dim_N = 7, ii[7] = {0, 0, 0, 0, 0, 0, 0}, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1) 
				for (l = 0; l < dim_N; l++) 
					if (! ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) return(0.);
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
			}
		if (matr[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
		double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
		for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matr[k][l0]; matr[k][l0] = 0.;
				for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
			}
	}
	return(matr[6][7]);
}
#undef  Message
