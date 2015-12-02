#include "stdafx.h"
#include "cgrid_qg.h"

////////////////////////////////////////////////////////////////////
//...������� ��������� ������ ���� ��������� �������� (���� � ����);
	const double CGrid_QG::cc[][2] = {
//	static const double c1[][2] = {
	{ 0.00000000000000e+000, 2.00000000000000e+000},
//	};
//	static const double c2[][2] = {
	{ 5.773502691896258e-001, 1.000000000000001e+000},
//	};
//	static const double c3[][2] = {
	{ 0.000000000000000e+000, 8.888888888888888e-001},
	{ 7.745966692414834e-001, 5.555555555555559e-001},
//	};
//	static const double c4[][2] = {
	{ 3.399810435848561e-001, 6.521451548625462e-001},
	{ 8.611363115940526e-001, 3.478548451374540e-001},
//	};
//	static const double c5[][2] = {
	{ 0.000000000000000e+000, 5.688888888888889e-001},
	{ 5.384693101056830e-001, 4.786286704993661e-001},
	{ 9.061798459386640e-001, 2.369268850561891e-001},
//	};
//	static const double c6[][2] = {
	{ 2.386191860831970e-001, 4.679139345726912e-001},
	{ 6.612093864662646e-001, 3.607615730481392e-001},
	{ 9.324695142031521e-001, 1.713244923791702e-001},
//	};
//	static const double c7[][2] = {
	{ 0.000000000000000e+000, 4.179591836734694e-001},
	{ 4.058451513773971e-001, 3.818300505051187e-001},
	{ 7.415311855993946e-001, 2.797053914892780e-001},
	{ 9.491079123427585e-001, 1.294849661688697e-001},
//	};
//	static const double c8[][2] = {
	{ 1.834346424956496e-001, 3.626837833783619e-001},
	{ 5.255324099163290e-001, 3.137066458778867e-001},
	{ 7.966664774136267e-001, 2.223810344533745e-001},
	{ 9.602898564975362e-001, 1.012285362903764e-001},
//	};
//	static const double c9[][2] = {
	{ 0.000000000000000e+000, 3.302393550012597e-001},
	{ 3.242534234038090e-001, 3.123470770400028e-001},
	{ 6.133714327005903e-001, 2.606106964029345e-001},
	{ 8.360311073266359e-001, 1.806481606948579e-001},
	{ 9.681602395076261e-001, 8.127438836157423e-002},
//	};
//	static const double c10[][2] = {
	{ 1.488743389816312e-001, 2.955242247147529e-001},
	{ 4.333953941292472e-001, 2.692667193099964e-001},
	{ 6.794095682990244e-001, 2.190863625159825e-001},
	{ 8.650633666889844e-001, 1.494513491505802e-001},
	{ 9.739065285171717e-001, 6.667134430868832e-002},
//	};
//	static const double c11[][2] = {
	{ 0.000000000000000e+000, 2.729250867779006e-001},
	{ 2.695431559523449e-001, 2.628045445102468e-001},
	{ 5.190961292068119e-001, 2.331937645919907e-001},
	{ 7.301520055740493e-001, 1.862902109277333e-001},
	{ 8.870625997680952e-001, 1.255803694649041e-001},
	{ 9.782286581460570e-001, 5.566856711617373e-002},
//	};
//	static const double c12[][2] = {
	{ 1.252334085114689e-001, 2.491470458134028e-001},
	{ 3.678314989981802e-001, 2.334925365383551e-001},
	{ 5.873179542866175e-001, 2.031674267230667e-001},
	{ 7.699026741943046e-001, 1.600783285433466e-001},
	{ 9.041172563704750e-001, 1.069393259953175e-001},
	{ 9.815606342467192e-001, 4.717533638651172e-002},
//	};
//	static const double c13[][2] = {
	{ 0.000000000000000e+000, 2.325515532308738e-001},
	{ 2.304583159551347e-001, 2.262831802628970e-001},
	{ 4.484927510364469e-001, 2.078160475368887e-001},
	{ 6.423493394403402e-001, 1.781459807619457e-001},
	{ 8.015780907333100e-001, 1.388735102197873e-001},
	{ 9.175983992229780e-001, 9.212149983772761e-002},
	{ 9.841830547185881e-001, 4.048400476531598e-002},
//	};
//	static const double c14[][2] = {
	{ 1.080549487073436e-001, 2.152638534631577e-001},
	{ 3.191123689278897e-001, 2.051984637212961e-001},
	{ 5.152486363581539e-001, 1.855383974779373e-001},
	{ 6.872929048116855e-001, 1.572031671581938e-001},
	{ 8.272013150697649e-001, 1.215185706879030e-001},
	{ 9.284348836635734e-001, 8.015808715975972e-002},
	{ 9.862838086968123e-001, 3.511946033175199e-002},
//	};
//	static const double c15[][2] = {
	{ 0.000000000000000e+000, 2.025782419255612e-001},
	{ 2.011940939974345e-001, 1.984314853271115e-001},
	{ 3.941513470775635e-001, 1.861610000155630e-001},
	{ 5.709721726085389e-001, 1.662692058169943e-001},
	{ 7.244177313601700e-001, 1.395706779261542e-001},
	{ 8.482065834104273e-001, 1.071592204671729e-001},
	{ 9.372733924007060e-001, 7.036604748810772e-002},
	{ 9.879925180204854e-001, 3.075324199611722e-002},
//	};
//	static const double c16[][2] = {
	{ 0.095012509837637440185, 0.189450610455068496285},
	{ 0.281603550779258913230, 0.182603415044923588867},
	{ 0.458016777657227386342, 0.169156519395002538189},
	{ 0.617876244402643748447, 0.149595988816576732081},
	{ 0.755404408355003033895, 0.124628971255533872052},
	{ 0.865631202387831743880, 0.095158511682492784810},
	{ 0.944575023073232576078, 0.062253523938647892863},
	{ 0.989400934991649932596, 0.027152459411754094852},
//	};
//	static const double c17[][2] = {
	{ 0.000000000000000e+000, 1.794464703562065e-001},
	{ 1.784841814958477e-001, 1.765627053669923e-001},
	{ 3.512317634538764e-001, 1.680041021564505e-001},
	{ 5.126905370864771e-001, 1.540457610768106e-001},
	{ 6.576711592166908e-001, 1.351363684685258e-001},
	{ 7.815140038968015e-001, 1.118838471934044e-001},
	{ 8.802391537269859e-001, 8.503614831717912e-002},
	{ 9.506755217687678e-001, 5.545952937398826e-002},
	{ 9.905754753144174e-001, 2.414830286854794e-002},
//	};
//	static const double c18[][2] = {
	{ 8.477501304173554e-002, 1.691423829631438e-001},
	{ 2.518862256915054e-001, 1.642764837458324e-001},
	{ 4.117511614628425e-001, 1.546846751262647e-001},
	{ 5.597708310739475e-001, 1.406429146706507e-001},
	{ 6.916870430603533e-001, 1.225552067114790e-001},
	{ 8.037049589725230e-001, 1.009420441062867e-001},
	{ 8.926024664975558e-001, 7.642573025488873e-002},
	{ 9.558239495713978e-001, 4.971454889497080e-002},
	{ 9.915651684209309e-001, 2.161601352648298e-002},
//	};
//	static const double c19[][2] = {
	{ 0.000000000000000e+000, 1.610544498487836e-001},
	{ 1.603586456402255e-001, 1.589688433939545e-001},
	{ 3.165640999636298e-001, 1.527660420658598e-001},
	{ 4.645707413759609e-001, 1.426067021736067e-001},
	{ 6.005453046616810e-001, 1.287539625393362e-001},
	{ 7.209661773352293e-001, 1.115666455473338e-001},
	{ 8.227146565371428e-001, 9.149002162245011e-002},
	{ 9.031559036148180e-001, 6.904454273764164e-002},
	{ 9.602081521348301e-001, 4.481422676569679e-002},
	{ 9.924068438435844e-001, 1.946178822972639e-002},
//	};
//	static const double c20[][2] = {
	{ 7.652652113349748e-002, 1.527533871307262e-001},
	{ 2.277858511416451e-001, 1.491729864726038e-001},
	{ 3.737060887154197e-001, 1.420961093183824e-001},
	{ 5.108670019508270e-001, 1.316886384491765e-001},
	{ 6.360536807265150e-001, 1.181945319615192e-001},
	{ 7.463319064601509e-001, 1.019301198172409e-001},
	{ 8.391169718222188e-001, 8.327674157670445e-002},
	{ 9.122344282513258e-001, 6.267204833410871e-002},
	{ 9.639719272779138e-001, 4.060142980038684e-002},
	{ 9.931285991850949e-001, 1.761400713915210e-002},
//	};
//	static const double c25[][2] = {
	{ 0.000000000000000e+000, 1.231760537267154e-001},
	{ 1.228646926107105e-001, 1.222424429903102e-001},
	{ 2.438668837209885e-001, 1.194557635357845e-001},
	{ 3.611723058093877e-001, 1.148582591457115e-001},
	{ 4.730027314457151e-001, 1.085196244742642e-001},
	{ 5.776629302412230e-001, 1.005359490670508e-001},
	{ 6.735663684734683e-001, 9.102826198296378e-002},
	{ 7.592592630373576e-001, 8.014070033500104e-002},
	{ 8.334426287608340e-001, 6.803833381235715e-002},
	{ 8.949919978782753e-001, 5.490469597583510e-002},
	{ 9.429745712289757e-001, 4.093915670130033e-002},
	{ 9.766639214595175e-001, 2.635498661503233e-002},
	{ 9.955569697904981e-001, 1.139379850102624e-002},
//	};
//	static const double c30[][2] = {
	{ 5.147184255531776e-002, 1.028526528935589e-001},
	{ 1.538699136085834e-001, 1.017623897484052e-001},
	{ 2.546369261678901e-001, 9.959342058679528e-002},
	{ 3.527047255308782e-001, 9.636873717464441e-002},
	{ 4.470337695380894e-001, 9.212252223778647e-002},
	{ 5.366241481420201e-001, 8.689978720108352e-002},
	{ 6.205261829892429e-001, 8.075589522942016e-002},
	{ 6.978504947933157e-001, 7.375597473770475e-002},
	{ 7.677774321048262e-001, 6.597422988218089e-002},
	{ 8.295657623827684e-001, 5.749315621761909e-002},
	{ 8.825605357920526e-001, 4.840267283059395e-002},
	{ 9.262000474292743e-001, 3.879919256962740e-002},
	{ 9.600218649683067e-001, 2.878470788333767e-002},
	{ 9.836681232797472e-001, 1.846646831108939e-002},
	{ 9.968934840746495e-001, 7.968192496166537e-003},
//	};
//	static const double c35[][2] = {
	{ 0.000000000000000e+000, 8.848679490710421e-002},
	{ 8.837134327565935e-002, 8.814053043027560e-002},
	{ 1.760510611659896e-001, 8.710444699718356e-002},
	{ 2.623529412092962e-001, 8.538665339209925e-002},
	{ 3.466015544308140e-001, 8.300059372885678e-002},
	{ 4.281375415178144e-001, 7.996494224232471e-002},
	{ 5.063227732414888e-001, 7.630345715544216e-002},
	{ 5.805453447497645e-001, 7.204479477255971e-002},
	{ 6.502243646658904e-001, 6.722228526908730e-002},
	{ 7.148145015566287e-001, 6.187367196607984e-002},
	{ 7.738102522869126e-001, 5.604081621237028e-002},
	{ 8.267498990922254e-001, 4.976937040135290e-002},
	{ 8.732191250252224e-001, 4.310842232617073e-002},
	{ 9.128542613593176e-001, 3.611011586346319e-002},
	{ 9.453451482078273e-001, 2.882926010889387e-002},
	{ 9.704376160392302e-001, 2.132297991149260e-002},
	{ 9.879357644438515e-001, 1.365082834836317e-002},
	{ 9.977065690996003e-001, 5.883433420443020e-003},
//	};
//	static const double c40[][2] = {
	{ 3.877241750605096e-002, 7.750594797842489e-002},
	{ 1.160840706752553e-001, 7.703981816424799e-002},
	{ 1.926975807013711e-001, 7.611036190062628e-002},
	{ 2.681521850072537e-001, 7.472316905796825e-002},
	{ 3.419940908257585e-001, 7.288658239580376e-002},
	{ 4.137792043716050e-001, 7.061164739128703e-002},
	{ 4.830758016861787e-001, 6.791204581523402e-002},
	{ 5.494671250951282e-001, 6.480401345660065e-002},
	{ 6.125538896679803e-001, 6.130624249292896e-002},
	{ 6.719566846141796e-001, 5.743976909939130e-002},
	{ 7.273182551899270e-001, 5.322784698393639e-002},
	{ 7.783056514265193e-001, 4.869580763507200e-002},
	{ 8.246122308333117e-001, 4.387090818567339e-002},
	{ 8.659595032122595e-001, 3.878216797447185e-002},
	{ 9.020988069688742e-001, 3.346019528254779e-002},
	{ 9.328128082786765e-001, 2.793700698002306e-002},
	{ 9.579168192137917e-001, 2.224584919416719e-002},
	{ 9.772599499837734e-001, 1.642105838195386e-002},
	{ 9.907262386994571e-001, 1.049828453115235e-002},
	{ 9.982377097105593e-001, 4.521277098533210e-003},
	};

/////////////////////////////////////////////
//...�������, ������������ ���������� ������;
void QG_generate(int N, long double cc[][2])
{
	long double * h  = (long double *)new_struct(1500*sizeof(long double)), 
					* an = (long double *)new_struct(30*sizeof(long double)),
					* bn = (long double *)new_struct(30*sizeof(long double)),
					* dn = (long double *)new_struct(30*sizeof(long double));

	int N1 = (N+1)/2, N2, N3, i, j, k, sg, 
		 nn, m, m1, m2, m3, m5, ms, ns;

	long double a, c1n, cn, cn1, d, dk, dks, fk0, fks,
					p, p1, p2, y2, t, s, s1, s2, r, r0, r1, r2, r3, rs, rc, r3c,
					rn, r1n, x0, x1, x2, eps, eps1 = 1e-15;

	if (N < 1 || ! cc) return;
   if (N < 2) {
			cc[0][0] = 0.;
			cc[0][1] = 2.;

			return;
	}
   if (N < 5)  j = 0; else
   if (N < 12) j = 1; else
   if (N < 21) j = 2; else j = 3;

	if (j) {
		c1n = sqr(cn = N+.5);

		for (i = 0; i < j; i++) {
			eps = eps1*c1n/(x0 = M_PI*(i+.75));
			ms  = 1;
lab18:
			bn[m1 = 0] = 1.;
			bn[++m1] = rn = sqr(r1n = x0*.5);

			while (bn[m1] > eps && m1 < N) {
				m2 = m1++;
				bn[m1] = rn*bn[m2]/sqr(m2+1.);
			}

			dn[m3 = 0] = r1n;
			do {
				m2 = m3++;
				dn[m3] = rn*dn[m2]/(m3*(m3+1.));
			}
			while (dn[m3] > eps && m3 < N-1);

			m5 = max(m1, m3);
			an[0] = 1.-sqr(.5/cn);

			for (k = 1; k <= m5; k++) 
				an[k] = an[k-1]*(1.- sqr((k+.5)/cn));

			for (k = 0; k < m1; k++) 
				bn[k+1] *= an[k];

			for (k = 0; k <= m3; k++) 
				dn[k] *= an[k];

			for (sg = 1, p = bn[0], k = 1; k <= m1; k++) 
				p += bn[k]*(sg = -sg);

			for (sg = 1, p1 = dn[0], k = 1; k <= m3; k++)
				p1 += dn[k]*(sg = -sg);

			x1 = x0+p/p1;
			if (fabsl(x1-x0) > eps) {
				ms++;
				x0 = x1;
				goto lab18;
			}

			cc[i][0] = 1.-sqr(x1)/(c1n*2.);

			x2 = x1*(cn-1.)/cn;
			cn1 = N-.5;

			bn[m1 = 0] = 1.;
			bn[++m1] = rn = sqr(r1n = x2*.5);

			while (bn[m1] > eps && m1 < N) {
				m2 = m1++;
				bn[m1] = rn*bn[m2]/sqr(m2+1.);
			}

			an[0] = 1.-sqr(.5/cn1);
			for (k = 1; k < m1; k++)
				an[k] = an[k-1]*(1.- sqr((k+.5)/cn1));

			for (k = 0; k < m1; k++) 
				bn[k+1] *= an[k];

			for (sg = 1, p2 = bn[0], k = 1; k <= m1; k++) 
				p2 += bn[k]*(sg = -sg);

			y2 = sqr(x1/cn);
			cc[i][1] = 2.*y2*(1.-y2*.25)/sqr(N*p2);
		}
   }

   N2 = N+1;
   N3 = 2*(N2/2);
	nn = N/2;

	if (N != N3) {
		nn = (N-1)/2;

		cc[N1-1][0] = 0.; 
		for (r = 1., i = 0; i < nn; i++) {
			r0 = 2.*(i+1);
			r *= r0/(r0-1.);
		}
      cc[N1-1][1] = 2.*sqr(r/N);
	}

	for (r = 1., i = 0; i < N2; i++) {
		r0 = 2.*(i+1);
		r *= r0/(r0+1.);
	}
   a = 4.*r/M_PI;

	for (k = j; k < nn; k++) {
		eps = eps1/sinl(fks = fk0 = M_PI*(k+.75)/(N+.5));
		ns  = 0;
		dk  = 0.;
lab10:
		ns++;
		h[m = 0] = 1./(8.*(N+1.5)*(rs = sinl(fks)));
		do {
			if (++m >= 1500) goto lab10;

			h[m] = h[m-1]*(m+.5)*(m+.5)/(2.*(m+1.)*(N+m+1.5)*rs);
		}
		while (h[m] > eps && h[m] < h[m-1]);

		r1 = (N+.5)*fks-M_PI*.25;
		rc = cosl(fks);
		t  = rc/rs;

		for (s = s1 = 0., i = 0; i < m; i++) {
			r2 = N+i+1.5;
			r3 = r2*fks-M_PI_2*(i+1.5);
			r3c = cosl(r3);

			s  += h[i]*r3c;
			s1 += h[i]*(r2*sinl(r3)+t*(i+1.)*r3c);
		}
		s  += cosl(r1);
		s1 += (N+.5)*sinl(r1);
		dks = dk+s/s1;
		fks = fk0+dks;

		if (fabsl(dks-dk) > eps) {
			dk = dks;
			goto lab10;
		}

		cc[k][0] = cosl(fks);
		rs       = sinl(fks);

		h[0] = 1.0/(8.0*(N2+1.5)*rs);
		s2   = h[0]*cosl((N2+1.5)*fks-M_PI*.75);

		for (i  = 1; i < m; i++) {
			h[i] = h[i-1]*(i+.5)*(i+.5)/(2.*(i+1.)*(N2+i+1.5)*rs);
			s2  += h[i]*cosl((N2+i+1.5)*fks-M_PI_2*(i+1.5));
		}

		s2 += cosl((N2+.5)*fks-M_PI*.25);
		r = a*N2*s2;
		d = sinl(fks);

		cc[k][1] = 4.*d*sqr(d/r);
	}

	delete [] h;
	delete [] an;
	delete [] bn;
	delete [] dn;

	for (i = 0; i < N1; i++) 
		for (k = i+1; k < N1; k++) 
			if (fabsl(cc[k][0]) < fabsl(cc[i][0])) {
				swap(cc[k][0], cc[i][0]);
				swap(cc[k][1], cc[i][1]);
			}
}

///////////////////////////////////////////////////////////
//...���������� ������ ��p���� 2M ��� ������ (���� � ����);
void CGrid_QG::QG2M_line(double * P, int M, const double c2M[][2], Num_State normal)
{
	double d1 [3] = {.5*(P[3]-P[0]), .5*(P[4]-P[1]), 0.}, 
			 xm1[3] = { d1[0]+P[0], d1[1]+P[1], 0. },
			 nX = -d1[1], 
			 nY =  d1[0], 
			 nZ =  d1[2], pp[1], F;
	int    i;

	F = QG_normal(nX, nY, nZ);
	if (normal == NULL_STATE) QG_normal(nX = d1[0], nY = d1[1], nZ = d1[2]);

	if (N_par < 1) add_params(1-N_par); 

	for (i = 0; i < M; i++) {
		pp[0] = c2M[i][1]*F;
		add_new_point(xm1[0]+c2M[i][0]*d1[0], xm1[1]+c2M[i][0]*d1[1], 0., nX, nY, nZ, pp); 

		pp[0] = c2M[i][1]*F;
		add_new_point(xm1[0]-c2M[i][0]*d1[0], xm1[1]-c2M[i][0]*d1[1], 0., nX, nY, nZ, pp); 
	}
}

///////////////////////////////////////////////////////////////////////
//...���������� ������ ��p���� 2M ��� ����������� ������ (���� � ����);
void CGrid_QG::QG2M_tria(double * P, int M, const double c2M[][2], Num_State normal)
{
	double d [3] = {.5*(P[6]-P[3]), .5*(P[7]-P[4]), .5*(P[8]-P[5])}, 
			 xm[3] = { d[0]+P[3], d[1]+P[4], d[2]+P[5] },
			 dm[3], d_i[3], xm_i[3], dm_i[3], pp[1], 
			 nX = (P[4]-P[1])*(P[8]-P[2])-(P[5]-P[2])*(P[7]-P[1]), 
			 nY = (P[5]-P[2])*(P[6]-P[0])-(P[3]-P[0])*(P[8]-P[2]), 
			 nZ = (P[3]-P[0])*(P[7]-P[1])-(P[4]-P[1])*(P[6]-P[0]), F;
	int    i, j;

	F = QG_normal(nX, nY, nZ)*.125;
	if (N_par < 1) add_params(1-N_par); 

	for (i = 0; i < M; i++) {
//		dm = c2M[i][0]*d;
		dm[0] = c2M[i][0]*d[0];
		dm[1] = c2M[i][0]*d[1];
		dm[2] = c2M[i][0]*d[2];

//		sum += c2M[i][1]*Func(xm+dm);
		d_i[0] = .5*(xm[0]+dm[0]-P[0]); xm_i[0] = d_i[0]+P[0];
		d_i[1] = .5*(xm[1]+dm[1]-P[1]); xm_i[1] = d_i[1]+P[1];
		d_i[2] = .5*(xm[2]+dm[2]-P[2]); xm_i[2] = d_i[2]+P[2];
		if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);
		for (j = 0; j < M; j++) {
			dm_i[0] = c2M[j][0]*d_i[0];
			dm_i[1] = c2M[j][0]*d_i[1];
			dm_i[2] = c2M[j][0]*d_i[2];

			pp[0] = c2M[i][1]*c2M[j][1]*(1.+c2M[j][0])*F;
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M[i][1]*c2M[j][1]*(1.-c2M[j][0])*F;
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}

//		sum += c2M[i][1]*Func(xm-dm);
		d_i[0] = .5*(xm[0]-dm[0]-P[0]); xm_i[0] = d_i[0]+P[0];
		d_i[1] = .5*(xm[1]-dm[1]-P[1]); xm_i[1] = d_i[1]+P[1];
		d_i[2] = .5*(xm[2]-dm[2]-P[2]); xm_i[2] = d_i[2]+P[2];
		if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);
		for (j = 0; j < M; j++) {
			dm_i[0] = c2M[j][0]*d_i[0];
			dm_i[1] = c2M[j][0]*d_i[1];
			dm_i[2] = c2M[j][0]*d_i[2];

			pp[0] = c2M[i][1]*c2M[j][1]*(1.+c2M[j][0])*F;
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M[i][1]*c2M[j][1]*(1.-c2M[j][0])*F;
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}
	}
}

///////////////////////////////////////////////////////////////////////////
//...���������� ������ ��p���� 2M ��� ��������������� ������ (���� � ����);
void CGrid_QG::QG2M_quad(double * P, int M, const double c2M[][2], Num_State normal)
{//...������� ���������� ����� ������� ��� � OpenGL (m_QUAD = 1);
	int    i, j, m_QUAD = 1;
	double d1 [3] = {.5*(P[3]-P[0]), .5*(P[4]-P[1]), .5*(P[5]-P[2])}, 
			 xm1[3] = { d1[0]+P[0], d1[1]+P[1], d1[2]+P[2] },
	       d2 [3] = {.5*(P[9]-P[6]), .5*(P[10]-P[7]), .5*(P[11]-P[8])}, 
			 xm2[3] = { d2[0]+P[6], d2[1]+P[7], d2[2]+P[8] },
			 dm1[3], dm2[3], d_i[3], xm_i[3], dm_i[3], pp[1], nX, nY, nZ, F1, F2, F3, F4;
	if (! m_QUAD) {
		d2[0] = -d2[0];
		d2[1] = -d2[1];
		d2[2] = -d2[2];
	}
	nX = (P[10]-P[4])*(P[11]-P[8])-(P[11]-P[5])*(P[10]-P[7]); 
	nY = (P[11]-P[5])*(P[9] -P[6])-(P[9] -P[3])*(P[11]-P[8]); 
	nZ = (P[9] -P[3])*(P[10]-P[7])-(P[10]-P[4])*(P[9] -P[6]);
	F4 = QG_normal(nX, nY, nZ)*.0625;

	nX = (P[7]-P[1])*(P[11]-P[8])-(P[8]-P[2])*(P[10]-P[7]); 
	nY = (P[8]-P[2])*(P[9] -P[6])-(P[6]-P[0])*(P[11]-P[8]); 
	nZ = (P[6]-P[0])*(P[10]-P[7])-(P[7]-P[1])*(P[9] -P[6]);
	F3 = QG_normal(nX, nY, nZ)*.0625;

	nX = (P[10]-P[4])*(P[5]-P[2])-(P[11]-P[5])*(P[4]-P[1]); 
	nY = (P[11]-P[5])*(P[3]-P[0])-(P[9] -P[3])*(P[5]-P[2]); 
	nZ = (P[9] -P[3])*(P[4]-P[1])-(P[10]-P[4])*(P[3]-P[0]);
	F2 = QG_normal(nX, nY, nZ)*.0625;

	nX = (P[7]-P[1])*(P[5]-P[2])-(P[8]-P[2])*(P[4]-P[1]); 
	nY = (P[8]-P[2])*(P[3]-P[0])-(P[6]-P[0])*(P[5]-P[2]); 
	nZ = (P[6]-P[0])*(P[4]-P[1])-(P[7]-P[1])*(P[3]-P[0]);
	F1 = QG_normal(nX, nY, nZ)*.0625;

	if (N_par < 1) add_params(1-N_par); 

	for (i = 0; i < M; i++) {
//		dm = c2M[i][0]*d;
		dm1[0] = c2M[i][0]*d1[0]; dm2[0] = c2M[i][0]*d2[0];
		dm1[1] = c2M[i][0]*d1[1]; dm2[1] = c2M[i][0]*d2[1];
		dm1[2] = c2M[i][0]*d1[2]; dm2[2] = c2M[i][0]*d2[2];

//		sum += c2M[i][1]*Func(xm+dm);
		d_i[0] = .5*(xm2[0]+dm2[0]-xm1[0]-dm1[0]); xm_i[0] = d_i[0]+xm1[0]+dm1[0];
		d_i[1] = .5*(xm2[1]+dm2[1]-xm1[1]-dm1[1]); xm_i[1] = d_i[1]+xm1[1]+dm1[1];
		d_i[2] = .5*(xm2[2]+dm2[2]-xm1[2]-dm1[2]); xm_i[2] = d_i[2]+xm1[2]+dm1[2];
		if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);
		for (j = 0; j < M; j++) {
			dm_i[0] = c2M[j][0]*d_i[0];
			dm_i[1] = c2M[j][0]*d_i[1];
			dm_i[2] = c2M[j][0]*d_i[2];

			pp[0] = c2M[i][1]*c2M[j][1]*(F1*(1.-c2M[i][0])*(1.-c2M[j][0])+F2*(1.+c2M[i][0])*(1.-c2M[j][0])+
												  F3*(1.-c2M[i][0])*(1.+c2M[j][0])+F4*(1.+c2M[i][0])*(1.+c2M[j][0]));
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M[i][1]*c2M[j][1]*(F1*(1.-c2M[i][0])*(1.+c2M[j][0])+F2*(1.+c2M[i][0])*(1.+c2M[j][0])+
												  F3*(1.-c2M[i][0])*(1.-c2M[j][0])+F4*(1.+c2M[i][0])*(1.-c2M[j][0]));
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}

//		sum += c2M[i][1]*Func(xm-dm);
		d_i[0] = .5*(xm2[0]-dm2[0]-xm1[0]+dm1[0]); xm_i[0] = d_i[0]+xm1[0]-dm1[0];
		d_i[1] = .5*(xm2[1]-dm2[1]-xm1[1]+dm1[1]); xm_i[1] = d_i[1]+xm1[1]-dm1[1];
		d_i[2] = .5*(xm2[2]-dm2[2]-xm1[2]+dm1[2]); xm_i[2] = d_i[2]+xm1[2]-dm1[2];
		if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);
		for (j = 0; j < M; j++) {
			dm_i[0] = c2M[j][0]*d_i[0];
			dm_i[1] = c2M[j][0]*d_i[1];
			dm_i[2] = c2M[j][0]*d_i[2];

			pp[0] = c2M[i][1]*c2M[j][1]*(F1*(1.+c2M[i][0])*(1.-c2M[j][0])+F2*(1.-c2M[i][0])*(1.-c2M[j][0])+
												  F3*(1.+c2M[i][0])*(1.+c2M[j][0])+F4*(1.-c2M[i][0])*(1.+c2M[j][0]));
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M[i][1]*c2M[j][1]*(F1*(1.+c2M[i][0])*(1.+c2M[j][0])+F2*(1.-c2M[i][0])*(1.+c2M[j][0])+
												  F3*(1.+c2M[i][0])*(1.-c2M[j][0])+F4*(1.-c2M[i][0])*(1.-c2M[j][0]));
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}
	}
}

/////////////////////////////////////////////////////////////
//...���������� ������ ��p���� 2M-1 ��� ������ (���� � ����);
void CGrid_QG::QG2M1_line(double * P, int M, const double c2M1[][2], Num_State normal)
{
	double d1 [3] = {.5*(P[3]-P[0]), .5*(P[4]-P[1]), 0.}, 
			 xm1[3] = { d1[0]+P[0], d1[1]+P[1], 0. },
			 nX = -d1[1], 
			 nY =  d1[0], 
			 nZ =  d1[2], pp[1], F;
	int    i;

	F = QG_normal(nX, nY, nZ);
	if (normal == NULL_STATE) QG_normal(nX = d1[0], nY = d1[1], nZ = d1[2]);

	if (N_par < 1) add_params(1-N_par); 

	pp[0] = c2M1[0][1]*F;
	add_new_point(xm1[0], xm1[1], 0., nX, nY, nZ, pp); 

	for (i = 1; i < M; i++) {
		pp[0] = c2M1[i][1]*F;
		add_new_point(xm1[0]+c2M1[i][0]*d1[0], xm1[1]+c2M1[i][0]*d1[1], 0., nX, nY, nZ, pp); 

		pp[0] = c2M1[i][1]*F;
		add_new_point(xm1[0]-c2M1[i][0]*d1[0], xm1[1]-c2M1[i][0]*d1[1], 0., nX, nY, nZ, pp); 
	}
}

/////////////////////////////////////////////////////////////////////////
//...���������� ������ ��p���� 2M-1 ��� ����������� ������ (���� � ����);
void CGrid_QG::QG2M1_tria(double * P, int M, const double c2M1[][2], Num_State normal)
{
	double d [3] = {.5*(P[6]-P[3]), .5*(P[7]-P[4]), .5*(P[8]-P[5])}, 
			 xm[3] = { d[0]+P[3], d[1]+P[4], d[2]+P[5] },
			 dm[3], d_i[3], xm_i[3], dm_i[3], pp[1], 
			 nX = (P[4]-P[1])*(P[8]-P[2])-(P[5]-P[2])*(P[7]-P[1]), 
			 nY = (P[5]-P[2])*(P[6]-P[0])-(P[3]-P[0])*(P[8]-P[2]), 
			 nZ = (P[3]-P[0])*(P[7]-P[1])-(P[4]-P[1])*(P[6]-P[0]), F;
	int    i, j;

	F = QG_normal(nX, nY, nZ)*.125;
	if (N_par < 1) add_params(1-N_par); 

//	sum += c2M1[0][1]*Func(xm);
	d_i[0] = .5*(xm[0]-P[0]); xm_i[0] = d_i[0]+P[0];
	d_i[1] = .5*(xm[1]-P[1]); xm_i[1] = d_i[1]+P[1];
	d_i[2] = .5*(xm[2]-P[2]); xm_i[2] = d_i[2]+P[2];
	if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);

	pp[0] = c2M1[0][1]*c2M1[0][1]*F;
	add_new_point(xm_i[0], xm_i[1], xm_i[2], nX, nY, nZ, pp); 

	for (j = 1; j < M; j++) {
		dm_i[0] = c2M1[j][0]*d_i[0];
		dm_i[1] = c2M1[j][0]*d_i[1];
		dm_i[2] = c2M1[j][0]*d_i[2];

		pp[0] = c2M1[0][1]*c2M1[j][1]*(1.+c2M1[j][0])*F;
		add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
		pp[0] = c2M1[0][1]*c2M1[j][1]*(1.-c2M1[j][0])*F;
		add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
	}
	for (i = 1; i < M; i++) {
//		dm = c2M1[i][0]*d;
		dm[0] = c2M1[i][0]*d[0];
		dm[1] = c2M1[i][0]*d[1];
		dm[2] = c2M1[i][0]*d[2];

//		sum += c2M1[i][1]*Func(xm+dm);
		d_i[0] = .5*(xm[0]+dm[0]-P[0]); xm_i[0] = d_i[0]+P[0];
		d_i[1] = .5*(xm[1]+dm[1]-P[1]); xm_i[1] = d_i[1]+P[1];
		d_i[2] = .5*(xm[2]+dm[2]-P[2]); xm_i[2] = d_i[2]+P[2];
		if (normal == NULL_STATE) QG_normal (nX = d_i[0], nY = d_i[1], nZ = d_i[2]);

		pp[0] = c2M1[i][1]*c2M1[0][1]*F;
		add_new_point(xm_i[0], xm_i[1], xm_i[2], nX, nY, nZ, pp); 

		for (j = 1; j < M; j++) {
			dm_i[0] = c2M1[j][0]*d_i[0];
			dm_i[1] = c2M1[j][0]*d_i[1];
			dm_i[2] = c2M1[j][0]*d_i[2];

			pp[0] = c2M1[i][1]*c2M1[j][1]*(1.+c2M1[j][0])*F;
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M1[i][1]*c2M1[j][1]*(1.-c2M1[j][0])*F;
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}

//		sum += c2M1[i][1]*Func(xm-dm);
		d_i[0] = .5*(xm[0]-dm[0]-P[0]); xm_i[0] = d_i[0]+P[0];
		d_i[1] = .5*(xm[1]-dm[1]-P[1]); xm_i[1] = d_i[1]+P[1];
		d_i[2] = .5*(xm[2]-dm[2]-P[2]); xm_i[2] = d_i[2]+P[2];
		if (normal == NULL_STATE) QG_normal (nX = d_i[0], nY = d_i[1], nZ = d_i[2]);

		pp[0] = c2M1[i][1]*c2M1[0][1]*F;
		add_new_point(xm_i[0], xm_i[1], xm_i[2], nX, nY, nZ, pp); 

		for (j = 1; j < M; j++) {
			dm_i[0] = c2M1[j][0]*d_i[0];
			dm_i[1] = c2M1[j][0]*d_i[1];
			dm_i[2] = c2M1[j][0]*d_i[2];

			pp[0] = c2M1[i][1]*c2M1[j][1]*(1.+c2M1[j][0])*F;
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M1[i][1]*c2M1[j][1]*(1.-c2M1[j][0])*F;
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
//...���������� ������ ��p���� 2M-1 ��� ��������������� ������ (���� � ����);
void CGrid_QG::QG2M1_quad(double * P, int M, const double c2M1[][2], Num_State normal)
{//...������� ���������� ����� ������� ��� � OpenGL (m_QUAD = 1);
	int    i, j, m_QUAD = 1;
	double d1 [3] = {.5*(P[3]-P[0]), .5*(P[4]-P[1]), .5*(P[5]-P[2])}, 
			 xm1[3] = { d1[0]+P[0], d1[1]+P[1], d1[2]+P[2] },
	       d2 [3] = {.5*(P[9]-P[6]), .5*(P[10]-P[7]), .5*(P[11]-P[8])}, 
			 xm2[3] = { d2[0]+P[6], d2[1]+P[7], d2[2]+P[8] },
			 dm1[3], dm2[3], d_i[3], xm_i[3], dm_i[3], pp[1], nX, nY, nZ, F1, F2, F3, F4;
	if (! m_QUAD) {
		d2[0] = -d2[0];
		d2[1] = -d2[1];
		d2[2] = -d2[2];
	}
	nX = (P[10]-P[4])*(P[11]-P[8])-(P[11]-P[5])*(P[10]-P[7]); 
	nY = (P[11]-P[5])*(P[9] -P[6])-(P[9] -P[3])*(P[11]-P[8]); 
	nZ = (P[9] -P[3])*(P[10]-P[7])-(P[10]-P[4])*(P[9] -P[6]);
	F4 = QG_normal(nX, nY, nZ)*.0625;

	nX = (P[7]-P[1])*(P[11]-P[8])-(P[8]-P[2])*(P[10]-P[7]); 
	nY = (P[8]-P[2])*(P[9] -P[6])-(P[6]-P[0])*(P[11]-P[8]); 
	nZ = (P[6]-P[0])*(P[10]-P[7])-(P[7]-P[1])*(P[9] -P[6]);
	F3 = QG_normal(nX, nY, nZ)*.0625;

	nX = (P[10]-P[4])*(P[5]-P[2])-(P[11]-P[5])*(P[4]-P[1]); 
	nY = (P[11]-P[5])*(P[3]-P[0])-(P[9] -P[3])*(P[5]-P[2]); 
	nZ = (P[9] -P[3])*(P[4]-P[1])-(P[10]-P[4])*(P[3]-P[0]);
	F2 = QG_normal(nX, nY, nZ)*.0625;

	nX = (P[7]-P[1])*(P[5]-P[2])-(P[8]-P[2])*(P[4]-P[1]); 
	nY = (P[8]-P[2])*(P[3]-P[0])-(P[6]-P[0])*(P[5]-P[2]); 
	nZ = (P[6]-P[0])*(P[4]-P[1])-(P[7]-P[1])*(P[3]-P[0]);
	F1 = QG_normal(nX, nY, nZ)*.0625;

	if (N_par < 1) add_params(1-N_par); 

//	sum += c2M1[0][1]*Func(xm);
	d_i[0] = .5*(xm2[0]-xm1[0]); xm_i[0] = d_i[0]+xm1[0];
	d_i[1] = .5*(xm2[1]-xm1[1]); xm_i[1] = d_i[1]+xm1[1];
	d_i[2] = .5*(xm2[2]-xm1[2]); xm_i[2] = d_i[2]+xm1[2];
	if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);

	pp[0] = c2M1[0][1]*c2M1[0][1]*(F1+F2+F3+F4);
	add_new_point(xm_i[0], xm_i[1], xm_i[2], nX, nY, nZ, pp); 

	for (j = 1; j < M; j++) {
		dm_i[0] = c2M1[j][0]*d_i[0];
		dm_i[1] = c2M1[j][0]*d_i[1];
		dm_i[2] = c2M1[j][0]*d_i[2];

		pp[0] = c2M1[0][1]*c2M1[j][1]*(F1*(1.-c2M1[j][0])+F2*(1.-c2M1[j][0])+F3*(1.+c2M1[j][0])+F4*(1.+c2M1[j][0]));
		add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp);
		pp[0] = c2M1[0][1]*c2M1[j][1]*(F1*(1.+c2M1[j][0])+F2*(1.+c2M1[j][0])+F3*(1.-c2M1[j][0])+F4*(1.-c2M1[j][0]));
		add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
	}
	for (i = 1; i < M; i++) {
//		dm = c2M1[i][0]*d;
		dm1[0] = c2M1[i][0]*d1[0]; dm2[0] = c2M1[i][0]*d2[0];
		dm1[1] = c2M1[i][0]*d1[1]; dm2[1] = c2M1[i][0]*d2[1];
		dm1[2] = c2M1[i][0]*d1[2]; dm2[2] = c2M1[i][0]*d2[2];

//		sum += c2M1[i][1]*Func(xm+dm);
		d_i[0] = .5*(xm2[0]+dm2[0]-xm1[0]-dm1[0]); xm_i[0] = d_i[0]+xm1[0]+dm1[0];
		d_i[1] = .5*(xm2[1]+dm2[1]-xm1[1]-dm1[1]); xm_i[1] = d_i[1]+xm1[1]+dm1[1];
		d_i[2] = .5*(xm2[2]+dm2[2]-xm1[2]-dm1[2]); xm_i[2] = d_i[2]+xm1[2]+dm1[2];
		if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);

		pp[0] = c2M1[i][1]*c2M1[0][1]*(F1*(1.-c2M1[i][0])+F2*(1.+c2M1[i][0])+F3*(1.-c2M1[i][0])+F4*(1.+c2M1[i][0]));
		add_new_point(xm_i[0], xm_i[1], xm_i[2], nX, nY, nZ, pp); 

		for (j = 1; j < M; j++) {
			dm_i[0] = c2M1[j][0]*d_i[0];
			dm_i[1] = c2M1[j][0]*d_i[1];
			dm_i[2] = c2M1[j][0]*d_i[2];

			pp[0] = c2M1[i][1]*c2M1[j][1]*(F1*(1.-c2M1[i][0])*(1.-c2M1[j][0])+F2*(1.+c2M1[i][0])*(1.-c2M1[j][0])+
													 F3*(1.-c2M1[i][0])*(1.+c2M1[j][0])+F4*(1.+c2M1[i][0])*(1.+c2M1[j][0]));
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M1[i][1]*c2M1[j][1]*(F1*(1.-c2M1[i][0])*(1.+c2M1[j][0])+F2*(1.+c2M1[i][0])*(1.+c2M1[j][0])+
													 F3*(1.-c2M1[i][0])*(1.-c2M1[j][0])+F4*(1.+c2M1[i][0])*(1.-c2M1[j][0]));
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}

//		sum += c2M1[i][1]*Func(xm-dm);
		d_i[0] = .5*(xm2[0]-dm2[0]-xm1[0]+dm1[0]); xm_i[0] = d_i[0]+xm1[0]-dm1[0];
		d_i[1] = .5*(xm2[1]-dm2[1]-xm1[1]+dm1[1]); xm_i[1] = d_i[1]+xm1[1]-dm1[1];
		d_i[2] = .5*(xm2[2]-dm2[2]-xm1[2]+dm1[2]); xm_i[2] = d_i[2]+xm1[2]-dm1[2];
		if (normal == NULL_STATE) QG_normal(nX = d_i[0], nY = d_i[1], nZ = d_i[2]);

		pp[0] = c2M1[i][1]*c2M1[0][1]*(F1*(1.+c2M1[i][0])+F2*(1.-c2M1[i][0])+F3*(1.+c2M1[i][0])+F4*(1.-c2M1[i][0]));
		add_new_point(xm_i[0], xm_i[1], xm_i[2], nX, nY, nZ, pp); 

		for (j = 1; j < M; j++) {
			dm_i[0] = c2M1[j][0]*d_i[0];
			dm_i[1] = c2M1[j][0]*d_i[1];
			dm_i[2] = c2M1[j][0]*d_i[2];

			pp[0] = c2M1[i][1]*c2M1[j][1]*(F1*(1.+c2M1[i][0])*(1.-c2M1[j][0])+F2*(1.-c2M1[i][0])*(1.-c2M1[j][0])+
													 F3*(1.+c2M1[i][0])*(1.+c2M1[j][0])+F4*(1.-c2M1[i][0])*(1.+c2M1[j][0]));
			add_new_point(xm_i[0]+dm_i[0], xm_i[1]+dm_i[1], xm_i[2]+dm_i[2], nX, nY, nZ, pp); 
			pp[0] = c2M1[i][1]*c2M1[j][1]*(F1*(1.+c2M1[i][0])*(1.+c2M1[j][0])+F2*(1.-c2M1[i][0])*(1.+c2M1[j][0])+
													 F3*(1.+c2M1[i][0])*(1.-c2M1[j][0])+F4*(1.-c2M1[i][0])*(1.-c2M1[j][0]));
			add_new_point(xm_i[0]-dm_i[0], xm_i[1]-dm_i[1], xm_i[2]-dm_i[2], nX, nY, nZ, pp); 
		}
	}
}

void CGrid_QG::QG_sheet(CMap * mp, int N_elem, int N_max)
{
	int k;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[k = size_of_map(2, NULL_GENUS)] == (CMap)SHEET_CELL
			 && fabs(mp[k+1]) > EE_ker && fabs(mp[k+2]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double pp[3] = {mp[1], mp[2], mp[3]}, temp_mpk1 = mp[k+1], temp_mpk2 = mp[k+2], A, B, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);
      mp[k+1] /= N_max;
      mp[k+2] /= N_max;

/////////////////////////////////
//...����������� ��������� �����;
		for (B = (mp[k+2]-temp_mpk2)*.5, m = 0; m < N_max; m++, B += mp[k+2])
		for (A = (mp[k+1]-temp_mpk1)*.5, l = 0; l < N_max; l++, A += mp[k+1]) {
			mp[1] = A;
			mp[2] = B;
			mp[3] = 0.; point_iso<double>(mp+1, pp, CZ, SZ, CY, SY, CX, SX);

			QG_sheet(mp, N_elem);
		}

//////////////////////////////
//...��������������� ��������;
		mp[1]   = pp[0]; 
		mp[2]   = pp[1]; 
		mp[3]   = pp[2]; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
	}
}

void CGrid_QG::QG_sheet(CMap * mp, int N_elem)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[k = size_of_map(2, NULL_GENUS)] == (CMap)SHEET_CELL 
			 && fabs(mp[k+1]) > EE_ker && fabs(mp[k+2]) > EE_ker) {
		double A = fabs(mp[k+1]), B = fabs(mp[k+2]), P[12], 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);

		P[0] = -A*.5; P[1] = -B*.5; P[2] = 0.;
		P[3] = -A*.5; P[4] =  B*.5; P[5] = 0.;
		P[6] =  A*.5; P[7] = -B*.5; P[8] = 0.;
		P[9] =  A*.5; P[10] = B*.5; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

		for (int i = N_ini; i < N; i++) {
			P[0] = X[i]; P[3] = 0.;
			P[1] = Y[i]; P[4] = 0.;
			P[2] = Z[i]; P[5] = 1.;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
		}
	}
}

void CGrid_QG::QG_ring_segment(CMap * mp, int N_elem, int N_max)
{
	int k;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[k = size_of_map(2, NULL_GENUS)] == (CMap)RING_SEGMENT 
			 && fabs(mp[k+1]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double temp_mp6 = mp[6], 
				temp_mpk1 = mp[k+1],
				temp_mpk2 = mp[k+2],
				temp_mpk3 = mp[k+3], d = (temp_mpk3-temp_mpk2)/N_max;
      mp[k+1] /= N_max;
      mp[k+3]  = temp_mpk2+d;

/////////////////////////////////
//...����������� ��������� �����;
		for (m = 0; m < N_max; m++, mp[k+2] = mp[k+3], mp[k+3] = temp_mpk2+d*(m+1))
		for (l = 0; l < N_max; l++) {
			mp[6] = temp_mp6+mp[k+1]*(l-(N_max-1.)/2.);
			QG_ring_segment(mp, N_elem);
		}
		mp[6]   = temp_mp6; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
		mp[k+3] = temp_mpk3;
	}
}

void CGrid_QG::QG_ring_segment(CMap * mp, int N_elem)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[k = size_of_map(2, NULL_GENUS)] == (CMap)RING_SEGMENT 
			 && fabs(mp[k+1]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], C = mp[k+3], RR, P[12],
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);

		P[0] = -A*.5; P[1]  = B; P[2] = 0.;
		P[3] = -A*.5; P[4]  = C; P[5] = 0.;
		P[6] =  A*.5; P[7]  = B; P[8] = 0.;
		P[9] =  A*.5; P[10] = C; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);
		P[3] = 0.; 
		P[4] = 0.; 
		P[5] = 1.; 
		point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);

		for (int i = N_ini; i < N; i++) {
			P[0] = Y[i]*cos(X[i]);
			P[1] = Y[i]*sin(X[i]); RR = Y[i];
			P[2] = 0.;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1]; nX[i] = P[3];
			Y[i] = P[1]+mp[2]; nY[i] = P[4];
			Z[i] = P[2]+mp[3]; nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*RR);
		}
	}
}

void CGrid_QG::QG_cyl_segment(CMap * mp, int N_elem, int N_max)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, CYL_GENUS) == mp[0] && mp[k = size_of_map(2, CYL_GENUS)] == (CMap)CYL_SEGMENT 
			 && fabs(mp[k+1]) > EE_ker && fabs(mp[k+2]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double temp_mp3 = mp[3], 
				 temp_mp6 = mp[6], 
				temp_mpk1 = mp[k+1],
				temp_mpk2 = mp[k+2], d = temp_mpk2/N_max;
      mp[k+1] /= N_max;
      mp[k+2]  = d;
		mp[3] += .5*(d-temp_mpk2);

/////////////////////////////////
//...����������� ��������� �����;
		for (m = 0; m < N_max; m++, mp[3] += d)
		for (l = 0; l < N_max; l++) {
			mp[6] = temp_mp6+mp[k+1]*(l-(N_max-1.)/2.);
			QG_cyl_segment(mp, N_elem);
		}
		mp[3]   = temp_mp3; 
		mp[6]   = temp_mp6; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
	}
}

void CGrid_QG::QG_cyl_segment(CMap * mp, int N_elem)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, CYL_GENUS) == mp[0] && mp[k = size_of_map(2, CYL_GENUS)] == (CMap)CYL_SEGMENT 
			 && fabs(mp[k+1]) > EE_ker && fabs(mp[k+2]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], RR = fabs(mp[7]), P[12], f = mp[7] < 0. ? -1. : 1.,
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);

		P[0] = -A*.5; P[1] = -B*.5; P[2] = 0.;
		P[3] = -A*.5; P[4] =  B*.5; P[5] = 0.;
		P[6] =  A*.5; P[7] = -B*.5; P[8] = 0.;
		P[9] =  A*.5; P[10] = B*.5; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

		for (int i = N_ini; i < N; i++) {
			P[3] = cos(X[i]); P[0] = RR*P[3]; P[3] *= f;
			P[4] = sin(X[i]); P[1] = RR*P[4]; P[4] *= f;
			P[5] = 0.;		   P[2] = Y[i];

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*RR);
		}
	}
}

void CGrid_QG::QG_sph_segment(CMap * mp, int N_elem, int N_max)
{
	int k;
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && mp[k = size_of_map(2, SPHERE_GENUS)] == (CMap)SPH_SEGMENT
			 && fabs(mp[k+1]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double temp_mp6 = mp[6], 
				temp_mpk1 = mp[k+1],
				temp_mpk2 = mp[k+2],
				temp_mpk3 = mp[k+3], d = (temp_mpk3-temp_mpk2)/N_max;
      mp[k+1] /= N_max;
      mp[k+3]  = temp_mpk2+d;

/////////////////////////////////
//...����������� ��������� �����;
		for (m = 0; m < N_max; m++, mp[k+2] = mp[k+3], mp[k+3] = temp_mpk2+d*(m+1))
		for (l = 0; l < N_max; l++) {
			mp[6] = temp_mp6+mp[k+1]*(l-(N_max-1.)/2.);
			QG_sph_segment(mp, N_elem);
		}
		mp[6]   = temp_mp6; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
		mp[k+3] = temp_mpk3;
	}
}

void CGrid_QG::QG_sph_segment(CMap * mp, int N_elem)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && mp[k = size_of_map(2, SPHERE_GENUS)] == (CMap)SPH_SEGMENT
			 && fabs(mp[k+1]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], C = mp[k+3], RR = fabs(mp[7]), P[12], f = (mp[7] < 0. ? -1. : 1.)/**(C < B ? -1. : 1.)*/, //...����������� ��� ������ � ����� �������;
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), g;

		P[0] = -A*.5; P[1]  = B; P[2] = 0.;
		P[3] = -A*.5; P[4]  = C; P[5] = 0.;
		P[6] =  A*.5; P[7]  = B; P[8] = 0.;
		P[9] =  A*.5; P[10] = C; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

		for (int i = N_ini; i < N; i++) {
			P[3] = cos(X[i])*(g = sin(Y[i])*f);	P[0] = RR*P[3]*f;
			P[4] = sin(X[i])* g;						P[1] = RR*P[4]*f;
			P[5] = cos(Y[i])* f;						P[2] = RR*P[5]*f;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*sqr(RR)*fabs(g));
		}
	}
}

void CGrid_QG::QG_spr_segment(CMap * mp, int N_elem, int N_max)
{
	int k;
	if (mp && ID_MAP(2, SPHEROID_GENUS) == mp[0] && mp[k = size_of_map(2, SPHEROID_GENUS)] == (CMap)SPR_SEGMENT
			 && fabs(mp[k+1]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double temp_mp6 = mp[6], 
				temp_mpk1 = mp[k+1],
				temp_mpk2 = mp[k+2],
				temp_mpk3 = mp[k+3], d = (temp_mpk3-temp_mpk2)/N_max;
      mp[k+1] /= N_max;
      mp[k+3]  = temp_mpk2+d;

/////////////////////////////////
//...����������� ��������� �����;
		for (m = 0; m < N_max; m++, mp[k+2] = mp[k+3], mp[k+3] = temp_mpk2+d*(m+1))
		for (l = 0; l < N_max; l++) {
			mp[6] = temp_mp6+mp[k+1]*l;
			QG_spr_segment(mp, N_elem);
		}
		mp[6]   = temp_mp6; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
		mp[k+3] = temp_mpk3;
	}
}

void CGrid_QG::QG_spr_segment(CMap * mp, int N_elem)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, SPHEROID_GENUS) == mp[0] && mp[k = size_of_map(2, SPHEROID_GENUS)] == (CMap)SPR_SEGMENT
			 && fabs(mp[k+1]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], C = mp[k+3], RR = fabs(mp[7]), rr = fabs(mp[8]), P[12], f = (mp[7] < 0. ? -1. : 1.)*(C < B ? -1. : 1.),
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), Co, Si, g, d, delt = rr/RR;

		P[0] = -A*.5; P[1]  = B; P[2] = 0.;
		P[3] = -A*.5; P[4]  = C; P[5] = 0.;
		P[6] =  A*.5; P[7]  = B; P[8] = 0.;
		P[9] =  A*.5; P[10] = C; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

		for (int i = N_ini; i < N; i++) {
			Co = cos(Y[i]); Si = sin(Y[i]); 
			d = sqrt(sqr(Si)+sqr(delt*Co));
			P[0] = cos(X[i])*Si; P[3] = P[0]*(g = f/d); P[0] *= rr;
			P[1] = sin(X[i])*Si; P[4] = P[1]* g; P[1] *= rr;				 
			P[2] = RR*Co;			P[5] = Co*delt*g;						 

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*sqr(RR)*Si*d);
		}
	}
}

void CGrid_QG::QG_sph_intrusion(CMap * mp, int N_elem, int N_max)
{
	int k;
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && mp[k = size_of_map(2, SPHERE_GENUS)] == (CMap)SPH_INTRUSION_CELL
			 && fabs(mp[k+1]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double temp_mp6 = mp[6], 
				temp_mpk1 = mp[k+1],
				temp_mpk2 = mp[k+2],
				temp_mpk3 = mp[k+3], d = (temp_mpk3-temp_mpk2)/N_max;
      mp[k+1] /= N_max;
      mp[k+3]  = temp_mpk2+d;

/////////////////////////////////
//...����������� ��������� �����;
		for (m = 0; m < N_max; m++, mp[k+2] = mp[k+3], mp[k+3] = temp_mpk2+d*(m+1))
		for (l = 0; l < N_max; l++) {
			mp[6] = temp_mp6+mp[k+1]*l;
			QG_sph_intrusion(mp, N_elem);
		}
		mp[6]   = temp_mp6; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
		mp[k+3] = temp_mpk3;
	}
}

void CGrid_QG::QG_sph_intrusion(CMap * mp, int N_elem)
{
	int i, k, N_ini = N;
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && mp[k = size_of_map(2, SPHERE_GENUS)] == (CMap)SPH_INTRUSION_CELL
			 && fabs(mp[k+1]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], C = mp[k+3], L = mp[k+4], RR = fabs(mp[7]), P[12], f = (mp[7] < 0. ? -1. : 1.)/**(C < B ? -1. : 1.)*/,
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), g, alpha, theta, fB, fC, fL, fY;

		P[0] = -A*.5; P[1]  = B; P[2] = 0.;
		P[3] = -A*.5; P[4]  = C; P[5] = 0.;
		P[6] =  A*.5; P[7]  = B; P[8] = 0.;
		P[9] =  A*.5; P[10] = C; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

//////////////////////////////////////////////////
//...��������� ���������� ��� ������� �����������;
		if (L < RR && RR > 0. && fabs(C-B) > 0.) {
			alpha = L/RR;
			theta = acos(alpha);
			fB = min(B, C);
			fC = max(B, C); B = fB; C = fC;
			for ( i = N_ini; i < N; i++) {
				fB = B;
				fC = C; 
				fY = 0.;
				fL = 1.;
				for (k = -4; k <= 4; k++)
				if (fabs(X[i]+mp[6]-k*M_PI_2) < theta) {
					fY = acos(alpha/cos(X[i]+mp[6]-k*M_PI_2)); break;
				}
				if ( Y[i] <= M_PI_2) {
					if (fB <  theta) fB = theta;
					if (fC <  theta) fB = fC = (hit[i] = 0)+theta;
					if (fB >= M_PI_2-fY) fB = fC = (hit[i] = 0)+M_PI_2-fY; 
					if (fC <= M_PI_2-fY) Y[i] = fB+(Y[i]-B)*(fL = (fC-fB)/(C-B)); else
					if (fC <= M_PI_2) Y[i] = fB+(Y[i]-B)*(fL = (M_PI_2-fY-fB)/(C-B)); else
											Y[i] = fB+(Y[i]-B)*(fL = (M_PI_2-fY-fB)/(M_PI_2-B));
				}
				else {
					if (fC >  M_PI-theta) fC = M_PI-theta;
					if (fB >  M_PI-theta) fC = fB = (hit[i] = 0)+M_PI-theta;
					if (fC <= M_PI_2+fY)  fC = fB = (hit[i] = 0)+M_PI_2+fY;
					if (fB >= M_PI_2+fY) Y[i] = fC-(C-Y[i])*(fL = (fC-fB)/(C-B)); else
					if (fB >= M_PI_2) Y[i] = fC-(C-Y[i])*(fL = (fC-M_PI_2-fY)/(C-B)); else
											Y[i] = fC-(C-Y[i])*(fL = (fC-M_PI_2-fY)/(C-M_PI_2));
				}
				set_param(0, i, get_param(0, i)*fL);
			}
		}

/////////////////////////////////////////////////////////////
//...������� ����� �������������� � ���������������� �������;
		for ( i = N_ini; i < N; i++) {
			P[3] = cos(X[i])*(g = sin(Y[i])*f);	P[0] = RR*P[3]*f;
			P[4] = sin(X[i])* g;						P[1] = RR*P[4]*f;
			P[5] = cos(Y[i])* f;						P[2] = RR*P[5]*f;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*sqr(RR)*fabs(g));
		}
	}
}

void CGrid_QG::QG_sht_intrusion(CMap * mp, int N_elem, int N_max)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[k = size_of_map(2, NULL_GENUS)] == (CMap)SHT_INTRUSION_CELL
			 && fabs(mp[k+1]) > EE_ker && fabs(mp[k+2]) > EE_ker) {

		double A = mp[k+1], B = mp[k+2], RR = mp[k+3], P[6],
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), fX, fY, fL, fR, RL, AB = min(A, B);

////////////////////////////////////////////////
//...���������� ���������� ��� ��������� ������;
		CMap * ring = get_map(2, NULL_GENUS, RING_SEGMENT);		
		if (ring) {
			ring[++(k = size_of_map(ring))] = (CMap)2.*M_PI;
			ring[++k] = (CMap)fabs(RR);
			ring[++k] = (CMap)min(fabs(A), fabs(B))*.5;

			QG_ring_segment(ring, N_elem, N_max);
			delete_struct	(ring);
		}

//////////////////////////////////////////////////
//...��������� ���������� ��� ������� �����������;
		for (int i = N_ini; i < N; i++) {
			RL = sqrt(X[i]*X[i]+Y[i]*Y[i]);
			fX = X[i]*(fR = RR/RL);
			fY = Y[i]* fR;
			if (A*fabs(Y[i]) < B*fabs(X[i])) fL = (sqrt(sqr(Y[i]/X[i])+1.)*A*.5-RR)/(AB*.5-RR);	else
			if (B*fabs(X[i]) < A*fabs(Y[i])) fL = (sqrt(sqr(X[i]/Y[i])+1.)*B*.5-RR)/(AB*.5-RR);
			P[0] = (X[i]-fX)*fL+fX; P[3] = nX[i];
			P[1] = (Y[i]-fY)*fL+fY; P[4] = nY[i];
			P[2] = 0.;				   P[5] = nZ[i];

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*fL*(fR+(1.-fR)*fL));
		}
	}
}

void CGrid_QG::QG_cone_segment(CMap * mp, int N_elem, int N_max)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, CONE_GENUS) == mp[0] && mp[k = size_of_map(2, CONE_GENUS)] == (CMap)CONE_SEGMENT 
			 && fabs(mp[k+1]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double temp_mp6 = mp[6], 
				temp_mpk1 = mp[k+1],
				temp_mpk2 = mp[k+2],
				temp_mpk3 = mp[k+3], d = (temp_mpk3-temp_mpk2)/N_max;
      mp[k+1] /= N_max;
      mp[k+3]  = temp_mpk2+d;

/////////////////////////////////
//...����������� ��������� �����;
		for (m = 0; m < N_max; m++, mp[k+2] = mp[k+3], mp[k+3] = temp_mpk2+d*(m+1))
		for (l = 0; l < N_max; l++) {
			mp[6] = temp_mp6+mp[k+1]*(l-(N_max-1.)/2.);
			QG_cone_segment(mp, N_elem);
		}
		mp[6]   = temp_mp6; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
		mp[k+3] = temp_mpk3;
	}
}

void CGrid_QG::QG_cone_segment(CMap * mp, int N_elem)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, CONE_GENUS) == mp[0] && mp[k = size_of_map(2, CONE_GENUS)] == (CMap)CONE_SEGMENT
			 && fabs(mp[k+1]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], C = mp[k+3], theta = fabs(mp[7]), P[12], f = (mp[7] < 0. ? -1. : 1.),
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), Co = cos(theta), Si = sin(theta), g;

		P[0] = -A*.5; P[1]  = B; P[2] = 0.;
		P[3] = -A*.5; P[4]  = C; P[5] = 0.;
		P[6] =  A*.5; P[7]  = B; P[8] = 0.;
		P[9] =  A*.5; P[10] = C; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

		for (int i = N_ini; i < N; i++) {
			P[3] = cos(X[i]); P[0] = (g = Y[i]*Si)*P[3]; P[3] *= Co*f; 
			P[4] = sin(X[i]); P[1] =  g*P[4];				P[4] *= Co*f;
			P[5] = -Si*f;		P[2] = Y[i]*Co;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*fabs(g));
		}
	}
}

void CGrid_QG::QG_torus_segment(CMap * mp, int N_elem, int N_max)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, TORUS_GENUS) == mp[0] && mp[k = size_of_map(2, TORUS_GENUS)] == (CMap)TORUS_SEGMENT 
			 && fabs(mp[k+1]) > EE_ker) {
		int l, m; N_max = max(1, N_max);
		double temp_mp6 = mp[6], 
				temp_mpk1 = mp[k+1],
				temp_mpk2 = mp[k+2],
				temp_mpk3 = mp[k+3], d = (temp_mpk3-temp_mpk2)/N_max;
      mp[k+1] /= N_max;
      mp[k+3]  = temp_mpk2+d;

/////////////////////////////////
//...����������� ��������� �����;
		for (m = 0; m < N_max; m++, mp[k+2] = mp[k+3], mp[k+3] = temp_mpk2+d*(m+1))
		for (l = 0; l < N_max; l++) {
			mp[6] = temp_mp6+mp[k+1]*(l-(N_max-1.)/2.);
			QG_torus_segment(mp, N_elem);
		}
		mp[6]   = temp_mp6; 
		mp[k+1] = temp_mpk1;
		mp[k+2] = temp_mpk2;
		mp[k+3] = temp_mpk3;
	}
}

void CGrid_QG::QG_torus_segment(CMap * mp, int N_elem)
{
	int k, N_ini = N;
	if (mp && ID_MAP(2, TORUS_GENUS) == mp[0] && mp[k = size_of_map(2, TORUS_GENUS)] == (CMap)TORUS_SEGMENT
			 && fabs(mp[k+1]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], C = mp[k+3], RR = fabs(mp[7]), rr = fabs(mp[8]), P[12], f = (mp[7] < 0. ? -1. : 1.),
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), g;

		P[0] = -A*.5; P[1]  = B; P[2] = 0.;
		P[3] = -A*.5; P[4]  = C; P[5] = 0.;
		P[6] =  A*.5; P[7]  = B; P[8] = 0.;
		P[9] =  A*.5; P[10] = C; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

		for (int i = N_ini; i < N; i++) {
			P[3] = cos(X[i]);	P[0] = P[3]*(RR+rr*(g = cos(Y[i]))); P[3] *= f*g;	
			P[4] = sin(X[i]);	P[1] = P[4]*(RR+rr* g); P[4] *= f*g;						
			P[5] = sin(Y[i]);	P[2] = rr*P[5];			P[5] *= f;						

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)*sqr(RR)*fabs(g));
		}
	}
}

void CGrid_QG::QG_ugolok_cell(CMap * mp, int N_elem)
{
	int k, N_ini = N, N_ini2;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[k = size_of_map(2, NULL_GENUS)] == (CMap)UGOLOK_CELL) {
		double A1 = mp[k+1], A2 = mp[k+2], B1 = mp[k+3], B2 = mp[k+4], 
				rad = mp[k+6], radc = mp[k+5], beta = mp[k+7]/180., P[12], 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), C0 = cos(M_PI_2*beta), S0 = sin(M_PI_2*beta), sgn = (beta > 1. ? -1. : 1.),
				t1 = 0., t2 = 0., C = cos(M_PI*beta), S = sin(M_PI*beta), T = fabs(tan(M_PI_2*(beta-1.))), R0 = radc*T;
		if (fabs(S) > EE_dop) {
			t1 = (B2+B1*C)/S; t2 = (B1+B2*C)/S;
		}
		CMap mp1[] = {ID_MAP(1, SPHERE_GENUS), R0*C0+radc*S0*sgn, R0*S0-radc*C0*sgn, 0., 0., 0., 0., radc, -1.},
			  mp2[] = {ID_MAP(1, SPHERE_GENUS), (t2+rad*T)*C0+B2*S0+rad*S0*sgn, (t2+rad*T)*S0-B2*C0-rad*C0*sgn, 0., 0., 0., 0., rad, -1.};

		P[0] = R0*C0; P[1] = R0*S0; P[2] = 0.;
		P[3] = A2*C0; P[4] = A2*S0; P[5] = 0.;
		P[6] = (t2+rad*T)*C0+B2*S0; P[7] = (t2+rad*T)*S0-B2*C0; P[8] = 0.;
		P[9] = A2*C0+B2*S0; P[10] = A2*S0-B2*C0; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE, NULL_STATE);

		P[0] = R0*C0; P[1] = -R0*S0; P[2] = 0.;
		P[3] = A1*C0; P[4] = -A1*S0; P[5] = 0.;
		P[6] = (t1+rad*T)*C0+B1*S0;  P[7] = -(t1+rad*T)*S0+B1*C0; P[8] = 0.;
		P[9] = A1*C0+B1*S0; P[10] = -A1*S0+B1*C0; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE, NULL_STATE);

		N_ini2 = N;
		if (fabs(radc) <= EE_ker) {
			P[0] = 0.; P[1] = 0.; P[2] = 0.;
			P[3] = (t1+rad*T)*C0+B1*S0; P[4] = -(t1+rad*T)*S0+B1*C0; P[5] = 0.;
			P[6] = (t2+rad*T)*C0+B2*S0; P[7] =  (t2+rad*T)*S0-B2*C0; P[8] = 0.;

			facet_QG(P, N_elem, NULL_STATE, NULL_STATE);
			QG_tria_circle(mp2, P, N_ini2);
		}
		else
		if (fabs(rad) <= EE_ker) {
			P[0] = t2*C0+B2*S0; P[1] = t2*S0-B2*C0; P[2] = 0.;
			P[3] = R0*C0; P[4] =  R0*S0; P[5] = 0.;
			P[6] = R0*C0; P[7] = -R0*S0; P[8] = 0.;
			
			facet_QG(P, N_elem, NULL_STATE, NULL_STATE);
			QG_tria_circle(mp1, P, N_ini2);
		}
		else {
			P[0] = R0*C0; P[1] =  R0*S0; P[2] = 0.;
			P[3] = R0*C0; P[4] = -R0*S0; P[5] = 0.;
			P[6] = (t2+rad*T)*C0+B2*S0;  P[7] =  (t2+rad*T)*S0-B2*C0;  P[8] = 0.;
			P[9] = (t1+rad*T)*C0+B1*S0; P[10] = -(t1+rad*T)*S0+B1*C0; P[11] = 0.;
	 
			facet_QG(P, N_elem, OK_STATE, NULL_STATE);
			QG_quad_circle(mp1, mp2, P, N_ini2);
		}
		P[3] = 0.;
		P[4] = 0.;
		P[5] = 1.;
		point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);

		for (int i = N_ini; i < N; i++) {
			P[0] = X[i]; 
			P[1] = Y[i]; 
			P[2] = Z[i]; 

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1]; nX[i] = P[3];
			Y[i] = P[1]+mp[2]; nY[i] = P[4];
			Z[i] = P[2]+mp[3]; nZ[i] = P[5];
		}
	}
}

/////////////////////////////////////////////////////////////////////
//...auxilliary functions for correction weigts and nodes for circle;
void CGrid_QG::QG_circle(CMap * mp)
{
	if (mp && ID_MAP(1, SPHERE_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && N_par > 0) {
		double R = fabs(mp[7]), R_inv = 1./R, P[6], r, t, S[6] = {0.,0.,0.,0.,0.,0.}, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);
		int i;
		for ( i = 0; i < N; i++) {
			P[0] = X[i]-mp[1];
			P[1] = Y[i]-mp[2];
			P[2] = Z[i]-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= R_inv;
			P[1] *= R_inv;
			P[2] = 0.;

			P[3] = nX[i]; S[0] += nX[i]*(r = get_param(0, i));
			P[4] = nY[i]; S[1] += nY[i]* r;
			P[5] = nZ[i]; S[2] += nZ[i]* r;
			point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[5] = 0.;

			r = P[0]*P[3]+P[1]*P[4];
			t = sqrt(fabs(1.+r*r-P[0]*P[0]-P[1]*P[1]))-fabs(r);
			if (r < 0.) {
				P[0] -= P[3]*t;
				P[1] -= P[4]*t;
			}
			else {
				P[0] += P[3]*t;
				P[1] += P[4]*t;
			}
			P[3] = P[0]*R;
			P[4] = P[1]*R;
			P[0] *= R;
			P[1] *= R;
			
			if ((r = sqrt(P[3]*P[3]+P[4]*P[4])) > EE) {
				P[3] *= (r = 1./r);
				P[4] *=  r;
			}
			else {
				P[3] = 1.;
				P[4] = 0.;
			}
			if (mp[7] < 0.) {
				P[3] = -P[3];
				P[4] = -P[4];
			}
			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];

			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			r = fabs(nX[i]*P[3]+nY[i]*P[4]+nZ[i]*P[5]);

			set_param(0, i, get_param(0, i)/r);
			nX[i] = P[3]; S[3] += nX[i]*(r = get_param(0, i));
			nY[i] = P[4]; S[4] += nY[i]* r;
			nZ[i] = P[5]; S[5] += nZ[i]* r;
		}

/////////////////////////////////////////////////////
//...��������� �������� (���������� � �� �� �������);
		t = abs_point(S, S+3); S[3] = -S[3]; S[4] = -S[4]; S[5] = -S[5];
		if ((r = abs_point(S, S+3)) < t)
		for (i = 0; i < N; i++) {
			nX[i] = -nX[i];
			nY[i] = -nY[i];
			nZ[i] = -nZ[i];
		}
	}
}

//////////////////////////////////////////////////////////////////////
//...auxilliary functions for correction weigts and nodes for ellipse;
void CGrid_QG::QG_ellipse(CMap * mp)
{
	if (mp && ID_MAP(1, CYL_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && fabs(mp[8]) > EE_ker && N_par > 0) {
		double A = fabs(mp[7]), B = fabs(mp[8]), A_inv = 1./A, B_inv = 1./B, P[6], r, t, S[6] = {0.,0.,0.,0.,0.,0.}, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);
		int i;
		for ( i = 0; i < N; i++) {
			P[0] = X[i]-mp[1];
			P[1] = Y[i]-mp[2];
			P[2] = Z[i]-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= A_inv;
			P[1] *= B_inv;
			P[2] = 0.;

			P[3] = nX[i]; S[0] += nX[i]*(r = get_param(0, i));
			P[4] = nY[i]; S[1] += nY[i]* r;
			P[5] = nZ[i]; S[2] += nZ[i]* r;
			point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[5] = 0.;

			r = P[0]*P[3]+P[1]*P[4];
			t = sqrt(fabs(1.+r*r-P[0]*P[0]-P[1]*P[1]))-fabs(r);
			if (r < 0.) {
				P[0] -= P[3]*t;
				P[1] -= P[4]*t;
			}
			else {
				P[0] += P[3]*t;
				P[1] += P[4]*t;
			}
			P[3] = P[0]*B;
			P[4] = P[1]*A;
			P[0] *= A;
			P[1] *= B;
			
			if ((r = sqrt(P[3]*P[3]+P[4]*P[4])) > EE) {
				P[3] *= (r = 1./r);
				P[4] *=  r;
			}
			else {
				P[3] = 1.;
				P[4] = 0.;
			}
			if (mp[7] < 0. || mp[8] < 0.) {
				P[3] = -P[3];
				P[4] = -P[4];
			}

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];

			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			r = fabs(nX[i]*P[3]+nY[i]*P[4]+nZ[i]*P[5]);

			set_param(0, i, get_param(0, i)/r);
			nX[i] = P[3]; S[3] += nX[i]*(r = get_param(0, i));
			nY[i] = P[4]; S[4] += nY[i]* r;
			nZ[i] = P[5]; S[5] += nZ[i]* r;
		}

/////////////////////////////////////////////////////
//...��������� �������� (���������� � �� �� �������);
		t = abs_point(S, S+3); S[3] = -S[3]; S[4] = -S[4]; S[5] = -S[5];
		if ((r = abs_point(S, S+3)) < t)
		for (i = 0; i < N; i++) {
			nX[i] = -nX[i];
			nY[i] = -nY[i];
			nZ[i] = -nZ[i];
		}
	}
}

void CGrid_QG::QG_tria_circle(CMap * mp, double * Po, int N_ini)
{
	if (mp && ID_MAP(1, SPHERE_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && N_par > 0) {
		double R = fabs(mp[7]), R_inv = 1./R, P[6], r, t, rr, RR, QQ, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]),
				tX = (Po[4]-Po[1])*(Po[8]-Po[2])-(Po[5]-Po[2])*(Po[7]-Po[1]), 
				tY = (Po[5]-Po[2])*(Po[6]-Po[0])-(Po[3]-Po[0])*(Po[8]-Po[2]), 
				tZ = (Po[3]-Po[0])*(Po[7]-Po[1])-(Po[4]-Po[1])*(Po[6]-Po[0]);
		QG_normal(tX, tY, tZ);

		for (int i = N_ini; i < N; i++) {
			rr = ((Po[3]-Po[0])*(Po[4]-Po[7])-(Po[4]-Po[1])*(Po[3]-Po[6]))/(nX[i]*(Po[4]-Po[7])-nY[i]*(Po[3]-Po[6]));

			P[0] = Po[0]+nX[i]*rr-mp[1];
			P[1] = Po[1]+nY[i]*rr-mp[2];
			P[2] = Po[2]+nZ[i]*rr-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= R_inv;
			P[1] *= R_inv;	P[2] = 0.;

			P[3] = nX[i];
			P[4] = nY[i];
			P[5] = nZ[i];
			point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);

			r = P[0]*P[3]+P[1]*P[4];
			t = sqrt(fabs(1.+r*r-P[0]*P[0]-P[1]*P[1]))-fabs(r);
			if (r < 0.) {
				P[0] -= P[3]*t;
				P[1] -= P[4]*t;
			}
			else {
				P[0] += P[3]*t;
				P[1] += P[4]*t;
			}
			P[0] *= R;
			P[1] *= R;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			P[0] += mp[1];
			P[1] += mp[2];
			P[2] += mp[3];

			RR = abs_point(P, Po);
			QQ = fabs(RR/rr);

			X[i] = Po[0]+(X[i]-Po[0])*QQ;
			Y[i] = Po[1]+(Y[i]-Po[1])*QQ;
			Z[i] = Po[2]+(Z[i]-Po[2])*QQ;
			nX[i] = tX;
			nY[i] = tY;
			nZ[i] = tZ;
			set_param(0, i, get_param(0, i)*sqr(QQ));
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////
//...auxilliary functions for correction weigts and nodes for curvilinear triangle;
//...�������: ���������� ����� ������� ��� ������, ��� ��� ����������;
void CGrid_QG::QG_tria_ellipse(CMap * mp, double * Po, int N_ini)
{
	if (mp && ID_MAP(1, CYL_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && fabs(mp[8]) > EE_ker && N_par > 0) {
		double A = fabs(mp[7]), B = fabs(mp[8]), A_inv = 1./A, B_inv = 1./B, P[6], r, t, rr, RR, QQ, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]),
				tX = (Po[4]-Po[1])*(Po[8]-Po[2])-(Po[5]-Po[2])*(Po[7]-Po[1]), 
				tY = (Po[5]-Po[2])*(Po[6]-Po[0])-(Po[3]-Po[0])*(Po[8]-Po[2]), 
				tZ = (Po[3]-Po[0])*(Po[7]-Po[1])-(Po[4]-Po[1])*(Po[6]-Po[0]);
		QG_normal(tX, tY, tZ);

		for (int i = N_ini; i < N; i++) {
			rr = ((Po[3]-Po[0])*(Po[4]-Po[7])-(Po[4]-Po[1])*(Po[3]-Po[6]))/(nX[i]*(Po[4]-Po[7])-nY[i]*(Po[3]-Po[6]));

			P[0] = Po[0]+nX[i]*rr-mp[1];
			P[1] = Po[1]+nY[i]*rr-mp[2];
			P[2] = Z[i]-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= A_inv;
			P[1] *= B_inv;
			P[2] = 0.;

			P[3] = nX[i];
			P[4] = nY[i];
			P[5] = nZ[i];
			point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[5] = 0.;

			P[3] = nX[i]*A_inv;
			P[4] = nY[i]*B_inv;

			if ((r = sqrt(P[3]*P[3]+P[4]*P[4])) > EE) {
				P[3] *= (r = 1./r);
				P[4] *=  r;
			}
			else {
				P[3] = 1.;
				P[4] = 0.;
			}

			r = P[0]*P[3]+P[1]*P[4];
			t = sqrt(fabs(1.+r*r-P[0]*P[0]-P[1]*P[1]))-fabs(r);
			if (r < 0.) {
				P[0] -= P[3]*t;
				P[1] -= P[4]*t;
			}
			else {
				P[0] += P[3]*t;
				P[1] += P[4]*t;
			}
			P[0] *= A;
			P[1] *= B;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			P[0] += mp[1];
			P[1] += mp[2];
			P[2] += mp[3];

			RR = sqrt(sqr(P[0]-Po[0])+sqr(P[1]-Po[1]));
			QQ = fabs(RR/rr);

			X[i] = Po[0]+(X[i]-Po[0])*QQ;
			Y[i] = Po[1]+(Y[i]-Po[1])*QQ;
			nX[i] = tX;
			nY[i] = tY;
			nZ[i] = tZ;
			set_param(0, i, get_param(0, i)*sqr(QQ));
		}
	}
}

void CGrid_QG::QG_quad_circle(CMap * mp1, CMap * mp2, double * Po, int N_ini)
{
	if ((! mp1 || ID_MAP(1, SPHERE_GENUS) == mp1[0] && fabs(mp1[7]) > EE_ker) &&
		 (! mp2 || ID_MAP(1, SPHERE_GENUS) == mp2[0] && fabs(mp2[7]) > EE_ker) && N_par > 0) {
		double R_inv, P[18], r, t, rr, RR, QQ, CZ1, SZ1, CY1, SY1, CX1, SX1, CZ2, SZ2, CY2, SY2, CX2, SX2, 
				tX, tY, tZ, F1, F2, F3, F4, tt, xx, GG;
		if (mp1) {
			CZ1 = cos(mp1[4]); SZ1 = sin(mp1[4]); 
			CY1 = cos(mp1[5]); SY1 = sin(mp1[5]);
			CX1 = cos(mp1[6]); SX1 = sin(mp1[6]);
		}
		if (mp2) {
			CZ2 = cos(mp2[4]); SZ2 = sin(mp2[4]); 
			CY2 = cos(mp2[5]); SY2 = sin(mp2[5]);
			CX2 = cos(mp2[6]); SX2 = sin(mp2[6]);
		}
		tX = (Po[10]-Po[4])*(Po[11]-Po[8])-(Po[11]-Po[5])*(Po[10]-Po[7]); 
		tY = (Po[11]-Po[5])*(Po[9] -Po[6])-(Po[9] -Po[3])*(Po[11]-Po[8]); 
		tZ = (Po[9] -Po[3])*(Po[10]-Po[7])-(Po[10]-Po[4])*(Po[9] -Po[6]);
		F4 = QG_normal(tX, tY, tZ)*.0625;

		tX = (Po[7]-Po[1])*(Po[11]-Po[8])-(Po[8]-Po[2])*(Po[10]-Po[7]); 
		tY = (Po[8]-Po[2])*(Po[9] -Po[6])-(Po[6]-Po[0])*(Po[11]-Po[8]); 
		tZ = (Po[6]-Po[0])*(Po[10]-Po[7])-(Po[7]-Po[1])*(Po[9] -Po[6]);
		F3 = QG_normal(tX, tY, tZ)*.0625;

		tX = (Po[10]-Po[4])*(Po[5]-Po[2])-(Po[11]-Po[5])*(Po[4]-Po[1]); 
		tY = (Po[11]-Po[5])*(Po[3]-Po[0])-(Po[9] -Po[3])*(Po[5]-Po[2]); 
		tZ = (Po[9] -Po[3])*(Po[4]-Po[1])-(Po[10]-Po[4])*(Po[3]-Po[0]);
		F2 = QG_normal(tX, tY, tZ)*.0625;

		tX = (Po[7]-Po[1])*(Po[5]-Po[2])-(Po[8]-Po[2])*(Po[4]-Po[1]); 
		tY = (Po[8]-Po[2])*(Po[3]-Po[0])-(Po[6]-Po[0])*(Po[5]-Po[2]); 
		tZ = (Po[6]-Po[0])*(Po[4]-Po[1])-(Po[7]-Po[1])*(Po[3]-Po[0]);
		F1 = QG_normal(tX, tY, tZ)*.0625;

		for (int i = N_ini; i < N; i++) {
			rr = ((Po[0]-X[i])*(Po[1]-Po[4])-(Po[1]-Y[i])*(Po[0]-Po[3]))/(nX[i]*(Po[1]-Po[4])-nY[i]*(Po[0]-Po[3]));
			P[9]  = -nX[i]*rr;
			P[10] = -nY[i]*rr;
			P[11] = -nZ[i]*rr;
			P[0] = P[12] = X[i]+nX[i]*rr;
			P[1] = P[13] = Y[i]+nY[i]*rr;
			P[2] = P[14] = Z[i]+nZ[i]*rr;

			rr = ((Po[6]-X[i])*(Po[7]-Po[10])-(Po[7]-Y[i])*(Po[6]-Po[9]))/(nX[i]*(Po[7]-Po[10])-nY[i]*(Po[6]-Po[9]));
			P[3] = P[15] = X[i]+nX[i]*rr;
			P[4] = P[16] = Y[i]+nY[i]*rr;
			P[5] = P[17] = Z[i]+nZ[i]*rr;

			xx = ((P[0]-Po[0])*(Po[3]-Po[0])+(P[1]-Po[1])*(Po[4]-Po[1])+(P[2]-Po[2])*(Po[5]-Po[2]))/((sqr(Po[3]-Po[0])+sqr(Po[4]-Po[1])+sqr(Po[5]-Po[2])));
			tt = ((X[i]-P[0])*(P[3]-P[0])+(Y[i]-P[1])*(P[4]-P[1])+(Z[i]-P[2])*(P[5]-P[2]))*(rr = 1./(sqr(P[3]-P[0])+sqr(P[4]-P[1])+sqr(P[5]-P[2])));
			
			if (mp1) {
				R_inv = 1./(RR = fabs(mp1[7])); 
				P[0] -= mp1[1];
				P[1] -= mp1[2];
				P[2] -= mp1[3];
				point_iso<double>(P, NULL, CX1, -SX1, CY1, -SY1, CZ1, -SZ1);
				P[0] *= R_inv;
				P[1] *= R_inv;	P[2] = 0.;

				P[6] = nX[i];
				P[7] = nY[i];
				P[8] = nZ[i];
				point_iso<double>(P+6, NULL, CX1, -SX1, CY1, -SY1, CZ1, -SZ1);

				r = P[0]*P[6]+P[1]*P[7];
				t = sqrt(fabs(1.+r*r-P[0]*P[0]-P[1]*P[1]))-fabs(r);
				if (r < 0.) {
					P[0] -= P[6]*t;
					P[1] -= P[7]*t;
				}
				else {
					P[0] += P[6]*t;
					P[1] += P[7]*t;
				}
				P[0] *= RR;
				P[1] *= RR;

				point_iso<double>(P, NULL, CZ1, SZ1, CY1, SY1, CX1, SX1);
				P[0] += mp1[1];
				P[1] += mp1[2];
				P[2] += mp1[3];
			}
			if (mp2) {
				R_inv = 1./(RR = fabs(mp2[7])); 
				P[3] -= mp2[1];
				P[4] -= mp2[2];
				P[5] -= mp2[3];
				point_iso<double>(P+3, NULL, CX2, -SX2, CY2, -SY2, CZ2, -SZ2);
				P[3] *= R_inv;
				P[4] *= R_inv;	P[5] = 0.;

				P[6] = nX[i];
				P[7] = nY[i];
				P[8] = nZ[i];
				point_iso<double>(P+6, NULL, CX2, -SX2, CY2, -SY2, CZ2, -SZ2);

				r = P[3]*P[6]+P[4]*P[7];
				t = sqrt(fabs(1.+r*r-P[3]*P[3]-P[4]*P[4]))-fabs(r);
				if (r < 0.) {
					P[3] -= P[6]*t;
					P[4] -= P[7]*t;
				}
				else {
					P[3] += P[6]*t;
					P[4] += P[7]*t;
				}
				P[3] *= RR;
				P[4] *= RR;

				point_iso<double>(P, NULL, CZ2, SZ2, CY2, SY2, CX2, SX2);
				P[3] += mp2[1];
				P[4] += mp2[2];
				P[5] += mp2[3];
			}
			RR = sqr(P[3]-P[0])+sqr(P[4]-P[1])+sqr(P[5]-P[2]);
			QQ = sqrt(RR*rr);

			X[i] = P[0]+P[9]*QQ;
			Y[i] = P[1]+P[10]*QQ;
			Z[i] = P[2]+P[11]*QQ;

			nX[i] = tX;
			nY[i] = tY;
			nZ[i] = tZ;

			rr *= ((X[i]-P[12])*(P[15]-P[12])+(Y[i]-P[13])*(P[16]-P[13])+(Z[i]-P[14])*(P[17]-P[14]));
			GG =((1.-xx)*(F1*(1.-rr)+F3*rr)+xx*(F2*(1.-rr)+F4*rr))/
				 ((1.-xx)*(F1*(1.-tt)+F3*tt)+xx*(F2*(1.-tt)+F4*tt));
			set_param(0, i, get_param(0, i)*QQ*GG);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//...auxilliary functions for correction weigts and nodes for curvilinear quadrangle;
void CGrid_QG::QG_quad_ellipse(CMap * mp1, CMap * mp2, double * Po, int N_ini)
{
	if ((! mp1 || ID_MAP(1, CYL_GENUS) == mp1[0] && fabs(mp1[7]) > EE_ker && fabs(mp1[8]) > EE_ker) &&
		 (! mp2 || ID_MAP(1, CYL_GENUS) == mp2[0] && fabs(mp2[7]) > EE_ker && fabs(mp2[8]) > EE_ker) && N_par > 0) {
		double A, B, A_inv, B_inv, P[18], r, t, rr, RR, QQ, CZ, SZ, CY, SY, CX, SX, 
				tX, tY, tZ, F1, F2, F3, F4, tt, xx, GG;

		tX = (Po[10]-Po[4])*(Po[11]-Po[8])-(Po[11]-Po[5])*(Po[10]-Po[7]); 
		tY = (Po[11]-Po[5])*(Po[9] -Po[6])-(Po[9] -Po[3])*(Po[11]-Po[8]); 
		tZ = (Po[9] -Po[3])*(Po[10]-Po[7])-(Po[10]-Po[4])*(Po[9] -Po[6]);
		F4 = QG_normal(tX, tY, tZ)*.0625;

		tX = (Po[7]-Po[1])*(Po[11]-Po[8])-(Po[8]-Po[2])*(Po[10]-Po[7]); 
		tY = (Po[8]-Po[2])*(Po[9] -Po[6])-(Po[6]-Po[0])*(Po[11]-Po[8]); 
		tZ = (Po[6]-Po[0])*(Po[10]-Po[7])-(Po[7]-Po[1])*(Po[9] -Po[6]);
		F3 = QG_normal(tX, tY, tZ)*.0625;

		tX = (Po[10]-Po[4])*(Po[5]-Po[2])-(Po[11]-Po[5])*(Po[4]-Po[1]); 
		tY = (Po[11]-Po[5])*(Po[3]-Po[0])-(Po[9] -Po[3])*(Po[5]-Po[2]); 
		tZ = (Po[9] -Po[3])*(Po[4]-Po[1])-(Po[10]-Po[4])*(Po[3]-Po[0]);
		F2 = QG_normal(tX, tY, tZ)*.0625;

		tX = (Po[7]-Po[1])*(Po[5]-Po[2])-(Po[8]-Po[2])*(Po[4]-Po[1]); 
		tY = (Po[8]-Po[2])*(Po[3]-Po[0])-(Po[6]-Po[0])*(Po[5]-Po[2]); 
		tZ = (Po[6]-Po[0])*(Po[4]-Po[1])-(Po[7]-Po[1])*(Po[3]-Po[0]);
		F1 = QG_normal(tX, tY, tZ)*.0625;

		for (int i = N_ini; i < N; i++) {
			rr = ((Po[0]-X[i])*(Po[1]-Po[4])-(Po[1]-Y[i])*(Po[0]-Po[3]))/(nX[i]*(Po[1]-Po[4])-nY[i]*(Po[0]-Po[3]));
			P[0] = P[12] = X[i]+nX[i]*rr;
			P[1] = P[13] = Y[i]+nY[i]*rr;
			P[2] = P[14] = Z[i];
			P[9] = -nX[i]*rr;
			P[10] = -nY[i]*rr;

			rr = ((Po[6]-X[i])*(Po[7]-Po[10])-(Po[7]-Y[i])*(Po[6]-Po[9]))/(nX[i]*(Po[7]-Po[10])-nY[i]*(Po[6]-Po[9]));
			P[3] = P[15] = X[i]+nX[i]*rr;
			P[4] = P[16] = Y[i]+nY[i]*rr;
			P[5] = P[17] = Z[i];

			xx = ((P[0]-Po[0])*(Po[3]-Po[0])+(P[1]-Po[1])*(Po[4]-Po[1]))/((sqr(Po[3]-Po[0])+sqr(Po[4]-Po[1])));
			tt = ((X[i]-P[0])*(P[3]-P[0])+(Y[i]-P[1])*(P[4]-P[1]))*(rr = 1./(sqr(P[3]-P[0])+sqr(P[4]-P[1])));
			if (mp1) {
				A_inv = 1./(A = fabs(mp1[7])); 
				B_inv = 1./(B = fabs(mp1[8])); 
				CZ = cos(mp1[4]); SZ = sin(mp1[4]); 
				CY = cos(mp1[5]); SY = sin(mp1[5]);
				CX = cos(mp1[6]); SX = sin(mp1[6]);

				P[0] -= mp1[1];
				P[1] -= mp1[2];
				P[2] -= mp1[3];
				point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
				P[0] *= A_inv;
				P[1] *= B_inv;
				P[2] = 0.;

				P[6] = nX[i];
				P[7] = nY[i];
				P[8] = nZ[i];
				point_iso<double>(P+6, NULL, CX, -SX, CY, -SY, CZ, -SZ);
				P[8] = 0.;

				P[6] = /*nX[i]*/P[6]*A_inv;
				P[7] = /*nY[i]*/P[7]*B_inv;

				if ((r = sqrt(P[6]*P[6]+P[7]*P[7])) > EE) {
					P[6] *= (r = 1./r);
					P[7] *=  r;
				}
				else {
					P[6] = 1.;
					P[7] = 0.;
				}

				r = P[0]*P[6]+P[1]*P[7];
				t = sqrt(fabs(1.+r*r-P[0]*P[0]-P[1]*P[1]))-fabs(r);
				if (r < 0.) {
					P[0] -= P[6]*t;
					P[1] -= P[7]*t;
				}
				else {
					P[0] += P[6]*t;
					P[1] += P[7]*t;
				}
				P[0] *= A;
				P[1] *= B;

				point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
				P[0] += mp1[1];
				P[1] += mp1[2];
				P[2] += mp1[3];
			}
			if (mp2) {
				A_inv = 1./(A = fabs(mp2[7])); 
				B_inv = 1./(B = fabs(mp2[8])); 
				CZ = cos(mp2[4]); SZ = sin(mp2[4]); 
				CY = cos(mp2[5]); SY = sin(mp2[5]);
				CX = cos(mp2[6]); SX = sin(mp2[6]);

				P[3] -= mp2[1];
				P[4] -= mp2[2];
				P[5] -= mp2[3];
				point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);
				P[3] *= A_inv;
				P[4] *= B_inv;
				P[5] = 0.;

				P[6] = nX[i];
				P[7] = nY[i];
				P[8] = nZ[i];
				point_iso<double>(P+6, NULL, CX, -SX, CY, -SY, CZ, -SZ);
				P[8] = 0.;

				P[6] = /*nX[i]*/P[6]*A_inv;
				P[7] = /*nY[i]*/P[7]*B_inv;

				if ((r = sqrt(P[6]*P[6]+P[7]*P[7])) > EE) {
					P[6] *= (r = 1./r);
					P[7] *=  r;
				}
				else {
					P[6] = 1.;
					P[7] = 0.;
				}

				r = P[3]*P[6]+P[4]*P[7];
				t = sqrt(fabs(1.+r*r-P[3]*P[3]-P[4]*P[4]))-fabs(r);
				if (r < 0.) {
					P[3] -= P[6]*t;
					P[4] -= P[7]*t;
				}
				else {
					P[3] += P[6]*t;
					P[4] += P[7]*t;
				}
				P[3] *= A;
				P[4] *= B;

				point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
				P[3] += mp2[1];
				P[4] += mp2[2];
				P[5] += mp2[3];
			}
			RR = (sqr(P[3]-P[0])+sqr(P[4]-P[1]));
			QQ = sqrt(RR*rr);

			X[i] = P[0]+P[9]*QQ;
			Y[i] = P[1]+P[10]*QQ;

			nX[i] = tX;
			nY[i] = tY;
			nZ[i] = tZ;

			rr *= ((X[i]-P[12])*(P[15]-P[12])+(Y[i]-P[13])*(P[16]-P[13]));
			GG =((1.-xx)*(F1*(1.-rr)+F3*rr)+xx*(F2*(1.-rr)+F4*rr))/
				 ((1.-xx)*(F1*(1.-tt)+F3*tt)+xx*(F2*(1.-tt)+F4*tt));

			set_param(0, i, get_param(0, i)*QQ*GG);
		}
	}
}

////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for sphere;
void CGrid_QG::QG_sphere(CMap * mp)
{
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && N_par > 0) {
		double R = fabs(mp[7]), P[6], r, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);
		for (int i = 0; i < N; i++) {
			P[0] = X[i]-mp[1];
			P[1] = Y[i]-mp[2];
			P[2] = Z[i]-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);

			P[3] = nX[i];
			P[4] = nY[i];
			P[5] = nZ[i];
			point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			
			if (( r = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2])) > EE) {
				P[3] = (P[0] *= (r = 1./r));
				P[4] = (P[1] *=  r);
				P[5] = (P[2] *=  r);
			}
			else {
				P[0] = P[3];
				P[1] = P[4];
				P[2] = P[5];
			}
			P[0] *= R;
			P[1] *= R;
			P[2] *= R;
			if (mp[7] < 0.) {
				P[3] = -P[3];
				P[4] = -P[4];
				P[5] = -P[5];
			}
			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];

			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			r = fabs(nX[i]*P[3]+nY[i]*P[4]+nZ[i]*P[5]);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)/r);
		}
	}
}

//////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for cylinder;
void CGrid_QG::QG_cylinder(CMap * mp)
{
	if (mp && ID_MAP(2, CYL_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && N_par > 0) {
		double R = fabs(mp[7]), P[6], r, 
				CZ = cos (mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos (mp[6]), SX = sin(mp[6]);
		for (int i = 0; i < N; i++) {
			P[0] = X[i]-mp[1];
			P[1] = Y[i]-mp[2];
			P[2] = Z[i]-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[3] = nX[i];
			P[4] = nY[i];
			P[5] = nZ[i];
			point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			
			P[3] = (P[0] *= (r = 1./sqrt(P[0]*P[0]+P[1]*P[1])));
			P[4] = (P[1] *= r);
			P[5] = 0.;
			P[0] *= R;
			P[1] *= R;
			if (mp[7] < 0.) {
				P[3] = -P[3];
				P[4] = -P[4];
			}
			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];
			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			swap(nX[i], P[3]);
			swap(nY[i], P[4]);
			swap(nZ[i], P[5]);
			set_param(0, i, get_param(0, i)/(r = fabs(nX[i]*P[3]+nY[i]*P[4]+nZ[i]*P[5])));
		}
	}
}

//////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for spheroid;
void CGrid_QG::QG_spheroid(CMap * mp)
{
	if (mp && ID_MAP(2, SPHEROID_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && fabs(mp[8]) > EE_ker && N_par > 0) {
		double A = fabs(mp[7]), B = fabs(mp[8]), A_inv = 1./A, B_inv = 1./B, P[6], r, t, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]);
		for (int i = 0; i < N; i++) {
			P[0] = X[i]-mp[1];
			P[1] = Y[i]-mp[2];
			P[2] = Z[i]-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= B_inv;
			P[1] *= B_inv;
			P[2] *= A_inv;

			P[3] = nX[i];
			P[4] = nY[i];
			P[5] = nZ[i];
			point_iso<double>(P+3, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			
			r = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
			t = sqrt(fabs(1.+r*r-P[0]*P[0]-P[1]*P[1]-P[2]*P[2]))-fabs(r);
			if (r < 0.) {
				P[0] -= P[3]*t;
				P[1] -= P[4]*t;
				P[2] -= P[5]*t;
			}
			else {
				P[0] += P[3]*t;
				P[1] += P[4]*t;
				P[2] += P[5]*t;
			}
			P[3] = (P[0] *= B);
			P[4] = (P[1] *= B);
			P[5] = (P[2] *= A);

			if ((r = sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5])) > EE) {
				P[3] *= (r = 1./r);
				P[4] *=  r;
				P[5] *=  r;
			}
			else {
				P[3] = 0.;
				P[4] = 0.;
				P[5] = 1.;
			}
			if (mp[7] < 0. || mp[8] < 0.) {
				P[3] = -P[3];
				P[4] = -P[4];
				P[5] = -P[5];
			}
			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];

			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			r = fabs(nX[i]*P[3]+nY[i]*P[4]+nZ[i]*P[5]);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];
			set_param(0, i, get_param(0, i)/r);
		}
	}
}

////////////////////////////////////////////////////////
//...��������� ������� ����������� ������ � �����������;
void CGrid_QG::QG_tria_sphere(CMap * mp, double * Po, int N_ini)
{
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && fabs(mp[7]) > EE_ker && N_par > 0) {
		double R = fabs(mp[7]), P[6], pp[6], r, t, rr, RR, QQ;
		for (int i = N_ini; i < N; i++) {
			pp[0] = (Po[4]-Po[1])*(Po[5]-Po[8])-(Po[5]-Po[2])*(Po[4]-Po[7]), 
			pp[1] = (Po[5]-Po[2])*(Po[3]-Po[6])-(Po[3]-Po[0])*(Po[5]-Po[8]);
			pp[2] = (Po[3]-Po[0])*(Po[4]-Po[7])-(Po[4]-Po[1])*(Po[3]-Po[6]);
			pp[3] = nY[i]*(Po[5]-Po[8])-nZ[i]*(Po[4]-Po[7]);
			pp[4] = nZ[i]*(Po[3]-Po[6])-nX[i]*(Po[5]-Po[8]);
			pp[5] = nX[i]*(Po[4]-Po[7])-nY[i]*(Po[3]-Po[6]);

			rr = sqrt((sqr(pp[0])+sqr(pp[1])+sqr(pp[2]))/(sqr(pp[3])+sqr(pp[4])+sqr(pp[5])));

			P[0] = Po[0]+(P[3] = nX[i])*rr-mp[1];
			P[1] = Po[1]+(P[4] = nY[i])*rr-mp[2];
			P[2] = Po[2]+(P[5] = nZ[i])*rr-mp[3];

			r = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
			t = sqrt(fabs(R*R-P[0]*P[0]-P[1]*P[1]-P[2]*P[2]+r*r))-fabs(r);

			if (r < 0.) { //...����������� ����� �� ��������� ����������� �����;
				P[0] -= P[3]*t;
				P[1] -= P[4]*t;
				P[2] -= P[5]*t;
			}
			else {
				P[0] += P[3]*t;
				P[1] += P[4]*t;
				P[2] += P[5]*t;
			}
			P[0] += mp[1];
			P[1] += mp[2];
			P[2] += mp[3];

			RR = sqrt(sqr(P[0]-Po[0])+sqr(P[1]-Po[1])+sqr(P[2]-Po[2]));
			QQ = fabs(RR/rr);

			X[i] = Po[0]+(X[i]-Po[0])*QQ;
			Y[i] = Po[1]+(Y[i]-Po[1])*QQ;
			Z[i] = Po[2]+(Z[i]-Po[2])*QQ;
			set_param(0, i, get_param(0, i)*sqr(QQ));
		}
	}
}

////////////////////////////////////////////////////////////
//...��������� ������� ��������������� ������ � �����������;
void CGrid_QG::QG_quad_sphere(CMap * mp1, CMap * mp2, double * Po, int N_ini)
{
	if ((! mp1 || ID_MAP(2, SPHERE_GENUS) == mp1[0] && fabs(mp1[7]) > EE_ker) &&
		 (! mp2 || ID_MAP(2, SPHERE_GENUS) == mp2[0] && fabs(mp2[7]) > EE_ker) && N_par > 0) {
		double P[12], pp[6], r, t, rr, RR, QQ;
		for (int i = N_ini; i < N; i++) {
			pp[0] = (Po[1]-Y[i])*(Po[5]-Po[2])-(Po[2]-Z[i])*(Po[4]-Po[1]), 
			pp[1] = (Po[2]-Z[i])*(Po[3]-Po[0])-(Po[0]-X[i])*(Po[5]-Po[2]);
			pp[2] = (Po[0]-X[i])*(Po[4]-Po[1])-(Po[1]-Y[i])*(Po[3]-Po[0]);
			pp[3] = nY[i]*(Po[5]-Po[2])-nZ[i]*(Po[4]-Po[1]), 
			pp[4] = nZ[i]*(Po[3]-Po[0])-nX[i]*(Po[5]-Po[2]);
			pp[5] = nX[i]*(Po[4]-Po[1])-nY[i]*(Po[3]-Po[0]);
			RR = sqrt((sqr(pp[0])+sqr(pp[1])+sqr(pp[2]))/(sqr(pp[3])+sqr(pp[4])+sqr(pp[5])));

			P[0] = X[i]-(P[9]  = (P[3] = nX[i])*RR);
			P[1] = Y[i]-(P[10] = (P[4] = nY[i])*RR);
			P[2] = Z[i]-(P[11] = (P[5] = nZ[i])*RR);

			if (mp1) {
				P[0] -= mp1[1];
				P[1] -= mp1[2];
				P[2] -= mp1[3];

				r = P[0]*P[3]+P[1]*P[4]+P[2]*P[5];
				t = sqrt(fabs(mp1[7]*mp1[7]-P[0]*P[0]-P[1]*P[1]-P[2]*P[2]+r*r))-fabs(r);

				if (r < 0.) { //...����������� ����� �� ��������� ����������� �����;
					P[0] -= P[3]*t;
					P[1] -= P[4]*t;
					P[2] -= P[5]*t;
				}
				else {
					P[0] += P[3]*t;
					P[1] += P[4]*t;
					P[2] += P[5]*t;
				}
				P[0] += mp1[1];
				P[1] += mp1[2];
				P[2] += mp1[3];
			}
			pp[0] = (Po[7]-Y[i])*(Po[11]-Po[8])-(Po[8]-Z[i])*(Po[10]-Po[7]), 
			pp[1] = (Po[8]-Z[i])*(Po[9] -Po[6])-(Po[6]-X[i])*(Po[11]-Po[8]);
			pp[2] = (Po[6]-X[i])*(Po[10]-Po[7])-(Po[7]-Y[i])*(Po[9] -Po[6]);
			pp[3] = nY[i]*(Po[11]-Po[8])-nZ[i]*(Po[10]-Po[7]), 
			pp[4] = nZ[i]*(Po[9] -Po[6])-nX[i]*(Po[11]-Po[8]);
			pp[5] = nX[i]*(Po[10]-Po[7])-nY[i]*(Po[9] -Po[6]);
			rr = sqrt((sqr(pp[0])+sqr(pp[1])+sqr(pp[2]))/(sqr(pp[3])+sqr(pp[4])+sqr(pp[5])));

			P[6] = X[i]+nX[i]*rr;
			P[7] = Y[i]+nY[i]*rr;
			P[8] = Z[i]+nZ[i]*rr;
			rr += RR;

			if (mp2) {
				P[6] -= mp2[1];
				P[7] -= mp2[2];
				P[8] -= mp2[3];

				r = P[6]*P[3]+P[7]*P[4]+P[8]*P[5];
				t = sqrt(fabs(mp2[7]*mp2[7]-P[6]*P[6]-P[7]*P[7]-P[8]*P[8]+r*r))-fabs(r);

				if (r < 0.) { //...����������� ����� �� ��������� ����������� �����;
					P[6] -= P[3]*t;
					P[7] -= P[4]*t;
					P[8] -= P[5]*t;
				}
				else {
					P[6] += P[3]*t;
					P[7] += P[4]*t;
					P[8] += P[5]*t;
				}
				P[6] += mp2[1];
				P[7] += mp2[2];
				P[8] += mp2[3];
			}
			RR = sqrt(sqr(P[6]-P[0])+sqr(P[7]-P[1])+sqr(P[8]-P[2]));
			QQ = fabs(RR/rr);

			X[i] = P[0]+P[9]*QQ;
			Y[i] = P[1]+P[10]*QQ;
			Z[i] = P[2]+P[11]*QQ;
			set_param(0, i, get_param(0, i)*QQ);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
//...���������� ������ ��� ������� ������ (����������� � ���������������);
void CGrid_QG::facet_QG(double * P, int N_elem, Num_State mode, Num_State normal)
{
	int shift = 0;
	if (! P) return;
	if (N_elem < 2)  N_elem = 1; else
	if (N_elem < 25) N_elem = min(N_elem,  20); else
	if (N_elem < 30) shift  = 85+(N_elem = 25); else
	if (N_elem < 35) shift  = 93+(N_elem = 30); else
	if (N_elem < 40) shift  =103+(N_elem = 35); else shift = 116+(N_elem = 40);
	if (N_elem % 2)
	switch (mode) {
		case    SPECIAL_STATE: QG2M1_line (P, (N_elem+1)/2, cc+(N_elem < 25 ? (N_elem*N_elem-1)/4 : shift), normal); break;
		case         OK_STATE: QG2M1_quad (P, (N_elem+1)/2, cc+(N_elem < 25 ? (N_elem*N_elem-1)/4 : shift), normal); break;
		default:               QG2M1_tria (P, (N_elem+1)/2, cc+(N_elem < 25 ? (N_elem*N_elem-1)/4 : shift), normal);
	}
	else
	switch (mode) {
		case    SPECIAL_STATE: QG2M_line(P, N_elem/2, cc+(N_elem < 25 ? N_elem*N_elem/4 : shift), normal); break;
		case         OK_STATE: QG2M_quad(P, N_elem/2, cc+(N_elem < 25 ? N_elem*N_elem/4 : shift), normal); break;
		default:               QG2M_tria(P, N_elem/2, cc+(N_elem < 25 ? N_elem*N_elem/4 : shift), normal);
	}
}


//////////////////////////////////////////////////////////////////////
//...auxilliary functions for weigts and nodes of dimensional elements;
void CGrid_QG::segms_QG(CMap * mp, int N_elem, int N_max)
{
	if (2 == map_dim(mp)) { //...����� ��������� ��� ��������� ���������;
		QG_sheet			 (mp, N_elem, N_max); //...�������;
		QG_ring_segment (mp, N_elem, N_max); //...�������;
		QG_cyl_segment  (mp, N_elem, N_max); //...�������;
		QG_sph_segment  (mp, N_elem, N_max); //...�������;
		QG_spr_segment  (mp, N_elem, N_max); //...�������;
		QG_sph_intrusion(mp, N_elem, N_max); //...�������;
		QG_sht_intrusion(mp, N_elem, N_max); //...�������;
		QG_cone_segment (mp, N_elem, N_max); //...�������;
		QG_torus_segment(mp, N_elem, N_max); //...�������;
		QG_ugolok_cell  (mp, N_elem); //...�������;
	}
}

///////////////////////////////////////////////////////
//...special quadrature for sphere penetration surface;
void CGrid_QG::sphere_intrusion_QG(CMap * mp, int N_elem, double L)
{
	int i, k, N_ini = N;
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && mp[k = size_of_map(2, SPHERE_GENUS)] == (CMap)SPH_SEGMENT
			 && fabs(mp[k+1]) > EE_ker) {
		double A = mp[k+1], B = mp[k+2], C = mp[k+3], RR = fabs(mp[7]), P[12], f = (mp[7] < 0. ? -1. : 1.)*(C < B ? -1. : 1.),
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), g, alpha, theta, fB, fC, fL, fY;

		P[0] = -A*.5; P[1]  = B; P[2] = 0.;
		P[3] = -A*.5; P[4]  = C; P[5] = 0.;
		P[6] =  A*.5; P[7]  = B; P[8] = 0.;
		P[9] =  A*.5; P[10] = C; P[11] = 0.;

		facet_QG(P, N_elem, OK_STATE);

//////////////////////////////////////////////////
//...��������� ���������� ��� ������� �����������;
		if (L < RR && RR > 0. && fabs(C-B) > 0.) {
			alpha = L/RR;
			theta = acos(alpha);
			fB = min(B, C);
			fC = max(B, C); B = fB; C = fC;
			for ( i = N_ini; i < N; i++) {
				fB = B;
				fC = C; 
				fY = 0.;
				fL = 1.;
				for (k = -4; k <= 4; k++)
				if (fabs(X[i]+mp[6]-k*M_PI_2) < theta) {
					fY = acos(alpha/cos(X[i]+mp[6]-k*M_PI_2)); break;
				}
				if ( Y[i] <= M_PI_2) {
					if (fB <  theta) fB = theta;
					if (fC <  theta) fB = fC = (hit[i] = 0)+theta;
					if (fB >= M_PI_2-fY) fB = fC = (hit[i] = 0)+M_PI_2-fY; 
					if (fC <= M_PI_2-fY) Y[i] = fB+(Y[i]-B)*(fL = (fC-fB)/(C-B)); else
					if (fC <= M_PI_2) Y[i] = fB+(Y[i]-B)*(fL = (M_PI_2-fY-fB)/(C-B)); else
											Y[i] = fB+(Y[i]-B)*(fL = (M_PI_2-fY-fB)/(M_PI_2-B));
				}
				else {
					if (fC >  M_PI-theta) fC = M_PI-theta;
					if (fB >  M_PI-theta) fC = fB = (hit[i] = 0)+M_PI-theta;
					if (fC <= M_PI_2+fY)  fC = fB = (hit[i] = 0)+M_PI_2+fY;
					if (fB >= M_PI_2+fY) Y[i] = fC-(C-Y[i])*(fL = (fC-fB)/(C-B)); else
					if (fB >= M_PI_2) Y[i] = fC-(C-Y[i])*(fL = (fC-M_PI_2-fY)/(C-B)); else
											Y[i] = fC-(C-Y[i])*(fL = (fC-M_PI_2-fY)/(C-M_PI_2));
				}
				set_param(0, i, get_param(0, i)*fL);
			}
		}

/////////////////////////////////////////////////////////////
//...������� ����� �������������� � ���������������� �������;
		for ( i = N_ini; i < N; i++) {
			P[3] = cos(X[i])*(g = sin(Y[i])*f);	P[0] = RR*P[3]*f;
			P[4] = sin(X[i])* g;						P[1] = RR*P[4]*f;
			P[5] = cos(Y[i])* f;						P[2] = RR*P[5]*f;

			point_iso<double>(P, NULL, CZ, SZ, CY, SY, CX, SX);
			X[i] = P[0]+mp[1];
			Y[i] = P[1]+mp[2];
			Z[i] = P[2]+mp[3];

			point_iso<double>(P+3, NULL, CZ, SZ, CY, SY, CX, SX);
			nX[i] = P[3];
			nY[i] = P[4];
			nZ[i] = P[5];

			set_param(0, i, get_param(0, i)*sqr(RR)*fabs(g));
		}
	}
}

////////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for all curves;
void CGrid_QG::QG_curve(CMap * mp)
{
	if (1 == map_dim(mp))
	if (mp[0] == ID_MAP(1, SPHERE_GENUS)) QG_circle (mp); else
	if (mp[0] == ID_MAP(1, CYL_GENUS))	  QG_ellipse(mp);
}

/////////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for all surface;
void CGrid_QG::QG_surface(CMap * mp)
{
	if (2 == map_dim(mp))
	if (mp[0] == ID_MAP(2,	 SPHERE_GENUS)) QG_sphere	(mp); else
	if (mp[0] == ID_MAP(2,		 CYL_GENUS)) QG_cylinder(mp); else
	if (mp[0] == ID_MAP(2, SPHEROID_GENUS)) QG_spheroid(mp);
}

///////////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for all curvilinear triangles;
void CGrid_QG::QG_tria_curve(CMap * mp, double * Po, int N_ini)
{
	if (1 == map_dim(mp))
	if (mp[0] == ID_MAP(1, SPHERE_GENUS)) QG_tria_circle (mp, Po, N_ini); else
	if (mp[0] == ID_MAP(1, CYL_GENUS))	  QG_tria_ellipse(mp, Po, N_ini);
}

/////////////////////////////////////////////////////////////////////////////////////
//...auxilliary functions for correction weigts and nodes for curvilinear quadrangle;
void CGrid_QG::QG_quad_curve(CMap * mp1, CMap * mp2, double * Po, int N_ini)
{
	if (1 == map_dim(mp1)) {
		if (mp1[0] == ID_MAP(1, SPHERE_GENUS)) QG_quad_circle (mp1, mp2, Po, N_ini); else
		if (mp1[0] == ID_MAP(1, CYL_GENUS))		QG_quad_ellipse(mp1, mp2, Po, N_ini);
	}
	else
	if (1 == map_dim(mp2)) {
		if (mp2[0] == ID_MAP(1, SPHERE_GENUS)) QG_quad_circle (mp1, mp2, Po, N_ini); else
		if (mp2[0] == ID_MAP(1, CYL_GENUS))		QG_quad_ellipse(mp1, mp2, Po, N_ini);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for all curvilinear triangles;
void CGrid_QG::QG_tria_surface(CMap * mp, double * Po, int N_ini)
{
	if (2 == map_dim(mp))
	if (mp[0] == ID_MAP(2, SPHERE_GENUS)) QG_tria_sphere(mp, Po, N_ini);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for correction weigts and nodes for all curvilinear quadrangles;
void CGrid_QG::QG_quad_surface(CMap * mp1, CMap * mp2, double * Po, int N_ini)
{
	if (2 == map_dim(mp1) || 2 == map_dim(mp2))
	if (mp1 && mp1[0] == ID_MAP(2, SPHERE_GENUS) ||
		 mp2 && mp2[0] == ID_MAP(2, SPHERE_GENUS)) QG_quad_sphere(mp1, mp2, Po, N_ini);
}

///////////////////////////////////////////////////////
//...������� ��� ���������� ���������� ��������� �����;
void CGrid_QG::QG_quad_bi(double * Po, int N_ini)
{
	if (N_par < 3) add_params(3-N_par);

	double nnX, nnY, nnZ;
	nnX = (Po[7]-Po[1])*(Po[5]-Po[2])-(Po[8]-Po[2])*(Po[4]-Po[1]); 
	nnY = (Po[8]-Po[2])*(Po[3]-Po[0])-(Po[6]-Po[0])*(Po[5]-Po[2]); 
	nnZ = (Po[6]-Po[0])*(Po[4]-Po[1])-(Po[7]-Po[1])*(Po[3]-Po[0]);
	QG_normal(nnX, nnY, nnZ);

	for (int i = N_ini; i < N; i++) {
		double tX, tY, tZ, F1, F2, F3, F4, rr, tt, xx, P[6];

		tX = nY[i]*(Z[i]-Po[2])-nZ[i]*(Y[i]-Po[1]);
		tY = nZ[i]*(X[i]-Po[0])-nX[i]*(Z[i]-Po[2]);
		tZ = nX[i]*(Y[i]-Po[1])-nY[i]*(X[i]-Po[0]);
		F1 = QG_normal(tX, tY, tZ);

		tX = nY[i]*(Z[i]-Po[8])-nZ[i]*(Y[i]-Po[7]);
		tY = nZ[i]*(X[i]-Po[6])-nX[i]*(Z[i]-Po[8]);
		tZ = nX[i]*(Y[i]-Po[7])-nY[i]*(X[i]-Po[6]);
		F2 = QG_normal(tX, tY, tZ);

		tX = (Po[4]-Po[1])*(Z[i]-Po[2])-(Po[5]-Po[2])*(Y[i]-Po[1]);
		tY = (Po[5]-Po[2])*(X[i]-Po[0])-(Po[3]-Po[0])*(Z[i]-Po[2]);
		tZ = (Po[3]-Po[0])*(Y[i]-Po[1])-(Po[4]-Po[1])*(X[i]-Po[0]);
		F3 = QG_normal(tX, tY, tZ);

		tX = (Po[10]-Po[7])*(Z[i]-Po[8])-(Po[11]-Po[8])*(Y[i]-Po[7]);
		tY = (Po[11]-Po[8])*(X[i]-Po[6])-(Po[9] -Po[6])*(Z[i]-Po[8]);
		tZ = (Po[9] -Po[6])*(Y[i]-Po[7])-(Po[10]-Po[7])*(X[i]-Po[6]);
		F4 = QG_normal(tX, tY, tZ);
		
		rr = (F2*F3)/(F1*F4);
		tt = rr/(1.+rr);

		P[0] = Po[0]+(Po[6]-Po[0])*tt;
		P[1] = Po[1]+(Po[7]-Po[1])*tt;
		P[2] = Po[2]+(Po[8]-Po[2])*tt;

		P[3] = Po[3]+(Po[9] -Po[3])*tt;
		P[4] = Po[4]+(Po[10]-Po[4])*tt;
		P[5] = Po[5]+(Po[11]-Po[5])*tt;

		xx = sqrt((sqr(X[i]-P[0])+sqr(Y[i]-P[1])+sqr(Z[i]-P[2]))/(sqr(P[3]-P[0])+sqr(P[4]-P[1])+sqr(P[5]-P[2])));

		set_param(1, i, tt); //...�������������� �������� (0 < tt < 1);
		set_param(2, i, xx); //...������������ ��������	(0 < xx < 1);
		nX[i] = nnX;			//...��������������� �������;
		nY[i] = nnY;
		nZ[i] = nnZ;
	}
}