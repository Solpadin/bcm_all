#include "stdafx.h"

#include "cgrid_el.h"
#include "cgrid_qg.h"

//////////////////////////////////////////////////////////////////////////
//                  ALL EXISTING TYPES OF GRID NODES                    //
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//...global function for construction of all existing types of grid nodes;
CGrid * CreateNodes(Num_Nodes id_NODES)
{
	switch (id_NODES) {
     case GRID_EL_NODES: return new CGrid_el;
     case GRID_QG_NODES: return new CGrid_QG;
	}
   return new CGrid;
}

/////////////////////////////////////////////////
//          ��������� ����� ��� ABAQUS         //
/////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//...���������� ����� ����� � ������ (���������� � ������ ������ ������� ������); <----------------!!!
void Inp_nodes_add(char * ch_NODES, CGrid * nd, int ID_node_set, int ID_element_set)
{
	const int STR_SIZE = 250;
	char * id_NODES = read_struct_ascii(ch_NODES);
	if (nd && id_NODES) {
		char  buff[2000], one_line[STR_SIZE+1], * pchar, temp = 0;
		unsigned long ppos_cur = 0, upper_limit, upper, upper_element, end_element;
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, "_add_nodes.inp");
		FILE * INP = fopen(buff, "w");
		if (INP) {
			int i, j, k, l, m = 0, NN = nd->N, elem;

/////////////////////////////////////////
//...�������� ���� � ��������� ���������;
			if (nd->hit && nd->cond && nd->cond_ptr && ID_node_set < 0) {
				for (k = 0; k < nd->cond[0]; k++)
				if (nd->cond[l = nd->cond_ptr[k]] == (int)GL_NODES && nd->cond[l+2] == ID_node_set) {
					for (j = 2; j <= nd->cond[l+1]; j++) if (nd->cond[l+1+j] < nd->N) 
						if (nd->add_new_point(nd->X[nd->cond[l+1+j]], nd->Y[nd->cond[l+1+j]], nd->Z[nd->cond[l+1+j]])) 
							 nd->hit[nd->N-1] = nd->cond[l+1+j]+1; //...��������� ���������;
				}
				else
				if (nd->cond[l] == (int)GL_NODES_GENERATE && nd->cond[l+2] == ID_node_set) {
					for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) if (j < nd->N) 
						if (nd->add_new_point(nd->X[nd->cond[l+1+j]], nd->Y[nd->cond[l+1+j]], nd->Z[nd->cond[l+1+j]])) 
							 nd->hit[nd->N-1] = nd->cond[l+1+j]+1;
				}
			}	
		
////////////////////////////////////////
//...������ �������������� ����� ������;
			FILE * INP_ADD = fopen("INP_ADD", "w");
			if (INP_ADD) {
				for (k = NN; k < nd->N; k++) //...�������������� ����;
				fprintf(INP_ADD, "%7i, %12g, %12g, %12g\n", nd->hit[k], nd->X[k], nd->Y[k], nd->Z[k]);
				fclose (INP_ADD);
			}

///////////////////////////////
//...��������� �������� ������;
			if (nd->cond && nd->cond_ptr && ID_element_set < 0) {

////////////////////////////////////////////
//...���� ������������� ��������� ���������;
				for (k = 0; ! m && k < nd->cond[0]; k++)
				if (nd->cond[l = nd->cond_ptr[k]] == (int)GL_ELEMENTS && nd->cond[l+2] == ID_element_set) m = 1;
				else
				if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE && nd->cond[l+2] == ID_element_set) m = 1;

////////////////////////////////////////////////
//...���� �������������� ���� � �������� ������;
				if (nd->geom && nd->geom_ptr && m)
				for (k = NN; k < nd->N; k++) { //...�������������� ����;
					for (j = nd->cond[l+1]; j > 1; j--) if (nd->cond[l+1+j] < nd->geom[0]) //...�������������� ��������; 
					for (m = nd->geom[(i = nd->geom_ptr[ nd->cond[l+1+j]])+1]; m > 2; m--) 
					if (nd->hit[k] == nd->geom[i+1+m]+1) nd->geom[i+1+m] = k; //...����� �� ���������� ���������;
				}
				for (k = 0; k < nd->N; k++) nd->hit[k] = k+1; //...��������������� ���������;
			}

//////////////////////////////////////////////////////////////////////////
//...������ ����������������� ������ � ������ *Part �������� ����� ������;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
			if ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");

////////////////////////////////////////////////////////////////////////////////
//...���� ��������� ����� � ������� *Part � �������� ��� �� ������ ����� � ����;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Node");
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE); swap(pchar[0], temp);
	
				fprintf(INP, "%s", id_NODES); swap(pchar[0], temp);
				fprintf(INP, "%s", one_line);
				for (k = 0; k < nd->N; k++)
				fprintf(INP, "%7i, %12g, %12g, %12g\n", nd->hit[k], nd->X[k], nd->Y[k], nd->Z[k]);

////////////////////////////////////////////////////////////////////
//...���� ��� ��������� ��������� � �������� �� � �������� ��������;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); end_element = (upper_element -= (upper_element < upper ? 1 : 0));
					if (! ::strncmp(id_NODES+ppos_cur, "B31",   3) || //...,���������� ��� ��������; 
						 ! ::strncmp(id_NODES+ppos_cur, "B32",	  3)) elem = GL_LINE_STRIP; else
					if (! ::strncmp(id_NODES+ppos_cur, "CPE3",  4) ||
						 ! ::strncmp(id_NODES+ppos_cur, "CPE6",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS3",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS6",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "S3",	  2)) elem = GL_TRIANGLES; else
					if (! ::strncmp(id_NODES+ppos_cur, "CPE4",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPE8",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS4",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "CPS8",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "S4R",	  3)) elem = GL_QUADS; else 
					if (! ::strncmp(id_NODES+ppos_cur, "C3D4",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "C3D10", 5)) elem = GL_TETRA; else 
					if (! ::strncmp(id_NODES+ppos_cur, "C3D8",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "C3D20", 5)) elem = GL_BOXS; else 
					if (! ::strncmp(id_NODES+ppos_cur, "C3D6",  4) || 
						 ! ::strncmp(id_NODES+ppos_cur, "C3D15", 5)) elem = GL_PENTA; else elem = ERR_STATE;
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
					fprintf(INP, "*Element, type=%s", one_line);

////////////////////////////////////////////
//...�������� ��� �������� ������ ���������;
					for (l = k = 0; k < nd->geom[0]; k++) {
						if (nd->geom[(++l)++] == elem) {
							fprintf(INP, "%7i", -nd->geom[l+1]);					
							for (j = 2; j < nd->geom[l]; j++)
							fprintf(INP, ",%7i", nd->geom[l+1+j]+1);					
							fprintf(INP, "\n");					
						}
						l += nd->geom[l];
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
				}

////////////////////////////////////////////////
//...�������� ���������� ������� �������� �����;
				fprintf(INP, "%s", id_NODES+end_element);
			}
		}
		fclose(INP);
	}
	delete_struct(id_NODES);
}

//////////////////////////////////////////////////////////////////////////////////
//...���������� ����� ��������� � ������ (��������� ��������� � ����� -- �������);
void Inp_elements_add(char * ch_NODES, CGrid * nd, int * ID_elements, int ID_part)
{
	const int STR_SIZE = 250;
	char * id_NODES = read_struct_ascii(ch_NODES);
	if (nd && id_NODES) {
		char  buff[2000], one_line[STR_SIZE+1], * pchar, temp = 0;
		unsigned long ppos_cur = 0, upper_limit, upper, upper_element, end_element, count_set, end_set, begin_set;
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, "_add_elements.inp");
		FILE * INP = fopen(buff, "w");
		if (INP && ID_part > 0 && nd->geom && nd->geom_ptr) {
			int i, j, j_max = 18, k, l, l0, m, set, elem, N_geom, NN_geom, N_set = 0, N_part_element_beg, N_part_element_end, max_element;

///////////////////////////////////////////////////
//...���� ������ � ����� ��������� ��������� *Part;
			for (j = 0; j < nd->geom[0]; j++) 
				if (nd->geom[nd->geom_ptr[j]+2]+ID_part == 0) break;
			N_part_element_beg = j;				

			for (j = N_part_element_beg; j < nd->geom[0]; j++) 
				if (nd->geom[nd->geom_ptr[j]+2]+ID_part != 0) break;
			N_part_element_end = j;

//////////////////////////////////////////////////////
//...���������� ������������ �������� ��������� *Part;
			for (max_element = -nd->geom[nd->geom_ptr[j = N_part_element_beg]+3]; j < N_part_element_end; j++) 
			if  (max_element < -nd->geom[nd->geom_ptr[j]+3]) max_element = -nd->geom[nd->geom_ptr[j]+3];
			NN_geom = nd->geom[0];

//////////////////////////////////
//...���������� �������� � ������;
			if (nd->cond && nd->cond_ptr && ID_elements) 
			for (set = 1; set <= ID_elements[0]; set++) {

//////////////////////////////////////////
//...���� ����������� ��������� ���������;
				for (m = k = 0; ! m && k < nd->cond[0]; k++)
				if ((nd->cond[l = nd->cond_ptr[k]] == (int)GL_ELEMENTS || nd->cond[l] == (int)GL_ELEMENTS_GENERATE) 
				  && nd->cond[l+2] == ID_elements[set]) m = 1;

//////////////////////////////////////////////////////////
//...���������, ��� ��� ��������� ������ � �������� *Part;
				if (m) {
					if (nd->cond[l] == (int)GL_ELEMENTS)
						for ( i = 2; m && i <= nd->cond[l+1]; i++) if (nd->geom[nd->geom_ptr[nd->cond[l+1+i]]+2]+ID_part != 0) m = 0;
					else
					if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE)
						for (i = nd->cond[l+3]; m && i <= nd->cond[l+4]; i += nd->cond[l+5]) if (nd->geom[nd->geom_ptr[i]+2]+ID_part != 0) m = 0;
				}

//////////////////////////////////////////////////////////////////////////
//...��������� ������� �� ���������� ��������� � �������� *Part ���������;
				if (m) {
					N_geom = nd->geom[0];

///////////////////////////////////////////////////////////////////////////////////////
//...������������ ��������� ��������� �� ������������ ��������� (���������� ���������);
					if (nd->cond[l] == (int)GL_ELEMENTS) {
						for ( j = 2; j <= nd->cond[l+1]; j++)
							nd->geom_ptr_add(nd->geom[nd->geom_ptr[nd->cond[l+1+j]]+1], N_geom);
					}
					else
					if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE) {
						for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5])
							nd->geom_ptr_add(nd->geom[nd->geom_ptr[j]+1], N_geom);
					}

//////////////////////////////////////////
//...���������������� ��������� ���������;
					int * geom_new = NULL;
					if ( (geom_new = (int *)new_struct((nd->geom_ptr[N_geom]+1)*sizeof(int))) != NULL) {
						memcpy(geom_new, nd->geom, (nd->geom_ptr[nd->geom[0]]+1)*sizeof(int)); delete_struct(nd->geom);
						nd->geom = geom_new;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...��������� ������ �� ��������� ������������ ��������� (��� ����) � ������������� ������� ����� ������ �� ������������;
						N_geom = nd->geom[0];
						if (nd->cond[l] == (int)GL_ELEMENTS) {
							for ( j = 2; j <= nd->cond[l+1]; j++) {
								l0 = nd->geom[nd->geom_ptr[nd->cond[l+1+j]]+1]+2;
								memcpy(nd->geom+(i = nd->geom_ptr[nd->geom[0]]), nd->geom+nd->geom_ptr[nd->cond[l+1+j]], l0*sizeof(int));
								nd->geom[i+3] = -(++nd->geom[0])+NN_geom-max_element;
							}
						}
						else
						if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE) {
							for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) { 
								l0 = nd->geom[nd->geom_ptr[j]+1]+2;
								memcpy(nd->geom+(i = nd->geom_ptr[nd->geom[0]]), nd->geom+nd->geom_ptr[j], l0*sizeof(int));
								nd->geom[i+3] = -(++nd->geom[0])+NN_geom-max_element;
							}
						}
					}
				}

/////////////////////////////////////////////////////////////////
//...������������ � ��������� �������������� ��������� ���������;
				if (m && N_geom < nd->geom[0]) {
					for (l = j = k = 0; k < nd->cond[0]; k++) {
						if (nd->cond[(++l)++] == GL_ELEMENTS || nd->cond[l-1] == GL_ELEMENTS_GENERATE) --j;
						l += nd->cond[l];
					}
					nd->cond_ptr_add(4, nd->cond[0]);

//////////////////////////////////////////////////////////////////////////
//...���������������� ������� � ������� ������ � �������������� ���������;
					int * cond_new = NULL;
					if ( (cond_new = (int *)new_struct((nd->cond_ptr[nd->cond[0]]+1)*sizeof(int))) != NULL) {
						memcpy(cond_new, nd->cond, (nd->cond_ptr[nd->cond[0]-1]+1)*sizeof(int)); delete_struct(nd->cond);
						nd->cond = cond_new;

//////////////////////
//...��������� ������;
						nd->cond[i = nd->cond_ptr[nd->cond[0]-1]] = GL_ELEMENTS_GENERATE;
						nd->cond[i+1] = 4;
						nd->cond[i+2] = --j;
						nd->cond[i+3] = N_geom;
						nd->cond[i+4] = nd->geom[0]-1;
						nd->cond[i+5] = 1;
						N_set++;
					}
				}
			}

////////////////////////////////////////////////////////////////////////////
//...������ ����������������� ������ � �������� *Part �������� ����� ������;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
			if ((upper = begin_set = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
				PPOS_CUR(id_NODES, pchar, begin_set, upper, "*Elset, elset=");
			}
			for ( j = 1; j < ID_part; j++) 
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
			if ((upper = count_set = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
				PPOS_CUR(id_NODES, pchar, count_set, upper, "*Elset, elset=");

/////////////////////////////////////////////////////////////////////////////////
//...���� ��������� ��������� � ������� *Part � �������� ��� �� ������ ���������;
				//pchar = id_NODES+50; swap(pchar[0], temp);
				//fprintf(INP, "%s", id_NODES); swap(pchar[0], temp); 
				//return;


				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type="); swap(pchar[0], temp);
				fprintf(INP, "%s", id_NODES); swap(pchar[0], temp);

////////////////////////////////////////////////////////////////////
//...���� ��� ��������� ��������� � �������� �� � �������� ��������;
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); end_element = (upper_element -= (upper_element < upper ? 1 : 0));
					
					nd->element_type(elem, id_NODES+ppos_cur);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
					
					fprintf(INP, "*Element, type=%s", one_line);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...�������� ��� �������� ������ ��������� �� ��������� *Part (������ ����� ��������� �� ������� ���������);
					for (l = nd->geom_ptr[k = N_part_element_beg]-1; k < N_part_element_end; k++) {
						if (nd->geom[(++l)++] == elem) {
							fprintf(INP, "%7i", -nd->geom[l+2]);					
							for (j = 3; j < nd->geom[l]; j++) {
								if (j != j_max) fprintf(INP, ",%7i", nd->hit[nd->geom[l+1+j]]);
								else			 fprintf(INP, ",\n%15i", nd->hit[nd->geom[l+1+j]]);
							}
							fprintf(INP, "\n");					
						}
						l += nd->geom[l];
					}
					for (l = nd->geom_ptr[k = NN_geom]-1; k < nd->geom[0]; k++) {
						if (nd->geom[(++l)++] == elem) {
							fprintf(INP, "%7i", -nd->geom[l+2]);					
							for (j = 3; j < nd->geom[l]; j++) {
								if (j != j_max) fprintf(INP, ",%7i", nd->hit[nd->geom[l+1+j]]);
								else			 fprintf(INP, ",\n%15i", nd->hit[nd->geom[l+1+j]]);
							}
							fprintf(INP, "\n");					
						}
						l += nd->geom[l];
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Element, type=");
				}

/////////////////////////////////////////////////////////////////////////////////////////////
//...������������� ��� ������������ ��������� � ���� �����, ���� �������� ����� ������������;
				ppos_cur = count_set;
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); end_set = (upper_element -= (upper_element < upper ? 1 : 0));
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
				}
				swap(id_NODES[end_set], temp); fprintf(INP, "%s", id_NODES+end_element);
				swap(id_NODES[end_set], temp);

///////////////////////////////////////////////////////////////////////
//...����� ����������� ��������� ��������� � ����������������� �������;
				if (nd->cond && nd->cond_ptr && ID_elements) 
				for (set = 1; set <= ID_elements[0]; set++) {

//////////////////////////////////////////
//...���� ����������� ��������� ���������;
					for (m = j = k = 0; ! m && k < nd->cond[0]-N_set; k++)
					if ((nd->cond[l = nd->cond_ptr[k]] == (int)GL_ELEMENTS || nd->cond[l] == (int)GL_ELEMENTS_GENERATE) && --j &&
						  nd->cond[l+2] == ID_elements[set]) m = 1;

//////////////////////////////////////////////////////////
//...���������, ��� ��� ��������� ������ � �������� *Part;
					if (m) {
						if (nd->cond[l] == (int)GL_ELEMENTS)
							for ( i = 2; m && i <= nd->cond[l+1]; i++) if (nd->geom[nd->geom_ptr[nd->cond[l+1+i]]+2]+ID_part != 0) m = 0;
						else
						if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE)
							for (i = nd->cond[l+3]; m && i <= nd->cond[l+4]; i += nd->cond[l+5]) if (nd->geom[nd->geom_ptr[i]+2]+ID_part != 0) m = 0;
					}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//...���� ��� ������������ ��������� ��������� � ����� �������� ��������� ��������� �� ������� ���������;
					if (m && N_geom < nd->geom[0]) {
						ppos_cur = begin_set;
						while ((upper_element = ppos_cur) < upper && ++j) {
							PPOS_CUR(id_NODES, pchar, upper_element, upper, "*");
							PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
						}
						if (! j) {
							ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
							if (pchar = strstr(one_line, ",")) {
								swap(pchar[0], temp); fprintf(INP, "*Elset, elset=%s_add, generate", one_line);
								swap(pchar[0], temp); fprintf(INP, "\xA");
							}
							else
							if (pchar = strstr(one_line, "\xA")) {
								swap(pchar[0], temp); fprintf(INP, "*Elset, elset=%s_add, generate", one_line);
								swap(pchar[0], temp); fprintf(INP, "%s", pchar);
							}
							fprintf(INP, "%s", pchar); /*!!! was "", now "%s" ???*/ i = nd->cond_ptr[nd->cond[0]-N_set];
							fprintf(INP, "%7i,%7i,%7i\n", -nd->geom[nd->geom_ptr[nd->cond[i+3]]+3], -nd->geom[nd->geom_ptr[nd->cond[i+4]]+3], nd->cond[i+5]);					
							N_set--;
						}
					}
				}

////////////////////////////////////////////////
//...�������� ���������� ������� �������� �����;
				fprintf(INP, "%s", id_NODES+end_set);
			}
		}
		fclose(INP);
	}
	delete_struct(id_NODES);
}

/////////////////////////////////////////////////////////////////////
//          INTERFACE FUNCTIONS FOR CONVERTING INP-FORMAT          //
/////////////////////////////////////////////////////////////////////
//...���������� ������ ���������� ������� (��������� �����);
unsigned long * Condit_list_nodes(char * id_NODES, unsigned long count, unsigned long upper_limit, int id_status = OK_STATE)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long ppos_cur = count, pposs, upper, upper_element, N_buf = 50, buf_incr = 20,
						* nodes_list = (unsigned long *)new_struct((N_buf*2+1)*sizeof(unsigned long));
		const int STR_SIZE = 250;
		int   elem, j = -1, j_elem = -1;
		char one_line[STR_SIZE+1], * pchar;

//////////////////////////////////////////////////////////////
//...���������� ����������� ��������c��� �� ������� *Assembly;
		if (id_status) {
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Assembly, name="); upper = ppos_cur;
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Assembly"); upper -= sizeof("*End Assembly");
			
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Instance, name="); upper_element = ppos_cur;
			while (upper_element < upper) {
				PPOS_CUR(id_NODES, pchar, ppos_cur,	upper, "*End Instance");
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*Instance, name=");
			}

//////////////////////////////////
//...��������� ������������ �����;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
			while ((upper_element = pposs = ppos_cur) < upper) {
				if (nodes_list[0] == N_buf) {
					unsigned long * new_nodes_list = nodes_list; nodes_list = (unsigned long *)new_struct(((N_buf += buf_incr)*2+1)*sizeof(unsigned long));
					memcpy(nodes_list, new_nodes_list, (new_nodes_list[0]*2+1)*sizeof(unsigned long)); delete_struct(new_nodes_list);
				}
				nodes_list[nodes_list[0]*2+1] = ppos_cur;
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

////////////////////////////////////////////////////////
//...,���������� ��� ���� (����� ��������� ���� ������);
				if (! ::strncmp(id_NODES+pposs, "CONTACT",  7) || 
						! ::strncmp(id_NODES+pposs, "LOADING",  7) || 
						! ::strncmp(id_NODES+pposs, "LATERAL_BOUND", 13) ||
						! ::strncmp(id_NODES+pposs, "BOTTOM_BOUND",  12)) elem = GL_NODES; else elem = ERR_STATE;
				if (elem == GL_NODES && strstr(one_line, "internal")) elem = ERR_STATE; 
				if (elem == GL_NODES && strstr(one_line, "generate")) elem = GL_NODES_GENERATE; 
	
				if (elem == GL_NODES || elem == GL_NODES_GENERATE)
				if ((pchar = strstr(one_line, ", generate")) != NULL || (pchar = strstr(one_line, "\xD")) != NULL || (pchar = strstr(one_line, "\xA")) != NULL) {
					pchar[0] = '\0';
					nodes_list[nodes_list[0]*2+2] = strlen(one_line);
					nodes_list[0]++;
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
			}
		}

/////////////////////////////////////
//...������������� ��� ������� *Part;
		ppos_cur = count;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");

//////////////////////////////////
//...��������� ������������ �����;
			while ((upper_element = ppos_cur) < upper) {
				if (nodes_list[0] == N_buf) {
					unsigned long * new_nodes_list = nodes_list; nodes_list = (unsigned long *)new_struct(((N_buf += buf_incr)*2+1)*sizeof(unsigned long));
					memcpy(nodes_list, new_nodes_list, (new_nodes_list[0]*2+1)*sizeof(unsigned long)); delete_struct(new_nodes_list);
				}
				nodes_list[nodes_list[0]*2+1] = ppos_cur;
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////
//...,���������� ��� ����;
				if (strstr(one_line, "internal")) elem = ERR_STATE; 
				if (strstr(one_line, "generate")) elem = GL_NODES_GENERATE; else elem = GL_NODES;

				if (elem == GL_NODES || elem == GL_NODES_GENERATE) 
				if ((pchar = strstr(one_line, ", generate")) != NULL || (pchar = strstr(one_line, "\xD")) != NULL || (pchar = strstr(one_line, "\xA")) != NULL) {
					pchar[0] = '\0';
					nodes_list[nodes_list[0]*2+2] = strlen(one_line);
					nodes_list[0]++;
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Nset, nset=");
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}
		return(nodes_list);
	}
	return(NULL); 
}

////////////////////////////////////////////////////////////////
//...���������� ������ ���������� ������� (��������� ���������);
unsigned long * Condit_list_elements(char * id_NODES, unsigned long count, unsigned long upper_limit)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long ppos_cur = count, upper, upper_element, N_buf = 50, buf_incr = 20,
						* elements_list = (unsigned long *)new_struct((N_buf*2+1)*sizeof(unsigned long));
		const int STR_SIZE = 250;
		int  elem, j = -1, j_elem = -1;
		char one_line[STR_SIZE+1], * pchar;

/////////////////////////////////////
//...������������� ��� ������� *Part;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");

//////////////////////////////////////
//...��������� ������������ ���������;
			while ((upper_element = ppos_cur) < upper) {
				if (elements_list[0] == N_buf) {
					unsigned long * new_elements_list = elements_list; elements_list = (unsigned long *)new_struct(((N_buf += buf_incr)*2+1)*sizeof(unsigned long));
					memcpy(elements_list, new_elements_list, (new_elements_list[0]*2+1)*sizeof(unsigned long)); delete_struct(new_elements_list);
				}
				elements_list[elements_list[0]*2+1] = ppos_cur;
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

//////////////////////////////
//...,���������� ��� ��������;
				if (strstr(one_line, "internal")) elem = ERR_STATE; 
				if (strstr(one_line, "generate")) elem = GL_ELEMENTS_GENERATE; else elem = GL_ELEMENTS; 

				if (elem == GL_ELEMENTS || elem == GL_ELEMENTS_GENERATE) 
				if ((pchar = strstr(one_line, ", internal")) != NULL || (pchar = strstr(one_line, ", generate")) != NULL || 
					 (pchar = strstr(one_line, "\xD")) != NULL || (pchar = strstr(one_line, "\xA")) != NULL) {
					pchar[0] = '\0';
					elements_list[elements_list[0]*2+2] = strlen(one_line);
					elements_list[0]++;
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Elset, elset=");
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}
		return(elements_list);
	}
	return(NULL); 
}

///////////////////////////////////////////////////////
//...���������� ������ ������������ � �������� *Part'e;
unsigned long * Part_surfaces_list(char * id_NODES, unsigned long count, unsigned long upper_limit, int ID_element_part)
{
	if (id_NODES &&  count < upper_limit) {
		unsigned long ppos_cur = count, upper, upper_element, N_buf = 50, buf_incr = 20,
						* surfaces_list = (unsigned long *)new_struct((N_buf*3+1)*sizeof(unsigned long));
		const int STR_SIZE = 250;
		int  k_side, N_part = ID_element_part;
		char one_line[STR_SIZE+1], * pchar, temp = 0;

/////////////////////////////////////
//...������������� ��� ������� *Part;
		PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		while ((upper = ppos_cur) < upper_limit) {
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part");
			if (--N_part == 0) {
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Surface, type=ELEMENT, name=");

//////////////////////////////////////
//...��������� ��������� ������������;
				while ((upper_element = ppos_cur) < upper) {
					PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
					ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

					while (ppos_cur < upper_element) {
						if (surfaces_list[0] == N_buf) {
							unsigned long * new_surfaces_list = surfaces_list; surfaces_list = (unsigned long *)new_struct(((N_buf += buf_incr)*3+1)*sizeof(unsigned long));
							memcpy(surfaces_list, new_surfaces_list, (new_surfaces_list[0]*3+1)*sizeof(unsigned long)); delete_struct(new_surfaces_list);
						}
						surfaces_list[surfaces_list[0]*3+1] = ppos_cur;
						ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);

						if (pchar = strstr(one_line, ",")) {
							swap(pchar[0], temp); surfaces_list[surfaces_list[0]*3+2] = strlen(one_line);	 k_side = atoi(pchar+3);
							swap(pchar[0], temp); surfaces_list[surfaces_list[0]*3+3] = k_side;	
							surfaces_list[0]++;
						}
					}
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Surface, type=ELEMENT, name=");
				}
			}
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
		}
		return(surfaces_list);
	}
	return(NULL); 
}

////////////////////////////////////////////////////////////
//...����������� ��������� �� 3D ������ Abaqus (���� *.gmt);
void Convert3D_gmt(char * ch_NODES, CGrid * nd, int ID_part)
{
	char * id_NODES = read_struct_ascii(ch_NODES), buff[2000];
	if (nd && id_NODES) {
		const int STR_SIZE = 250, ID_hexa_elem = 300;
		unsigned long ppos_cur = 0, upper_limit, upper, upper_element, count, count_part, count_length, * elements_list = NULL;
		char one_line[STR_SIZE+1], * pchar, temp = '\0';
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, ".gmt");
		FILE * GMT = fopen(buff, "w");
		if (GMT) {
			int j, k, l, m, N_part = ID_part, N_part_nodes_beg, N_part_nodes_end, N_part_element_beg, N_part_element_end, dim = 3;

//////////////////////////
//...������ ����� �������;
			char job_name[STR_SIZE+1], model_name[STR_SIZE+1];
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "** Job name: ");
			ONE_LINE(id_NODES, pchar, ppos_cur, upper_limit, job_name, STR_SIZE); count = 0; upper = (unsigned long)strlen(job_name);
			PPOS_CUR(job_name, pchar, count, upper, " Model name: "); job_name[upper-2] = 0;

			if (pchar) pchar[0] = 0; strcpy(model_name, job_name+count); 	
			fprintf(GMT, "[MESH NAME]\n%s\n%s\n<END>\n%i   // ����������� ������������\n", job_name, model_name, dim);
			fprintf(GMT, "//------------------------------------------------------------------------------\n");

///////////////////////////////////////////////
//...���� ������ � ����� ����� ��������� *Part;
			for (N_part_nodes_beg = j = 0; j < nd->N; j++) if (! j || nd->hit[j] < nd->hit[j-1]) {//...������ ���� part'�;
				if (N_part--) N_part_nodes_beg = j;				
				else break;
			}
			N_part_nodes_end = j;

///////////////////////////////////////////////////
//...���� ������ � ����� ��������� ��������� *Part;
			for (j = 0; j < nd->geom[0]; j++) 
				if (nd->geom[nd->geom_ptr[j]+2]+ID_part == 0) break;
			N_part_element_beg = j;				

			for (j = N_part_element_beg; j < nd->geom[0]; j++) 
				if (nd->geom[nd->geom_ptr[j]+2]+ID_part != 0) break;
			N_part_element_end = j;

/////////////////////////////////
//...����� ������ ������� part'�;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name="); count_part = ppos_cur;
			if ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); 
			}
			for ( j = 1; j < ID_part; j++) {
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
				if ((upper = ppos_cur) < upper_limit) {
					PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); count_part = ppos_cur; 
				}
			}

////////////////////////////////////////////////
//...������ ���������� ��������� (������������);
			unsigned long Color[4] = {0x32C8C8C8, 0x00808080, 0x0080FFFF, 0xC80000FF}, index = 0, index_max = 4, count_light = 0, count_max = 7;
			elements_list = Condit_list_elements(id_NODES, 0, upper_limit); ppos_cur = count_part;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");

			fprintf(GMT, "[MATERIAL COLORS]\n");
			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
			
				if ( (pchar = strstr(one_line, ", material")) != NULL) {
					pchar[0] = temp; count_length = strlen(one_line);

					for (count = 0; count < elements_list[0]; count++)  
						if (elements_list[count*2+2] == count_length && ! strncmp(one_line, id_NODES+elements_list[count*2+1], count_length)) break;
					
					if (count < elements_list[0]) {
						count_light++;
						float blue  = (float)(((count_light/1)%2) ? ((Color[index] & 0x000000ff) >>  0) : (Color[index] & 0xff000000) >> 24)/255.f,
								green = (float)(((count_light/2)%2) ? ((Color[index] & 0x0000ff00) >>  8) : (Color[index] & 0xff000000) >> 24)/255.f,
								red   = (float)(((count_light/4)%2) ? ((Color[index] & 0x00ff0000) >> 16) : (Color[index] & 0xff000000) >> 24)/255.f;
						fprintf(GMT, "%7i   %3.2g   %3.2g   %3.2g\n", (int)count+1, red, green, blue);
						fprintf(GMT, "�������� ID=%i   // ������������� � ��������: %s\n<END>\n", (int)count+1, one_line);
						if (count_light == count_max) {
							 count_light = 0; index++;
							 if (index == index_max-1) {
								  count_max = 4; count_light = 2;
							 }
							 if (index == index_max) {
								  index = 0; count_max = 7;
								  Color[0] -= 0x10202020;
								  Color[1] -= 0x00101010;
								  Color[2] -= 0x00102020;
								  Color[3] -= 0x20000020;
							 }
						}
					}
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");
			}
			fprintf(GMT, "0   // ������� ����� ������\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");
			delete_struct(elements_list);

//////////////////////////////////
//...������ ����� (������ ������);
			fprintf(GMT, "[FIGURES]\n0   // ������� ����� ������\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");

/////////////////////////
//...������ ����� ������;
			fprintf(GMT, "[NODES]\n// �������������� ���������:\n//__________X_______Y_______Z:\n           1.0     1.0     1.0\n");
			fprintf(GMT, "// ���������� �� �����:\n//_ID_______X_______Y_______Z:\n");
			for (k = N_part_nodes_beg; k < N_part_nodes_end; k++) //...������� ����;
			fprintf(GMT, "%7i   %12g    %12g    %12g\n", nd->hit[k], nd->X[k], nd->Y[k], nd->Z[k]);
			fprintf(GMT, "0   // ������� ����� ������\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");

////////////////////////////
//...������ �������� ������;
			fprintf(GMT, "[ELEMENTS]\n// ���������� �� ���������:\n");
			fprintf(GMT, "//_ID__mat__type_____1_____2_____3_____4_____5_____6_____7_____8    ������� ����\n");
			fprintf(GMT, "// �������������� ���� (�� 10)\n//              _____9____10____11____12____13____14____15____16____17____18____19\n//              ____20____21 ... ������ �����\n");
			fprintf(GMT, "//���-�� ���������� ����� �� �����:   __1__2__3__4__5:\n");
			if (nd->geom && nd->geom_ptr && nd->hit) {
				for (l = nd->geom_ptr[k = N_part_element_beg];	k < N_part_element_end; l = nd->geom_ptr[++k]) { //...������� ��������;
					if (nd->geom[l] == (int)GL_BOXS)	 m = ID_hexa_elem; else m = ERR_STATE;
					if (m != ERR_STATE) {
						fprintf(GMT, "%7i   %7i   %7i", -nd->geom[l+3], -nd->geom[l+4], m);
						for (j = 3; j < nd->geom[l+1]; j++) if (nd->geom[l+j+2] < nd->N)
						fprintf(GMT, "   %7i", nd->hit[nd->geom[l+j+2]]);
						fprintf(GMT, "\n");
					}
				}
			}
			fprintf(GMT, "0   // ������� ����� ������\n");
			fprintf(GMT, "//------------------------------------------------------------------------------\n");
			fclose (GMT);
		}
	}
	delete_struct(id_NODES);
}

/////////////////////////////////
//...��������� ����������� �����;
void SetLoading_nodes(CGrid * nd, int elem_int, int k_side)
{
	int i, m1, m2, m3, m4, m5, m6, m7, m8;
	switch (nd->geom[i = nd->geom_ptr[elem_int]]) {
		case GL_BOXS: if (nd->geom[i+1] == 11) {
			m1 = nd->geom[i+5];
			m2 = nd->geom[i+6];
			m3 = nd->geom[i+7];
			m4 = nd->geom[i+8];
			m5 = nd->geom[i+9];
			m6 = nd->geom[i+10];
			m7 = nd->geom[i+11];
			m8 = nd->geom[i+12];
			if (k_side == 6) {//...face lateral S6;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m5] = -abs(nd->hit[m5]);
				nd->hit[m8] = -abs(nd->hit[m8]);
				nd->hit[m4] = -abs(nd->hit[m4]);
			}
			if (k_side == 4) {//...face lateral S4;
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m7] = -abs(nd->hit[m7]);
				nd->hit[m6] = -abs(nd->hit[m6]);
			}
			if (k_side == 3) {//...face lateral S3;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m6] = -abs(nd->hit[m6]);
				nd->hit[m5] = -abs(nd->hit[m5]);
			}
			if (k_side == 5) {//...face lateral S5;
				nd->hit[m4] = -abs(nd->hit[m4]);
				nd->hit[m8] = -abs(nd->hit[m8]);
				nd->hit[m7] = -abs(nd->hit[m7]);
				nd->hit[m3] = -abs(nd->hit[m3]);
			}
			if (k_side == 1) {//...face basis S1;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m4] = -abs(nd->hit[m4]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m2] = -abs(nd->hit[m2]);
			}
			if (k_side == 2) {//...face basis S2;
				nd->hit[m5] = -abs(nd->hit[m5]);
				nd->hit[m6] = -abs(nd->hit[m6]);
				nd->hit[m7] = -abs(nd->hit[m7]);
				nd->hit[m8] = -abs(nd->hit[m8]);
			}
		}  break;
		case GL_PENTA: if (nd->geom[i+1] == 9) {
			m1 = nd->geom[i+5];
			m2 = nd->geom[i+6];
			m3 = nd->geom[i+7];
			m4 = nd->geom[i+8];
			m5 = nd->geom[i+9];
			m6 = nd->geom[i+10];
			if (k_side == 3) {//...face lateral S3;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m5] = -abs(nd->hit[m5]);
				nd->hit[m4] = -abs(nd->hit[m4]);
			}
			if (k_side == 4) {//...face lateral S4;
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m6] = -abs(nd->hit[m6]);
				nd->hit[m5] = -abs(nd->hit[m5]);
			}
			if (k_side == 5) {//...face lateral S5;
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m4] = -abs(nd->hit[m4]);
				nd->hit[m6] = -abs(nd->hit[m6]);
			}
			if (k_side == 1) {//...face basis S1;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m2] = -abs(nd->hit[m2]);
			}
			if (k_side == 2) {//...face basis S2;
				nd->hit[m4] = -abs(nd->hit[m4]);
				nd->hit[m5] = -abs(nd->hit[m5]);
				nd->hit[m6] = -abs(nd->hit[m6]);
			}
		}  break;
		case GL_TETRA: if (nd->geom[i+1] == 7) {
			m1 = nd->geom[i+5];
			m2 = nd->geom[i+6];
			m3 = nd->geom[i+7];
			m4 = nd->geom[i+8];
			if (k_side == 1) {//...face basis S1;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m2] = -abs(nd->hit[m2]);
			}
			if (k_side == 3) {//...face lateral S3;
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m4] = -abs(nd->hit[m4]);
			}
			if (k_side == 4) {//...face lateral S4;
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m4] = -abs(nd->hit[m4]);
			}
			if (k_side == 2) {//...face lateral S2;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m4] = -abs(nd->hit[m4]);
			}  
		}  break;
		case GL_PYRAMID: if (nd->geom[i+1] == 8) {
			m1 = nd->geom[i+5];
			m2 = nd->geom[i+6];
			m3 = nd->geom[i+7];
			m4 = nd->geom[i+8];
			m5 = nd->geom[i+9];
			if (k_side == 2) {//...face lateral S2;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m5] = -abs(nd->hit[m5]);
			}
			if (k_side == 3) {//...face lateral S3;
				nd->hit[m2] = -abs(nd->hit[m2]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m5] = -abs(nd->hit[m5]);
			}
			if (k_side == 4) {//...face lateral S4;
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m4] = -abs(nd->hit[m4]);
				nd->hit[m5] = -abs(nd->hit[m5]);
			}
			if (k_side == 5) {//...face lateral S5;
				nd->hit[m4] = -abs(nd->hit[m4]);
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m5] = -abs(nd->hit[m5]);
			}  
			if (k_side == 1) {//...fase basis S1;
				nd->hit[m1] = -abs(nd->hit[m1]);
				nd->hit[m4] = -abs(nd->hit[m4]);
				nd->hit[m3] = -abs(nd->hit[m3]);
				nd->hit[m2] = -abs(nd->hit[m2]);
			}
		}  break;
	}
}

/////////////////////////////////////////////////////////////////////////
//...����������� �������� � ����������� �� 3D ������ Abaqus (���� *.prb);
void Convert3D_prb(char * ch_NODES, CGrid * nd, int ID_part)
{
	char * id_NODES = read_struct_ascii(ch_NODES), buff[2000];
	if (nd && id_NODES) {
		const int STR_SIZE = 250; char one_line[STR_SIZE+1], material_line[STR_SIZE+1], * pchar, * ppos, * material_name, temp = '\0';
		unsigned long ppos_cur = 0, ppos_material, ppos_save, upper_limit, upper, upper_element, upper_material, k_side, 
						  count, count_part, count_beg, count_length, count_el, * nodes_list = NULL, * elements_list = NULL, * surfaces_list = NULL;
		user_Count (id_NODES, 0, upper_limit, '\x0');
		
		strcpy(buff, ch_NODES);	pchar = ::strrchr(buff, '.'); if (pchar) pchar[0] = 0; strcat(buff, ".prb");
		FILE * PRB = fopen(buff, "w");
		if (PRB) {
			int k, l, j, j_node, j_max = 10, n_steps = 0, buf_steps = 20, buf_incr = 10, m[3];
			double * steps = (double *)new_struct(buf_steps*sizeof(double)), sum_steps = 0., step;

//////////////////////////
//...������ ����� �������;
			char job_name[STR_SIZE+1], model_name[STR_SIZE+1];
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "** Job name: "); count_beg = ppos_cur;
			ONE_LINE(id_NODES, pchar, ppos_cur, upper_limit, job_name, STR_SIZE); count = 0; upper = (unsigned long)strlen(job_name);
			PPOS_CUR(job_name, pchar, count, upper, " Model name: "); job_name[upper-2] = 0;

			if (pchar) pchar[0] = 0; strcpy(model_name, job_name+count); 	
			fprintf(PRB, "[PROBLEM NAME]\n%s\n%s\n<END>\n", job_name, model_name);
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

//////////////////////////////////////////////////
//...������ ���������� ������ (�����, ��� ������);
			fprintf(PRB, "[START]\n0   // 0 - ���� � \"0\", ����� - �������� ����\n// ���� � ����������� ��� ��������:\nD:\\Progs\\UWAY\\DATA\\Mike_Test\\Name_RESTART.dat\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

			fprintf(PRB, "[PROBLEM]\n1 3   // ��� ������ � ��� ����������� ���������\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

////////////////////////////////////////////////
//...������� � ������ ���������� ������� (����);
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Step, name=");
			while ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Step"); upper -= sizeof("*End Step");
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				if (n_steps == buf_steps) {
					double * prev_steps = steps; steps = (double *)new_struct((buf_steps += buf_incr)*sizeof(double));
					memcpy(steps, prev_steps, n_steps*sizeof(double));
					delete_struct(prev_steps);
				}
				if ((ppos = strstr(one_line, ",")) != NULL) steps[n_steps] = user_strtod(ppos+1);
				else													  steps[n_steps] = 1.;
				sum_steps += steps[n_steps]; steps[n_steps++] = sum_steps;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Step, name=");
			}
			fprintf(PRB, "[TIMES]\n0.0	// ����� ������    (T0)\n");
			fprintf(PRB, "%2.1f	// ����� ��������� (Tn)\n",   sum_steps);
			fprintf(PRB, "%2.1f	// ��� �� �������  (Step)\n", sum_steps/n_steps);
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

///////////////////////
//...����� �����������;
			fprintf(PRB, "[RESULT OUTS]\n0		// \"0\" - ����� ����������� � ������� ����������\n");
			fprintf(PRB, "2		// ��� ������\n// 0 - ���� + ����� �������������� (��������)\n// 1 - ���� + ���� (��������)\n// 2 - ���� + ������ ������� (��������)\n// 3 - ��� ����\n");
			fprintf(PRB, "0		// ��� ������� ������ �����������\n// 0 - ������ UWay  (������������)\n// 1 - ������ FEMAP (����������������)\n");
			fprintf(PRB, "// ������� �������, � ������� ����� ������������� ����� ����������� ��������:\n");
			fprintf(PRB, "//____1_______2_______3_______4_______5_______6_______7_______8_______9______10:");
			for (j = 0; j < n_steps; j++) {
				if (j%j_max) fprintf(PRB, "%8.1f", steps[j]);
				else			 fprintf(PRB, "\n%8.1f", steps[j]);
			}
			if (j%j_max) fprintf(PRB, "      0      //<<< \"0\" - the end of data\n");
			else			 fprintf(PRB, "\n      0      //<<< \"0\" - the end of data\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");
			
///////////////////////
//...��������� �������;
			fprintf(PRB, "[STATE 0]\n0   // 0 - ��� ��������� ������� = 0.0;\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

/////////////////////////////////////////////
//...������ ������� ����������� ����� ������;
			fprintf(PRB, "[FIXED NODES]\n// ������ ��������� ����������� ����� � �����������:\n0.0   ����������� � t=0.0\n");
			fprintf(PRB, "//_N/F__ID__X__Y__Z:   ����������� (N/F: \"0\"-����; \"1\"-������; ID-����/������)\n");

//////////////////////////////////////////////////////////////////////////////
//...����� � ������ ������� ����������� �� ������� step'� (��������� �������);
			nodes_list = Condit_list_nodes(id_NODES, 0, upper_limit); ppos_cur = upper = count_beg;
			PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*Step, name=");
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Boundary"); m[0] = m[1] = m[2] = 0;

			while (ppos_cur < upper) {
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);

				if ((ppos = strstr(one_line, ".")) != NULL)
				if ((pchar = strstr(ppos, ", ")) != NULL) {
					swap(temp, pchar[0]); count_length = strlen(ppos+1);
					for (count = 0; count < nodes_list[0]; count++)  
						if (nodes_list[count*2+2] == count_length && ! strncmp(ppos+1, id_NODES+nodes_list[count*2+1], count_length)) break;
					swap(temp, pchar[0]); k = (int)count; 
					if ((ppos = strstr(pchar, ", ")) != NULL) {
						ppos[0] = temp; j = atoi(pchar+1);
						if (0 < j && j < 4) m[j-1] = 1;
					}
					ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
					if (one_line[0] != '*' && (pchar = strstr(one_line, ",")) != NULL) {
						if ((ppos = strstr(pchar, ", ")) != NULL) {
							ppos[0] = temp; j = atoi(pchar+1);
							if (0 < j && j < 4) m[j-1] = 1;
						}
						ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
						if (one_line[0] != '*' && (pchar = strstr(one_line, ",")) != NULL) {
							if ((ppos = strstr(pchar, ", ")) != NULL) {
								ppos[0] = temp; j = atoi(pchar+1);
								if (0 < j && j < 4) m[j-1] = 1;
							}
						}
					}
				}
///////////////////////////////////
//...������ ���������� �����������;
				if (nd->cond_ptr && nd->cond && nd->hit) {
					l = nd->cond_ptr[k];
					if (nd->cond[l] == (int)GL_NODES) {
						for (j = 2; j <= nd->cond[l+1]; j++) if (nd->cond[l+1+j] < nd->N)
							fprintf(PRB, "    0   %7i   %7i   %7i   %7i\n", nd->hit[nd->cond[l+1+j]], m[0], m[1], m[2]);
					}
					else
					if (nd->cond[l] == (int)GL_NODES_GENERATE) {
						for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) if (j < nd->N)
							fprintf(PRB, "    0   %7i   %7i   %7i   %7i\n", nd->hit[j], m[0], m[1], m[2]);
					}
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Boundary"); m[0] = m[1] = m[2] = 0;
			}
			fprintf(PRB, "<END>\n");
			fprintf(PRB, "<END>   //<<< the end of data\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");
			delete_struct(steps);

/////////////////////////////////
//...����� ������ ������� part'�;
			ppos_cur = 0;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name="); count_part = ppos_cur;
			if ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); 
			}
			for ( j = 1; j < ID_part; j++) {
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Part, name=");
				if ((upper = ppos_cur) < upper_limit) {
					PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Part"); upper -= sizeof("*End Part"); count_part = ppos_cur; 
				}
			}

/////////////////////////////////////////////////////////////////////////////
//...������ ���������� ��������� (���������� ��������� ���������� � �������);
			int material_first = 1;
			elements_list = Condit_list_elements(id_NODES, 0, upper_limit); ppos_cur = count_part;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");
			
			fprintf(PRB, "[MATERIAL DATA]\n");
			while ((upper_element = ppos_cur) < upper) {
				PPOS_CUR(id_NODES, pchar, upper_element, upper, "*"); upper_element -= (upper_element < upper ? 1 : 0);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper_element, one_line, STR_SIZE);
			
				if ( (pchar = strstr(one_line, ", material=")) != NULL) {
					swap(pchar[0], temp); count_length = strlen(one_line);

					for (count = 0; count < elements_list[0]; count++)  
						if (elements_list[count*2+2] == count_length && ! strncmp(one_line, id_NODES+elements_list[count*2+1], count_length)) break;
					
					swap(pchar[0], temp); material_name = pchar+strlen(", material=");
					if ( (ppos = strstr(material_name, "\xD")) != NULL || (ppos = strstr(material_name, "\xA")) != NULL) ppos[0] = temp;
					count_length = strlen(material_name);

					if (count < elements_list[0]) {
						if (! material_first)
						fprintf(PRB, "//................................................\n");

//////////////////////////////////////////////
//...���������� �������� � ���������� ��� ���;
						ppos_material = ppos_cur;
						PPOS_CUR(id_NODES, pchar, ppos_material, upper_limit, "*Material, name=");

						while ((upper_material = ppos_material) < upper_limit) {
							PPOS_CUR(id_NODES, pchar, upper_material, upper_limit, "*Material, name=");  upper_material -= strlen("*Material, name=");
							ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE); ppos_save = ppos_material;

							if ((ppos = strstr(material_line, "\xD")) != NULL || (ppos = strstr(material_line, "\xA")) != NULL) ppos[0] = temp;
							if (strlen(material_line) == count_length && ! strncmp(material_name, material_line, count_length)) {
								double Gamma = 0.;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Density");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										Gamma = user_strtod(material_line); 
									}
								}
								double E = 0., nju = 0.; int  model = 0, aniso = 0;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Elastic");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										E = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										nju = user_strtod(pchar+1); 
									}
									model = 1; 
								}
								double E1 = 0., E2 = 0., E3 = 0., nju1 = 0., nju2 = 0., nju3 = 0., G1 = 0., G2 = 0., G3 = 0.;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Elastic, type=ENGINEERING CONSTANTS");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										E1 = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										E2 = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										E3 = user_strtod(ppos+1); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										nju1 = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										nju2 = user_strtod(ppos+1); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										nju3 = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										G1 = user_strtod(ppos+1); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										G2 = user_strtod(pchar+1); 
									}
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										G3 = user_strtod(material_line); 
									}
									aniso = 2; 
								}
								double angle_frict = 0., flow_ratio = 0., angle_dilat = 0.;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Drucker Prager");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										angle_frict = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										flow_ratio = user_strtod(pchar+1); 
									}
									if ((pchar = strstr(ppos+1, ",")) != NULL || (pchar = strstr(ppos+1, "\xD")) != NULL || (pchar = strstr(ppos+1, "\xA")) != NULL) {
										angle_dilat = user_strtod(ppos+1); 
									}
								   model = 9;
								}
								double yield_stress = 0., plastic_strain = 0.;
								ppos_material = ppos_save;
								PPOS_CUR(id_NODES, pchar, ppos_material, upper_material, "*Drucker Prager Hardening");
								if (ppos_material < upper_material) {
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									ONE_LINE(id_NODES, pchar, ppos_material, upper_material, material_line, STR_SIZE);
									if ((pchar = strstr(material_line, ",")) != NULL || (pchar = strstr(material_line, "\xD")) != NULL || (pchar = strstr(material_line, "\xA")) != NULL) {
										yield_stress = user_strtod(material_line); 
									}
									if ((ppos = strstr(pchar+1, ",")) != NULL || (ppos = strstr(pchar+1, "\xD")) != NULL || (ppos = strstr(pchar+1, "\xA")) != NULL) {
										plastic_strain = user_strtod(pchar+1); 
									}
								}

//////////////////////////////
//...������ ������ ����������;
								fprintf(PRB, "%i	%i  %i     // ������: ���������, ������, �����������\n", (int)count+1, model, aniso);
								switch (model) {
								case 1: switch (aniso) {
									case 0: fprintf(PRB, "// E      v       Gamma:\n");
									fprintf(PRB, "%g     %g     %g\n", E, nju, Gamma);
									break;
									case 2: fprintf(PRB, "// E1      E2      E3      v1      v2      v3      G1      G2      G3       Gamma:\n");
									fprintf(PRB, "%g     %g     %g     %g     %g     %g     %g     %g     %g     %g\n", E1, E2, E3, nju1, nju2, nju3, 
										G1, G2, G3, Gamma);
									}
								break;
								case 9:  switch (aniso) {
									case 0: fprintf(PRB, "// E      v      C     fi    Gamma\n");
									fprintf(PRB, "%g     %g     %g     %g     %g\n", E, nju, yield_stress, angle_frict, Gamma);
									break;
									case 2: fprintf(PRB, "// E1      E2      E3      v1      v2      v3      G1      G2      G3      C     fi       Gamma:\n");
									fprintf(PRB, "%g     %g     %g     %g     %g     %g     %g     %g     %g     %g     %g     %g\n", E1, E2, E3, nju1, nju2, nju3, 
										G1, G2, G3, yield_stress, angle_frict, Gamma);
									}
								}
							}
							ppos_material = upper_material+strlen("*Material, name=");
						}
						fprintf(PRB, "%s		// �������� ���������\n<END>\n", material_name);
						material_first = 0;
					}
				}
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Solid Section, elset=");
			}
			fprintf(PRB, "0   //<<< ������� ����� ������\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

//////////////////////////////////////////
//...����� ���� ��������� (������ ������);
			fprintf(PRB, "[VARIABLE MATERIALS]\n");
			fprintf(PRB, "0   // ����/��� ���������� ��������� (\"0\" - ���)\n");
			fprintf(PRB, "0   //<<< ������� ����� ������\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");
				 
/////////////////////////////////////
//...����� � ������ �������� � �����;
			surfaces_list = Part_surfaces_list(id_NODES, 0, upper_limit, ID_part);
			fprintf(PRB, "[LOADS]\n"); sum_steps = 0.;
			PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Step, name=");
			while ((upper = ppos_cur) < upper_limit) {
				PPOS_CUR(id_NODES, pchar, upper,	upper_limit, "*End Step"); upper -= sizeof("*End Step");

////////////////////
//...���������� ���;
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
				ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);

				if ((ppos = strstr(one_line, ",")) != NULL) step = user_strtod(ppos+1);
				else													  step = 1.;
				count_beg = ppos_cur;

/////////////////////////////////////////////////
//...���� ��� ��������� �������� �������� � ����;
				while (ppos_cur < upper) {
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Dload");
					if (pchar != NULL) {
							ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
							ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
////////////////////////////////////////////////////
//...���������� ��� �������� �������� � �� ��������;
							pchar = strstr(one_line, ","); swap(temp, pchar[0]); 
							if (! ::strncmp (pchar+1, " GRAV", 5)) {
								double grav = 0., nX = 0., nY = -1., nZ = 0.;
								if ((ppos = strstr(pchar+1, ",")) != NULL) {
									grav = user_strtod(ppos+1);
									if ((ppos = strstr(ppos+1, ",")) != NULL) {
										nX	= user_strtod(ppos+1);
										if ((ppos = strstr(ppos+1, ",")) != NULL) {
											nY	= user_strtod(ppos+1);
											if ((ppos = strstr(ppos+1, ",")) != NULL) nZ = user_strtod(ppos+1);
										}
									}
								}
//////////////////////////////////
//...������ �������� - ����������;
								fprintf(PRB, "103   // ��� �������� - ����������� ���\n<���� �������>   // �������� ��������\n");
								fprintf(PRB, "0   // ������� �������� � �����: 0-RAM; 1-ASCII; 2-Binary\n");
								fprintf(PRB, "1   // �-��� ���������� �� ������� - f(t)\n// 1  �������-�������� �-���   (n+1), n - ���-�� �������\n// t - ��������� �������\n");
								fprintf(PRB, " %g    %g\n", sum_steps, sum_steps+step);
								fprintf(PRB, "<END>\n// ��������� �-��� ���������� (���������� ��������) �� �������, ������� ����������\n// ������ (����) �� ������������ �� ����������� �������� �������������� ��������\n");
								fprintf(PRB, " 0.0    1.0\n");
								fprintf(PRB, "<END>\n");
								fprintf(PRB, "1   // �-��� ���������� �� ���������������� ��������� - f(x)\n// 1  ���������� �������������  (1)\n// p - ��������� ��������\n");
								fprintf(PRB, " %g\n", grav);
								fprintf(PRB, "<END>\n");
								fprintf(PRB, " %g %g %g			// �����������\n", nX, nY, nZ);
								fprintf(PRB, "1.0				// ���-�� �������� �� ��������� ������\n");
								fprintf(PRB, "0   // \"+1\" - ��� �������� ..., \"-1\" - �������� ������� ..., \"0\" - ��� ���� ��������\n");
								fprintf(PRB, "//..............................................................................\n");
							}
							swap(temp, pchar[0]); 
					}
				}
				ppos_cur = count_beg;

//////////////////////////////////////////////////////
//...���� ��� ��������� ������������� �������� � ����;
				while (ppos_cur < upper) {
					PPOS_CUR(id_NODES, pchar, ppos_cur, upper, "*Dsload");
					if (pchar != NULL) {
							ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
							ONE_LINE(id_NODES, pchar, ppos_cur, upper, one_line, STR_SIZE);
/////////////////////////////////////////////////////////
//...���������� ��� ������������� �������� � �� ��������;
							pchar = strstr(one_line, ","); swap(temp, pchar[0]); 
							if (! ::strncmp (pchar+1, " TRVEC", 6) || ! ::strncmp (pchar+1, " TRSHR", 6)) {
								double value = 0., nX = 0., nY = -1., nZ = 0.;
								if ((ppos = strstr(pchar+1, ",")) != NULL) {
									value = user_strtod(ppos+1);
									if ((ppos = strstr(ppos+1, ",")) != NULL) {
										nX	= user_strtod(ppos+1);
										if ((ppos = strstr(ppos+1, ",")) != NULL) {
											nY	= user_strtod(ppos+1);
											if ((ppos = strstr(ppos+1, ",")) != NULL) nZ = user_strtod(ppos+1);
										}
									}
								}
//////////////////////////////////////////
//...������ �������� - ������������� ����;
								fprintf(PRB, "107 // ��� �������� - �������������� �� ����������� (�����) ��������\n<������������  �������������� �� ����������� ��������>   // �������� ��������\n");
								fprintf(PRB, "0   // ������� �������� � �����: 0-RAM; 1-ASCII; 2-Binary\n");
								fprintf(PRB, "1   // �-��� ���������� �� ������� - f(t)\n// 1  �������-�������� �-���   (n+1), n - ���-�� �������\n// t - ��������� �������\n");
								fprintf(PRB, " %g    %g\n", sum_steps, sum_steps+step);
								fprintf(PRB, "<END>\n// ��������� �-��� ���������� (���������� ��������) �� �������, ������� ����������\n// ������ (����) �� ������������ �� ����������� �������� �������������� ��������\n");
								fprintf(PRB, " 0.0    1.0\n");
								fprintf(PRB, "<END>\n");
								fprintf(PRB, "1   // �-��� ���������� �� ���������������� ��������� - f(x)\n// 1  ����������                 (1)  p = const\n// p - ��������� ��������\n");
								fprintf(PRB, " %g\n", value);
								fprintf(PRB, "<END>\n");
								fprintf(PRB, " %g %g %g			// �����������\n", nX, nY, nZ);
								fprintf(PRB, "1.0				// ���-�� �������� �� ��������� ������\n");
								fprintf(PRB, "1   // \"+1\" - ��� �������� ..., \"-1\" - �������� ������� ..., \"0\" - ��� ���� ��������\n");
								fprintf(PRB, "// �� ����������� ��������:\n");
								fprintf(PRB, "//_____1_____2_____3_____4_____5_____6_____7_____8_____9____10:");
////////////////////////////////////////////////////////////////////////////
//...��������� ����� ��������� � �������� ������ ����� ���������� ���������;
								ppos = ::strrchr(one_line, '.')+1; // - ��� �����������, �� ������� ������ ��������� �������;
								count_length = strlen(ppos);

								for (count = 0; count < surfaces_list[0]; count++)  
								if ( count_length < surfaces_list[count*3+2] && ! strncmp(ppos, id_NODES+surfaces_list[count*3+1]+1, count_length)) {
									k_side = surfaces_list[count*3+3];
////////////////////////////////////////////////
//...������������ ����������� ����� �����������;
									for (count_el = 0; count_el < elements_list[0]; count_el++)  
									if ( surfaces_list[count*3+2] == elements_list[count_el*2+2] && ! strncmp(id_NODES+surfaces_list[count*3+1], id_NODES+elements_list[count_el*2+1], surfaces_list[count*3+2])) {
										k = (int)(count_el+nodes_list[0]);
////////////////////////////////////////////////////////////////////////////////////////////////////
//...������������� ��� ��������, ������������ �����������, � �������� ���� ��� ������ (������ ����);
										if (nd->cond_ptr && nd->cond) {
											l = nd->cond_ptr[k];
											if (nd->cond[l] == (int)GL_ELEMENTS) {
												for (j = 2; j <= nd->cond[l+1]; j++)
													SetLoading_nodes(nd, nd->cond[l+1+j], k_side);
											}
											else
											if (nd->cond[l] == (int)GL_ELEMENTS_GENERATE) {
												for (j = nd->cond[l+3]; j <= nd->cond[l+4]; j += nd->cond[l+5]) 
													SetLoading_nodes(nd, j, k_side);
											}
										}
									}
								}
///////////////////////////////////////
//...�������� ������� ���������� �����;
								for (j_node = 0, k = 0; k < nd->N; k++)
								if (nd->hit[k] < 0) {
									if (! (j_node % j_max))	fprintf(PRB, "\n ");
									fprintf(PRB, " %7i", -nd->hit[k]); j_node++;
								}
								fprintf(PRB, "\n");
								fprintf(PRB, "//..............................................................................\n");
								for (j_node = 0, k = 0; k < nd->N; k++) nd->hit[k] = abs(nd->hit[k]); //...���������� �������;
							}
							swap(temp, pchar[0]); 
					}
				}
				ppos_cur = count_beg;

////////////////////////////////////////////////////////////////
//...���� ��� ��������� ������ ������������ ��� �������� � ����;
//...

/////////////////////////////////
//...��������� � ���������� ����;
				sum_steps += step;
				PPOS_CUR(id_NODES, pchar, ppos_cur, upper_limit, "*Step, name=");
			}
			fprintf(PRB, "0		//<<<<<<<<<<<<<<<<<<<<<<< the end of data\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");

/////////////////
//...����� �����;
			fprintf(PRB, "[MESH CHANGE]\n");
			fprintf(PRB, "0   // 0 - ����� ���������\n");
			fprintf(PRB, "//------------------------------------------------------------------------------\n");
			fclose (PRB);
		}
		delete_struct(nodes_list);
		delete_struct(elements_list);
		delete_struct(surfaces_list);
	}
	delete_struct(id_NODES);
}
