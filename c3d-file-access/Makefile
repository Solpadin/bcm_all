f       = f
c       = cpp
o       = o
a       = a
r       =

.SUFFIXES: .$f .$c .$o

CC      = mpicxx -O3 -c 
CF      = mpif77 -O3 -c 
LINK    =

LIB     = c3d-file-access.$a
LIBOBJS = \
	./Read_C3D_Channel.$o \
	./Read_C3D_Data.$o \
	./Read_C3D_Header.$o \
	./Read_C3D_Marker.$o \
#	./Read_C3D_Parameters.$o \

$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
	ranlib $(LIB)

.$f.$o:
	$(CF) $<

.$c.$o:
	$(CC) $<

clean:
	-rm *.o


