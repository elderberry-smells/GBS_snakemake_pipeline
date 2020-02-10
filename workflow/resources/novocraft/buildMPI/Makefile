
#ifndef OSTYPE
OSTYPE=$(shell uname -s)
#endif
OS=$(shell cat /etc/issue | cut -f 1 "--delimiter= " | head -1)

HDF5=-lhdf5 -lhdf5_cpp 
HDF5INC=/usr/include/hdf5/serial
HDF5LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial

HDF5= -lhdf5_cpp -lhdf5
HDF5INC=/usr/include/hdf5
HDF5LIB=/usr/lib/x86_64-linux-gnu/hdf5
ifneq ($(wildcard /usr/include/hdf5/serial), )
	HDF5INC=/usr/include/hdf5/serial
	HDF5LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial
endif


ifdef CRAY_MPICH2
	INCLUDE=-I ${MPICH_DIR}/include
	CFLAGS=-m64 -c -O3 -DCRAY_MPICH2
	MPICC=g++
	LDLIBSOPTIONS=-pthread -lcrypto -lz -lbz2 -ldl -Wl,-rpath,/opt/gcc/4.8.1/snos/lib64 -l mpich -L ${MPICH_DIR}/lib
	MAP=-Wl,-Map=
else
	MPICC=mpicxx
	ifeq ($(OSTYPE),Darwin)
		MACHTYPE=$(shell echo $MACHTYPE)
		ifeq ($(MACHTYPE), x86_64-apple-darwin10.0)
			CFLAGS=-arch x86_64 -m64 -c -O3 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5 -I. -I/usr/local/include
			LNFLAGS=-m64 -mmacosx-version-min=10.5  -isysroot /Developer/SDKs/MacOSX10.5.sdk 
		else
		    CFLAGS=-arch x86_64 -m64 -c -O3 -I. -I/usr/local/include
		    LNFLAGS=-m64
		endif
		MAP=-Wl,-map 
		LDLIBSOPTIONS= -lpthread -lcrypto -lz -lbz2
	else
		ifeq ($(OSTYPE),SunOS)
			CFLAGS=-m64 -c -O3 -I/usr/sfw/include
			LNFLAGS=-m64
			MAP=-Wl,-m >
			LDLIBSOPTIONS= -lpthread -lrt -lz -lbz2 -L /usr/sfw/lib/amd64 -R /usr/sfw/lib/amd64 -lcrypto -ldl
		else
			LNFLAGS=-m64 -Wl,--eh-frame-hdr
			MAP=-Wl,-Map=
			LDLIBSOPTIONS= ./lib/libtcmalloc.a ./lib/libunwind.a -pthread -lcrypto -lz -lbz2 -ldl -llzma
			ifeq ($(OS),Scientific)
				CFLAGS=-m64 -c -O3 -I./unwind -I/usr/lib64 
			else
				CFLAGS=-m64 -c -O3 -I./unwind
			endif
		endif
	endif
endif


# Build Targets
all: novoalignMPI novoalignCSMPI

novoalignMPI:	./lib/novoalignMPI.a mpidriver.o trace.o
	echo >novoalignMPI.map # for Mac ld as map file must exist!!!
	${MPICC}  ${LNFLAGS} -o novoalignMPI mpidriver.o trace.o ./lib/novoalignMPI.a ${LDLIBSOPTIONS} ${MAP}novoalignMPI.map 

xsqreader.o: xsqreader.cc 
	${MPICC}  ${CFLAGS} ${INCLUDE} -I${HDF5INC} -o $@ $?

%.o: %.cc 
	${MPICC}  ${CFLAGS} ${INCLUDE} -o $@ $?

novoalignCSMPI: ./lib/novoalignCSMPI.a mpidriver.o trace.o xsqreader.o
	echo >novoalignCSMPI.map # for Mac ld as map file must exist!!!
	${MPICC}  ${LNFLAGS} -o novoalignCSMPI  mpidriver.o trace.o xsqreader.o ./lib/novoalignCSMPI.a ${HDF5} -L${HDF5LIB} ${LDLIBSOPTIONS} ${MAP}novoalignCSMPI.map 

# Clean Targets
clean:
	rm novoalignMPI novoalignCSMPI mpidriver.o trace.o xsqreader.o
