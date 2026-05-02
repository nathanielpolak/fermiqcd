#MYPATH=$(HOME)/TP2/TP_Linux
MYPATH=$(HOME)/home/Teaching/DEA-Orsay/TP2/TP_Linux

FERMIQCDLIB = $(MYPATH)/FermiQCD/FermiQCD/Version_3.3/Libraries

GRACELIB = $(MYPATH)/Grace

MYLIB = $(MYPATH)/UseFermiQCD/Include



all: Autocorrel.o grace_np.o
	g++ -o a.out Autocorrel.o grace_np.o

poly: PolyakovStudy.o grace_np.o
	g++ -o a.out PolyakovStudy.o grace_np.o

grace_np.o: $(GRACELIB)/grace_np.c
	gcc -c -o grace_np.o -I $(GRACELIB) $(GRACELIB)/grace_np.c

Autocorrel.o: Autocorrel.cpp
	g++ -O3  -c -o Autocorrel.o -I $(FERMIQCDLIB) -I $(MYLIB) -I $(GRACELIB) Autocorrel.cpp

PolyakovStudy.o: PolyakovStudy.cpp
	g++ -O3  -c -o PolyakovStudy.o -I $(FERMIQCDLIB) -I $(MYLIB) -I $(GRACELIB)  PolyakovStudy.cpp

clean:
	rm  a.out Autocorrel.o grace_np.o
