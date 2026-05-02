MYPATH=$(HOME)/TP2/TP_Linux
#MYPATH=$(HOME)/home/Teaching/DEA-Orsay/TP2/TP_Linux


FERMIQCDLIB = $(MYPATH)/FermiQCD/FermiQCD/Version_3.3/Libraries

GRACELIB = $(MYPATH)/Grace

MYLIB = $(MYPATH)/UseFermiQCD/Include


all: Thermalise.o grace_np.o
	g++ -o a.out Thermalise.o grace_np.o

grace_np.o: $(GRACELIB)/grace_np.c
	gcc -c -o grace_np.o -I $(GRACELIB) $(GRACELIB)/grace_np.c

Thermalise.o: Thermalise.cpp
	g++ -O3  -c -o Thermalise.o -I $(FERMIQCDLIB) -I $(MYLIB) -I $(GRACELIB)  Thermalise.cpp

clean:
	rm  a.out Thermalise.o grace_np.o
