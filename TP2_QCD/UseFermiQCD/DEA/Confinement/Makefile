
MYPATH=$(HOME)/TP2/TP_Linux

FERMIQCDLIB = $(MYPATH)/FermiQCD/FermiQCD/Version_3.3/Libraries

GRACELIB = $(MYPATH)/Grace

MYLIB = $(MYPATH)/UseFermiQCD/Include


all: Confine.o grace_np.o CPULoop.o
	g++  -o a.out Confine.o grace_np.o CPULoop.o

nograce: NoGraceConfine.o grace_np.o CPULoop.o
	g++  -o a.out NoGraceConfine.o grace_np.o CPULoop.o

CPULoop.o:  CPULoop.c
	 g++  -c -o CPULoop.o -I include CPULoop.c

grace_np.o: $(GRACELIB)/grace_np.c
	gcc -c -o grace_np.o -I $(GRACELIB) $(GRACELIB)/grace_np.c

Confine.o: Confine.cpp
	g++ -O3  -c -o Confine.o -I $(FERMIQCDLIB) -I $(MYLIB) -I include -I $(GRACELIB) Confine.cpp

NoGraceConfine.o: NoGraceConfine.cpp
	g++ -O3  -c -o NoGraceConfine.o -I $(FERMIQCDLIB) -I $(MYLIB) -I include -I $(GRACELIB) NoGraceConfine.cpp

clean:
	rm  a.out Confine.o CPULoop.o grace_np.o NoGraceConfine.o