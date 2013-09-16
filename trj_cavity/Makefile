CC = g++

OBJS =		trj_cavity.o

TARGET =	trj_cavity

GROMACS =	/home/tp334/software/gromacs-4.5.3

CFLAGS =	-Wall -I$(GROMACS)/include/gromacs -g -O1 -B/usr/bin

LDFLAGS=	-L$(GROMACS)/lib -lgmx

SOURCE =	src/

LIBS =		-lm -pthread -ldl

USRLIBS =	$(SOURCE)Atom.cpp $(SOURCE)AtomReader.cpp $(SOURCE)AtomWriter.cpp $(SOURCE)Coordinates.cpp $(SOURCE)Grid.cpp $(SOURCE)ForceField.cpp


$(TARGET):	$(OBJS)
	$(CC) $(USRLIBS) -o $(TARGET) $(OBJS) $(CFLAGS) $(LIBS) $(LDFLAGS) 
	

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
