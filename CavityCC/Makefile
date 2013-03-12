CC = g++

OBJS =		CavityCC.o

TARGET =	CavityCC

CFLAGS =	-Wall -I/home/tp334/software/gromacs-4.5.3/include/gromacs -g -O0 -B/usr/bin

LDFLAGS=	-L/home/tp334/software/gromacs-4.5.3/lib -lgmx

SOURCE =	src/

LIBS =		-lm -pthread -ldl

USRLIBS =	$(SOURCE)Atom.cpp $(SOURCE)AtomReader.cpp $(SOURCE)AtomWriter.cpp $(SOURCE)Coordinates.cpp $(SOURCE)Grid.cpp $(SOURCE)Log.cpp $(SOURCE)ForceField.cpp $(SOURCE)Cavity.cpp


$(TARGET):	$(OBJS)
	$(CC) $(USRLIBS) -o $(TARGET) $(OBJS) $(CFLAGS) $(LIBS) $(LDFLAGS) 
	

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
