CC = g++

OBJS =		trj_cavity.o

TARGET =	trj_cavity

CFLAGS =	-Wall -I/home/tp334/software/gromacs-4.5.3/include/gromacs -g -O1 -B/usr/bin

LDFLAGS=	-L/home/tp334/software/gromacs-4.5.3/lib -lgmx

SOURCE =	src/

LIBS =		-lm -pthread -ldl

USRLIBS =	$(SOURCE)Atom.cpp $(SOURCE)AtomReader.cpp $(SOURCE)AtomWriter.cpp $(SOURCE)Coordinates.cpp $(SOURCE)Grid.cpp


$(TARGET):	$(OBJS)
	$(CC) $(USRLIBS) -o $(TARGET) $(OBJS) $(CFLAGS) $(LIBS) $(LDFLAGS) 
	

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
