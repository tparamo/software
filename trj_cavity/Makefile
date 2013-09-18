CC = g++

OBJS =		trj_cavity.o

TARGET =	trj_cavity

#Version 4.5.x
#LDFLAGS=	-L/home/mary/software/gromacs-4.5.3/lib -lgmx
#CFLAGS =	-Wall -I/home/mary/software/gromacs-4.5.3/include/gromacs -g -O1 -B/usr/bin

#Version 4.6.x
#LDFLAGS=	-L/home/mary/software/gromacs-4.6.3/build/src/gmxlib -lgmx
#CFLAGS =	-Wall -I/mary/tp334/software/gromacs-4.6.3/build/include/gromacs -g -O1 -B/usr/bin -w

#Universal
LDFLAGS=	-L$(GMXLDLIB) -lgmx
CFLAGS =	-Wall -w -I$(GMXLDLIB)/../include/gromacs -g -O1 -B/usr/bin

SOURCE =	src/

LIBS =		-lm -pthread -ldl

USRLIBS =	$(SOURCE)Atom.cpp $(SOURCE)AtomReader.cpp $(SOURCE)AtomWriter.cpp $(SOURCE)Coordinates.cpp $(SOURCE)Grid.cpp $(SOURCE)ForceField.cpp


$(TARGET):	$(OBJS)
	$(CC) $(USRLIBS) -o $(TARGET) $(OBJS) $(CFLAGS) $(LIBS) $(LDFLAGS) 
	
	

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
