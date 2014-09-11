#
# Linux/UNIX Makefile for pofad version 1.07
#

# Objects
OBJS = information.o nexusdata.o organismsdata.o GenFunctions.o ReadNexusFile.o pofad.o

# Flags
CFLAGS = -O3

# Compiler
CC = g++ $(CFLAGS)

#Folder containning binairies of the program
EXEDIR = ../bin

a1b2c3d4:
	@echo " "
	@echo " "
	@echo "  To use the POFAD version 1.07 Makefile, type"
	@echo " "
	@echo "     make install       to compile POFAD vers. 1.06 and then install the "
	@echo "                          executables in $(EXEDIR), then remove the object"
	@echo "                          files to save space"
	@echo "     make all           to compile the whole package and leaving the executables"
	@echo "                          in the current directoty and leaving the object files"
	@echo "     make put           to move the executables into $(EXEDIR)"
	@echo "     make clean         to remove all object files and executables from the"
	@echo "                          current directory"
	@echo "     make pofad         to compile pofad and leaving the executables in the"
	@echo "                          current directoty and leaving the object files"
	@echo " "

introduce:
	@echo ""
	@echo "Building POFAD version 1.07"

all:        introduce pofad
	@echo ""
	@echo "Finished compiling"

install:        all put clean
	@echo "Done compiling POFAD."
	@echo ""

put:
	@echo "Moving POFAD version 1.07 binaries in $(EXEDIR)"
	@mkdir -p $(EXEDIR)
	@cp pofad $(EXEDIR)

clean:
	@echo "Removing object files to save space"
	@rm -f *.o
	@rm -f pofad


pofad: $(OBJS)
	$(CC) $(OBJS) -o pofad -lm

pofad.o: pofad.h 
	$(CC) -c pofad.cpp

ReadNexusFile.o: ReadNexusFile.h
	$(CC) -c ReadNexusFile.cpp

organismsdata.o: organismsdata.h
	$(CC) -c organismsdata.cpp

nexusdata.o: nexusdata.h
	$(CC) -c nexusdata.cpp

information.o: information.h
	$(CC) -c information.cpp

GenFunctions.o: GenFunctions.h
	$(CC) -c GenFunctions.cpp

