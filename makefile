EXE=cratersim
SOURCES=Clarkson-Delaunay.c cratersim.c
OBJECTS=$(SOURCES:.c=.o)

# Main target
all: $(EXE)


# CFLG
# LIBS
# CLEAN

#  Msys/MinGW
ifeq "$(OS)" "Windows_NT"
CFLG=-O3 -Wall -DUSEGLEW
LIBS=-lfreeglut -lglew32 -lglu32 -lopengl32 -lm
CLEAN=rm -f *.exe *.o *.a
else
#  OSX
ifeq "$(shell uname)" "Darwin"
RES=$(shell uname -r|sed -E 's/(.).*/\1/'|tr 12 21)
CFLG=-O3 -Wall -Wno-deprecated-declarations -DRES=$(RES)
LIBS=-framework GLUT -framework OpenGL
#  Linux/Unix/Solaris
else
CFLG=-O3 -Wall
LIBS=-lglut -lGLU -lGL -lm
endif
#  OSX/Linux/Unix/Solaris
CLEAN=rm -f $(EXE) *.o *.a
endif




# $@ = all
# $< = first prerequisite. 
# $^ = all prerequisites. 



# Compile rules
.c.o: # build .o files from .c files
	gcc -c $(CFLG)  $^

.cpp.o:
	g++ -c $(CFLG)  $^

#  Link
$(EXE): $(OBJECTS)
	gcc $(CFLG) -o $@ $^ $(LIBS)



#  Clean
clean:
	$(CLEAN)
