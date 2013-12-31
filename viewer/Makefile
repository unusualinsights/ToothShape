CC=g++

all: tview


SRCS = tview.cpp plylib.cpp viewer.cpp modelseg.cpp
OBJS = $(SRCS:.cpp=.o)

INCLUDE = -I/usr/X11R6/include

XLIBS = -L/lusr/X11/lib -lXi -lXmu -lX11 -lm
CFLAGS = -O2

GLUT_LIB_LOCATION=/usr/X11R6/lib
GLUT_INC_LOCATION=/usr/X11R6/include
GL_LIBS=-L${GLUT_LIB_LOCATION} -lglut -lGL -lGLU


tview: $(OBJS)
	$(CC) $(CFLAGS) -I$(INCLUDE) $(OBJS) $(GL_LIBS) $(XLIBS) -o tview 

tview.o: tview.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) tview.cpp
plylib.o: plylib.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) plylib.cpp
viewer.o: viewer.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) viewer.cpp
modelseg.o: modelseg.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) modelseg.cpp


clean:
	rm -f *.o core.*
	rm -f tview

