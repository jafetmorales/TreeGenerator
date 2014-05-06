CC= g++
EXEC= trees
OBJS= main.o treegenerator.o treerenderer.o 
LIBS+= -L/usr/local/lib/ -lpng -ljpeg -ltiff -lGLEW -lz -framework OpenGL -framework GLUT
INCL+= -I/usr/local/include/
FLAGS += -g -Wno-deprecated-declarations

$(EXEC) : $(OBJS)
	$(CC) -g $(OBJS) $(LIBS) -o $(EXEC)

%.o : %.cpp
	$(CC) $(FLAGS) $(INCL) -c $<