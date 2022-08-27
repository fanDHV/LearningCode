Prerequirement packages:
Grphics:
  libX11-dev  (might be installed)
  libXi-dev
  libXmu-dev
  libgmp-dev  (might be installed)
  libmpfr-dev (might be installed)

CGAL:
  libcgal-dev
  libcgal8
  libcgal-demo
  liboost_system

OpenGL:
  freeglut3-dev
  libgl-dev
  libgl1
  libgl1-mesa-dev
  libqt4-opengl
  libqt4-opengl-dev	

1. run:
./sgp inputs/25.off 0.1 inputs/25.in 

2.OpenGL direction (keys):
  z+w (or only w)
  z+a (or only a)
  z+s (or only s)
  z+d (or only d)
  
3. select soure and destination face, p,q then are selected on the
the faces.

##############################################################################
- 2015/08/25: back to work on it again for comparing with a global solution.
This is note for using Qt Creator.
  1. create a new project 
  2. add necessary files into the project
  3. edit the .pro file, especially for required libraries, add the following
     ESOURCES +=
     QMAKE_CXXFLAGS += -frounding-math
     LIBS += -lCGAL -lCGAL_Core\
        -lX11 -lXi -lXmu -lglut -lGL -lGLU\
        -lgmp -lmpfr\
        -lboost_system

     OTHER_FILES += \
     Makefile
  4. in termial, type 'qmake' to create the makefile
  5. finally is compile and run

######################################################################
Le
1. qmake -project           # create "project".pro
2. edit the .pro file, especially for required libraries, add the following
     ESOURCES +=
     QMAKE_CXXFLAGS += -frounding-math
     LIBS += -lCGAL -lCGAL_Core\
        -lX11 -lXi -lXmu -lglut -lGL -lGLU\
        -lgmp -lmpfr\
        -lboost_system

     OTHER_FILES += \
     Makefile
3. qmake -o Makefile sgp.pro 
4. make
5. ./sgp inputs/*.off 0.1 inputs/*.in 
EX:  ./sgp data-ok/terrain-sphere100.off 0.05 data-ok/terrain-sphere100.in
 ./sgp data-ok/terrain-sphere100.off 0.05 data-ok/terrain-sphere100_1.in


###################
Note: To run and take  
File main.cpp: l102:
bool f_super_test = true; 
files .off, .in put in the folder containing main.cpp



