app: foo.cpp Vector.h Bounds.h Bounds.cpp Atom.h SquareOffsetStruct.h AtomGrid.h AtomGrid.cpp SphHarmTools.h SphHarmTools.cpp Grid.h ContourTools.h ContourTools.cpp PythonTools.cpp PythonTools.h Snap.h DataTools.h DataTools.cpp AtomTools.h AtomTools.cpp
	g++ -std=c++0x -ansi -Wall -I/usr/include/python2.7 foo.cpp Atom.h Vector.h Bounds.h Bounds.cpp SquareOffsetStruct.h AtomGrid.h AtomGrid.cpp SphHarmTools.h SphHarmTools.cpp Grid.h ContourTools.h ContourTools.cpp PythonTools.cpp PythonTools.h Snap.h DataTools.h DataTools.cpp AtomTools.h AtomTools.cpp -lpython2.7