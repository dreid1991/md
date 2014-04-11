FLAGS = -std=c++11

H = $(wildcard *.h)
C = $(wildcard *.cpp)
app: ../blank.cpp $(H) $(C) 
	g++ $(FLAGS) -Wall -I/usr/include/python2.7 ../main.cpp $(H) $(C) -lpython2.7


