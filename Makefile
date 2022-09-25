.PHONY : all
CXXFLAGS = -O3 --std=c++17 -static-libgcc -static-libstdc++ -idirafter /usr/include -idirafter /usr/include/eigen3
LDFLAGS = -Wl,-allow-multiple-definition -lglew -lglfw3 -lgdi32 -lopengl32
CXX = x86_64-w64-mingw32-g++

all: freeflow.exe
freeflow.exe: freeflow.cc
	$(CXX) $(CXXFLAGS) freeflow.cc $(LDFLAGS) -o freeflow.exe