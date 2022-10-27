grid: freeflow.cc SignedDistanceFieldRenderer.hpp GridRenderer.hpp shader.hpp VectorFieldRenderer.hpp
	x86_64-w64-mingw32-g++ -O3 --std=c++17 -o freeflow.exe \
	-idirafter /usr/include freeflow.cc \
	-static-libgcc -static-libstdc++ \
	-Wl,-allow-multiple-definition -lglew -lglfw3 -lgdi32 -lopengl32