grid: freeflow.cc SignedDistanceFieldRenderer.hpp GridRenderer.hpp shader.hpp VectorFieldRenderer.hpp
	clang++-mp-14 -fopenmp -O3 --std=c++20 -o grid \
	-I/opt/local/include \
	-I/opt/local/include/opencv4 \
	freeflow.cc \
	 /opt/local/lib/libglfw.dylib \
	 /opt/local/lib/libGLEW.dylib \
	 -framework OpenGL \
	 -framework OpenCL