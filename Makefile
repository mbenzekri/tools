all: check_graph build_graph
	echo "Build complete."

check_graph: check_graph.cpp
	g++ -O3 -std=c++17 check_graph.cpp -o check_graph
	
build_graph: build_graph.cpp
	g++ -O3 -std=c++17 build_graph.cpp -lz -lbz2 -o build_graph

clean:
	rm -f check_graph build_graph
	rm -f *.o

