GDAL_CFLAGS := $(shell gdal-config --cflags)
GDAL_LIBS   := $(shell gdal-config --libs)

all: check_graph build_graph
	echo "Build complete."

check_graph: check_graph.cpp
	g++ -O3 -std=c++17 check_graph.cpp -o check_graph
	
build_graph: build_graph.cpp
	g++ -O3 -std=c++17 build_graph.cpp $(GDAL_CFLAGS) -lz -lbz2 $(GDAL_LIBS) -o build_graph

clean:
	rm -f check_graph build_graph
	rm -f *.o

test: check_graph build_graph
	rm -f ../out/alsace.graph ../out/alsace-highway.osm.pbf ../out/alsace_edge.fgb ../out/alsace_node.fgb ../out/alsace_restriction.fgb
	osmium tags-filter ../in/alsace-latest.osm.pbf w/highway r/type=restriction -o ../out/alsace-highway.osm.pbf
	./build_graph ../out/alsace-highway.osm.pbf ../out/alsace.graph --fgb
	./check_graph ../out/alsace.graph
	echo "Test complete."