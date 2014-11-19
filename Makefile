


all:
	g++ *.cpp -std=c++11 -pthread -O3

test-run:
	./a.out --experiment=test --replicate=2 --LOD=lod --genome=genome --generations=50
