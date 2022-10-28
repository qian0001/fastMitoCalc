#BaseCoverage:
#	g++ BaseCoverage.cpp -o BaseCoverage

BaseCoverage: BaseCoverage.o
	g++ BaseCoverage.o -o BaseCoverage

BaseCoverage.o: BaseCoverage.cpp
	g++ -c BaseCoverage.cpp -o BaseCoverage.o
clean:
	rm -f BaseCoverage.o BaseCoverage
