##################################################
# 2Dvector is also essential for s(q) calculation#
##################################################
### issure: make RDF.exe  to compile #########
CC = g++ -O3

clean:
	/bin/rm RDF.ext

oldfunction2.o: oldfunction.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c oldfunction.cxx -o oldfunction2.o

oldfunction3.o: oldfunction.cxx MapCorrelation.h
	$(CC) -DDimension=3 -c oldfunction.cxx -o oldfunction3.o

shearevent2.o: shearevent.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c shearevent.cxx -o shearevent2.o

shearevent3.o: shearevent.cxx MapCorrelation.h
	$(CC) -DDimension=3 -c shearevent.cxx -o shearevent3.o

grid2.o: grid.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c grid.cxx -o grid2.o

grid3.o: grid.cxx MapCorrelation.h
	$(CC) -DDimension=3 -c grid.cxx -o grid3.o

MapCorrelationIO2.o: MapCorrelationIO.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c MapCorrelationIO.cxx -o MapCorrelationIO2.o

MapCorrelationIO3.o: MapCorrelationIO.cxx MapCorrelation.h
	$(CC) -DDimension=3 -c MapCorrelationIO.cxx -o MapCorrelationIO3.o

localstructure2.o: localstructure.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c localstructure.cxx -o localstructure2.o

localstructure3.o: localstructure.cxx MapCorrelation.h
	$(CC) -DDimension=3 -c localstructure.cxx -o localstructure3.o

cluster2.o: cluster.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c cluster.cxx -o cluster2.o

cluster3.o: cluster.cxx MapCorrelation.h
	$(CC) -DDimension=3 -c cluster.cxx -o cluster3.o

cell3.o:	cell.cxx MapCorrelation.h
	$(CC) -DDimension=3 -c cell.cxx -o cell3.o

MapCorrelation3.o:	MapCorrelation.cxx  MapCorrelation.h
	$(CC) -DDimension=3 -c MapCorrelation.cxx -o MapCorrelation3.o

RDF.exe:       md.analysis.wCell_latest_version.cxx	
	$(CC) -o RDF.exe md.analysis.wCell_latest_version.cxx

MapCorrelation3:	cell3.o MapCorrelation3.o potential.h cluster3.o MapCorrelation.h localstructure3.o grid3.o shearevent3.o MapCorrelationIO3.o ran2.o oldfunction3.o
	$(CC) -DDimension=3 -lm -o MapCorrelation3 ran2.cxx cell3.o cluster3.o MapCorrelation3.o localstructure3.o grid3.o shearevent3.o MapCorrelationIO3.o oldfunction3.o

cell2.o:	cell.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c cell.cxx -o cell2.o

ran2.o: ran2.cxx
	$(CC) -c ran2.cxx -o ran2.o

MapCorrelation2.o:	MapCorrelation.cxx MapCorrelation.h
	$(CC) -DDimension=2 -c MapCorrelation.cxx -o MapCorrelation2.o

MapCorrelation2:	cell2.o MapCorrelation2.o potential.h cluster2.o MapCorrelation.h localstructure2.o MapCorrelationIO2.o grid2.o shearevent2.o ran2.o oldfunction2.o
	$(CC) -DDimension=2 -o MapCorrelation2 ran2.o cell2.o MapCorrelation2.o cluster2.o localstructure2.o MapCorrelationIO2.o grid2.o shearevent2.o oldfunction2.o -lm

