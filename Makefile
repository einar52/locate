
CFLAGS=-g

testL : locS
	locS -i 018.event 
testF : fitTT
	fitTT -m 6 -i mdh3.4.data
f=fitTT.o golubc.o travelTime.c

fitTT : $f
	cc -g -o fitTT $f -lm 

l=locate.o golubc.o travelTime.o sgeom.o
locS : $l
	cc -g -o locS $l -lm
install : locS
	cp locS ~eik/bin

clean :
	rm -f core* *.o
