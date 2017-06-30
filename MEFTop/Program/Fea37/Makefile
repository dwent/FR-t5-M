# Makefile for Generating Feature


OBJS = \
	main.o	\
	ReadInfor.o	\
	GenSSE.o	\
	CalFeature.o	\
	IOprotein.o	\
	MathTool.o	\
	NCACO.o	\


GG = g++ -c -O3 -g -rdynamic

GetFea : $(OBJS)
	g++ -O3 -o GenFea $(OBJS) -g -rdynamic
	
all : $(OBJS)
	g++ -O3 -o GenFea $(OBJS) -g -rdynamic

main.o : main.cpp NCACO.hpp CalFeature.hpp
	$(GG) main.cpp

GenSSE.o : GenSSE.cpp GenSSE.hpp
	$(GG) GenSSE.cpp

ReadInfor.o : ReadInfor.cpp ReadInfor.hpp GenSSE.hpp
	$(GG) ReadInfor.cpp

CalFeature.o : CalFeature.cpp CalFeature.hpp ReadInfor.hpp
	$(GG) CalFeature.cpp

IOprotein.o : IOprotein.cpp IOprotein.hpp
	$(GG) IOprotein.cpp

MathTool.o : MathTool.cpp MathTool.hpp
	$(GG) MathTool.cpp

NCACO.o : NCACO.cpp NCACO.hpp MathTool.hpp IOprotein.hpp
	$(GG) NCACO.cpp



.PHONY clean:
clean :
	-rm GenFea $(OBJS)

inobj :
	mv *.o obj

outobj :
	mv obj/*.o .


