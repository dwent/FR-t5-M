# Makefile for Assembly Protocol


OBJS = \
	Main.o	\
	HelixAxis.o	\
	MathTool.o	\
	HydroCluster.o	\
	HydroFace.o	\
	IOprotein.o	\
	CoreStruct.o	\
		
	

GG = g++ -c -O2 -g

SHCore : $(OBJS)
	g++ -O2 -o HCsta $(OBJS) -g

Main.o : Main.cpp
	$(GG) Main.cpp

HelixAxis.o : HelixAxis.cpp HelixAxis.hpp MathTool.hpp
	$(GG) HelixAxis.cpp
	
MathTool.o : MathTool.cpp MathTool.hpp
	$(GG) MathTool.cpp

HydroCluster.o :HydroCluster.cpp HydroCluster.hpp MathTool.hpp HydroFace.hpp IOprotein.hpp
	$(GG) HydroCluster.cpp

HydroFace.o :HydroFace.cpp HydroFace.hpp HelixAxis.hpp
	$(GG) HydroFace.cpp

IOprotein.o : IOprotein.cpp IOprotein.hpp
	$(GG) IOprotein.cpp

CoreStruct.o : CoreStruct.cpp CoreStruct.hpp HydroFace.hpp HydroCluster.hpp
	$(GG) CoreStruct.cpp



.PHONY clean:
clean :
	-rm HCsta $(OBJS)

inobj :
	mv *.o obj

outobj :
	mv obj/*.o .

