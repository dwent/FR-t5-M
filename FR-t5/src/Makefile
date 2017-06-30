# Makefile of Sequence-Structure-Threading-Alignment
#				/**********************************************************
#				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
#				All rights reserved.
#				NAME:		Makefile
#				ABSTRACT:	MAKEFILE COMPILING
#				VERSION:	0.01
#				AUTHOR:	Miao Zhichao
#				DATE:		2012.03.19 Mon
#				CONTACT:	chichaumiau@gmail.com
#
#				NOTICE: This is free software and the source code is freely
#				available. You are free to redistribute or modify under the
#				conditions that (1) this notice is not removed or modified
#				in any way and (2) any modified versions of the program are
#				also available for free.
#						** Absolutely no Warranty **
#				***********************************************************/
#
#USAGE: (COMPLIE)./make
#	 (CLEAN)  ./make clean
#
#THIS NAME CAN BE CHOSEN BETWEEN MIAOZHICHAO AND HUYUN
#USE THE MIAOZHICAO FOR NAME IT WILL COMPILE THE NEW JTHREADER WORKFLOW
#CHANGE IT TO HUYUN, IT'D COMPILE THE OLD FR-T5 WORKFLOW
NAME=MIAOZHICHAO
#OBJECT FILES
OBJS = \
	main.o	\
	Benchmark.o	\
	Smith-Waterman.o	\
	BasicDef.o	\
	IOprotein.o	\
	MathTool.o	\
	Matrices.o	\
	Options.o	\
	PrintArray.o	\
	ScoreMatrix.o	\
	TraceArrow.o	\
	Trace.o	\
	GetSelfSM.o	\
	FiveBeadEnergy.o	\
#COMPILE PARAMETERS
GG = g++ -c -g  -O3 -D $(NAME) -Wall
#EXECUTIVE, LINK ALL OBJECT FILES
binary : $(OBJS)
	g++ -o ../bin/JTHREADER_$(NAME) $(OBJS)

#DEPENDENCIES:
main.o : main.cpp Benchmark.hpp Smith-Waterman.hpp  FiveBeadEnergy.hpp
	$(GG) main.cpp

Benchmark.o : Benchmark.cpp Benchmark.hpp Smith-Waterman.hpp 
	$(GG) Benchmark.cpp

Smith-Waterman.o : Smith-Waterman.cpp Smith-Waterman.hpp Matrices.hpp Trace.hpp TraceArrow.hpp \
			Options.hpp GetSelfSM.hpp FiveBeadEnergy.hpp MathTool.hpp
	$(GG) Smith-Waterman.cpp

BasicDef.o : BasicDef.cpp BasicDef.hpp
	$(GG) BasicDef.cpp

GetSelfSM.o : GetSelfSM.cpp GetSelfSM.hpp MathTool.hpp IOprotein.hpp FiveBeadEnergy.hpp
	$(GG) GetSelfSM.cpp

IOprotein.o : IOprotein.cpp IOprotein.hpp
	$(GG) IOprotein.cpp

MathTool.o : MathTool.cpp MathTool.hpp
	$(GG) MathTool.cpp

Matrices.o : Matrices.cpp Matrices.hpp BasicDef.hpp Matrix.hpp PrintArray.hpp Sequence.hpp
	$(GG) Matrices.cpp

Options.o : Options.cpp Options.hpp
	$(GG) Options.cpp

PrintArray.o : PrintArray.cpp PrintArray.hpp BasicDef.hpp Matrix.hpp
	$(GG) PrintArray.cpp

ScoreMatrix.o : ScoreMatrix.cpp ScoreMatrix.hpp Matrix.hpp BasicDef.hpp 
	$(GG) ScoreMatrix.cpp

TraceArrow.o : TraceArrow.cpp TraceArrow.hpp BasicDef.hpp PrintArray.hpp Trace.hpp ScoreMatrix.hpp
	$(GG) TraceArrow.cpp

Trace.o : Trace.cpp Trace.hpp Matrix.hpp BasicDef.hpp PrintArray.hpp Sequence.hpp Matrices.hpp
	$(GG) Trace.cpp

FiveBeadEnergy.o : FiveBeadEnergy.cpp FiveBeadEnergy.hpp MathTool.hpp IOprotein.hpp 
	$(GG) FiveBeadEnergy.cpp

#HOW TO CLEAN
.PHONY clean:
clean :
	-rm $(OBJS) ../bin/JTHREADER*

inobj :
	mv *.o obj

outobj :
	mv obj/*.o .

