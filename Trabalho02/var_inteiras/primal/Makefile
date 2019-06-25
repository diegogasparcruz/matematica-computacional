CC = g++
CFLAGS1 = -c -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -DILOSTRICTPOD
CFLAGS2 =    -m64 -O -O0 -fPIC -fexceptions -DNDEBUG -DIL_STD -DILOSTRICTPOD -Wall
LIBS =  -lm 
INCLUDES = -I/opt/ibm/ILOG/CPLEX_Studio_Community129/cplex/include/ -I/opt/ibm/ILOG/CPLEX_Studio_Community129/concert/include/
OBJS = Primal.o
SRCS = Primal.cpp
LIBS = -L/opt/ibm/ILOG/CPLEX_Studio_Community129/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L/opt/ibm/ILOG/CPLEX_Studio_Community129/concert/lib/x86-64_linux/static_pic -lconcert -lm -lpthread -ldl
HDRS =  
EXEC = Primal

compila: Primal.cpp
#para compilar no lia no cplex versão 11
	${CC} ${CFLAGS1} ${INCLUDES} ${SRCS}
	${CC} ${CFLAGS2} ${INCLUDES} ${OBJS} -o ${EXEC} ${LIBS}

executa: compila
	./Primal ./instancias/dados
