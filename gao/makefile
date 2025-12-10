CC = g++

gaomain: iofile.o gao.o gaiodata.o BasicES.o
	$(CC) main.cpp -Wall -Wextra -I ../ga/ -I ../ -L ../ga/ -l ga  -o main.exe

iofile.o:
	$(CC) -c -Wall -Wextra ../iofile/iofile.cpp

gao.o: gaiodata.o
	$(CC) -c -Wall -Wextra -I ../ga/ -I ../ gao.cpp

gaiodata.o: BasicES.o
	$(CC) -Wall -Wextra -I ../ga/ -I ../ gaiodata.cpp ../ga/GAAllele.o

BasicES.o:
	$(CC) -Wall -Wextra BasicES.cpp

#g++ main.cpp iofile.o gao.o gaiodata.o BasicES.o -Wall -Wextra -I ../ga/ -I ../ -L ../ga/ -l ga -o main.exe