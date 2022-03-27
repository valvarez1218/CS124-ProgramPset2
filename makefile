all: strassen.cpp
	g++ -std=c++17 -O2 -Wall -Wextra strassen.cpp -o strassen -lm -lpthread

clean:
	rm strassen