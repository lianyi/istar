CC=clang++

all: bin/encode bin/usr

bin/encode: src/encode.cpp
	${CC} -o $@ $^ -std=c++11 -O2

bin/usr: src/usr.cpp
	${CC} -o $@ $^ -std=c++11 -O2 -mavx

clean:
	rm -f bin/*
