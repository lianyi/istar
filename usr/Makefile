CC=clang++ -O2

all: bin/score bin/encode bin/usr

bin/score: obj/score.o
	${CC} -o $@ $^

bin/encode: obj/encode.o
	${CC} -o $@ $^ -L${OPENBABEL_ROOT}/lib -lopenbabel

bin/usr: obj/main.o
	${CC} -o $@ $^ -pthread -L${OPENBABEL_ROOT}/lib -lopenbabel -L${BOOST_ROOT}/lib -lboost_system -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_date_time -L${POCO_ROOT}/lib -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT} -lmongoclient

obj/score.o: src/score.cpp
	${CC} -o $@ $< -c -std=c++11 -DNDEBUG -Wall

obj/encode.o: src/encode.cpp
	${CC} -o $@ $< -c -std=c++11 -DNDEBUG -Wall -I${OPENBABEL_ROOT}/include

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++11 -DNDEBUG -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${OPENBABEL_ROOT}/include -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/include

clean:
	rm -f bin/* obj/*
