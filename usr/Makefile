CC=clang++

all: bin/encode bin/usrcat

bin/encode: src/encode.cpp
	${CC} -o $@ $^ -std=c++11 -O2 -I${OPENBABEL_ROOT}/include -L${OPENBABEL_ROOT}/lib -lopenbabel

bin/usrcat: src/usrcat.cpp
	${CC} -o $@ $^ -std=c++11 -O2 -pthread -Wno-deprecated-declarations -Wno-deprecated-register -I${OPENBABEL_ROOT}/include -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/Foundation/include -I${POCO_ROOT}/Net/include -L${OPENBABEL_ROOT}/lib -lopenbabel -L${BOOST_ROOT}/lib -lboost_system -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_date_time -L${POCO_ROOT}/lib/Linux/x86_64 -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT}/build/linux2/64/release/use-system-all -lmongoclient

clean:
	rm -f bin/*
