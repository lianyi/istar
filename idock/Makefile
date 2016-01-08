CC=g++ -O2 -flto

bin/idock: obj/scoring_function.o obj/box.o obj/quaternion.o obj/io_service_pool.o obj/safe_counter.o obj/receptor.o obj/ligand.o obj/grid_map_task.o obj/monte_carlo_task.o obj/random_forest_test.o obj/main.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_thread -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -L${MONGODBCXXDRIVER_ROOT}/sharedclient -lmongoclient -L${CURL_ROOT}/lib -lcurl

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++14 -DNDEBUG -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/include -I${CURL_ROOT}/include

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++14 -DNDEBUG -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT}

clean:
	rm -f bin/idock obj/*.o
