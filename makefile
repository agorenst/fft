driver: driver.cpp fft.o
	/s/gcc-4.8.0/bin/g++ -Wall -std=c++0x driver.cpp fft.o -o driver

fft.o: fft.h fft.cpp
	/s/gcc-4.8.0/bin/g++ -Wall -std=c++0x fft.cpp -o fft.o -c

clean:
	rm -f driver *~ *.o
