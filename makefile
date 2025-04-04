Heisenberg: Heisenberg.cpp
	g++ Heisenberg.cpp -o Heisenberg.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Heisenberg.x
	python3 Heisenberg.py

photon_counting: photon_counting.cpp
	g++ photon_counting.cpp -o photon_counting.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./photon_counting.x
	#python3 photon_counting.py