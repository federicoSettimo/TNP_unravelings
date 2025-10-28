Heisenberg: Heisenberg.cpp
	g++ Heisenberg.cpp -o Heisenberg.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Heisenberg.x
	python3 Heisenberg.py

photon_counting_MCWF: photon_counting_MCWF.cpp
	g++ photon_counting_MCWF.cpp -o photon_counting_MCWF.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./photon_counting_MCWF.x
	python3 photon_counting_MCWF.py

check_CP_div_Heisenberg: check_CP_div_Heisenberg.cpp
	g++ check_CP_div_Heisenberg.cpp -o check_CP_div_Heisenberg.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./check_CP_div_Heisenberg.x
	python3 rates.py