g++ -std=c++11 -o pathtracer.exe pathtracer.cpp -fopenmp && pathtracer.exe > output.ppm

g++ -o pathtracer_transparent pathtracer_transparent.cpp -std=c++11 -O3 && pathtracer_transparent > output.ppm

