gaussian:
	g++ -std=c++17 -Werror -pedantic main.cpp -o main

.PHONY: clean

clean:
	rm main
