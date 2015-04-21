clear
astyle "*.cpp" --indent=force-tab --style=java / -A2 --recursive
astyle "*.c" --indent=force-tab --style=java / -A2 --recursive
astyle "*.h" --indent=force-tab --style=java / -A2 --recursive
find . -name "*.orig" -type f -delete
rm RubiksCubeGA.exe
mpicxx RubiksCubeGA.cpp -o RubiksCubeGA.exe
mpirun -np 16 ./RubiksCubeGA.exe

