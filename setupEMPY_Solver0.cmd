call my_env\Scripts\activate
mkdir build
cmake -S ../EMPY_Solver -B build -G "Visual Studio 17 2022"
cmake --build build --config release 
copy build\EMPY_Solver\Release\EMPY_Solver.cp310-win_amd64.pyd bin\Release\EMPY_Solver.pyd