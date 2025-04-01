call my_env\Scripts\activate
git clone https://github.com/kamearia/EMPY_Solver.git ../GIT_EMPY_Solver
mkdir build
cmake -S ../GIT_EMPY_Solver -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON -Dpybind11_DIR="C:\EMSolution\EMSolpy5\EMPY_Analysis\my_env/Lib/site-packages/pybind11/share/cmake/pybind11"
cmake --build build --config release 
copy build\EMPY_Solver\Release\EMPY_Solver.cp310-win_amd64.pyd bin\Release\EMPY_Solver.pyd
rd /s /q build
rd /s /q ..\GIT_EMPY_Solver