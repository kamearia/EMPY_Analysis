# call my_env\Scripts\activate
mkdir build
# Conda環境のPythonを参照するため、PYTHON_EXECUTABLEなどを設定
cmake -S ../EMPY_Solver -B build \
      -G "Visual Studio 17 2022" \
      -DPYTHON_EXECUTABLE="C:\Users\kamea\anaconda3\envs\ngsolve_env_new\python.exe" \
      -DPYTHON_INCLUDE_DIR="C:\Users\kamea\anaconda3\envs\ngsolve_env_new\include" \
      -DNGS_LIB_DIR="C:\Users\kamea\anaconda3\envs\ngsolve_env_new\Library\lib"
cmake --build build --config release 
#  copy build\EMPY_Solver\Release\EMPY_Solver.cp310-win_amd64.pyd bin\Release\EMPY_Solver.pyd
copy build\EMPY_Solver\Release\EMPY_Solver.cp311-win_amd64.pyd bin\Release\EMPY_Solver.pyd