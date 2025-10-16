set CMAKE_PREFIX_PATH=C:\EMSolution\EMSolpy5\EMPY_Analysis\em_env\Lib\site-packages
set pybind11_DIR=
mkdir build
# 2. CMakeの実行
# PYTHON_EXECUTABLEを仮想環境内のものに強制指定
cmake -S ../EMPY_Field/EMPY_Field -B build -G "Visual Studio 17 2022"

# cmake -S ../EMPY_Field/EMPY_Field -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON
cmake --build build --config release 
copy build\Release\EMPY_Field.pyd bin\Release\EMPY_Field.pyd
rd /s /q build

