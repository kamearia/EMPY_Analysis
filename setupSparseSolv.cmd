git clone https://github.com/JP-MARs/SparseSolv.git
mkdir build
cmake -S SparseSolv -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON -Dpybind11_DIR="C:\EMSolution\EMSolpy5\EMPY_Analysis\my_env/Lib/site-packages/pybind11/share/cmake/pybind11"
cmake --build build --config release 
copy SparseSolv\Release\SparseSolvPy.cp310-win_amd64.pyd bin\Release\SparseSolvPy.pyd
rd /s /q build
rd /s /q SparseSolv