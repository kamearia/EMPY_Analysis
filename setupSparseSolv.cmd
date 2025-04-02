git clone https://github.com/JP-MARs/SparseSolv.git
copy SparseSolveCMakeLists.txt SparseSolv\CMakeLists.txt
mkdir build
cmake -S SparseSolv -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON
cmake --build build --config release 
copy SparseSolv\Release\SparseSolvPy.cp310-win_amd64.pyd bin\Release\SparseSolvPy.pyd
rd /s /q build
rd /s /q SparseSolv