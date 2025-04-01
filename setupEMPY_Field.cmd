call my_env\Scripts\activate
git clone https://github.com/kamearia/EMPY_Field.git  ../GIT_EMPY_Field
mkdir build
cmake -S ../GIT_EMPY_Field/EMPY_Field -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON
cmake --build build --config release 
copy build\Release\EMPY_Field.pyd bin\Release\EMPY_Field.pyd
rd /s /q build
rd /s /q ..\GIT_EMPY_Field

