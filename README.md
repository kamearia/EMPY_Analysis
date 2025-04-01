インストール

EMPY_Analysisのインストール

1．install directoryを作成する
>mkdir EMPY
>cd EMPY

2.EMPY_Analysis directoryを読み込む
>git clone https://github.com/kamearia/EMPY_Analysis.git
>cd EMPY_Analysis

3. python環境作成
>python -m venv my_env
>my_env\Scripts\activate
>pip install -r requirements.txt

４．JP_MARｓ/SpaerseSolve
>call my_env\Scripts\activate
>git clone https://github.com/JP-MARs/SparseSolv.git
>mkdir build
>cmake -S SparseSolv -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON -Dpybind11_DIR=" Install directory\EMPY_Analysis\my_env/Lib/site-packages/pybind11/share/cmake/pybind11"
>cmake --build build --config release 
>copy SparseSolv\Release\SparseSolvPy.cp310-win_amd64.pyd bin\Release\SparseSolvPy.pyd
>rd /s /q build
>rd /s /q SparseSolv

5. 
＞setupEMPY_Field

6. 
＞setupEMPY_Solver

７．実行test
>jupyter notebook EddyCurrent/A-2_Phi_Potential_BathPlate_with_Reg.ipynb

