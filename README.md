EMPY_Analysisのインストール  

1 install directoryを作成する  
C:\...> mkdir EMPY  
C:\...> cd EMPY  

2. EMPY_Analysis directoryを読み込む  
C:\...\EMPY> git clone https://github.com/kamearia/EMPY_Analysis.git  
C:\...\EMPY> cd EMPY_Analysis  

3. python環境作成  
C:\...\EMPY\EMPY_Analysis> python -m venv my_env  
C:\...\EMPY\EMPY_Analysis> my_env\Scripts\activate  
(my_env)C:\...\EMPY\EMPY_Analysis> pip install -r requirements.txt  

4. モジュールインストール  
(my_env)C:\...\EMPY\EMPY_Analysis> setup  
以下のファイルが作成されたことを確認  
bin\Relsese\SparseSolvPy.pyd    
bin\Relsese\EMPY_Solver.pyd  
bin\Relsese\EMPY_Field.pyd  

5. jupyter実行 (EMPY_Solver使用）  
(my_env)C:\...\EMPY\EMPY_Analysis> jupyter notebook "EddyCurrent\A-2_Phi_Potential_BathPlate_with_Reg.ipynb"  

Run->Run All Cells  
保存された図と一致する結果を確認される。

6. Jupyter Labの実行
(my_env)C:\...\EMPY\EMPY_Analysis>jupyter lab

Linuxでのインストール

1. install directoryを作成し、EMPY_Analysisを読み込む  
$ mkdir EMPY  
$ cd EMPY  
$ git clone https://github.com/kamearia/EMPY_Analysis.git  
$ cd EMPY_Analysis  

2. python環境作成（Python 3.12）  
$ python3.12 -m venv my_env  
$ source my_env/bin/activate  
(my_env)$ pip install ngsolve numpy scipy matplotlib sympy mpmath jupyter  
requirements.txt には Windows 専用パッケージ（pywin32、pywinpty）が含まれるため、Linux ではそのまま使えない。

3. モジュールについて  
bin/Release/EMPY_Solver.cpython-312-x86_64-linux-gnu.so が含まれており、Python 3.12（x86_64）でそのまま使える。  
他の Python バージョンを使う場合は、EMPY_Solver リポジトリ（https://github.com/kamearia/EMPY_Solver）からビルドする。手順は同リポジトリの EMPY_Solver/LINUX.md を参照。  
SparseSolvPy と EMPY_Field は Linux 版がまだない。include/MatrixSolver.py は SparseSolvPy を import するため、Linux ではその行をコメントアウトするか try/except で囲むこと（cpp_solver="EMPY" の利用では SparseSolvPy は使われない）。

4. jupyter実行  
Windows の場合と同じ。
