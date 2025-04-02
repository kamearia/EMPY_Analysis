インストール

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
C:\...\EMPY\EMPY_Analysis> pip install -r requirements.txt  

4. モジュールインストール  
C:\...\EMPY\EMPY_Analysis> setup
以下のファイルが作成されたことを確認
bin\Relsese\SparseSolvPy.pyd  
bin\Relsese\EMPY_Solver.pyd
bin\Relsese\EMPY_Field.pyd

5. jupyter実行 (EMPY_Solver使用）  
C:\...\EMPY\EMPY_Analysis> jupyter notebook "A-2_Phi_Potential_BathPlate_with_Reg.ipynb"

Run->Run All Cells  

6. jupyter実行 (JP_MARｓ/SparseSolve使用）  
A-2 Phi_Potential_BathPlate_with_Reg.ipynにおいて、cpp_solver="EMPY"を、cpp_solver="JP_MARs"に変更 
cpp_solver="EMPY"
#cpp_solver="JP_MARs"
-->
cpp_solver="EMPY"  
#cpp_solver="JP_MARs"
Kernel->Restart Kernel and Run All Cells... 
JP_MARｓ/SparseSolvでは収束しないことが確認される。  
