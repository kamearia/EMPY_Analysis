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
