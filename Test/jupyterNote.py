import os
import sys

# ----------------------------------------------------------------------
# ！！重要！！ 特定した DLL のフォルダパスを正確に設定してください
# ----------------------------------------------------------------------
dll_paths_to_add = [
    # 以前見つかった libiomp5md.dll のパス
    r"C:\EMSolution\EMSolpy5\EMPY_Analysis\my_env\bin", 
    
    # ここに、Step 1 で特定した tbb.dll が入っているフォルダのパスを追加
    # 例: r"C:\Users\username\Anaconda3\Library\bin"
    # 例: r"C:\Program Files\Intel\MKL\...\bin"
    # 例: r"C:\Windows\System32" 
    # 例: r"C:\Users\...\Anaconda3\envs\base\Library\bin" # Anaconda環境の場合
    r"C:\Users\kamea\AppData\Roaming\Zoom\bin\OpenVINO_A",

    # 見つけたパスを追記してください 
    # r"TBB_DLL_PATH_HERE",
    # r"MKL_DLL_PATH_HERE",
    r"C:\EMSolution\EMSolpy5\EMPY_Analysis\my_env\Library\bin" 
]

# パスを結合
new_path_entries = os.pathsep.join(dll_paths_to_add)

# 環境変数 PATH の先頭に DLL 探索パスを強制的に追加（最も重要）
# これで Jupyter カーネルがシステム DLL より先に、これらのパスを参照します
os.environ['PATH'] = new_path_entries + os.pathsep + os.environ.get('PATH', '')
print(f"PATHにDLL探索パスを追加: {new_path_entries}")

# sys.path にも追加
for p in dll_paths_to_add:
    if p not in sys.path:
        sys.path.append(p)

print("環境変数を設定しました。ngsolveをインポートします。")
# ngsolve のインポートを試みる
from ngsolve import *
print("NGSolveのインポートに成功しました！")

# EMPY_Fieldのインポートに進む
from EMPY_Field import *
print("EMPY_Fieldのインポートも成功しました！")