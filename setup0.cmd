@echo off
setlocal enabledelayedexpansion

:: ----------------------------------------------------
:: 0. 初期セットアップ（仮想環境の作成と依存パッケージのインストール）
:: ----------------------------------------------------
:: 0.1 仮想環境 em_env が存在しない場合のみ作成
if not exist em_env (
    echo.
    echo === 0.1 仮想環境 em_env を作成 ===
    python -m venv em_env
    if errorlevel 1 goto :cleanup_and_exit
)

:: 0.2 仮想環境のアクティベート
call em_env\Scripts\activate
if errorlevel 1 goto :cleanup_and_exit

:: 0.3 依存パッケージのインストール
if not exist requirements.txt (
    echo [ERROR] requirements.txt が存在しません。依存パッケージのインストールをスキップします。
    echo C++モジュールのビルドにはngsolveなどの依存パッケージが必要です。
) else (
    echo.
    echo === 0.3 依存パッケージをインストール ===
    pip install -r requirements.txt
    if errorlevel 1 goto :cleanup_and_exit
)

:: ----------------------------------------------------
:: 1. CMAKE環境変数の設定 (ポータブル版)
:: ----------------------------------------------------
set CMAKE_PREFIX_PATH=%CD%\em_env\Lib\site-packages
set pybind11_DIR=

:: ----------------------------------------------------
:: 2. bin\Release のバックアップ処理（安定版）
:: ----------------------------------------------------
echo.
echo === 2. bin\Release のバックアップ確認と作成 ===
set "TARGET_DIR=bin\Release"
set "TEMP_TS_FILE=.\temp_timestamp.txt"

if exist %TARGET_DIR% (
    :: PowerShellでタイムスタンプを一時ファイルに書き出し（最もクリーンな方法）
    powershell -command "(Get-Date -format yyyyMMdd_HHmmss) | Out-File -FilePath %TEMP_TS_FILE% -Encoding ASCII"
    if errorlevel 1 goto :cleanup_and_exit

    :: ファイルからタイムスタンプを読み込み（余分な空白・改行を完全に排除）
    for /f %%i in (%TEMP_TS_FILE%) do set "TIMESTAMP_STR=%%i"
    del %TEMP_TS_FILE%
    
    set "BACKUP_DIR=!TARGET_DIR!_!TIMESTAMP_STR!"

    echo 古い %TARGET_DIR% を !BACKUP_DIR! にバックアップします。
    :: move コマンドを使用 (安定性が高い)
    move %TARGET_DIR% !BACKUP_DIR!
    if errorlevel 1 goto :cleanup_and_exit
)

:: bin フォルダが存在しない場合もあるため、binとReleaseを両方作成
if not exist bin mkdir bin
mkdir %TARGET_DIR%
if errorlevel 1 goto :cleanup_and_exit


:: ----------------------------------------------------
:: 3. EMPY_Field のビルド
:: ----------------------------------------------------
echo.
echo === 3. EMPY_Field のビルド開始 ===
mkdir build
if errorlevel 1 goto :cleanup_and_exit
cmake -S ../EMPY_Field/EMPY_Field -B build -G "Visual Studio 17 2022"
if errorlevel 1 goto :cleanup_and_exit
cmake --build build --config Release
if errorlevel 1 goto :cleanup_and_exit
copy build\Release\EMPY_Field.pyd %TARGET_DIR%\EMPY_Field.pyd
if errorlevel 1 goto :cleanup_and_exit
rd /s /q build

:: ----------------------------------------------------
:: 4. EMPY_Solver のビルド
:: ----------------------------------------------------
echo.
echo === 4. EMPY_Solver のビルド開始 ===
mkdir build
if errorlevel 1 goto :cleanup_and_exit
cmake -S ../EMPY_Solver -B build -G "Visual Studio 17 2022"
if errorlevel 1 goto :cleanup_and_exit
cmake --build build --config Release
if errorlevel 1 goto :cleanup_and_exit
:: 確定したパスでコピー
copy build\EMPY_Solver\Release\EMPY_Solver.cp312-win_amd64.pyd %TARGET_DIR%\EMPY_Solver.pyd
if errorlevel 1 goto :cleanup_and_exit
rd /s /q build

:: ----------------------------------------------------
:: 5. SparseSolvPy のビルド
:: ----------------------------------------------------
echo.
echo === 5. SparseSolvPy のビルド開始 ===
git clone https://github.com/JP-MARs/SparseSolv.git
if errorlevel 1 goto :cleanup_and_exit
copy SparseSolveCMakeLists.txt SparseSolv\CMakeLists.txt
if errorlevel 1 goto :cleanup_and_exit

:: SparseSolvをビルド
mkdir build
if errorlevel 1 goto :cleanup_and_exit
cmake -S SparseSolv -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON
if errorlevel 1 goto :cleanup_and_exit
cmake --build build --config Release
if errorlevel 1 goto :cleanup_and_exit

:: 確定したコピー元パスでコピー
copy SparseSolv\Release\SparseSolvPy.cp312-win_amd64.pyd %TARGET_DIR%\SparseSolvPy.pyd
if errorlevel 1 goto :cleanup_and_exit

rd /s /q build
rd /s /q SparseSolv
if errorlevel 1 goto :cleanup_and_exit

:: ----------------------------------------------------
:: 6. 正常終了
:: ----------------------------------------------------
echo.
echo === 全てのセットアップとビルドが正常に完了しました ===
goto :EOF


:: ----------------------------------------------------
:: エラー処理ルーチン
:: ----------------------------------------------------
:cleanup_and_exit
echo.
echo [FATAL ERROR] コマンドの実行に失敗しました。一時ファイルをクリーンアップします。
:: 存在チェックを入れて、エラーを抑止しながら強制削除
if exist build (rd /s /q build)
if exist SparseSolv (rd /s /q SparseSolv)
if exist %TEMP_TS_FILE% (del %TEMP_TS_FILE%)
endlocal
exit /b 1