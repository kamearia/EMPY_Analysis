@echo off
setlocal enabledelayedexpansion

:: ----------------------------------------------------
:: 0. �����Z�b�g�A�b�v�i���z���̍쐬�ƈˑ��p�b�P�[�W�̃C���X�g�[���j
:: ----------------------------------------------------
:: 0.1 ���z�� em_env �����݂��Ȃ��ꍇ�̂ݍ쐬
if not exist em_env (
    echo.
    echo === 0.1 ���z�� em_env ���쐬 ===
    python -m venv em_env
    if errorlevel 1 goto :cleanup_and_exit
)

:: 0.2 ���z���̃A�N�e�B�x�[�g
call em_env\Scripts\activate
if errorlevel 1 goto :cleanup_and_exit

:: 0.3 �ˑ��p�b�P�[�W�̃C���X�g�[��
if not exist requirements.txt (
    echo [ERROR] requirements.txt �����݂��܂���B�ˑ��p�b�P�[�W�̃C���X�g�[�����X�L�b�v���܂��B
    echo C++���W���[���̃r���h�ɂ�ngsolve�Ȃǂ̈ˑ��p�b�P�[�W���K�v�ł��B
) else (
    echo.
    echo === 0.3 �ˑ��p�b�P�[�W���C���X�g�[�� ===
    pip install -r requirements.txt
    if errorlevel 1 goto :cleanup_and_exit
)

:: ----------------------------------------------------
:: 1. CMAKE���ϐ��̐ݒ� (�|�[�^�u����)
:: ----------------------------------------------------
set CMAKE_PREFIX_PATH=%CD%\em_env\Lib\site-packages
set pybind11_DIR=

:: ----------------------------------------------------
:: 2. bin\Release �̃o�b�N�A�b�v�����i����Łj
:: ----------------------------------------------------
echo.
echo === 2. bin\Release �̃o�b�N�A�b�v�m�F�ƍ쐬 ===
set "TARGET_DIR=bin\Release"
set "TEMP_TS_FILE=.\temp_timestamp.txt"

if exist %TARGET_DIR% (
    :: PowerShell�Ń^�C���X�^���v���ꎞ�t�@�C���ɏ����o���i�ł��N���[���ȕ��@�j
    powershell -command "(Get-Date -format yyyyMMdd_HHmmss) | Out-File -FilePath %TEMP_TS_FILE% -Encoding ASCII"
    if errorlevel 1 goto :cleanup_and_exit

    :: �t�@�C������^�C���X�^���v��ǂݍ��݁i�]���ȋ󔒁E���s�����S�ɔr���j
    for /f %%i in (%TEMP_TS_FILE%) do set "TIMESTAMP_STR=%%i"
    del %TEMP_TS_FILE%
    
    set "BACKUP_DIR=!TARGET_DIR!_!TIMESTAMP_STR!"

    echo �Â� %TARGET_DIR% �� !BACKUP_DIR! �Ƀo�b�N�A�b�v���܂��B
    :: move �R�}���h���g�p (���萫������)
    move %TARGET_DIR% !BACKUP_DIR!
    if errorlevel 1 goto :cleanup_and_exit
)

:: bin �t�H���_�����݂��Ȃ��ꍇ�����邽�߁Abin��Release�𗼕��쐬
if not exist bin mkdir bin
mkdir %TARGET_DIR%
if errorlevel 1 goto :cleanup_and_exit


:: ----------------------------------------------------
:: 3. EMPY_Field �̃r���h
:: ----------------------------------------------------
echo.
echo === 3. EMPY_Field �̃r���h�J�n ===
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
:: 4. EMPY_Solver �̃r���h
:: ----------------------------------------------------
echo.
echo === 4. EMPY_Solver �̃r���h�J�n ===
mkdir build
if errorlevel 1 goto :cleanup_and_exit
cmake -S ../EMPY_Solver -B build -G "Visual Studio 17 2022"
if errorlevel 1 goto :cleanup_and_exit
cmake --build build --config Release
if errorlevel 1 goto :cleanup_and_exit
:: �m�肵���p�X�ŃR�s�[
copy build\EMPY_Solver\Release\EMPY_Solver.cp312-win_amd64.pyd %TARGET_DIR%\EMPY_Solver.pyd
if errorlevel 1 goto :cleanup_and_exit
rd /s /q build

:: ----------------------------------------------------
:: 5. SparseSolvPy �̃r���h
:: ----------------------------------------------------
echo.
echo === 5. SparseSolvPy �̃r���h�J�n ===
git clone https://github.com/JP-MARs/SparseSolv.git
if errorlevel 1 goto :cleanup_and_exit
copy SparseSolveCMakeLists.txt SparseSolv\CMakeLists.txt
if errorlevel 1 goto :cleanup_and_exit

:: SparseSolv���r���h
mkdir build
if errorlevel 1 goto :cleanup_and_exit
cmake -S SparseSolv -B build -G "Visual Studio 17 2022" -DPYBIND_EXPORT=ON
if errorlevel 1 goto :cleanup_and_exit
cmake --build build --config Release
if errorlevel 1 goto :cleanup_and_exit

:: �m�肵���R�s�[���p�X�ŃR�s�[
copy SparseSolv\Release\SparseSolvPy.cp312-win_amd64.pyd %TARGET_DIR%\SparseSolvPy.pyd
if errorlevel 1 goto :cleanup_and_exit

rd /s /q build
rd /s /q SparseSolv
if errorlevel 1 goto :cleanup_and_exit

:: ----------------------------------------------------
:: 6. ����I��
:: ----------------------------------------------------
echo.
echo === �S�ẴZ�b�g�A�b�v�ƃr���h������Ɋ������܂��� ===
goto :EOF


:: ----------------------------------------------------
:: �G���[�������[�`��
:: ----------------------------------------------------
:cleanup_and_exit
echo.
echo [FATAL ERROR] �R�}���h�̎��s�Ɏ��s���܂����B�ꎞ�t�@�C�����N���[���A�b�v���܂��B
:: ���݃`�F�b�N�����āA�G���[��}�~���Ȃ��狭���폜
if exist build (rd /s /q build)
if exist SparseSolv (rd /s /q SparseSolv)
if exist %TEMP_TS_FILE% (del %TEMP_TS_FILE%)
endlocal
exit /b 1