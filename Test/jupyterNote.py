import os
import sys

# ----------------------------------------------------------------------
# �I�I�d�v�I�I ���肵�� DLL �̃t�H���_�p�X�𐳊m�ɐݒ肵�Ă�������
# ----------------------------------------------------------------------
dll_paths_to_add = [
    # �ȑO�������� libiomp5md.dll �̃p�X
    r"C:\EMSolution\EMSolpy5\EMPY_Analysis\my_env\bin", 
    
    # �����ɁAStep 1 �œ��肵�� tbb.dll �������Ă���t�H���_�̃p�X��ǉ�
    # ��: r"C:\Users\username\Anaconda3\Library\bin"
    # ��: r"C:\Program Files\Intel\MKL\...\bin"
    # ��: r"C:\Windows\System32" 
    # ��: r"C:\Users\...\Anaconda3\envs\base\Library\bin" # Anaconda���̏ꍇ
    r"C:\Users\kamea\AppData\Roaming\Zoom\bin\OpenVINO_A",

    # �������p�X��ǋL���Ă������� 
    # r"TBB_DLL_PATH_HERE",
    # r"MKL_DLL_PATH_HERE",
    r"C:\EMSolution\EMSolpy5\EMPY_Analysis\my_env\Library\bin" 
]

# �p�X������
new_path_entries = os.pathsep.join(dll_paths_to_add)

# ���ϐ� PATH �̐擪�� DLL �T���p�X�������I�ɒǉ��i�ł��d�v�j
# ����� Jupyter �J�[�l�����V�X�e�� DLL ����ɁA�����̃p�X���Q�Ƃ��܂�
os.environ['PATH'] = new_path_entries + os.pathsep + os.environ.get('PATH', '')
print(f"PATH��DLL�T���p�X��ǉ�: {new_path_entries}")

# sys.path �ɂ��ǉ�
for p in dll_paths_to_add:
    if p not in sys.path:
        sys.path.append(p)

print("���ϐ���ݒ肵�܂����Bngsolve���C���|�[�g���܂��B")
# ngsolve �̃C���|�[�g�����݂�
from ngsolve import *
print("NGSolve�̃C���|�[�g�ɐ������܂����I")

# EMPY_Field�̃C���|�[�g�ɐi��
from EMPY_Field import *
print("EMPY_Field�̃C���|�[�g���������܂����I")