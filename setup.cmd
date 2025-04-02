rem call setupEnv
mkdir bin
cd bin
mkdir Release
cd ..
call setupSparseSolv
call setupEMPY_Solver
call setupEMPY_Field
jupyter notebook "EddyCurrent/A-2 Phi_Potential_BathPlate_with_Reg.ipynb"