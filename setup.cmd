call setupEnv
mkdir bin
cd bin
mkdir Release
cd ..
call setupSparseSolv
call setupEMPY_Solver
call setupEMPY_Field
jupyter notebook "EddyCurrent/A-2 Phi Potential BathPlate with Reg.ipynb"