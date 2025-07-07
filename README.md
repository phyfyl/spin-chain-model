# Spin Chain Model
## Description
### (1) This is a code obtaining the gound state of a spin-chain model under a series of external magnetic fields.
### (2) The gound state is obtained by simulated annealing algorithm.
### (3) There are two versions: python and fortran. The python version can be directly run without compile and the parameters can be changed within the source code. The fortran version is need to compile before using, the compile command is:
mpiifort magmom.f90 -o magmom
### if you use the mpiifort compiler. After compiling, you need a input file named "input.txt" (parameters) to use the program "magmom", a sample input file has been put in the folder.
### (5) The computing speed of python version is lower than fortran version, and the fortran version is more convenient to use beacause of the existence of the "input.txt" file. I recommand the fortran version if you use.

## The description of the parameters of "input.txt":
### Nls: the number of spins of the spin chain.
### Bj: the strength of exchange coupling between adjacent spins: >0: AFM, <0: FM
### Bk: the uniaxial magnetic anisoropic strength of the spin chain.
### Jt(Jb): the relative strength of the exchange coupling between the Nls(2) th and Nls-1(1) th spins considering the possible boundary effects.
### Kt(Kb): the relative strength of the uniaxial magnetic anisoropic strength of the Nls(1) th spins considering the possible boundary effects.
### Bmin(Bmax): the minimum(maximum) value of the magnetic field.
### Nbs: the number of magnetic fields betwenn Bmin and Bmax.
### T_init(T_min): the initial(final) temperature of the simulated annealing algorithm.
### alpha: the decreasing rate of the temperature.
### iterations: the number of ietrations for every value of temperature.
