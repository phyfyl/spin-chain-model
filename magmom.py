from matplotlib import pyplot
import numpy as np
import math
from random import *
import scipy.linalg as sla
import tinyarray
import time
from math import *
import cmath
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mpi4py import MPI

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

tau_x = tinyarray.array([[0, 1], [1, 0]])
tau_z = tinyarray.array([[1, 0], [0, -1]])
tau_y = tinyarray.array([[0, -1j], [1j, 0]])
tau_0 = tinyarray.array([[1, 0], [0, 1]])
tau_plus = tau_x + 1j*tau_y
tau_minus = tau_x -1j*tau_y
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
sigma_plus = (sigma_x + 1.0j*sigma_y)/2
sigma_minus = (sigma_x - 1.0j*sigma_y)/2
kron = np.kron
sin=np.sin
cos=np.cos
exp=np.exp
pi=np.pi
muB=5.78838*10**(-5)
nl=7
Bj=5.14
Bk=0.3*Bj
Jt=0.05
Kt=5
Jb=0.05
Kb=5
Bmin=0.5
Bmax=9
Nb=300
T_init = 10 # 初始温度
T_min = 1e-8   # 最小温度
alpha = 0.999   # 温度衰减因子
iterations = 10000  # 每个温度下的迭代次数

def Em(Bz,t):
    
    zeeman_u=0
    afm_u=0
    pma_u=0
    
    for i in np.arange(2,nl-3):
        zeeman_u=zeeman_u-Bz*cos(t[i])
        afm_u=afm_u+Bj/2*cos(t[i+1]-t[i])
        pma_u=pma_u-Bk/2*cos(t[i])**2

    zeeman_u=zeeman_u-Bz*cos(t[0])-Bz*cos(t[nl-1])
    afm_u=afm_u+Bj/2*(Jb*cos(t[1]-t[0])+Jt*cos(t[nl-1]-t[nl-2]))
    pma_u=pma_u-Bk/2*(Kb*cos(t[0])**2+Kt*cos(t[nl-1])**2)
            
    u=zeeman_u+afm_u+pma_u
    
    return u

# 随机生成初始解 [-pi, pi] 的范围
def random_solution():
    return np.random.uniform(-pi, pi, size=7)

# 生成邻域解，通过对当前解稍微扰动得到
def neighbor_solution(t):
    dtmax=0.05
    dtmin=0.005
    perturbation = np.random.uniform(-dtmax, dtmax, size=7)
    for i in range(len(perturbation)):
        dt = perturbation[i]
        if -dtmin < dt < 0:
            dt =-dtmin 
        elif 0 < dt < dtmin:
            dt = dtmin
        perturbation[i] = dt
    # perturbation = np.clip(perturbation, -dtmin, dtmin)
    new_solution = t + perturbation
    # 将新解的角度限制在 [-pi, pi] 的范围内
    for i in range(len(new_solution)):
        tn = new_solution[i]
        if tn < -pi:
            tn += 2*pi
        elif tn > pi:
            tn -= 2*pi
        new_solution[i] = tn
    return new_solution

# 退火算法
# 退火算法
def simulated_annealing(Bz):
    # 初始化当前解和能量
    current_solution = random_solution()
    current_energy = Em(Bz,current_solution)
    best_solution = current_solution
    best_energy = current_energy
    
    T = T_init  # 初始化温度
    
    while T > T_min:
        for _ in range(iterations):
            # 生成新解并计算其能量
            new_solution = neighbor_solution(current_solution)
            new_energy = Em(Bz,new_solution)
            
            # 判断是否接受新解
            if new_energy < current_energy:
                current_solution = new_solution
                current_energy = new_energy
            else:
                # 按一定概率接受更差的解，概率与温度 T 相关
                delta_energy = new_energy - current_energy
                acceptance_probability = np.exp(-delta_energy / T)
                if np.random.rand() < acceptance_probability:
                    current_solution = new_solution
                    current_energy = new_energy
            
            # 更新最佳解
            if current_energy < best_energy:
                best_solution = current_solution
                best_energy = current_energy
        
        # 降低温度
        T *= alpha
    print(best_solution*180/pi)
    print(best_energy)
    return best_solution*180/pi,best_energy

thetas=np.zeros([nl,Nb])
Es=np.zeros(Nb)
Bs=np.zeros(Nb)
for i in np.arange(rank,Nb,size):
    B=Bmin+i*(Bmax-Bmin)/(Nb-1)
    Bs[i]=B
    theta,E=simulated_annealing(B)
    thetas[:,i]=theta[:]
    Es[i]=E

thetas=comm.reduce(thetas,op=MPI.SUM,root=0)
Es=comm.reduce(Es,op=MPI.SUM,root=0)
Bs=comm.reduce(Bs,op=MPI.SUM,root=0)

if rank==0:
    #把thetas的结果储存成txt文件
    np.savetxt('thetas.txt',thetas)
    np.savetxt('Es.txt',Es)
    np.savetxt('Bs.txt',Bs)

