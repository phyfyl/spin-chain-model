MODULE VariablesModule
    IMPLICIT NONE
    INTEGER :: Nls, NBs,iterations
    REAL :: Bj, Bk, Jt,Kt,Bmin, Bmax,T_init, T_min, alpha
END MODULE VariablesModule

SUBROUTINE ReadVariables
    USE VariablesModule
    IMPLICIT NONE
    INTEGER :: unit_number, iunit,pos
    CHARACTER(LEN=100) :: line
 !   CHARACTER(LEN=20) :: var_name
 !   REAL :: var_value

    ! 定义文件的逻辑单元号
    unit_number = 10

    ! 打开文件 input.txt
    open(unit=unit_number, file='./input.txt', status='old', action='read')

    ! 逐行读取文件
    do while (.true.)
        read(unit_number, '(A)', iostat=iunit) line
        if (iunit /= 0) exit  ! 读取结束或出错

        ! 查找并解析 Nls
        pos = index(line, 'Nls=')
        if (pos > 0) then
            read(line(pos+4:), *, iostat=iunit) Nls
        end if

        pos = index(line, 'Bj=')
        if (pos  > 0) then
            read(line(pos+3:), *, iostat=iunit) Bj
        end if

        pos = index(line, 'Bk=')
        if (pos  > 0) then
            read(line(pos+3:), *, iostat=iunit) Bk
        end if

        pos = index(line, 'Jt=')
        if (pos  > 0) then
            read(line(pos+3:), *, iostat=iunit) Jt
        end if

        pos = index(line, 'Kt=')
        if (pos  > 0) then
            read(line(pos+3:), *, iostat=iunit) Kt
        end if

        pos = index(line, 'Bmin=')
        if (pos  > 0) then
            read(line(pos+5:), *, iostat=iunit) Bmin
        end if

        pos = index(line, 'Bmax=')
        if (pos  > 0) then
            read(line(pos+5:), *, iostat=iunit) Bmax
        end if

        pos = index(line, 'NBs=')
        if (pos  > 0) then
            read(line(pos+4:), *, iostat=iunit) NBs
        end if

        pos = index(line, 'T_init=')
        if (pos  > 0) then
            read(line(pos+7:), *, iostat=iunit) T_init
        end if

        pos = index(line, 'T_min=')
        if (pos  > 0) then
            read(line(pos+6:), *, iostat=iunit) T_min
        end if

        pos = index(line, 'alpha=')
        if (pos  > 0) then
            read(line(pos+6:), *, iostat=iunit) alpha
        end if

        pos = index(line, 'iterations=')
        if (pos  > 0) then
            read(line(pos+11:), *, iostat=iunit) iterations
        end if

    end do

    ! 关闭文件
    close(unit_number)

    ! 从文件中逐行读取内容
 !   DO
  !      READ(unit_number, '(A)', IOSTAT=ios) line
   !     IF (ios /= 0) EXIT  ! 如果读取结束或出现错误，退出循环

        ! 解析每一行的变量名和值
  !      IF (INDEX(line, '=') /= 0) THEN
   !         READ(line, '(A20, "=", F10.5)', IOSTAT=ios) var_name, var_value
!
            ! 根据变量名存储对应的值
 !           SELECT CASE (TRIM(ADJUSTL(var_name)))
  !              CASE ('Nls')
  !                  Nls = INT(var_value)
   !             CASE ('Bj')
   !                 Bj = var_value
   !             CASE ('Bk')
   !                 Bk = var_value
   !             CASE ('Jt')
   !                 Jt = var_value
   !             CASE ('Kt')
   !                 Kt = var_value
   !             CASE ('Bmin')
   !                 Bmin = var_value
   !             CASE ('Bmax')
   !                 Bmax = var_value
   !             CASE ('NBs')
   !                 NBs = INT(var_value)
   !             CASE ('T_init')
  !                  T_init = var_value
  !              CASE ('T_min')
  !                  T_min = var_value
  !              CASE ('alpha')
  !                  alpha = var_value
  !              CASE ('iterations')
  !                  iterations = INT(var_value)
  !          END SELECT
  !      END IF
  !  END DO

    ! 关闭文件
 !   CLOSE(unit_number)
END SUBROUTINE ReadVariables

module constants
  implicit none
  real(8), parameter :: pi = 3.141592653589793d0
end module constants

function Em(Bz, t) result(u)
  use constants
  use VariablesModule
  implicit none
  integer :: i
  real(8) :: Bz, u
  real(8) :: t(Nls),afm_u,zeeman_u,pma_u

  afm_u=0d0
  zeeman_u=0d0
  pma_u=0d0
  do i=1,Nls-2
     afm_u=afm_u+Bj/2*cos(t(i)-t(i+1)) 
  end do

  afm_u=afm_u+Bj/2*Jt*cos(t(Nls-1)-t(Nls)) 

  do i=1,Nls-1
     zeeman_u=zeeman_u-Bz*cos(t(i))
     pma_u=pma_u-Bk/2*cos(t(i))**2
  end do

  zeeman_u=zeeman_u-Bz*cos(t(Nls))
  pma_u=pma_u-Bk/2*Kt*cos(t(Nls))**2

  u=afm_u+pma_u+zeeman_u

end function Em

subroutine simulated_annealing(Bz, best_solution, best_energy)
  use constants
  use VariablesModule
  implicit none
  real(8) :: Bz, current_energy, new_energy, best_energy, T,delta_energy
  real(8) :: current_solution(Nls),best_solution(Nls), new_solution(Nls), perturbation(Nls),acceptance_probability
  integer :: i
  real(8),external :: Em
  real(8) :: rand_val

  call random_solution(current_solution)
  current_energy = Em(Bz, current_solution)
  best_energy = current_energy
  best_solution = current_solution

  T = T_init
  do while (T > T_min)
    do i = 1, iterations
      call neighbor_solution(current_solution, new_solution)
      new_energy = Em(Bz, new_solution)
      
      if (new_energy < current_energy) then
        current_solution = new_solution
        current_energy = new_energy
      else
        delta_energy = new_energy - current_energy
        acceptance_probability = exp(-delta_energy / T)
        call random_number(rand_val)
        if (rand_val < acceptance_probability) then
          current_solution = new_solution
          current_energy = new_energy
        endif
      endif
      
      if (current_energy < best_energy) then
        best_solution = current_solution
        best_energy = current_energy
      endif
    end do
    T = T * alpha
  end do
end subroutine simulated_annealing

subroutine random_solution(solution)
  use constants
  use VariablesModule
  implicit none
  real(8) :: solution(Nls)
  real(8) :: random_vals(Nls)
  integer :: i

  call random_number(random_vals)
  do i = 1, Nls
    solution(i) = -pi + 2.0d0 * pi * random_vals(i)
  end do
end subroutine random_solution

subroutine neighbor_solution(current_solution, new_solution)
  use constants
  use VariablesModule
  implicit none
  real(8) :: current_solution(Nls), new_solution(Nls)
  real(8) :: perturbation(Nls), dt
  integer :: i
  real(8), parameter :: dtmax = 0.05d0, dtmin = 0.005d0

  ! 生成扰动
  call random_number(perturbation)
  perturbation = -dtmax + (2.0d0 * dtmax) * perturbation

  do i = 1, Nls
    dt = perturbation(i)
    if (abs(dt) < dtmin) then
      if (dt < 0) then
        dt = -dtmin
      else
        dt = dtmin
      endif
    endif
    new_solution(i) = current_solution(i) + dt

    ! 确保新解在 [-pi, pi] 之间
    if (new_solution(i) < -pi) then
      new_solution(i) = new_solution(i) + 2.0d0 * pi
    elseif (new_solution(i) > pi) then
      new_solution(i) = new_solution(i) - 2.0d0 * pi
    endif
  end do
end subroutine neighbor_solution

program annealing_mpi
  use mpi
  use constants
  use VariablesModule
  implicit none
  
  integer :: ierr, rank, num_cpu, i,j
  real(8), allocatable :: thetas(:,:),thetas_mpi(:,:)
  real(8), allocatable :: Es(:), Bs(:),Es_mpi(:),Bs_mpi(:),best_solution(:)
  real(8) :: B, B_min,B_max,dB, best_energy
  integer :: num_B,num_theta

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_cpu, ierr)
  
  call ReadVariables

  if (rank==0) then
    PRINT *, 'In PrintVariables Subroutine:'
    PRINT *, 'Nls  = ', Nls
    PRINT *, 'Bj   = ', Bj
    PRINT *, 'Bk   = ', Bk
    PRINT *, 'Jt   = ', Jt
    PRINT *, 'Kt   = ', Kt
    PRINT *, 'Bmin = ', Bmin
    PRINT *, 'Bmax = ', Bmax
    PRINT *, 'NBs  = ', NBs
    PRINT *, 'T_init  = ', T_init
    PRINT *, 'T_min  = ', T_min
    PRINT *, 'alpha  = ', alpha
    PRINT *, 'iterations  = ', iterations
  endif

  num_B=NBs
  num_theta=Nls
  B_min=Bmin
  B_max=Bmax
  allocate(thetas_mpi(num_theta,num_B))
  allocate(Bs_mpi(num_B),Es_mpi(num_B))
  allocate(thetas(num_theta,num_B))
  allocate(Bs(num_B),Es(num_B),best_solution(num_theta))

  dB=(B_max-B_min)/(num_B-1)
  Bs_mpi=0d0; Es_mpi=0d0; thetas_mpi=0d0
  Bs=0d0; Es=0d0; thetas=0d0
  ! 计算模拟
  do i = 1+rank, num_B, num_cpu
    B = B_min + (i-1)*dB
    Bs_mpi(i) = B
    call simulated_annealing(B, best_solution, best_energy)
    thetas_mpi(:, i) = best_solution*180/pi
    Es_mpi(i) = best_energy
  end do

  ! 收集数据
  call MPI_ALLREDUCE(thetas_mpi, thetas, size(thetas), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(Es_mpi, Es, num_B, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(Bs_mpi, Bs, num_B, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! 主进程保存结果
  if (rank == 0) then

    open(unit=20, file='thetas.txt', status='unknown')
    do i = 1, num_B
      do j = 1, Nls
      if (j.eq.Nls)then
      write(20, '(1(1PE20.12))') thetas(j, i)
      else
      write(20, '(1(1PE20.12))',ADVANCE='NO') thetas(j, i)
      endif
    end do
    end do
    close(20)

    open(unit=40, file='Bs.txt', status='unknown')
    do i =1 ,num_B
       write(40, '(1(1PE20.12))') Bs(i)
    end do
    close(40)
  endif

  call MPI_FINALIZE(ierr)
  deallocate(Bs_mpi,Es_mpi,thetas_mpi,Bs,Es,thetas,best_solution)
end program annealing_mpi

