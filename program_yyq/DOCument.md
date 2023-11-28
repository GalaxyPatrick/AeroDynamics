# FluidSolver_23012

**翼型**：NACA23012

**NACA23012.m**：根据公式求解翼型并对翼型进行离散化

- m, k1, t 为 NACA23012 固定参数
- cpoint 为前后区域分割点位置
- gap_smallx 为翼型前端离散点分布间隔
- gap_largex 为翼型后端离散点分布间隔

**solver_NACA23012.m**：求解关于 $\gamma$ 的线性方程组，得到各控制点的涡强，同时给出 Zhukovsky 公式计算得到的升力系数
$$
\begin{cases}
\displaystyle{\left[ -\nabla_i \ \sum_{j \neq i}\frac{\gamma_j}{2\pi} \int_{s_j} \arctan \left( \frac{y_i - y_j}{x_i - x_j} \right) \mathrm{d}s_j + \boldsymbol{v}_{\mathrm{\infty}} \right] \cdot \boldsymbol{n}_i = 0} &\quad,\quad (i=1,\cdots,N-1) \\[12pt]
\gamma_1 + \gamma_N = 0 &\quad,\quad (\text{Kutta's condition})
\end{cases}
$$

$$
\text{Zhukovsky's formula:\quad\quad}\varGamma = \sum_{i} \gamma_i s_i \quad,\quad c_L = \frac{2\varGamma}{|\boldsymbol{v}_{\mathrm{\infty}}|}
$$

*（认为涡顺时针转为正值，故而诱导的速度是顺时针旋转90度，面元排列顺序也为顺时针排列，环量积分也是沿顺时针进行）*

- unit 为控制点位置
- number 为控制点（面元）数目
- tangent 为每个面元的切向量（未归一化）
- deltas 为每个面元的长度
- normal 为每个面元的法向量
- alpha 为攻角（单位：度）
- v_infty 为来流速度

**cp_calculator.m**：根据各控制点的涡强计算各控制点的压力系数，并通过压力积分的方式给出升力系数
$$
v_{ti} = \left[ \frac{\gamma_i}{2}-\nabla_i \ \sum_{j \neq i}\frac{\gamma_j}{2\pi} \int_{s_j} \arctan \left( \frac{y_i - y_j}{x_i - x_j} \right) \mathrm{d}s_j + \boldsymbol{v}_{\mathrm{\infty}} \right] \cdot \boldsymbol{\tau}_i
$$

$$
c_{pi} = 1 - \left( \frac{v_{ti}}{|\boldsymbol{v}_{\mathrm{\infty}}|} \right)^2 \quad,\quad c_L = \sum_i c_{pi} s_i (\boldsymbol{\tau}_i \cdot \boldsymbol{e}_y)
$$

- gamma 为各控制点的压力系数
- unit，tangent，deltas，normal，number，v_infty，alpha 与 solver_NACA23012.m 中意义相同
- x_d_c 为控制点在弦上对应的位置

**mainlab.m**：启动文件

**coord**：存储翼型控制点位置、控制点涡强值

- x_u_c.csv，x_l_c.csv，y_u_c.csv，y_l_c.csv 记录控制点坐标
- x_d_c.csv 记录控制点在弦上对应的位置
- unit.csv，deltas.csv，tangent.csv，normal.csv，gamma.csv 记录内容与程序中同名全局变量相同
- cl.txt 记录攻角与对应的升力系数

___

现有参数下，程序在正攻角较高（接近10度）和负攻角时表现不是很好。在这两种情况下的问题主要是：

1. 两种方法计算得到的差别不小
2. 至少有一种方法的结果与 $\mathrm{Re = 10^6}$ 的实验数据契合不是很好

另一方面，面元选得过密似乎会出现解算出错的问题，解算出的值会偏离，目前无法采用过密的网格

___

