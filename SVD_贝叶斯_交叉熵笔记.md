# 三大主题笔记：SVD、贝叶斯过程、交叉熵

> 以下三部分内容均为手写推导 + 自编数值算例，矩阵用 $\LaTeX$ 表示，便于直接复习使用。

---

## 一、奇异值分解（SVD）

### 1.1 定理陈述

任意矩阵 $A \in \mathbb{R}^{m \times n}$（不要求方阵、不要求可逆），都可以分解为

$$
A = U\Sigma V^{T}
$$

其中：

| 矩阵 | 维度 | 性质 |
|---|---|---|
| $U$ | $m \times m$ | 左奇异向量，$U^{T}U = I_m$，列为 $AA^{T}$ 的特征向量 |
| $\Sigma$ | $m \times n$ | 对角矩阵，对角元 $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_r > 0$，$r = \mathrm{rank}(A)$ |
| $V$ | $n \times n$ | 右奇异向量，$V^{T}V = I_n$，列为 $A^{T}A$ 的特征向量 |

### 1.2 分解步骤（算法流程）

**Step 1：计算 $A^{T}A$**

$A^{T}A$ 是 $n \times n$ 对称半正定矩阵，其特征值 $\lambda_i \geq 0$。

**Step 2：求 $A^{T}A$ 的特征值与特征向量，得到 $V$ 与 $\Sigma$**

$$
\sigma_i = \sqrt{\lambda_i}, \qquad v_i \text{ 为 } A^{T}A \text{ 对应 } \lambda_i \text{ 的单位特征向量}
$$

将 $V = [\,v_1,\ v_2,\ \cdots,\ v_n\,]$ 按 $\sigma_i$ 从大到小排列。

**Step 3：计算 $U$ 的前 $r$ 列**

$$
u_i = \frac{Av_i}{\sigma_i}, \qquad i = 1, 2, \dots, r
$$

**Step 4：补全 $U$（若 $r < m$）**

用 $AA^{T}$ 的零特征值对应的特征向量（或对 $\{u_1,\dots,u_r\}$ 做 Gram–Schmidt 正交补）填满 $U$ 的剩余列。

**Step 5：按** $A = U\Sigma V^{T}$ **拼出分解并回代验证。**

> 几何意义：$V^T$ 旋转/反射输入空间 → $\Sigma$ 沿各正交方向拉伸 $\sigma_i$ 倍 → $U$ 旋转/反射到输出空间。即"任何线性变换 = 旋转 + 伸缩 + 再旋转"。

### 1.3 自编算例

设

$$
A = \begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{bmatrix}_{3\times 2}
$$

**(1) 计算 $A^{T}A$：**

$$
A^{T}A = \begin{bmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \end{bmatrix}\begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{bmatrix} = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}
$$

**(2) 求特征值：** 特征方程 $\det(A^TA - \lambda I) = (2-\lambda)^2 - 1 = 0$，得

$$
\lambda_1 = 3, \qquad \lambda_2 = 1 \;\Rightarrow\; \sigma_1 = \sqrt{3},\ \sigma_2 = 1
$$

**(3) 求 $V$：**

- $\lambda_1 = 3$：解 $(A^TA - 3I)v = 0$，得 $v_1 = \dfrac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ 1 \end{bmatrix}$
- $\lambda_2 = 1$：解 $(A^TA - I)v = 0$，得 $v_2 = \dfrac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ -1 \end{bmatrix}$

$$
V = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 & 1 \\ 1 & -1 \end{bmatrix}
$$

**(4) 求 $U$ 的前两列：**

$$
u_1 = \frac{Av_1}{\sigma_1} = \frac{1}{\sqrt{3}}\cdot\frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ 1 \\ 2 \end{bmatrix} = \frac{1}{\sqrt{6}}\begin{bmatrix} 1 \\ 1 \\ 2 \end{bmatrix}, \qquad
u_2 = \frac{Av_2}{\sigma_2} = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ -1 \\ 0 \end{bmatrix}
$$

**(5) 补第三列 $u_3$：** $AA^{T}$ 的特征值含 $0$，对应零空间方向满足 $x_1 + x_3 = 0,\ x_2 = 0$，取

$$
u_3 = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 \\ 0 \\ -1 \end{bmatrix}
$$

**(6) 最终分解：**

$$
\boxed{
A = \underbrace{\begin{bmatrix} \tfrac{1}{\sqrt{6}} & \tfrac{1}{\sqrt{2}} & \tfrac{1}{\sqrt{2}} \\[4pt] \tfrac{1}{\sqrt{6}} & -\tfrac{1}{\sqrt{2}} & 0 \\[4pt] \tfrac{2}{\sqrt{6}} & 0 & -\tfrac{1}{\sqrt{2}} \end{bmatrix}}_{U}
\underbrace{\begin{bmatrix} \sqrt{3} & 0 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}}_{\Sigma}
\underbrace{\frac{1}{\sqrt{2}}\begin{bmatrix} 1 & 1 \\ 1 & -1 \end{bmatrix}}_{V^{T}}
}
$$

**(7) 回代验证（第一列）：**

$$
\sigma_1 u_1 v_{11} + \sigma_2 u_2 v_{21} = \tfrac{1}{2}\begin{bmatrix}1\\1\\2\end{bmatrix} + \tfrac{1}{2}\begin{bmatrix}1\\-1\\0\end{bmatrix} = \begin{bmatrix}1\\0\\1\end{bmatrix} \checkmark
$$

第二列同理可得 $(0,1,1)^T$，与 $A$ 完全一致（已用 NumPy `linalg.svd` 数值复核，符号约定不同但分解等价）。

---

## 二、贝叶斯过程：从建模到后验推断的完整案例

### 2.1 贝叶斯推断的核心流程

贝叶斯推断的全部过程可概括为一条链：

$$
\underbrace{p(\theta)}_{\text{先验}} \;\xrightarrow{\;\text{乘以似然}\;}\; \underbrace{p(\theta \mid D)}_{\text{后验}} \;\propto\; \underbrace{p(D \mid \theta)}_{\text{似然}}\, \underbrace{p(\theta)}_{\text{先验}}
$$

**五步建模法：**

1. **定参数（Parameterization）**：明确未知量 $\theta$ 的含义；
2. **定先验（Prior）**：根据历史经验或领域知识给出 $p(\theta)$；
3. **定似然（Likelihood）**：写出数据生成机制 $p(D \mid \theta)$；
4. **求后验（Posterior）**：贝叶斯定理 + （共轭/数值）推断；
5. **做决策**：用后验均值、后验区间、后验预测分布输出结论。

### 2.2 实际案例：App 新功能按钮的点击率（CTR）评估

**项目背景：** 我在项目中负责评估一个新按钮的点击率。直接上线全量用户风险高，先对小流量桶做实验：投放 $n = 30$ 个用户，观察到 $k = 12$ 次点击。但该产品历史 CTR 约为 $20\%$（老版本长期数据），样本量小，直接用 $\hat\theta = 12/30 = 0.4$ 做决策波动太大。用贝叶斯方法把"历史经验"作为先验、把"小样本实验"作为似然，融合得到后验。

#### 2.2.1 第一步：参数化

未知参数 $\theta \in (0,1)$：新按钮的真实点击率。

#### 2.2.2 第二步：先验建模——Beta 分布

历史 CTR 约 $20\%$，且我们对这个先验的"置信强度"相当于已观察过约 $10$ 个用户（其中 $2$ 次点击）。选 Beta 分布：

$$
p(\theta) = \mathrm{Beta}(\theta; \alpha_0, \beta_0) = \frac{\theta^{\alpha_0 - 1}(1-\theta)^{\beta_0 - 1}}{B(\alpha_0, \beta_0)}, \qquad \alpha_0 = 2,\ \beta_0 = 8
$$

先验均值 $E[\theta] = \dfrac{\alpha_0}{\alpha_0+\beta_0} = 0.2$，与历史一致。

#### 2.2.3 第三步：似然建模——二项分布

30 个用户相互独立，每人点击概率为 $\theta$，故

$$
p(D \mid \theta) = \binom{n}{k}\theta^{k}(1-\theta)^{n-k} = \binom{30}{12}\theta^{12}(1-\theta)^{18}
$$

#### 2.2.4 第四步：求后验——共轭更新

由贝叶斯定理：

$$
p(\theta \mid D) \propto p(D \mid \theta)\,p(\theta) \propto \theta^{12}(1-\theta)^{18}\cdot \theta^{1}(1-\theta)^{7} = \theta^{13}(1-\theta)^{25}
$$

这正是 $\mathrm{Beta}(14, 26)$ 的核。归一化后：

$$
\boxed{\,p(\theta \mid D) = \mathrm{Beta}(\theta;\ 14,\ 26)\,}
$$

> **Beta–Binomial 共轭规律：** $\mathrm{Beta}(\alpha_0,\beta_0) + (k, n-k)$ 次二项试验 $\Rightarrow$ $\mathrm{Beta}(\alpha_0 + k,\ \beta_0 + n - k)$。可记忆为"先验伪计数 + 实验计数"。

#### 2.2.5 第五步：决策输出

| 估计量 | 数值 | 说明 |
|---|---|---|
| 先验均值 | $0.200$ | 纯历史经验 |
| MLE（频率学派） | $0.400$ | 只用 30 个样本，方差大 |
| **后验均值** | $\dfrac{14}{14+26} = \mathbf{0.350}$ | 向先验收缩，更稳健 |
| 后验标准差 | $\sqrt{\dfrac{14\cdot 26}{(40)^2 \cdot 41}} \approx 0.075$ | 不确定性量化 |
| $P(\theta > 0.2 \mid D)$ | $> 0.99$ | 新按钮显著优于旧版的概率 |

**结论：** 后验表明新按钮 CTR 大概率优于旧版 $20\%$，但后验均值 $0.35$ 比 MLE 的 $0.4$ 更保守；可据此做灰度放量决策，并继续用新数据迭代后验——这正体现了贝叶斯"序贯更新"的思想：今天实验的后验，就是明天实验的先验。

---

## 三、交叉熵：定义、原理与推导

### 3.1 定义

对离散分布 $P$（真实分布）与 $Q$（模型分布），**交叉熵**定义为

$$
H(P, Q) = -\sum_{x} P(x)\,\log Q(x)
$$

在分类问题中标签常为 one-hot，对单样本退化为 $-\log q_{y}$（$y$ 为真实类别），对整个训练集取平均即**交叉熵损失**：

$$
\mathcal{L} = -\frac{1}{N}\sum_{i=1}^{N}\Big[\,y_i \log p_i + (1-y_i)\log(1-p_i)\,\Big] \qquad (\text{二分类})
$$

### 3.2 原理：交叉熵从何而来？

交叉熵不是凭空设计的损失函数，它来自**最大似然估计**。以二分类（逻辑回归）为例：

#### Step 1：写出数据似然

设 $y_i \in \{0,1\}$，模型输出 $p_i = \sigma(w^{T}x_i) \in (0,1)$ 为 $P(y_i=1 \mid x_i)$。样本独立同分布，似然为

$$
L(w) = \prod_{i=1}^{N} p_i^{\,y_i}(1-p_i)^{\,1-y_i}
$$

#### Step 2：取对数似然

$$
\log L(w) = \sum_{i=1}^{N}\Big[\,y_i \log p_i + (1-y_i)\log(1-p_i)\,\Big]
$$

#### Step 3：最大化似然 $\Leftrightarrow$ 最小化负对数似然

$$
\hat w = \arg\max_w \log L(w) = \arg\min_w \underbrace{\Big(-\frac{1}{N}\log L(w)\Big)}_{\text{交叉熵损失 } \mathcal{L}}
$$

**结论：交叉熵 = 负对数似然取平均。** 因此最小化交叉熵在统计上等价于最大似然估计，这是它有坚实理论根基的原因。

#### 补充：与 KL 散度的关系

$$
H(P,Q) = \underbrace{H(P)}_{\text{熵（与 }Q\text{ 无关）}} + \underbrace{D_{KL}(P\,\|\,Q)}_{\text{非负，}Q=P\text{ 时为 }0}
$$

训练时 $P$（数据）固定，故最小化交叉熵 = 最小化 KL 散度 = 让模型分布逼近真实分布。

### 3.3 数值推导算例（3 个样本）

设三个样本的真实标签与模型预测概率为：

| 样本 | $y_i$ | $p_i$ | 单样本损失 $-\big[y_i\log p_i + (1-y_i)\log(1-p_i)\big]$ |
|---|---|---|---|
| 1 | $1$ | $0.9$ | $-\ln 0.9 = 0.1054$ |
| 2 | $0$ | $0.1$ | $-\ln(1-0.1) = -\ln 0.9 = 0.1054$ |
| 3 | $1$ | $0.7$ | $-\ln 0.7 = 0.3567$ |

$$
\mathcal{L} = \frac{0.1054 + 0.1054 + 0.3567}{3} = \frac{0.5675}{3} \approx \boxed{0.1891}
$$

**逐条解读：**

- 样本 1、2 预测"很自信且正确"，损失小（$0.105$）；
- 样本 3 预测 $p=0.7$ 但真实为 $1$，置信度不够，损失升至 $0.357$；
- 若某样本预测完全错误（如 $y=1, p \to 0$），则 $-\ln p \to +\infty$——交叉熵对"自信地犯错"施加**指数级惩罚**，这是它相比均方误差（MSE）更适合分类的核心原因。

再对比 MSE 在同一批数据上的表现：$\mathrm{MSE} = \frac{1}{3}\big[(1-0.9)^2 + (0-0.1)^2 + (1-0.7)^2\big] = \frac{0.11}{3} \approx 0.0367$。MSE 梯度在 $\sigma$ 饱和区趋近于零（梯度消失），而交叉熵的梯度为 $\frac{\partial \mathcal{L}}{\partial z} = p - y$（$z$ 为 logits），形式简洁且不饱和，训练效率更高。

### 3.4 应用场景小结

| 场景 | 用法 |
|---|---|
| 二分类 / 多分类 | softmax + 交叉熵，深度学习的默认分类损失 |
| 多标签分类 | 对每个标签独立用二值交叉熵（BCE） |
| 模型蒸馏 | 用教师模型输出的软标签计算"软交叉熵" |
| 生成模型 | VAE 中的重构损失即交叉熵项 |
| 与信息论的联系 | 等价于最小化 KL 散度，衡量"编码真实分布所需的多余比特数" |

---

## 结语

三者看似独立，实则统一于一条主线：**SVD 是矩阵层面的最优低秩近似（几何视角），贝叶斯推断是参数层面的最优信念更新（概率视角），交叉熵是分布层面的最优匹配准则（信息论视角）**。它们共同构成了机器学习"建模—推断—优化"三大环节的数学基础。
