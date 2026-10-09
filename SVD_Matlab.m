%% ========================================================================
%  奇异值分解（SVD）原理详解
%  核心思想：A = U * S * V'
%     - V' ：输入空间中的旋转（正交变换）
%     - S  ：沿坐标轴的缩放（奇异值）
%     - U  ：输出空间中的旋转（正交变换）
%  几何意义：任意线性变换 = 旋转 + 缩放 + 旋转
% ========================================================================
clear; clc; close all;
rng(42);  % 固定随机种子，便于复现

%% ========== 第一部分：手动实现 SVD（基于特征值分解原理）==========
fprintf('==================== 第一部分：手动实现 SVD ====================\n');

% 构造一个示例矩阵（非方阵，更能体现 SVD 的普适性）
A = [3, 1, 1;
     1, 3, 1;
     -1, 1, 3;
     2, 0, 1];
[m, n] = size(A);
fprintf('矩阵 A 尺寸: %d x %d\n', m, n);
disp(A);

%% --- 步骤 1：计算 A'*A 的特征值分解 ---
% 原理：A'A 是对称半正定矩阵，可正交对角化：A'A = V * Λ * V'
%       其中 Λ 的特征值 = 奇异值的平方，V 的列向量 = 右奇异向量
fprintf('\n--- 步骤1：对 A''A 做特征值分解得到 V 和奇异值 ---\n');

AtA = A' * A;
[V_manual, D_manual] = eig(AtA);

% eig 返回的特征值未排序，需要按降序排列
[eigenvalues, idx] = sort(diag(D_manual), 'descend');
V_manual = V_manual(:, idx);

% 奇异值 = 特征值的平方根（数值上要防止负的舍入误差）
singular_values = sqrt(max(eigenvalues, 0));

fprintf('A''A 的特征值（降序）: ');
fprintf('%.6f  ', eigenvalues);
fprintf('\n对应的奇异值:          ');
fprintf('%.6f  ', singular_values);
fprintf('\n');

%% --- 步骤 2：由 A*V = U*S 反推左奇异向量 U ---
% 原理：A * v_i = σ_i * u_i  =>  u_i = A * v_i / σ_i
fprintf('\n--- 步骤2：由 u_i = A*v_i/σ_i 计算左奇异向量 U ---\n');

U_manual = zeros(m, n);
for i = 1:n
    if singular_values(i) > 1e-12   % 避免除以零
        U_manual(:, i) = A * V_manual(:, i) / singular_values(i);
    else
        U_manual(:, i) = 0;         % 零奇异值对应的方向
    end
end

% 若 m > n，还需补全 U 使其成为 m×m 正交矩阵
% （通过 Gram-Schmidt 对零空间补充正交向量）
if m > n
    U_full = zeros(m, m);
    U_full(:, 1:n) = U_manual;
    % 对剩余的列用随机向量做 Gram-Schmidt 正交化
    for j = n+1:m
        v = randn(m, 1);
        % 减去已在前面的投影
        for k = 1:j-1
            v = v - (U_full(:,k)' * v) * U_full(:,k);
        end
        v = v / norm(v);
        U_full(:, j) = v;
    end
    U_manual = U_full;
end

% 构造对角矩阵 S
S_manual = zeros(m, n);
for i = 1:min(m,n)
    S_manual(i, i) = singular_values(i);
end

fprintf('手动计算的 U:\n'); disp(U_manual);
fprintf('手动计算的 S:\n'); disp(S_manual);
fprintf('手动计算的 V:\n'); disp(V_manual);

%% --- 步骤 3：验证 A = U * S * V' ---
fprintf('\n--- 步骤3：验证 A = U*S*V'' ---\n');
A_recon = U_manual * S_manual * V_manual';
fprintf('重构矩阵:\n'); disp(A_recon);
fprintf('重构误差 ||A - U*S*V''||_F = %.4e\n', norm(A - A_recon, 'fro'));

%% --- 步骤 4：验证正交性 ---
fprintf('\n--- 步骤4：验证 U、V 的正交性 ---\n');
fprintf('||U''U - I||_F = %.4e\n', norm(U_manual'*U_manual - eye(m), 'fro'));
fprintf('||V''V - I||_F = %.4e\n', norm(V_manual'*V_manual - eye(n), 'fro'));

%% ========== 第二部分：与 MATLAB 内置 svd 对比 ==========
fprintf('\n==================== 第二部分：与内置 svd 对比 ====================\n');
[U_builtin, S_builtin, V_builtin] = svd(A);

fprintf('内置 svd 的奇异值: ');
fprintf('%.6f  ', diag(S_builtin)');
fprintf('\n手动 的奇异值:     ');
fprintf('%.6f  ', singular_values);
fprintf('\n');

% 注意：奇异向量可能相差一个符号（因为 ±u 都满足 SVD）
% 这里用绝对值比较
fprintf('U 与内置 U 的列方向差异（应全为 ±1）:\n');
sign_diff_U = sign(diag(U_manual(:, 1:min(m,n))' * U_builtin(:, 1:min(m,n))));
disp(sign_diff_U');

%% ========== 第三部分：SVD 的几何意义可视化 ==========
fprintf('\n==================== 第三部分：SVD 几何意义演示 ====================\n');
fprintf('演示：单位圆 --V''--> --S--> --U--> 椭圆\n');

% 用 2x2 矩阵，方便可视化
A2 = [2, 1;
      0.5, 1.5];
[U2, S2, V2] = svd(A2);

% 生成单位圆上的点
theta = linspace(0, 2*pi, 200);
circle = [cos(theta); sin(theta)];

% 逐步变换
step1 = V2' * circle;          % 第一步：右奇异向量旋转
step2 = S2 * step1;            % 第二步：奇异值缩放（圆变椭圆）
step3 = U2 * step2;            % 第三步：左奇异向量旋转
direct = A2 * circle;          % 直接应用 A

figure('Position', [100 100 1200 300]);

% 原始单位圆
subplot(1,4,1);
plot(circle(1,:), circle(2,:), 'b-', 'LineWidth', 2); hold on;
quiver(0,0,V2(1,1),V2(2,1),'r','LineWidth',2,'MaxHeadSize',0.5);
quiver(0,0,V2(1,2),V2(2,2),'g','LineWidth',2,'MaxHeadSize',0.5);
axis equal; grid on; axis([-3 3 -3 3]);
title('① 单位圆 + V 的列向量');
xlabel('x'); ylabel('y');

% 第一步：V' 旋转
subplot(1,4,2);
plot(step1(1,:), step1(2,:), 'b-', 'LineWidth', 2); hold on;
axis equal; grid on; axis([-3 3 -3 3]);
title('② 经 V'' 旋转');
xlabel('x'); ylabel('y');

% 第二步：S 缩放
subplot(1,4,3);
plot(step2(1,:), step2(2,:), 'm-', 'LineWidth', 2); hold on;
quiver(0,0,S2(1,1),0,'r','LineWidth',2,'MaxHeadSize',0.5);
quiver(0,0,0,S2(2,2),'g','LineWidth',2,'MaxHeadSize',0.5);
axis equal; grid on; axis([-3 3 -3 3]);
title(sprintf('③ 经 S 缩放 (σ₁=%.2f, σ₂=%.2f)', S2(1,1), S2(2,2)));
xlabel('x'); ylabel('y');

% 第三步：U 旋转（等价于 A * circle）
subplot(1,4,4);
plot(step3(1,:), step3(2,:), 'r-', 'LineWidth', 2); hold on;
plot(direct(1,:), direct(2,:), 'k--', 'LineWidth', 1.5);  % 验证一致
axis equal; grid on; axis([-3 3 -3 3]);
title('④ 经 U 旋转 = A·circle');
legend('U*S*V''*circle', 'A*circle', 'Location', 'best');
xlabel('x'); ylabel('y');

sgtitle('SVD 几何意义：A = U \cdot S \cdot V^T  = 旋转 → 缩放 → 旋转');

%% ========== 第四部分：SVD 的应用——低秩近似 ==========
fprintf('\n==================== 第四部分：低秩近似 ====================\n');

% 生成一个低秩 + 噪声的矩阵
m2 = 20; n2 = 15; true_rank = 3;
U_true = orth(randn(m2, true_rank));
V_true = orth(randn(n2, true_rank));
Sigma_true = diag([10, 5, 2]);
A_noisy = U_true * Sigma_true * V_true' + 0.1 * randn(m2, n2);

[U_n, S_n, V_n] = svd(A_noisy);
sv = diag(S_n);

fprintf('奇异值（前 8 个）: ');
fprintf('%.4f  ', sv(1:min(8,end)));
fprintf('\n');

% 用前 k 个奇异值做近似，观察误差
figure('Position', [100 100 900 350]);
subplot(1,2,1);
bar(sv); grid on;
title('奇异值分布'); xlabel('索引 k'); ylabel('σ_k');

subplot(1,2,2);
errors = zeros(length(sv), 1);
for k = 1:length(sv)
    A_k = U_n(:,1:k) * S_n(1:k,1:k) * V_n(:,1:k)';
    errors(k) = norm(A_noisy - A_k, 'fro') / norm(A_noisy, 'fro');
end
plot(errors, 'o-', 'LineWidth', 2); grid on;
title('秩-k 近似相对误差'); xlabel('k'); ylabel('相对误差');
sgtitle('SVD 低秩近似：奇异值衰减越快，压缩效果越好');

fprintf('秩-3 近似相对误差: %.4f\n', errors(3));
fprintf('秩-5 近似相对误差: %.4f\n', errors(5));

%% ========== 第五部分：SVD 与特征值分解的关系 ==========
fprintf('\n==================== 第五部分：SVD 与特征值分解的关系 ====================\n');
fprintf('A''A 的非零特征值 = A 的奇异值的平方\n');
fprintf('AA'' 的非零特征值 = A 的奇异值的平方\n\n');

eig_AtA = sort(eig(A'*A), 'descend');
eig_AAt = sort(eig(A*A'), 'descend');
sv_sq   = singular_values.^2;

fprintf('A''A  的特征值: ');
fprintf('%.6f  ', eig_AtA');
fprintf('\nAA''  的特征值: ');
fprintf('%.6f  ', eig_AAt');
fprintf('\n奇异值的平方:   ');
fprintf('%.6f  ', sv_sq);
fprintf('\n');

% 关系总结
fprintf('\n【核心结论】\n');
fprintf('  1) A 的非零奇异值 = sqrt(eig(A''A)) = sqrt(eig(AA''))\n');
fprintf('  2) V 的列 = A''A 的特征向量（右奇异向量）\n');
fprintf('  3) U 的列 = AA'' 的特征向量（左奇异向量）\n');
fprintf('  4) 秩(A) = 非零奇异值个数\n');
fprintf('  5) ||A||_2 = σ_max,  ||A||_F = sqrt(sum(σ_i^2))\n');
fprintf('  6) cond(A) = σ_max / σ_min\n');
fprintf('  本矩阵: 秩=%d, 谱范数=%.4f, 条件数=%.4f\n', ...
        sum(singular_values > 1e-10), singular_values(1), ...
        singular_values(1)/singular_values(end));