clc;
clear;

% 螺距
p = 0.55;
% 圈数
n = 16;

% 计算螺线的参数范围
theta_end = n * 2 * pi;
theta = linspace(0, theta_end, 1000);

% 计算螺线的极坐标
r = p * theta / (2 * pi);

% 将极坐标转换为笛卡尔坐标
x = r .* cos(theta);
y = r .* sin(theta);

% 绘制阿基米德螺线
figure;
plot(x, y, 'b');
hold on;

% 计算螺线的长度
dx_dtheta = @(theta) p / (2 * pi) * cos(theta) - (p * theta / (2 * pi)) .* sin(theta);
dy_dtheta = @(theta) p / (2 * pi) * sin(theta) + (p * theta / (2 * pi)) .* cos(theta);

integrand = @(theta) sqrt(dx_dtheta(theta).^2 + dy_dtheta(theta).^2);
L = integral(integrand, 0, theta_end);
fprintf('阿基米德螺线的长度为: %.2f 米\n', L);

% 运动速度
v = 1; % 米每秒

% 总时间
total_time = 300; % 秒

% 预先分配存储位置的数组
x_positions = zeros(total_time + 1, 1);
y_positions = zeros(total_time + 1, 1);
r_positions = zeros(total_time + 1, 1);
theta_positions = zeros(total_time + 1, 1);

% 计算每一秒的位置
for t = 0:total_time - 1
    distance_traveled = v * t;
    % 计算对应的角度 theta
    theta_func = @(s) integral(integrand, s, theta_end) - distance_traveled;
    theta_t = fzero(theta_func, theta_end / 2); % 估计初值为 theta_end 的一半

    % 计算在该时刻的极坐标
    r_t = p * theta_t / (2 * pi);

    % 计算在该时刻的笛卡尔坐标
    x_t = r_t * cos(theta_t);
    y_t = r_t * sin(theta_t);

    % 存储位置
    x_positions(t + 1) = x_t;
    y_positions(t + 1) = y_t;
    r_positions(t + 1) = r_t;
    theta_positions(t + 1) = theta_t;
end

% 打印每一秒的位置的极坐标
fprintf('时间 (秒)\t r (米)\t\t θ (弧度)\n');
for t = 0:total_time
    fprintf('%d\t\t %.2f\t\t %.2f\n', t, r_positions(t + 1), theta_positions(t + 1));
end

% 在图中标出每一秒的位置
plot(x_positions, y_positions, 'ro', 'MarkerSize', 5, 'LineWidth', 1);
for t = 1:total_time
    if mod(t, 50) == 0 % 每50秒标注一次时间
        text(x_positions(t+1), y_positions(t+1), sprintf('  %d s', t), 'FontSize', 8);
    end
end

% 标注图例和标题
legend('阿基米德螺线', '物体位置');
title('阿基米德螺线上的运动物体 (每秒标注位置)');
xlabel('X坐标');
ylabel('Y坐标');
axis equal;
grid on;






