clc;
clear;
% 假设交点坐标为 (x_inter, y_inter)
x_inter = 0;
y_inter = 0;

% 假设线段端点坐标存储在两个数组中，x1, y1为线段起点，x2, y2为线段终点
x1 = [2.7103 , 1.9957]; % 线段起点x坐标
y1 = [0.5618 , 1.0202];    % 线段起点y坐标
x2 = [2.1806 , -0.4254];   % 线段终点x坐标
y2 = [1.7719 , 2.3146];    % 线段终点y坐标

% 线段数量
n_segments = length(x1);

% 创建新图像
figure;
hold on; % 保持图像，以便在同一幅图上画多条线

for i = 1:n_segments
    % 计算线段的方向向量
    dx = x2(i) - x1(i);
    dy = y2(i) - y1(i);

    % 计算线段的长度
    length_segment = sqrt(dx^2 + dy^2);

    % 计算单位方向向量
    dx_unit = dx / length_segment;
    dy_unit = dy / length_segment;

    % 增长线段，每边增长27.5cm
    x1_extended = x1(i) - 0.275 * dx_unit;
    y1_extended = y1(i) - 0.275 * dy_unit;
    x2_extended = x2(i) + 0.275 * dx_unit;
    y2_extended = y2(i) + 0.275 * dy_unit;

    % 靠近点15cm的线段端点
    x1_close = x1_extended - 0.15 * dy_unit
    y1_close = y1_extended + 0.15 * dx_unit
    x2_close = x2_extended - 0.15 * dy_unit
    y2_close = y2_extended + 0.15 * dx_unit

    % 远离点15cm的线段端点
    x1_far = x1_extended + 0.15 * dy_unit
    y1_far = y1_extended - 0.15 * dx_unit
    x2_far = x2_extended + 0.15 * dy_unit
    y2_far = y2_extended - 0.15 * dx_unit

    % 画出增长后的线段
    line([x1_extended, x2_extended], [y1_extended, y2_extended], 'Color', 'r');

    % 画出靠近点15cm的线段
    line([x1_close, x2_close], [y1_close, y2_close], 'Color', 'b');

    % 画出远离点15cm的线段
    line([x1_far, x2_far], [y1_far, y2_far], 'Color', 'g');


    
    %text(x1_far,y1_far, '(x1_far,y1_far)');
    text(x1_far,y1_far,['(',num2str(x1_far),',',num2str(y1_far),')'])
    text(x2_far,y2_far,['(',num2str(x2_far),',',num2str(y2_far),')'])
    text(x1_close,y1_close,['(',num2str(x1_close),',',num2str(y1_close),')'])
    text(x2_close,y2_close,['(',num2str(x2_close),',',num2str(y2_close),')'])


    % 画出矩形
    % 连接靠近点15cm的线段端点到远离点15cm的线段端点
    line([x1_close, x1_far], [y1_close, y1_far], 'Color', 'k');
    line([x1_far, x2_far], [y1_far, y2_far], 'Color', 'k');
    line([x2_far, x2_close], [y2_far, y2_close], 'Color', 'k');
    line([x2_close, x1_close], [y2_close, y1_close], 'Color', 'k');
end

% 设置图像属性
axis equal;
grid on;
xlabel('X');
ylabel('Y');
title('t=414s时龙头发生碰撞的临界情况图示');
hold off;




