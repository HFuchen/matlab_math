clc;
clear;
% 参数设置
a = 0;  % 阿基米德螺线初始半径
b = 0.55/(2*pi);  % 螺线间距
theta_1 = linspace(0*pi, 32*pi, 4000);  % 螺线的角度范围

% 阿基米德螺线方程
r_1 = a + b * theta_1;

% 将螺线的极坐标转换为直角坐标
x_1 = r_1 .* cos(theta_1);
y_1 = r_1 .* sin(theta_1);

% 圆的参数
R = 2.86;  % 圆的半径
theta_2 = linspace(0, 2*pi, 100);  % 圆的角度范围

% 偏移后的圆心坐标
x_2 = 4.4 + R * cos(theta_2);  % 圆心沿 x 轴偏移
y_2 = 2.4 + R * sin(theta_2);  % 圆心在 y 轴方向无偏


x_4 = zeros(1, 222);
y_4 = zeros(1, 222);
x_4h = zeros(1, 223);
x_4h = zeros(1, 223);

figure;

%%
% 每个时间点
a = 0;  % 阿基米德螺线初始半径
b = 0.55/(2*pi);  % 螺线间距
for t = 423:426

     cla;

    theta = 0;
    theta_min = 1000;
    for i = 0:4000000
        theta = theta + 0.0001;
        S = 0.0876*(1/2 * theta * sqrt(1 + theta^2) + 1/2 * log(abs(theta + sqrt(1 + theta^2))));
        if abs(S - 442.91 + t) < 0.001
            if abs(S - 442.91 + t) < theta_min
                theta_min = theta;
            end
        end
    end
    display(theta_min)

    x_3(1) = (a + b * theta_min) .* cos(theta_min);
    y_3(1) = (a + b * theta_min) .* sin(theta_min);

    table(1,1) = x_3(1);
    table(1,2) = y_3(1);



x_4(1) = x_3;
y_4(1) = y_3;




    %%
    %判定龙头位置
    a = 0;  % 阿基米德螺线初始半径
    b = 0.55/(2*pi);  % 螺线间距
    
    theta_00 =  theta_min;
    
    for n = 1:1
    
        theta_now = theta_00;
        theta_22 = theta_00;
        theta_11 = theta_00;
        theta_aa = 1000;
        
        for j = 0:3000
        
            theta_11 = theta_11 + 0.01;
        
            for i = 0:700
        
                theta_22 = theta_22 + 0.01;
            
                aa = 0;  
                bb = 0.55/(2*pi); 
                
                r_11 = aa + bb * theta_11;
                x_11 = r_11 .* cos(theta_11);
                y_11 = r_11 .* sin(theta_11);
                
                R = 2.86; 
                
                x_22 = x_3 + R * cos(theta_22);
                y_22 = y_3 + R * sin(theta_22);
            
                if abs(x_11 - x_22) < 0.1 && abs(y_11 - y_22) < 0.1
                    if abs(theta_11 - theta_now) < theta_aa && theta_11 - theta_now > 0
                        theta_aa = abs(theta_11 - theta_now);
                        theta_00 = theta_11;
                    end
                end
            
            end
        end
    
        display(theta_00)
    
        x_3 = (a + b * theta_00) .* cos(theta_00);
        y_3 = (a + b * theta_00) .* sin(theta_00);
        table(2,1) = x_3;
        table(2,2) = y_3;
    end
    
    
x_4(2) = x_3;
y_4(2) = y_3;
x_4h(1) = x_3;
y_4h(1) = y_3;
    
    
    
    %%
    %判定每段位置
    c = 0;
    for n = 2:223
    
        theta_now = theta_00;
        theta_22 = theta_00;
        theta_11 = theta_00;
        theta_aa = 1000;
        
        for j = 0:3000
        
            theta_11 = theta_11 + 0.001;
        
            for i = 0:700
        
                theta_22 = theta_22 + 0.01;
            
                aa = 0;  
                bb = 0.55/(2*pi); 
                
                r_11 = aa + bb * theta_11;
                x_11 = r_11 .* cos(theta_11);
                y_11 = r_11 .* sin(theta_11);
                
                R = 1.65; 
                
                x_22 = x_3 + R * cos(theta_22);
                y_22 = y_3 + R * sin(theta_22);
            
                if abs(x_11 - x_22) < 0.25 && abs(y_11 - y_22) < 0.25
                    if abs(theta_11 - theta_now) < theta_aa && theta_11 - theta_now > 0
                        theta_aa = abs(theta_11 - theta_now);
                        theta_00 = theta_11;
                    end
                end
            
            end
        end
    
        display(theta_00)
    
        x_3 = (a + b * theta_00) .* cos(theta_00);
        y_3 = (a + b * theta_00) .* sin(theta_00);
        
        

        table(n+1,1) = x_3;
        table(n+1,2) = y_3;
        if (n+1) < 223
            x_4(n+1) = x_3;
            y_4(n+1) = y_3;
        end
            x_4h(n) = x_3;
            y_4h(n) = y_3;

    end






x_5 = x_4;
y_5 = y_4;





%%
% 绘制阿基米德螺线
%figure; 
%plot(x_1, y_1, 'b', 'LineWidth', 1.5);  % 螺线图像
%hold on;

% 绘制偏移后的圆
%plot(x_2, y_2, 'r', 'LineWidth', 1.5);  % 偏移圆的图像

%x_3 = (a + b * theta_00) .* cos(theta_00);
%y_3 = (a + b * theta_00) .* sin(theta_00);

% 绘制交点
%plot(x_5, y_5, '-' , 'LineWidth', 20);
%plot(x_4, y_4, 'o' , 'LineWidth', 1.2);
%plot(x_3, y_3, 'go' , 'MarkerSize', 7);



%    t1 = 23.3599
%    t2 = 24.1747
%    x_n = (a + b * t1) .* cos(t1);
%    y_n = (a + b * t1) .* sin(t1);
%    x_m = (a + b * t2) .* cos(t2);
%    y_m = (a + b * t2) .* sin(t2);
%    plot(x_m, y_m, 'o' , 'LineWidth', 20);
%    plot(x_n, y_n, 'o' , 'LineWidth', 20);










%filename = 'result1.xlsx';                    % 文件名
%writematrix(table, filename);   % 将 A 写入 Excel 文件



%%



% 假设交点坐标为 (x_inter, y_inter)
x_inter = 0;
y_inter = 0;

% 假设线段端点坐标存储在两个数组中，x1, y1为线段起点，x2, y2为线段终点
x1 = x_4; % 线段起点x坐标
y1 = y_4;    % 线段起点y坐标
x2 = x_4h;   % 线段终点x坐标
y2 = y_4h;    % 线段终点y坐标

% 线段数量
n_segments = length(x1);

% 创建新图像

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
    x1_close = x1_extended - 0.15 * dy_unit;
    y1_close = y1_extended + 0.15 * dx_unit;
    x2_close = x2_extended - 0.15 * dy_unit;
    y2_close = y2_extended + 0.15 * dx_unit;

    % 远离点15cm的线段端点
    x1_far = x1_extended + 0.15 * dy_unit;
    y1_far = y1_extended - 0.15 * dx_unit;
    x2_far = x2_extended + 0.15 * dy_unit;
    y2_far = y2_extended - 0.15 * dx_unit;


    %x_n = (a + b * t1) .* cos(t1);
    %y_n = (a + b * t1) .* sin(t1);
    %x_m = (a + b * t2) .* cos(t2);
    %y_m = (a + b * t2) .* sin(t2);
    %plot(x_m, y_m, 'o' , 'LineWidth', 20);
    %plot(x_n, y_n, 'o' , 'LineWidth', 20);

    % 画出增长后的线段
    %line([x1_extended, x2_extended], [y1_extended, y2_extended], 'Color', 'r');

    plot(x_1, y_1, 'b', 'LineWidth', 0.02);  % 螺线图像
    hold on;

    % 画出靠近点15cm的线段
    line([x1_close, x2_close], [y1_close, y2_close], 'Color', 'b');

    % 画出远离点15cm的线段
    line([x1_far, x2_far], [y1_far, y2_far], 'Color', 'g');

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

drawnow;

end





