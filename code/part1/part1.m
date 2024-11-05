clc;
clear;
% 参数设置
a = 0;  % 阿基米德螺线初始半径
b = 0.55/(2*pi);  % 螺线间距
theta_1 = linspace(0*pi, 32*pi, 400);  % 螺线的角度范围

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


x_4 = zeros(1, 223);
y_4 = zeros(1, 223);



%%
% 每个时间点
a = 0;  % 阿基米德螺线初始半径
b = 0.55/(2*pi);  % 螺线间距
for t = 300:300
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

    table(1,t+1) = x_3(1);
    table(2,t+1) = y_3(1);



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
        table(3,t+1) = x_3;
        table(4,t+1) = y_3;
    end
    
    
x_4(2) = x_3;
y_4(2) = y_3;
    
    
    
    %%
    %判定每段位置
    c = 0;
    for n = 2:5
    
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
        
        

        table(2*n+1,t+1) = x_3;
        table(2*n+2,t+1) = y_3;

        x_4(n+1) = x_3;
        y_4(n+1) = y_3;
    end


    end



x_5 = x_4
y_5 = y_4





%%
% 绘制阿基米德螺线
figure; 
plot(x_1, y_1, 'b', 'LineWidth', 1.5);  % 螺线图像
hold on;

% 绘制偏移后的圆
plot(x_2, y_2, 'r', 'LineWidth', 1.5);  % 偏移圆的图像

x_3 = (a + b * theta_00) .* cos(theta_00);
y_3 = (a + b * theta_00) .* sin(theta_00);

% 绘制交点
plot(x_5, y_5, '--' , 'LineWidth', 1.5);
plot(x_4, y_4, 'o' , 'LineWidth', 2);
plot(x_3, y_3, 'go' , 'MarkerSize', 7);

% 调整图像比例
axis equal;
title('Intersection of Archimedean Spiral and Offset Circle');
xlabel('x');
ylabel('y');
grid on;

% 显示图例
legend('Archimedean Spiral', 'Offset Circle', 'Intersection Points');

hold off;


filename = 'result1.xlsx';                    % 文件名
writematrix(table, filename);   % 将 A 写入 Excel 文件



