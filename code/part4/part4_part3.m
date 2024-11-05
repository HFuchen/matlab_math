clc;
clear;
% 参数设置

%建立每个板凳的比较算子
com = 1;

R = 4.3235;
%更改龙头位置
theta_min = 5.6 * pi

a = 0;  % 阿基米德螺线初始半径
b = 1.7/(2*pi);  % 螺线间距
theta_1 = linspace(R/b, 32*pi, 4000);  % 螺线的角度范围


% 阿基米德螺线方程
r_1l = a + b * theta_1;
% 将螺线的极坐标转换为直角坐标
x_1l = r_1l .* cos(theta_1);
y_1l = r_1l .* sin(theta_1);

% 阿基米德螺线方程
r_1r = a - b * theta_1;
% 将螺线的极坐标转换为直角坐标
x_1r = r_1r .* cos(theta_1);
y_1r = r_1r .* sin(theta_1);

% 圆的参数1
theta_2 = linspace(0, 2*pi, 100);  % 圆的角度范围

% 偏移后的圆心坐标1
x_2 = R * cos(theta_2);  % 圆心沿 x 轴偏移
y_2 = R * sin(theta_2);  % 圆心在 y 轴方向无偏

% 圆弧参数
%切入点
x_p1 = (a + b * R/b) .* cos(R/b);
y_p1 = (a + b * R/b) .* sin(R/b);
%切出点
x_p2 = (a - b * R/b) .* cos(R/b);
y_p2 = (a - b * R/b) .* sin(R/b);

x_q1 = x_p2 + (x_p1 - x_p2)/6;
y_q1 = y_p2 + (y_p1 - y_p2)/6;
x_q2 = x_p1 - (x_p1 - x_p2)/3;
y_q2 = y_p1 - (y_p1 - y_p2)/3;

circle_radius = R;
arc_radius1 = R/3; % 小圆弧半径
arc_radius2 = 2 * arc_radius1; % 大圆弧半径
arc_center1 = [x_q1, y_q1]; % 小圆弧圆心
arc_center2 = [x_q2, y_q2]; % 大圆弧圆心

p = R ;
theta_p = p/b

x_bei = theta_p;
y_chu = 2 * pi;
theta_pi = mod(x_bei, y_chu)
type = tan(theta_pi)

% 定义圆弧的角度范围
theta_arc1 = linspace(pi + theta_pi + pi, 2 * pi + theta_pi+ pi, 100);
theta_arc2 = linspace(0 + theta_pi+ pi, pi + theta_pi+ pi, 100);

% 计算圆弧的坐标
x_arc1 = arc_center1(1) + arc_radius1 * cos(theta_arc1);
y_arc1 = arc_center1(2) + arc_radius1 * sin(theta_arc1);
x_arc2 = arc_center2(1) + arc_radius2 * cos(theta_arc2);
y_arc2 = arc_center2(2) + arc_radius2 * sin(theta_arc2);

%theta_arc3 = linspace(0 , 2 * pi, 100);
%theta_arc4 = linspace(0 , 2 * pi, 100);

% 计算圆弧的坐标
%x_arc1 = arc_center1(1) + arc_radius1 * cos(theta_arc3);
%y_arc1 = arc_center1(2) + arc_radius1 * sin(theta_arc3);
%x_arc2 = arc_center2(1) + arc_radius2 * cos(theta_arc4);
%y_arc2 = arc_center2(2) + arc_radius2 * sin(theta_arc4);




%%
% 每个时间点

for iiii = 1:88




t = 14 - 13.5827 + iiii;
    theta = 0;
    theta_min = 1000;
    for i = 0:4000000

        theta = theta + 0.0001;
        S = b*(1/2 * theta * sqrt(1 + theta^2) + 1/2 * log(abs(theta + sqrt(1 + theta^2))));
        if abs(S - 35.080245087274040 - t) < 0.001
            if abs(S - 35.080245087274040 - t) < theta_min
                theta_min = theta;
            end
        end
    end
    display(theta_min)










x_3 = (a - b * theta_min) .* cos(theta_min)
y_3 = (a - b * theta_min) .* sin(theta_min)
x_4(1) = x_3;
y_4(1) = y_3;
table111(1,iiii) = x_3;
table111(2,iiii) = y_3;
com = com + 1;

    %%
    %判定龙头位置
    
    theta_00 =  theta_min;
    
    for n = 1:1
    
        theta_now = theta_00;
        theta_22 = theta_00;
        theta_11 = theta_00;
        theta_aa = 1000;
        
        for j = 0:50000
        
            theta_11 = theta_11 - 0.001;
        
            for i = 0:800
        
                theta_22 = theta_22 - 0.01;
            
                aa = 0;  
                bb = b; 
                
                r_11 = aa - bb * theta_11;
                x_11 = r_11 .* cos(theta_11);
                y_11 = r_11 .* sin(theta_11);
                
                Rnn = 2.86; 
                
                x_22 = x_3 + Rnn * cos(theta_22);
                y_22 = y_3 + Rnn * sin(theta_22);

                if abs(x_11 - x_22) < 0.01 && abs(y_11 - y_22) < 0.01
                    q_tge = theta_11;
                    if abs(theta_11 - theta_now) < theta_aa && (theta_11 - theta_now) < 0
                        theta_aa = abs(theta_11 - theta_now);
                        theta_00 = theta_11;
                    end
                end
            
            end
        end
    
        display(theta_00)
    
        x_3 = (a - b * theta_00) .* cos(theta_00);
        y_3 = (a - b * theta_00) .* sin(theta_00);
    end


    zhi_tou = x_3^2 + y_3^2
    if zhi_tou > R^2
    
        x_4(2) = x_3;
        y_4(2) = y_3;
        x_4h(1) = x_3;
        y_4h(1) = y_3;


    else


        %%
        %判定龙头首部位置
        theta_tou1 =  theta_min;
        % 阿基米德螺线方程
        r_tou1 = a - b * theta_tou1;
        % 将螺线的极坐标转换为直角坐标
        x_tou1 =  r_tou1 .* cos(theta_tou1);
        y_tou1 =  r_tou1 .* sin(theta_tou1);

        table
        
        %x_tou1 =  -0.878478;
        %y_tou1 =  3.48184;
        
        
        %判定龙头后部位置
        
        %圆弧小
        x1x = x_q1; y1x = y_q1; r1x = arc_radius1;
        %圆弧大
        x1d = x_q2; y1d = y_q2; r1d = arc_radius2;
        
        %相交圆
        x2 = x_tou1; y2 = y_tou1; r2 = 2.86;
        
        
        
        intersect_points1 = circle_intersection(x1x, y1x, r1x, x2, y2, r2);
        intersect_points2 = circle_intersection(x1d, y1d, r1d, x2, y2, r2);
        
        %判断小圆弧
        if isempty(intersect_points1)
            disp('两个圆不相交');
        else
            fprintf('交点坐标: (%.2f, %.2f) 和 (%.2f, %.2f)\n', intersect_points1(1, 1), intersect_points1(1, 2), intersect_points1(2, 1), intersect_points1(2, 2));
        x_tou2x = intersect_points1(1, 1);
        y_tou2x = intersect_points1(1, 2);
        x_3 = x_tou2x;
        y_3 = y_tou2x;
        end
        
        
        
        %判断大圆弧
        if isempty(intersect_points2)
            disp('两个圆不相交');
        else
            fprintf('交点坐标: (%.2f, %.2f) 和 (%.2f, %.2f)\n', intersect_points2(1, 1), intersect_points2(1, 2), intersect_points2(2, 1), intersect_points2(2, 2));
        x_tou2d = intersect_points2(2, 1);
        y_tou2d = intersect_points2(2, 2);
        x_3 = x_tou2d;
        y_3n = y_tou2d;
        end

        x_4(2) = x_3
        y_4(2) = y_3
        x_4h(1) = x_3
        y_4h(1) = y_3
    
    end

    table111(3,iiii) = x_3;
    table111(4,iiii) = y_3;
    com = com + 1;
    
    %%
    %判定龙身位置
    time1 = 0
    for n = 2:20
    
        theta_now = theta_00;
        theta_22 = theta_00;
        theta_11 = theta_00;
        theta_aa = 1000;
        
        for j = 0:5000
        
            theta_11 = theta_11 - 0.001;
        
            for i = 0:700
        
                theta_22 = theta_22 - 0.01;
            
                aa = 0;  
                bb = b; 
                
                r_11 = aa - bb * theta_11;
                x_11 = r_11 .* cos(theta_11);
                y_11 = r_11 .* sin(theta_11);
                
                Rkk = 1.65; 
                
                x_22 = x_3 + Rkk * cos(theta_22);
                y_22 = y_3 + Rkk * sin(theta_22);

                if abs(x_11 - x_22) < 0.05 && abs(y_11 - y_22) < 0.05
                    if abs(theta_11 - theta_now) < theta_aa && (theta_11 - theta_now) < 0
                        theta_aa = abs(theta_11 - theta_now);
                        theta_00 = theta_11;
                    end
                end
            
            end
        end
    
        display(theta_00)
    
        x_3 = (a - b * theta_00) .* cos(theta_00)
        y_3 = (a - b * theta_00) .* sin(theta_00)

        R2 = R^2
        zhi = (x_3^2 + y_3^2)

        x_4h(n) = x_3;
        y_4h(n) = y_3;

        if zhi < R^2
            break;
        end


        x_4(n+1) = x_3;
        y_4(n+1) = y_3;
        table111(2 * com - 1,iiii) = x_3;
        table111(2 * com,iiii) = y_3;
        com = com + 1;


    end
        x_4(n) = 0;
        y_4(n) = 0;
        x_4h(n) = 0;
        y_4h(n) = 0;


%%
%判定龙身位置1

x_shen = x_4h(n - 1);
y_shen = y_4h(n - 1);

k = zeros(1, 20);

x3_s1 = zeros(1, 222);
y3_s1 = zeros(1, 222);
x3_s2 = zeros(1, 223);
y3_s2 = zeros(1, 223);

x3_s1(1) = x_shen;
y3_s1(1) = y_shen;
x3_s = x_shen; 
y3_s = y_shen;

for i = 0:10
        
    
    
    %圆弧小
    x1x = x_q1; y1x = y_q1; r1x = arc_radius1;
    %圆弧大
    x1d = x_q2; y1d = y_q2; r1d = arc_radius2;
    
    %相交圆(不断更改)
    x3 = x3_s; y3 = y3_s; r3 = 1.65;
    
    
    
    intersect_points3 = circle_intersection(x1x, y1x, r1x, x3, y3, r3);
    intersect_points4 = circle_intersection(x1d, y1d, r1d, x3, y3, r3);
    
    %判断小圆弧
    if isempty(intersect_points3)
        %disp('两个圆不相交');
    else
        %fprintf('交点坐标: (%.2f, %.2f) 和 (%.2f, %.2f)\n', intersect_points3(1, 1), intersect_points3(1, 2), intersect_points3(2, 1), intersect_points3(2, 2));
        x_shen2x = intersect_points3(1, 1);
        y_shen2x = intersect_points3(1, 2);
        x_shen1 = x_shen2x;
        y_shen1 = y_shen2x;
    
    end
    
    %判断大圆弧
    if isempty(intersect_points4)
        %disp('两个圆不相交');
    else
        %fprintf('交点坐标: (%.2f, %.2f) 和 (%.2f, %.2f)\n', intersect_points4(1, 1), intersect_points4(1, 2), intersect_points4(2, 1), intersect_points4(2, 2));
        x_shen2d = intersect_points4(2, 1);
        y_shen2d = intersect_points4(2, 2);
        x_shen1 = x_shen2d;
        y_shen1 = y_shen2d;
    end


    %if (i+1) < 10
        x3_s1(i+2) = x_shen1;
        y3_s1(i+2) = y_shen1;
    %end
    x3_s2(i+1) = x_shen1;
    y3_s2(i+1) = y_shen1;

    x3_s = x_shen1;
    y3_s = y_shen1;

    table111(2 * com - 1,iiii) = x_shen1;
    table111(2 * com,iiii) = y_shen1;
    com = com + 1;

    lizi = i;

    q_shen = (y_shen1 - y_p1)^2 + (x_shen1 - x_p1)^2
    1.65^2
        if ((y_shen1 - y_p1)^2 + (x_shen1 - x_p1)^2) < 1.65^2
           (y_shen1 - y_p1)^2 + (x_shen1 - x_p1)^2
        break;
        end


end
    
    
    
















%%
%判定龙身位置2
x_cen = x3_s;
y_cen = y3_s;

for n = 1:1

    theta_now = 0;
    theta_22 = 0;
    theta_11 = 0;
    l_min = 1000;
    theta_aa = 1000;
    
    for j = 0:3000
    
        theta_11 = theta_11 + 0.01;
    
        for i = 0:700
    
            theta_22 = theta_22 + 0.01;
            
            %阿基米德
            aa = 0;  
            bb = b; 
            
            r_11 = aa + bb * theta_11;
            x_11 = r_11 .* cos(theta_11);
            y_11 = r_11 .* sin(theta_11);
            
            %圆
            R0 = 1.65; 
            
            x_22 = x_cen + R0 * cos(theta_22);
            y_22 = y_cen + R0 * sin(theta_22);
        
            if abs(x_11 - x_22) < 0.1 && abs(y_11 - y_22) < 0.1
                if (x_11^2 + y_11^2) > R^2
                    if (x_11^2 + y_11^2) < l_min 
                    theta_aa = abs(theta_11 - theta_now);
                    theta_00 = theta_11;
                    l_min = x_11^2 + y_11^2;
                    end

                end
            end
        
        end
    end

    display(theta_00)

    x_3 = (a + b * theta_00) .* cos(theta_00);
    y_3 = (a + b * theta_00) .* sin(theta_00);
    %table(3,t+1) = x_3;
    %table(4,t+1) = y_3;
end

    i = lizi;
    x3_s1(i+2) = x_cen;
    y3_s1(i+2) = y_cen;

    x3_s2(i+2) = x_3;
    y3_s2(i+2) = y_3;
    lizi_1 = i+2


    x3_s1(lizi_1 + 1) = x_3;
    y3_s1(lizi_1 + 1) = y_3;

    table111(2 * com - 1,iiii) = x_3;
    table111(2 * com,iiii) = y_3;
    com = com + 1;


    
%%
%判定龙身位置3

c = 0;
for n = 1:222

    theta_now = theta_00;
    theta_22 = theta_00;
    theta_11 = theta_00;
    theta_aa = 1000;
    
    for j = 0:3000
    
        theta_11 = theta_11 + 0.001;
    
        for i = 0:700
    
            theta_22 = theta_22 + 0.01;
        
            aa = 0;  
            bb = b; 
            
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
    
    

    %table(2*n+1,t+1) = x_3;
    %table(2*n+2,t+1) = y_3;
    if (n+1) < 223
        x3_s1(lizi_1 + n + 1) = x_3;
        y3_s1(lizi_1 + n + 1) = y_3;
    end
        x3_s2(lizi_1 + n) = x_3;
        y3_s2(lizi_1 + n) = y_3;

    table111(2 * com - 1,iiii) = x_3;
    table111(2 * com,iiii) = y_3;
    com = com + 1;

end






end






%%
figure; 


%%
%绘制龙


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



    % 画出靠近点15cm的线段
    line([x1_close, x2_close], [y1_close, y2_close], 'Color', 'b');

    % 画出远离点15cm的线段
    line([x1_far, x2_far], [y1_far, y2_far], 'Color', 'g');

    % 画出矩形
    % 连接靠近点15cm的线段端点到远离点15cm的线段端点
    line([x1_close, x1_far], [y1_close, y1_far], 'Color', 'k', 'LineWidth', 1);
    line([x1_far, x2_far], [y1_far, y2_far], 'Color', 'k', 'LineWidth', 1);
    line([x2_far, x2_close], [y2_far, y2_close], 'Color', 'k', 'LineWidth', 1);
    line([x2_close, x1_close], [y2_close, y1_close], 'Color', 'k', 'LineWidth', 1);
end
    
    
%%
%画龙身矩形1


% 假设交点坐标为 (x_inter, y_inter)
x_inter = 0;
y_inter = 0;

% 假设线段端点坐标存储在两个数组中，x1, y1为线段起点，x2, y2为线段终点
x11 = x3_s1; % 线段起点x坐标
y11 = y3_s1; % 线段起点y坐标

%判断小圆弧
x21 = x3_s2;% 线段终点x坐标
y21 = y3_s2;% 线段终点y坐标




% 线段数量
n_segments = length(x11);

% 创建新图像

hold on; % 保持图像，以便在同一幅图上画多条线

for i = 1:n_segments
    % 计算线段的方向向量
    dx = x21(i) - x11(i);
    dy = y21(i) - y11(i);

    % 计算线段的长度
    length_segment = sqrt(dx^2 + dy^2);

    % 计算单位方向向量
    dx_unit = dx / length_segment;
    dy_unit = dy / length_segment;

    % 增长线段，每边增长27.5cm
    x1_extended = x11(i) - 0.275 * dx_unit;
    y1_extended = y11(i) - 0.275 * dy_unit;
    x2_extended = x21(i) + 0.275 * dx_unit;
    y2_extended = y21(i) + 0.275 * dy_unit;

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

    hold on;

    % 画出靠近点15cm的线段
    line([x1_close, x2_close], [y1_close, y2_close], 'Color', 'b');

    % 画出远离点15cm的线段
    line([x1_far, x2_far], [y1_far, y2_far], 'Color', 'g');

    % 画出矩形y
    % 连接靠近点15cm的线段端点到远离点15cm的线段端点
    line([x1_close, x1_far], [y1_close, y1_far], 'Color', 'k', 'LineWidth', 1);
    line([x1_far, x2_far], [y1_far, y2_far], 'Color', 'k', 'LineWidth', 1);
    line([x2_far, x2_close], [y2_far, y2_close], 'Color', 'k', 'LineWidth', 1);
    line([x2_close, x1_close], [y2_close, y1_close], 'Color', 'k', 'LineWidth', 1);
end
        









%%


% 绘制阿基米德螺线
plot(x_1l, y_1l,'g', x_1r, y_1r, 'b' , 'LineWidth', 0.2);  % 螺线图像
hold on;

% 绘制点
%plot(x_cen, y_cen,  'o' , 'LineWidth', 1);
%plot(x_shen, y_shen,  'o' , 'LineWidth', 1);
%plot(x_shen1, y_shen1,  'o' , 'LineWidth', 1);

%plot(x_q1 , y_q1,  'o' , 'LineWidth', 0.2);
%plot(x_q2 , y_q2,  'o' , 'LineWidth', 0.2);

% 绘制偏移后的圆
plot(x_2, y_2, 'r', 'LineWidth', 0.2);  % 偏移圆的图像

% 绘制圆弧
plot(x_arc1, y_arc1, 'k-', 'LineWidth', 0.2);
plot(x_arc2, y_arc2, 'k-', 'LineWidth', 0.2);

% 调整图像比例
axis equal;
title('Intersection of Archimedean Spiral and Offset Circle');
xlabel('x');
ylabel('y');
grid on;

% 显示图例
legend('Archimedean Spiral', 'Offset Circle');

hold off;















%%
function [intersect_points] = circle_intersection(x1, y1, r1, x2, y2, r2)
    % 计算两圆心之间的距离
    d = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    
    % 判断圆是否相交
    if d > r1 + r2 || d < abs(r1 - r2)
        intersect_points = []; % 返回空数组
        return;
    end
    
    % 计算交点距离
    a = (r1^2 - r2^2 + d^2) / (2 * d);
    h = sqrt(r1^2 - a^2);
    
    % 计算交点基准点坐标
    P2_x = x1 + a * (x2 - x1) / d;
    P2_y = y1 + a * (y2 - y1) / d;
    
    % 计算两个交点坐标
    x3_1 = P2_x + h * (y2 - y1) / d;
    y3_1 = P2_y - h * (x2 - x1) / d;
    
    x3_2 = P2_x - h * (y2 - y1) / d;
    y3_2 = P2_y + h * (x2 - x1) / d;
    
    intersect_points = [x3_1, y3_1; x3_2, y3_2];
end









