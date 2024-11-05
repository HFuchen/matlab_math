x1 = 1.93291; y1 = 1.96367;
x2 = 2.68316; y2 = 0.249729;

x3 = 2.30894; y3 = 1.02282;


% 调用函数
distance = point_to_line(x1, y1, x2, y2, x3, y3);
fprintf('The shortest distance from (%d, %d) to the line is: %f\n', x3, y3, distance);


function distance = point_to_line(x1, y1, x2, y2, x3, y3)
    % 计算直线的系数A, B, C
    A = y2 - y1;
    B = x1 - x2;
    C = x2 * y1 - x1 * y2;

    % 计算点到直线的距离
    distance = abs(A * x3 + B * y3 + C) / sqrt(A^2 + B^2);
end