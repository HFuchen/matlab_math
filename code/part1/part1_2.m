clc;
clear;
a = 0;  % 阿基米德螺线初始半径
b = 0.55/(2*pi);  % 螺线间距
for t = 0:300
    theta = 0;
    theta_min = 1000;
    for i = 0:400000
        theta = theta + 0.001;
        S = 0.0876*(1/2 * theta * sqrt(1 + theta^2) + 1/2 * log(abs(theta + sqrt(1 + theta^2))));
        if abs(S - 442.91 + t) < 0.01
            if abs(S - 442.91 + t) < theta_min
                theta_min = theta;
            end
        end
    end
    display(theta_min)
    x_3 = (0.55/(2*pi) * theta_min) .* cos(theta_min)
    y_3 = (0.55/(2*pi) * theta_min) .* sin(theta_min)
end



