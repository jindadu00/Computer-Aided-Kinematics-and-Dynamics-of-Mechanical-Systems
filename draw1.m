function draw1

l1 = 10;
l2 = 6; % 杆长

lrect = 8;
wrect = 4;% 矩形长宽

xtri = 0;
ytri = -1.5;
ltri = 2.5;
thetatri = 0;% 三角形边长/中心坐标/转角

pause_time = 0.1;
dt = 0.01;
axis_range =  [-20 20 -20 20];
t=0:dt:10;

theta1 = 4*t;
theta2 = cos(t);
theta3 = 10*t;

x1 = 0*t;
y1 = 0*t;

x2 = l1*cos(theta1);
y2 = l1*sin(theta1);

x3 = x2+l2*cos(theta2);
y3 = y2+l2*sin(theta2);


fig = plot(x1,y1); % 初始化图形

%设置杆件、矩形、三角形的线性和颜色等参数：
line_handle = line('LineStyle','-', ...
                 'Marker','o', ...
                 'LineWidth',6, ...
                 'Color','b', ...
                 'MarkerSize',5, ...
                 'MarkerEdgeColor','c');
             
rec_handle = line('LineStyle','-', ...
                 'LineWidth',2, ...
                 'Color','g');
             
tri_handle = line('LineStyle','-', ...
                 'LineWidth',2, ...
                 'Color','g');

%绘制三角形：             
[X Y] = drawtri(xtri,ytri,ltri,thetatri);
set(tri_handle,'xdata',X, ...
    'ydata',Y); 
% fill(X,Y,'r');
i=1;
n=length(t);

while ishandle(fig) % 绘画循环

    set(line_handle,'xdata',[x1(i) x2(i) x3(i)], ...
        'ydata',[y1(i) y2(i) y3(i)]);
    
    [X,Y] = drawrect(x3(i),y3(i),lrect,wrect,theta3(i));        
    set(rec_handle,'xdata',X, ...
        'ydata',Y);
     
    axis(axis_range);
    axis square
    
    drawnow
    
    i=i+1;
    if i>n
        i=1;
    end
    
    pause(pause_time);
end

function [X Y]=drawrect(x0,y0,l,w,theta)
% 求解矩形四角坐标
%x0 y0 中心坐标
%l w 长宽
%theta 转角
T=[cos(theta) sin(theta);-sin(theta) cos(theta)];
coord = [x0 y0]+[l/2 w/2]*T;
x1=coord(1);y1=coord(2);
coord = [x0 y0]+[l/2 -w/2]*T;
x2=coord(1);y2=coord(2);
coord = [x0 y0]+[-l/2 -w/2]*T;
x3=coord(1);y3=coord(2);
coord = [x0 y0]+[-l/2 +w/2]*T;
x4=coord(1);y4=coord(2);
X = [x1 x2 x3 x4 x1];
Y = [y1 y2 y3 y4 y1];

function [X Y]=drawtri(x0,y0,l,theta)
% 求解三角形三点坐标
%x0 y0 中心坐标
%l  边长长
%theta 转角
l = l/2/sin(pi/3);
T=[cos(theta+pi/2) sin(theta+pi/2);-sin(theta+pi/2) cos(theta+pi/2)];
coord = [x0 y0]+[l 0]*T;
x1=coord(1);y1=coord(2);
T=[cos(theta+7*pi/6) sin(theta+7*pi/6);-sin(theta+7*pi/6) cos(theta+7*pi/6)];
coord = [x0 y0]+[l 0]*T;
x2=coord(1);y2=coord(2);
T=[cos(theta+11*pi/6) sin(theta+11*pi/6);-sin(theta+7*pi/6) cos(theta+7*pi/6)];
coord = [x0 y0]+[l 0]*T;
x3=coord(1);y3=coord(2);

X = [x1 x2 x3 x1];
Y = [y1 y2 y3 y1];






             
             

