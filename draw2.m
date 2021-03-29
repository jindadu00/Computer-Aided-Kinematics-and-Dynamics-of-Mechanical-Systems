function draw2

R = 5;
x0 = -10;
y0 = 0;

pause_time = 0.1;
dt = 0.1;
axis_range =  [-20 30 -20 30];
t=0:dt:10;

x1 = x0+3*t;
y1 = y0+0*t+R;

x2 = x1+R*cos(x1/R);
y2 = y1-R*sin(x1/R);




fig = plot(x1,y1-R); % 初始化图形 画地面
set(fig,'LineWidth',3,'Color','b');

line_handle = line('LineStyle','-', ...
                 'LineWidth',3, ...
                 'Color','c');
                         
circle_handle = rectangle('Curvature',[1,1], ...
                          'EdgeColor','g', ...
                          'LineWidth',3);

i=1;
n=length(t);

while ishandle(fig) % 绘画循环
 
    set(circle_handle,'Position',[x1(i)-R,y1(i)-R,2*R,2*R]);
    
    set(line_handle,'xdata',[x1(i) x2(i)], ...
        'ydata',[y1(i) y2(i)]);
    
    axis(axis_range);
    axis square
    
    drawnow
    
    i=i+1;
    if i>n
        i=1;
    end
    
    pause(pause_time);
end








             
             

