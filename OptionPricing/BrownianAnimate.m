% Specify a number of points, like
n = 20

% and a temperature or velocity, such as

s = .02

% The best values for these two parameters depend upon 
% the speed of your particular computer. Generate n random 
% points with (x,y) coordinates between -1/2 and +1/2

x = rand(n,1)-0.5;
y = rand(n,1)-0.5;

% Plot the points in a square with sides at -1 and +1. 
% Save the handle for the vector of points and set its 
% EraseMode to xor. This tells the MATLAB graphics system
% not to redraw the entire plot when the coordinates of one
% point are changed, but to restore the background color in
% the vicinity of the point using an "exclusive or" operation. 

h = plot(x,y,'.');
axis([-1 1 -1 1])
axis square
grid off
set(h,'EraseMode','xor','MarkerSize',18)

% Now begin the animation. Here is an infinite while loop, which you 
% will eventually break out of by typing <ctrl>-c. Each time through 
% the loop, add a small amount of normally distributed random noise to 
% the coordinates of the points. Then, instead of creating an entirely 
% new plot, simply change the XData and YData properties of the original plot. 

while 1
   drawnow
   x = x + s*randn(n,1);
   y = y + s*randn(n,1);
   set(h,'XData',x,'YData',y)
end


