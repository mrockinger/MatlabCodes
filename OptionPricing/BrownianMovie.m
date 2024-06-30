% If you increase the number of points in the Brownian motion
% example to something like n = 300 and s = .02, the motion is
% no longer very fluid; it takes too much time to draw each time
% step. It becomes more effective to save a predetermined number
% of frames as bitmaps and to play them back as a movie.

% First, decide on the number of frames, say

nframes = 50;
n=50;
s=0.2;

% Next, set up the first plot as before, except do not use EraseMode.

x = rand(n,1)-0.5;
y = rand(n,1)-0.5;
h = plot(x,y,'.')
set(h,'MarkerSize',18);
axis([-1 1 -1 1])
axis square
grid off

% Now, allocate enough memory to save the complete movie,

M = moviein(nframes);

% This sets aside a large matrix with nframes columns. Each
% column is long enough to save one frame. The total amount of
% memory required is proportional to the number of frames, and
% to the area of the current axes; it is independent of the
% complexity of the particular plot. For 50 frames and the
% default axes, over 7.5 megabytes of memory is required. This
% example is using square axes, which is slightly smaller, so
% only about 6 megabytes is required.

% Generate the movie and use getframe to capture each frame.

for k = 1:nframes
   x = x + s*randn(n,1);
   y = y + s*randn(n,1);
   set(h,'XData',x,'YData',y)
   M(:,k) = getframe;
end

% Finally, play the movie 30 times.

movie(M,30)


