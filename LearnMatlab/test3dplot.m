function test3dplot()
Z = peaks(20)
figure(1)
subplot(2,1,2)
h = surf(Z)
return
colormap hot
shading interp
set(h,'EdgeColor','k')
light('Position',[-2,2,20])
lighting 
phongmaterial([0.4,0.6,0.5,30])
set(h,'FaceColor',[0.7 0.7 0],...      
    'BackFaceLighting','lit')
view([30,25])
set(gca,'CameraViewAngleMode','Manual')
axis([5 15 5 15 -8 8])
set(gca,'ZTickLabel','Negative||Positive')
set(gca,'PlotBoxAspectRatio',[2.5 2.5 1])
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Function Value')
title('Peaks')
