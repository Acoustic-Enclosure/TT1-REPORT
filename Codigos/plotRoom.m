function plotRoom(roomDimensions,receiverCoord,sourceCoord,figHandle)
% PLOTROOM Helper function to plot 3D room with receiver/transmitter points
figure(figHandle)
X = [0;roomDimensions(1);roomDimensions(1);0;0];
Y = [0;0;roomDimensions(2);roomDimensions(2);0];
Z = [0;0;0;0;0];
figure;
hold on;
plot3(X,Y,Z,"k",LineWidth=1.5);   
plot3(X,Y,Z+roomDimensions(3),"k",LineWidth=1.5); 
set(gca,"View",[-28,35]); 
for k=1:length(X)-1
    plot3([X(k);X(k)],[Y(k);Y(k)],[0;roomDimensions(3)],"k",LineWidth=1.5);
end
grid on
xlabel("X (m)")
ylabel("Y (m)")
zlabel("Z (m)")
plot3(sourceCoord(1),sourceCoord(2),sourceCoord(3),"bx",LineWidth=2)
plot3(receiverCoord(1),receiverCoord(2),receiverCoord(3),"ro",LineWidth=2)
end