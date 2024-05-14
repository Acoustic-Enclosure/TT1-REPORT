function [] = PlotSurfaces(RoomSurface)
    figure
    hold on
    for i = 1:length(RoomSurface)       
        if strcmp(RoomSurface(i).Orientation,'XY')
            [YSize,XSize] = size(RoomSurface(i).AbsMap);
            [X,Y] = meshgrid(1:XSize,1:YSize);
            surf(X,Y,100*RoomSurface(i).FirstVertexPos(3)*ones(YSize,XSize),RoomSurface(i).AbsMap, ...
                "EdgeColor","none","FaceAlpha",0.3)
        elseif strcmp(RoomSurface(i).Orientation,'XZ')
            [YSize,XSize] = size(RoomSurface(i).AbsMap);
            [X,Y] = meshgrid(1:XSize,1:YSize);
            surf(X,100*RoomSurface(i).FirstVertexPos(2)*ones(YSize,XSize),Y,RoomSurface(i).AbsMap, ...
                "EdgeColor","none","FaceAlpha",0.3)
        elseif strcmp(RoomSurface(i).Orientation,'YZ')
            [YSize,XSize] = size(RoomSurface(i).AbsMap);
            [X,Y] = meshgrid(1:XSize,1:YSize);
            surf(100*RoomSurface(i).FirstVertexPos(1)*ones(YSize,XSize),X,Y,RoomSurface(i).AbsMap, ...
                "EdgeColor","none","FaceAlpha",0.3)
        else
            warning('Not a valid orientation')          
        end
        colorbar
    end
    view(-30,30)
    c = colorbar;
    c.Label.String = 'Absorption Coefficient (0 to 1)';
    xlabel('EJE X'); ylabel('EJE Y'); zlabel('EJE Z');
    hold off
end