clear; close all; clc;

%% Config
roomDim = [6 4 3.4]; %[X,Y,Z]
V = prod(roomDim);
ConcreteAbsCoeff = 0.02;
PanelAbs = 0.8;
DoorAbs = 0.30;
PanelDim = [0.3 1.33];

%% Create room
RoomSurface(1) = Panel('XY',[0,0,0],[roomDim(1) roomDim(2),0], ConcreteAbsCoeff); % Floor
RoomSurface(2) = Panel('XY',[0,0,roomDim(3)],[roomDim(1) roomDim(2) roomDim(3)], ConcreteAbsCoeff); %Ceiling
RoomSurface(3) = Panel('XZ',[roomDim(1),0,0],[0 0 roomDim(3)], 0.02); %First Wall
RoomSurface(4) = Panel('YZ',[roomDim(1) roomDim(2) 0],[roomDim(1),0,roomDim(3)], ConcreteAbsCoeff); %Second Wall
RoomSurface(5) = Panel('XZ',[0,roomDim(2),0],[roomDim(1) roomDim(2) roomDim(3)], ConcreteAbsCoeff); %Third Wall
RoomSurface(6) = Panel('YZ',[0 0 0],[0,roomDim(2),roomDim(3)], ConcreteAbsCoeff); %Fourth Wall
RoomSurface = RoomSurface.addPanel_onTop(4,[0,1.5],[2.2,2.5],DoorAbs); %Door panel

AbsorptionBefore = GetAbsorption(RoomSurface);
    T60Before = 0.161*(V/AbsorptionBefore)

%% Calculate surface required
RequiredT60 = 0.5

ReqAbs = 0.161*(V/RequiredT60);
AbsError = ReqAbs - AbsorptionBefore;
AbsPerPanel = (prod(PanelDim)*(PanelAbs-ConcreteAbsCoeff));
PanelsRequired = round(AbsError/AbsPerPanel,0)


%% Add panels to surface
PosXZ = 0.5;
PosYZ = 0.5;
PosXY = 0.5;
while PanelsRequired > 0
    if PosXZ < roomDim(1)-0.5-PanelDim(1)
        RoomSurface = RoomSurface.addPanel_onTop(3,[1,PosXZ],[1+PanelDim(2),PosXZ+PanelDim(1)],PanelAbs); %Acoustic panel
        RoomSurface = RoomSurface.addPanel_onTop(5,[1,PosXZ],[1+PanelDim(2),PosXZ+PanelDim(1)],PanelAbs); %Acoustic panel
        PanelsRequired = PanelsRequired-2;
        PosXZ = PosXZ+PanelDim(1)+0.05;
    else
        if PosYZ < roomDim(2)-0.5-PanelDim(1)
            RoomSurface = RoomSurface.addPanel_onTop(6,[1,PosYZ],[1+PanelDim(2),PosYZ+PanelDim(1)],PanelAbs); %Acoustic panel
            PosYZ = PosYZ+PanelDim(1)+0.05;
            PanelsRequired = PanelsRequired-1;
        else
            if PosXY < roomDim(1)-0.5-PanelDim(1)
                RoomSurface = RoomSurface.addPanel_onTop(2,[1.33,PosXY],[1+PanelDim(2),PosXY+PanelDim(1)],PanelAbs); %Acoustic panel
                PosXY = PosXY+PanelDim(1)+0.05;
                PanelsRequired = PanelsRequired-1;
            else
                warning('Not posible with current panel absorption, no more space to alocate panels')
                fprintf('%i panels not alocated',PanelsRequired)
                PanelsRequired = 0;
            end
        end
    end    
end

%% Calculate New Absorption
AbsorptionAfter = GetAbsorption(RoomSurface);
T60After = 0.161*(V/AbsorptionAfter)
PlotSurfaces(RoomSurface)