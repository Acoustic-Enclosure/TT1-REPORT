function Absorption = GetAbsorption(RoomSurface)
    Absorption = 0;
    for i = 1:length(RoomSurface)    
        Absorption = Absorption + RoomSurface(i).Absorption;
    end
end