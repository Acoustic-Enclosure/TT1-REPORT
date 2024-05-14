classdef Panel
    properties
        Orientation
        FirstVertexPos
        SecondVertexPos
        Surface
        BackgroundAbsCoeff
        Absorption
        AbsMap
    end

    methods
        function obj = Panel(Orientation,FirstVertexPosition,SecondVertexPosition,AbsorptionCoefficient)
            % Constructor of Panel
            obj.BackgroundAbsCoeff = AbsorptionCoefficient; 
            obj.FirstVertexPos = FirstVertexPosition;
            obj.SecondVertexPos = SecondVertexPosition;
            if strcmp('XY',Orientation)
                SizeX = round(abs(SecondVertexPosition(2)-FirstVertexPosition(2))*100,0);
                SizeY = round(abs(SecondVertexPosition(1)-FirstVertexPosition(1))*100,0);
                tempMap = zeros(SizeX,SizeY);
                obj.AbsMap = AbsorptionCoefficient+tempMap;
                FirstVertexPosition(3) = [];
                SecondVertexPosition(3) = [];
            elseif strcmp('YZ',Orientation)       
                SizeX = round(abs(SecondVertexPosition(3)-FirstVertexPosition(3))*100,0);
                SizeY = round(abs(SecondVertexPosition(2)-FirstVertexPosition(2))*100,0);
                tempMap = zeros(SizeX,SizeY);
                obj.AbsMap = AbsorptionCoefficient+tempMap;
                FirstVertexPosition(1) = [];
                SecondVertexPosition(1) = [];
            elseif strcmp('XZ',Orientation) 
                SizeX = round(abs(SecondVertexPosition(3)-FirstVertexPosition(3))*100,0);
                SizeY = round(abs(SecondVertexPosition(1)-FirstVertexPosition(1))*100,0);
                tempMap = zeros(SizeX,SizeY);
                obj.AbsMap = AbsorptionCoefficient+tempMap;
                FirstVertexPosition(2) = [];
                SecondVertexPosition(2) = [];
            else
                warning('Not a valid Orientation');
            end
            obj.Surface = abs(prod(FirstVertexPosition - SecondVertexPosition));
            obj.Orientation = Orientation;
            obj.Absorption = AbsorptionCoefficient*obj.Surface;  
            
        end

        function obj = addPanel_onTop(obj,Idx,FirstVertexPosition,SecondVertexPosition,AbsorptionCoefficient)          
            S = prod(FirstVertexPosition-SecondVertexPosition);
            obj(Idx).Absorption = obj(Idx).Absorption + S*(AbsorptionCoefficient - obj(Idx).BackgroundAbsCoeff);
            obj(Idx).AbsMap(round(1+(FirstVertexPosition(1)*100),0):round(SecondVertexPosition(1)*100,0), ...
                round(1+(FirstVertexPosition(2)*100),0):round(SecondVertexPosition(2)*100,0)) = AbsorptionCoefficient;
        end
    end
end