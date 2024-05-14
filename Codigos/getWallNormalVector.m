function N = getWallNormalVector(surfaceofimpact)
% GETWALLNORMALVECTOR Get the normal vector of a surface
switch surfaceofimpact
    case 1
        N = [1 0 0];
    case 2
        N = [-1 0 0];
    case 3
        N = [0 1 0];
    case 4
        N = [0 -1 0];
    case 5
        N = [0 0 1];
    case 6
        N = [0 0 -1];
end

end