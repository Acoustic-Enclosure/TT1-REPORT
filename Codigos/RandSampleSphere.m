
function X=RandSampleSphere(N)
% RANDSAMPLESPHERE Return random ray directions

% Sample the unfolded right cylinder
z = 2*rand(N,1)-1
lon = 2*pi*rand(N,1);

% Convert z to latitude
z(z<-1) = -1;
z(z>1) = 1;
lat = acos(z);

% Convert spherical to rectangular co-ords
s = sin(lat);
x = cos(lon).*s;
y = sin(lon).*s;

X = [x y z];
end