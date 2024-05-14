clear; close all; clc;

%% SetUp
SetUpStruct.room = [10 8 4];
SetUpStruct.src_pos = [2 2 2];
SetUpStruct.mic_pos = [5 5 1.8];
SetUpStruct.mic_radius = 0.0875;
impResTime = 2;

plotRoom(SetUpStruct.room,SetUpStruct.mic_pos,SetUpStruct.src_pos,1)

%% Generate Rays
N = 5000;
rng(0)
rays = RandSampleSphere(N);

%% Reflections and Scattering Coefficients
FVect = [125 250 500 1000 2000 4000];

A = [0.08 0.09 0.12 0.16 0.22 0.24;
    0.08 0.09 0.12 0.16 0.22 0.24;
    0.08 0.09 0.12 0.16 0.22 0.24;
    0.08 0.09 0.12 0.16 0.22 0.24;
    0.08 0.09 0.12 0.16 0.22 0.24;
    0.08 0.09 0.12 0.16 0.22 0.24].';

R = sqrt(1-A);

D = [0.05 0.3 0.7 0.9 0.92 0.94;
     0.05 0.3 0.7 0.9 0.92 0.94; 
     0.05 0.3 0.7 0.9 0.92 0.94; 
     0.05 0.3 0.7 0.9 0.92 0.94; 
     0.01 0.05 0.1 0.2 0.3 0.5;
     0.01 0.05 0.1 0.2 0.3 0.5];

%% Energy Histogram
histTimeStep = 0.0010;
nTBins = round(impResTime/histTimeStep);
nFBins = length(FVect);
TFHist = zeros(nTBins,nFBins);

%% Ray Tracing
for iBand = 1:nFBins
    % Perform ray tracing independently for each frequency band.
    for iRay = 1:size(rays,1)
        % Select ray direction
        ray = rays(iRay,:);
        % All rays start at the source/transmitter
        ray_xyz = SetUpStruct.mic_pos;
        % Set initial ray direction. This direction changes as the ray is
        % reflected off surfaces.
        ray_dxyz = ray;
        % Initialize ray travel time. Ray tracing is terminated when the
        % travel time exceeds the impulse response length.
        ray_time = 0;
        % Initialize the ray energy to a normalized value of 1.     Energy
        % decreases when the ray hits a surface.
        ray_energy = 1;

        while (ray_time <= impResTime)

            % Determine the surface that the ray encounters
            [surfaceofimpact,displacement] = getImpactWall(ray_xyz,...
                                             ray_dxyz,SetUpStruct.room);
            
            % Determine the distance traveled by the ray
            distance = sqrt(sum(displacement.^2));

            % Determine the coordinates of the impact point
            impactCoord = ray_xyz+displacement;

            % Update ray location/source
            ray_xyz = impactCoord;

            % Update cumulative ray travel time
            c = 343; % speed of light (m/s)
            ray_time = ray_time+distance/c;

            % Apply surface reflection to ray's energy
            % This is the amount of energy that is not lost through
            % absorption.
            ray_energy = ray_energy*R(surfaceofimpact,iBand);

            % Apply diffuse reflection to ray energy
            % This is the fraction of energy used to determine what is
            % detected at the receiver
            rayrecv_energy = ray_energy*D(surfaceofimpact,iBand);

            % Determine impact point-to-receiver direction.
            rayrecvvector = SetUpStruct.mic_pos-impactCoord;

            % Determine the ray's time of arrival at receiver.
            distance = sqrt(sum(rayrecvvector.*rayrecvvector));
            recv_timeofarrival = ray_time+distance/c;

            if recv_timeofarrival>impResTime
                break
            end

            % Determine amount of diffuse energy that reaches the receiver.
            % See (5.20) in [2].

            % Compute received energy
            N = getWallNormalVector(surfaceofimpact);
            cosTheta = sum(rayrecvvector.*N)/(sqrt(sum(rayrecvvector.^2)));
            cosAlpha = sqrt(sum(rayrecvvector.^2)-SetUpStruct.mic_radius^2)/sum(rayrecvvector.^2);
            E = (1-cosAlpha)*2*cosTheta*rayrecv_energy;

            % Update energy histogram
            tbin = floor(recv_timeofarrival/histTimeStep + 0.5);
            TFHist(tbin,iBand) = TFHist(tbin,iBand) + E;

            % Compute a new direction for the ray.
            % Pick a random direction that is in the hemisphere of the
            % normal to the impact surface.
            d = rand(1,3);
            d = d/norm(d);
            if sum(d.*N)<0
                d = -d;
            end

            % Derive the specular reflection with respect to the incident
            % wall
            ref = ray_dxyz-2*(sum(ray_dxyz.*N))*N;

            % Combine the specular and random components
            d = d/norm(d);
            ref = ref/norm(ref);
            ray_dxyz = D(surfaceofimpact,iBand)*d+(1-D(surfaceofimpact,iBand))*ref;
            ray_dxyz = ray_dxyz/norm(ray_dxyz);
        end
    end
end

%% View
figure(1)
bar(histTimeStep*(0:size(TFHist)-1),TFHist)
grid on
xlabel("Time (s)")
legend(["125 Hz","250 Hz","500 Hz","1000 Hz","2000 Hz","4000 Hz"])

%% Generate Impulse Response
fs = 44100;
V = prod(SetUpStruct.room);
t0 = ((2*V*log(2))/(4*pi*c^3))^(1/3); % eq 5.45 in [2]
poissonProcess = [];
timeValues = [];
t = t0;
while (t<impResTime)
    timeValues = [timeValues t]; %#ok
    % Determine polarity.
    if (round(t*fs)-t*fs) < 0 
        poissonProcess = [poissonProcess 1]; %#ok
    else
        poissonProcess = [poissonProcess -1];%#ok
    end
    % Determine the mean event occurence (eq 5.44 in [2])
    mu = min(1e4,4*pi*c^3*t^2/V); 
    % Determine the interval size (eq. 5.44 in [2])
    deltaTA = (1/mu)*log(1/rand); % eq. 5.43 in [2])
    t = t+deltaTA;
end
randSeq = zeros(ceil(impResTime*fs),1);
for index=1:length(timeValues)
    randSeq(round(timeValues(index)*fs)) = poissonProcess(index);
end
flow = [115 225 450 900 1800 3600];
fhigh = [135 275 550 1100 2200 4400];
NFFT = 8192;
win = hann(882,"symmetric");
sfft = dsp.STFT(Window = win,OverlapLength=441,FFTLength=NFFT,FrequencyRange="onesided");
isfft = dsp.ISTFT(Window=win,OverlapLength=441,FrequencyRange="onesided");
F = sfft.getFrequencyVector(fs);
RCF = zeros(length(FVect),length(F));
for index0 = 1:length(FVect)
    for index=1:length(F)
        f = F(index);
        if f<FVect(index0) && f>=flow(index0)
            RCF(index0,index) = .5*(1+cos(2*pi*f/FVect(index0)));
        end
        if f<fhigh(index0) && f>=FVect(index0)
            RCF(index0,index) = .5*(1-cos(2*pi*f/(FVect(index0)+1)));
        end
    end
end
frameLength = 441;
numFrames = length(randSeq)/frameLength;
y = zeros(length(randSeq),6);
for index=1:numFrames
    x = randSeq((index-1)*frameLength+1:index*frameLength);
    X = sfft(x);
    X = X.*RCF.';
    y((index-1)*frameLength+1:index*frameLength,:) = isfft(X);
end
impTimes = (1/fs)*(0:size(y,1)-1);
hisTimes = histTimeStep/2 + histTimeStep*(0:nTBins);
W = zeros(size(impTimes,2),numel(FVect));
BW = fhigh-flow;
for k=1:size(TFHist,1)
    gk0 = floor((k-1)*fs*histTimeStep)+1;
    gk1 = floor(k*fs*histTimeStep);
    yy = y(gk0:gk1,:).^2;
    val = sqrt(TFHist(k,:)./sum(yy,1)).*sqrt(BW/(fs/2));
    for iRay=gk0:gk1
        W(iRay,:)= val;
    end
end
%% Create Impulse Response
y = y.*W;
ip = sum(y,2);
ip = ip./max(abs(ip));
vectorTiempo = (1/fs)*(0:numel(ip)-1);
figure
plot(vectorTiempo,ip)
grid on
xlabel("Time (s)")
ylabel("Impulse Response")

schroeder_cumsum = cumsum(flipud(ip.^2));
schroeder_normalized = schroeder_cumsum / max(schroeder_cumsum);
L = flipud(10*log10(schroeder_normalized));
[SignalSize, ~] = size(ip); 

vectorTiempo = (1:SignalSize)/fs;
[LF_EDT,LF_T20,LF_T30,LF_T60] = GetLinearFits(L,vectorTiempo);

[~, EDT_Idx] = min(abs(LF_EDT+60));
T60delEDT = vectorTiempo(EDT_Idx)
[~, T20_Idx] = min(abs(LF_T20+60));
T60delT20 = vectorTiempo(T20_Idx)
[~, T30_Idx] = min(abs(LF_T30+60));
T60delT30 = vectorTiempo(T30_Idx)
[~, T60_Idx] = min(abs(LF_T60+60));
T60delT60 = vectorTiempo(T60_Idx)

figure
plot(vectorTiempo,L,'LineWidth',2)
title('Schroeder Integration')
xlabel('Time (s)')
ylabel('Decay (Db)')
hold on
plot(vectorTiempo,LF_EDT)
plot(vectorTiempo,LF_T20)
plot(vectorTiempo,LF_T30)
plot(vectorTiempo,LF_T60)
ylim([-100, 0])
legend('Energy Decay Cruve','Linear fit for EDT','Linear fit for T20','Linear fit for T30', ...
     'Linear fit for T60','Location','southwest')
title('Linear fits')
xlabel('Time (s)')
ylabel('Decay (Db)')
hold off





