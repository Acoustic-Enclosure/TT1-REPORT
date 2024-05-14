clc; clear; close all;

%% Method Selection
    Method = "Sweep Sine";
    %Method = "Gunshot";

if strcmp(Method,"Sweep Sine")
    %Method config
    SampleRate = 44100;			% Sample rate (Hz)
    samplesPerFrame = 1024; 	% Samples per frame
    durationPerRun = 3; 	    % Duration per Run (s)
    outputLevel = -6;		    % Excitation Level (dBFS)
    nbWarmUps = 0;			    % Number of warm-up runs
    sweepRange = [10 22000];    % Sweep start/stop frequency (Hz)
    sweepDur = 1;			    % Sweep Duration (s)
    %Impulse Response Adquisition
    Results = Sweep(SampleRate,samplesPerFrame,durationPerRun, ...
        outputLevel,nbWarmUps,sweepRange,sweepDur,'ImpulseResponse.wav');
else
    %Method config
    SamplesPerFrame = 1024; 
    fs = 44100; % Sampling frequency (Hz)
    duration = 1; % Duration of the sound (seconds)
    attack_time = 0.01; % Duration of the attack (seconds)
    decay_time = 0.05; % Duration of the decay (seconds)
    %Impulse Response Adquisition
    Results = Gunshot(fs,SamplesPerFrame,duration,attack_time, ...
        decay_time,'ImpulseResponse.wav');
end
%% Calculate T using function from .wav
T_func = iosr.acoustics.irStats('SimulatedRIR.wav')

%% Read response from .wav
[impulseResponse_wav, fs_wav] = audioread('SimulatedRIR.wav');
if size(impulseResponse_wav,2) > 1
    impulseResponse_wav = mean(impulseResponse_wav, 2); %Create one channel
end
sound_duration = length(impulseResponse_wav)/fs_wav;    %Sound_duration
vectorTiempo = (0:1/fs_wav:sound_duration-1/fs_wav)';   %Time vector
hsq = impulseResponse_wav.^2;   %Squared Impulse Response

figure(1)
subplot(2,2,1)
plot(vectorTiempo,hsq)
title('Squared Impulse Response')
ylabel('Amplitude')
xlabel('Time (s)')

%% Energy Decay Curve
IRf = ACweight_time_filter(3,hsq);
IRfh = abs(hilbert(IRf));
IRfhn = (IRfh-min(IRfh))./(max(IRfh)-min(IRfh));
windowWidth = 300; 
kernel = ones(1,windowWidth) / windowWidth;
IRfhnf = filter(kernel,1, IRfhn);
E_t = 20*log10(IRfhnf/max(IRfhnf));

subplot(2,2,2)
plot(vectorTiempo,E_t)
title('Energy Decay Curve')
xlabel('Time (s)')
ylabel('Energy (Db)')

%% Schroeder Integration
schroeder_cumsum = cumsum(flipud(hsq));
schroeder_normalized = schroeder_cumsum / max(schroeder_cumsum);
L = flipud(10*log10(schroeder_normalized));

subplot(2,2,3)
plot(vectorTiempo,L)
title('Schroeder Integration')
xlabel('Time (s)')
ylabel('Decay (Db)')

%% Linear Fits
[LF_EDT,LF_T20,LF_T30,LF_T60] = GetLinearFits(L,vectorTiempo);

subplot(2,2,4)
hold on
plot(vectorTiempo,LF_EDT)
plot(vectorTiempo,LF_T20)
plot(vectorTiempo,LF_T30)
plot(vectorTiempo,LF_T60)
ylim([-100, 0])
legend('Linear fit for EDT','Linear fit for T20','Linear fit for T30', ...
     'Linear fit for T60','Location','southwest')
title('Linear fits')
xlabel('Time (s)')
ylabel('Decay (Db)')
hold off

%% Reverberance time calculation
[~, EDT_Idx] = min(abs(LF_EDT+60));
T60delEDT = vectorTiempo(EDT_Idx)
[~, T20_Idx] = min(abs(LF_T20+60));
T60delT20 = vectorTiempo(T20_Idx)
[~, T30_Idx] = min(abs(LF_T30+60));
T60delT30 = vectorTiempo(T30_Idx)
[~, T60_Idx] = min(abs(LF_T60+60));
T60delT60 = vectorTiempo(T60_Idx)

%% Acoustic parameters
V = 6*4*3.5; %Volume of the room
S = 120; %Surface of the room

[IdealSRIR, Absorption, DecayRate, Clarity, ExpectedSoundStrength] = CalcAcousticParams(V,S,T60delT60,vectorTiempo,hsq);

fprintf("Absorption = %f \n",Absorption)
fprintf("Decay Rate = %f \n",DecayRate)
fprintf("Clarity = %f \n",Clarity)
fprintf("ExpectedSoundStrength = %f \n",ExpectedSoundStrength)

figure(2)
plot(vectorTiempo,IdealSRIR)
title('Predicted Squared Impulse Response')
xlabel('Time (s)')
ylabel('Amplitude')