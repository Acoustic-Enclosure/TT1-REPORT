clear; close all; clc;

%% Config
SetUpStruct.Fs = 16000;
SetUpStruct.room = [6 4 3.4];
SetUpStruct.mic_pos = [3 2 1.7;3.5 2.5 1.7;2.7 1.7 1.7];
SetUpStruct.T60 = 0.85;
SetUpStruct.abs_weights = [1,1,1,1,0.5,0.5];
SetUpStruct.src_pos = [0.3 0.3 0.3];
[SetUpStruct.AbsCoeffs, SetUpStruct.OkFlag] = ISM_AbsCoeff('t60', SetUpStruct.T60, SetUpStruct.room, ...
    SetUpStruct.abs_weights, 'LehmannJohansson');
SetUpStruct.RefCoeffs = sqrt(1-SetUpStruct.AbsCoeffs);

fprintf('Needed Absorption Coefficientes = [%f,%f,%f,%f,%f,%f] \n', SetUpStruct.AbsCoeffs);

plot3(SetUpStruct.src_pos(:,1),SetUpStruct.src_pos(:,2),SetUpStruct.src_pos(:,3),'ro-','markersize',4);
hold on
plot3(SetUpStruct.mic_pos(:,1),SetUpStruct.mic_pos(:,2),SetUpStruct.mic_pos(:,3),'ko','markerfacecolor',ones(1,3)*.6);
axis equal; axis([0 SetUpStruct.room(1) 0 SetUpStruct.room(2) 0 SetUpStruct.room(3)]);
box on; xlabel('x-axis (m)'); ylabel('y-axis (m)'); zlabel('z-axis (m)'); grid on;
title('Room and layout')

%% Simulate using Image Source Method
[NumberOfMics,~] = size(SetUpStruct.mic_pos);
for i = 1:NumberOfMics
    [temp_RIRt, okf] = fast_ISM_RoomResp(SetUpStruct.Fs,SetUpStruct.RefCoeffs,'t60',SetUpStruct.T60,SetUpStruct.src_pos, ...
        SetUpStruct.mic_pos(i,:),SetUpStruct.room,'Delta_dB',70);
    
    if okf == 1
        [Length_temp_RIRt, ~] = size(temp_RIRt);
        if i ==1
            RIRt(:,i) = temp_RIRt;     
        else
            [LengthRIRt, ~] = size(RIRt);
            if Length_temp_RIRt > LengthRIRt 
                RIRt(LengthRIRt+1:Length_temp_RIRt,i:-1:1) = 0;
            else
                temp_RIRt(end:LengthRIRt) = 0;
                RIRt(:,i) = temp_RIRt;
            end
        end
    end
end
RIRts = mean(RIRt,2).^2;
%% Create time vector
[SignalSize,~] = size(RIRts);
vectorTiempo = (1:SignalSize)/SetUpStruct.Fs;

%% Plot RIR
figure
subplot(1,2,1)
plot(vectorTiempo,RIRts)
xlabel('Time(s)'); ylabel('Amplitude'); title('Squared RIR');
% audiowrite('SimulatedRIR.wav',RIRt,SetUpStruct.Fs);

%% Calculate RT and Linear Fit
schroeder_cumsum = cumsum(flipud(RIRts));
schroeder_normalized = schroeder_cumsum / max(schroeder_cumsum);
L = flipud(10*log10(schroeder_normalized));

[LF_EDT,LF_T20,LF_T30,LF_T60] = GetLinearFits(L,vectorTiempo);

%% T60 calc
[~, EDT_Idx] = min(abs(LF_EDT+60));
T60delEDT = vectorTiempo(EDT_Idx)
[~, T20_Idx] = min(abs(LF_T20+60));
T60delT20 = vectorTiempo(T20_Idx)
[~, T30_Idx] = min(abs(LF_T30+60));
T60delT30 = vectorTiempo(T30_Idx)
[~, T60_Idx] = min(abs(LF_T60+60));
T60delT60 = vectorTiempo(T60_Idx)

%% Plot RT and fits
subplot(1,2,2)
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
