function [IdealSRIR, Absorption, DecayRate, Clarity, ExpectedSoundStrength] = CalcAcousticParams(Volume,Surface,ReverberanceTime,vectorTiempo,hsq)
    Absorption = 0.161*Volume/Surface*ReverberanceTime;

    %RT = 0.161*V/sum(s*a);
    
    DecayRate = 60/ReverberanceTime; %Decay rate from reverberation time
    IdealSRIR = DecayRate*exp(-13.8*vectorTiempo/ReverberanceTime); %Distribution in time of impulse response squared, assuming exponential decay
    
    [~, ClarIdx] = min(abs(vectorTiempo-0.8));
    t_0a80 = vectorTiempo(1:ClarIdx); 
    t_80aEnd = vectorTiempo(ClarIdx+1:end); %Time intervals for integration
    hsq_0a80 = hsq(1:ClarIdx); 
    hsq_80aEnd = hsq(ClarIdx+1:end); %Amplitude intervals for integration
    
    Clarity = 10*log10(trapz(t_0a80,hsq_0a80)/trapz(t_80aEnd,hsq_80aEnd)); %Clarity calculation
    ExpectedSoundStrength = 10*log10(ReverberanceTime/Volume)+45; %Expected sound strength
end