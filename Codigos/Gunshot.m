function capture = Gunshot(fs,L,duration, attack_time, decay_time, NameofTheFile)
    device = "ASIO4ALL v2"; % Audio device name
    
    recChMap = 1;			% Recorder channel mapping
    playChMap = 1;			% Player channel mapping
    nbPlayCh = 1;			% Total number of playback channels
    nbRecCh = 1;			% Total number of recorder channels
    %% Create the excitation signal
    
    % Generate time vector
    t = (0:1/fs:duration)';
    % Generate the envelope
    envelope = ones(size(t));
    envelope(t <= attack_time) = linspace(0, 1, length(t(t <= attack_time)));
    envelope(t > attack_time & t <= attack_time + decay_time) = linspace(1, 0.1, length(t(t > attack_time & t <= attack_time + decay_time)));
    envelope(t > attack_time + decay_time) = 0.1 * exp(-(t(t > attack_time + decay_time) - (attack_time + decay_time)) / (decay_time / 5));
    gunshot = envelope .* (randn(size(t))); % Adding noise to simulate gunshot
    % Normalize the sound
    exc = (gunshot / max(abs(gunshot)));
    %% 
    
    % Allocate the input/output buffers
    sequenceLength = size(exc,1);
    bufExc = dsp.AsyncBuffer(sequenceLength+L);
    bufRec = dsp.AsyncBuffer(sequenceLength+2*L);
    
    % Copy the excitation to the output buffer (including one extra
    % frame of silence to account for minimum latency of one frame)
    write(bufExc,exc);
    write(bufExc,zeros(L,nbPlayCh));
    
    %% Play and capture using the selected device
    disp("Recording...")
    
    % Setup the capture device
    apr = audioPlayerRecorder( ...
        fs,Device=device, ...
        PlayerChannelMapping=playChMap, ...
        RecorderChannelMapping=recChMap);
    setup(apr,zeros(L,nbPlayCh));
    
    % Playback and capture loop
    while bufExc.NumUnreadSamples > 0
        x = read(bufExc,L);
        [y,under,over] = apr(x);
        write(bufRec,y);
        if under>0 || over>0
            error("Underrun or overrun occurred, terminating measurement")
        end
    end
    
    % Release the audio device
    release(apr);
    
    %% Compute the results
    disp("Computing results...")
    read(bufRec,L);
    capture = read(bufRec,sequenceLength);
    % subplot(2,1,1)
    % plot(t,exc)
    % title('Excitation Signal');
    % subplot(2,1,2)
    % plot(t,ir)
    % title('Impulse Response');
    audiowrite(NameofTheFile,capture,fs);
end