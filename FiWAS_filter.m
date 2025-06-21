% FiWAS system filter
clc;clear; close all;
% Prerequest: signal processing toolbox, wavelet toolbox
% raw data should be in ".mat" format, as a variable named "voltage_array"
% such .mat file should contain sampling rate information in the file name,
% e.g. "heart raw data_5k.mat" (for sampling rate 5kHz)

%% Change accordingly
file_name = "lung raw data 5k";
% the threshold to dectect the peak as a heart beat
heart_beat_threshold = 0.3;
% cutoff this period of time at the beginning and end of the signal
cut_off_front_t = 0; % unit: s 
cut_off_end_t = 0;
% 0 for not writing audio, 1 for writing
audio_write_flag = 0;
% write raw data to audio
raw_data_audio_flag = 0;
% plot setting
title_font_size = 35;
label_font_size = 25;
xytick_font_size = 20;

%% Initialization
disp(file_name);
filename_char = convertStringsToChars(file_name);
mat_file = file_name + ".mat";
load(file_name); 
disp('[INFO] Successfully open raw data file');
% dbr = deep breathing; mdbr = midium deep breathing;
% sbr = shallow breathing; ndbr = nose deep breathing
if contains(filename_char,'dbr') || contains(filename_char,'sbr')...
        ||contains(filename_char,'lung')
    lung_filter = 1; % run lung filter
    disp('[INFO] Recognized as lung data');
else
    lung_filter = 0;
    disp('[INFO] Recognized as heart data');
end
if contains(filename_char,'25k')
    fs = 25000;
elseif contains(filename_char,'50k')
    fs = 50000;
elseif contains(filename_char,'10k')
    fs = 10000;
elseif contains(filename_char,'20k')
    fs = 20000;
elseif contains(filename_char,'5k')
    fs = 5000;
elseif contains(filename_char,'2k')
    fs = 2000;
elseif contains(filename_char,'1k')
    fs = 1000;
else
    fs = round(length(voltage_array)/Record_time);
end
Sample1 = voltage_array; % voltage_array is the raw data stored in mat file
raw_signal = voltage_array;
run_time = length(Sample1)/fs; % acquire data duration
sample_count = 1:length(Sample1);
sample_time = sample_count/fs;

%% Raw data plot
raw_data_plot = figure(1);
set(raw_data_plot,'Position',[50,50,1200,950])
subplot(2,1,1)
plot(sample_time, Sample1);
set(gca,'FontSize',xytick_font_size);
xlim([0 run_time])
grid on;
title("Raw Data",'FontSize',title_font_size)
xlabel("Time/s",'FontSize',label_font_size)
ylabel("Voltage/mV",'FontSize',label_font_size)
% raw data FFT plot
L = length(Sample1); 
sample_fft = fft(Sample1);
P2 = abs(sample_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
subplot(2,1,2)
plot(f,P1);
set(gca,'FontSize',xytick_font_size);
if lung_filter == 1
    xlim([20 500])
else
    xlim([2 200])
end
title("FFT of Raw Signal",'FontSize',title_font_size)
xlabel("f (Hz)",'FontSize',label_font_size)
grid on;
ylabel("|P1(f)|",'FontSize',label_font_size)

%% Filter

if lung_filter == 0
    disp('[INFO] Lunch heart filter for fiber sensor');
    disp('[PROCESS] Running high/low pass filter');
    Sample1 = high_pass(Sample1, fs, 40);
    Sample1 = low_pass(Sample1, fs, 150);
    disp('[PROCESS] Running notch filter');
    fprintf("Noise freq:");
    for i = 1:3
        multicancel = 1;
        for j = 1:multicancel
            L = length(Sample1);
            sample_fft = fft(Sample1);
            P2 = abs(sample_fft/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = fs*(0:(L/2))/L;
            % custimized notch filter can automatically find sharp peak
            % frequencies from the frequency range setted by FFT
            Sample1 = notch(50*i-40,50*i+1,fs,f,P1,Sample1,1);
        end
    end
    fprintf("\n");

    disp('[PROCESS] Cutting off abnormal part')
    Sample1 = cutoff_front(fs,cut_off_front_t,length(voltage_array),Sample1);
    Sample1 = cutoff_end(fs,cut_off_end_t,length(Sample1),Sample1);
    raw_signal = cutoff_front(fs,cut_off_front_t,length(voltage_array),raw_signal);
    raw_signal = cutoff_end(fs,cut_off_end_t,length(raw_signal),raw_signal);
    run_time = length(Sample1)/fs;
    sample_count = 1:length(Sample1);
    sample_time = sample_count/fs;
    disp('[PROCESS] Normalize signal')
    amp_rate = 1/max(abs(Sample1));
    Sample1 = amp_rate*Sample1;
    disp('[RESULT] Heart rate (beats/s):')
    heart_rate = heart_rate_calculation(Sample1,heart_beat_threshold,0.2,fs,run_time);
    disp(heart_rate)
    disp('[PROCESS] Normalize signal')
    amp_rate = 1/max(abs(Sample1));
    Sample1 = amp_rate*Sample1;
    amp_rate = 1/max(abs(raw_signal));
    raw_signal = amp_rate*raw_signal;
    heart_snr = snr(Sample1, raw_signal-Sample1);
    disp('[Result] Heart sound Signal-to-Noise Ratio is:');
    disp(heart_snr);
else
    disp('[INFO] Lunch Lung filter for fiber sensor');
    disp('[PROCESS] Running high/low pass filter');
    Sample1 = high_pass(Sample1, fs, 150);
    Sample1 = low_pass(Sample1, fs, 500);

    disp('[PROCESS] Running notch filter');
    for i = 1:20
        multicancel = 2;
        for j = 1:multicancel
            L = length(Sample1);
            sample_fft = fft(Sample1);
            P2 = abs(sample_fft/L);
            P1 = P2(1:floor(L/2+1));
            P1(2:end-1) = 2*P1(2:end-1);
            f = fs*(0:(L/2))/L;
            Sample1 = notch(50*i-49,50*i+2,fs,f,P1,Sample1,1);
        end
    end


    disp('[PROCESS] Filtering with wavelet transform');
    Sample1=wdenoise(Sample1,8,Wavelet='dmey');

    disp('[PROCESS] Cutting off abnormal part')
    Sample1 = cutoff_front(fs,cut_off_front_t,length(voltage_array),Sample1);
    Sample1 = cutoff_end(fs,cut_off_end_t,length(Sample1),Sample1);
    raw_signal = cutoff_front(fs,cut_off_front_t,length(voltage_array),raw_signal);
    raw_signal = cutoff_end(fs,cut_off_end_t,length(raw_signal),raw_signal);
    run_time = length(Sample1)/fs;
    sample_count = 1:length(Sample1);
    sample_time = sample_count/fs;

    disp('[PROCESS] Normalize signal')
    amp_rate = 1/max(abs(Sample1));
    Sample1 = amp_rate*Sample1;
    amp_rate = 1/max(abs(raw_signal));
    raw_signal = amp_rate*raw_signal;
    lung_snr = snr(Sample1, raw_signal-Sample1);
    disp('[Result] Lung sound Signal-to-Noise Ratio is:');
    disp(lung_snr);

end

%% Plot Filtered Data
filtered_plot = figure(2);
set(filtered_plot,'Position',[50,50,1200,950])
subplot(2,1,1)
plot(sample_time,Sample1);
set(gca,'FontSize',xytick_font_size);
xlim([0 run_time])
title("Filtered Data",'FontSize',title_font_size)
xlabel("Time/s",'FontSize',label_font_size)
% ylabel("Amplitude",'FontSize',label_font_size)
ylabel("Voltage/mV",'FontSize',label_font_size)
grid on;
ylim([-1.2 1.2])

L = length(Sample1);
sample_fft = fft(Sample1);
P2 = abs(sample_fft/L);
P1 = P2(1:floor(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
subplot(2,1,2)
plot(f,P1);
set(gca,'FontSize',xytick_font_size);
if lung_filter == 1
    xlim([20 500])
else
    xlim([2 200])
end
title("FFT of Signal after filter",'FontSize',title_font_size)
xlabel("f (Hz)",'FontSize',label_font_size)
grid on;
ylabel("|P1(f)|",'FontSize',label_font_size)

%% Save the audio signal as a WAV file
if audio_write_flag == 1
    disp("[PROCESS] Writing filtered signal audio");
    amp_rate = 0.9/max(abs(Sample1)); % normalize signal for audio conversion
    Sample1 = amp_rate*Sample1;
    wav_file = file_name + "_SF.wav";
    audiowrite(wav_file, Sample1, fs,"BitsPerSample",24);
end
if raw_data_audio_flag == 1
    disp("[PROCESS] Writing raw signal audio");
    amp_rate = 0.9/max(abs(voltage_array));  % normalize signal for audio conversion
    Sample1 = amp_rate*voltage_array;
    wav_file = file_name + "_Raw.wav";
    audiowrite(wav_file, Sample1, fs,"BitsPerSample",24);
end

%% Functions
function filtered = notch(low_limit, high_limit,fs,freq_array,freq_amp,raw_data,varargin)
    f = freq_array;
    P1 = freq_amp;
    Sample1 = raw_data;
    peak_range = find(f>low_limit & f<high_limit);
    if isempty(peak_range)
        fprintf("\n[ERROR] freq range around %.2fHz too narrow\n",mean(low_limit,high_limit));
    else
        [~, peak_index] = max(P1(peak_range));
        peak_index = peak_index+peak_range(1)-1;
        peak_freq = round(f(peak_index),2);
        % Notch filter 
        notch_freq = peak_freq; 
        if nargin == 6
            bw = 5;
        else
            bw = varargin{1}; %  Bandwidth for the notch filter
        end
        [bn, an] = iirnotch(notch_freq/(fs/2), bw/(fs/2));
        Sample1 = filter(bn, an, Sample1); 
        fprintf("%d ",peak_freq);
    end
    filtered = Sample1;
end

function filtered = low_pass(raw_data, fs, freq)
    %Low-pass filter
    Fc_lp = freq; % Cut-off frequency of the filter in Hz, adjust as needed
    [b_lp, a_lp] = butter(4, Fc_lp/(fs/2)); % 4th order Butterworth low-pass filter
    filtered = filter(b_lp, a_lp, raw_data);
end

function filtered = high_pass(raw_data, fs, freq)
    % High-pass filter
    Fc_hp = freq; % High-pass cut-off frequency in Hz, adjust as needed
    [b_hp, a_hp] = butter(4, Fc_hp/(fs/2), 'high'); % 4nd order Butterworth high-pass filter
    filtered = filter(b_hp, a_hp,raw_data);
end

function cutoffed = cutoff_front(fs,cutoff_time,orignial_len,Sample1)
    cutoff_time = cutoff_time*fs;
    if cutoff_time == 0
        cutoff_time = 1;
    end
    cutoffed = Sample1(round(cutoff_time):orignial_len);
end

function cutoffed = cutoff_end(fs,cutoff_time,orignial_len,Sample1)
    cutoff_time = cutoff_time*fs;
    % make sure the length of data is an even number
    if cutoff_time == 0
        cutoffed = Sample1;
    else
        if mod(orignial_len-2*cutoff_time,2) == 0
            cutoffed = Sample1(1:round(orignial_len-cutoff_time+1));
        else
            cutoffed = Sample1(1:round(orignial_len-cutoff_time));
        end
    end
end

function heart_rate = heart_rate_calculation(raw_data,heart_threshold,min_peak_distance,fs,period)
    Sample1 = raw_data;
    [~,peak_loc] = findpeaks(abs(Sample1),'MinPeakProminence', heart_threshold,...
        'MinPeakDistance', min_peak_distance*fs);
%     figure(4);
%     peak_loc = peak_loc/fs;
%     findpeaks(abs(Sample1),'MinPeakProminence', heart_threshold,...
%         'MinPeakDistance', min_peak_distance*fs);
    heart_rate = round(length(peak_loc)/2/period*60);
    
end








