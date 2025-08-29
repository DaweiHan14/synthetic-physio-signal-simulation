% Robotic Arm Control Signal Generation and Visualization System
% This program simulates the generation of EEG, EMG, and magnetic sensor data for robotic arm control research.
% The data generation code and results can be found in the supplementary materials.

clear; close all; clc;

%% Parameter Settings
num_participants = 30; 
sampling_rate = 1000; % Signal sampling rate (Hz)
duration = 60; % Data acquisition duration (seconds)
time_points = sampling_rate * duration; % Total number of data points

% Initialize data storage matrices
eeg_data = zeros(num_participants, time_points, 32); % 32-channel EEG data
emg_data = zeros(num_participants, time_points, 8);  % 8-channel EMG data
magnetic_data = zeros(num_participants, time_points, 6); % 6-axis magnetic/force-torque data

%% Generate EEG Data (Simulated EEG signals)
eeg_frequency_ranges = [1, 4; 4, 8; 8, 13; 13, 30; 30, 100]; % EEG bands: delta, theta, alpha, beta, gamma
fprintf('Generating EEG data...\n');

for p = 1:num_participants
    for ch = 1:32
        % Generate signals for each frequency band and superimpose them
        eeg_signal = zeros(1, time_points);
        for band = 1:size(eeg_frequency_ranges, 1)
            freq_range = eeg_frequency_ranges(band, :);
            center_freq = mean(freq_range);
            bandwidth = diff(freq_range);
            
            % Generate oscillatory signal for this frequency band
            band_signal = 0.1 * (band/5) * sin(2 * pi * center_freq * (0:time_points-1)/sampling_rate + rand()*2*pi);
            
            % Add amplitude modulation to simulate EEG rhythms
            modulation = 0.5 + 0.5 * sin(2 * pi * 0.1 * (0:time_points-1)/sampling_rate);
            band_signal = band_signal .* modulation;
            
            eeg_signal = eeg_signal + band_signal;
        end
        
        % Add noise and inter-individual variability
        eeg_data(p, :, ch) = eeg_signal + 0.02 * randn(1, time_points) + 0.05 * rand();
    end
end

%% Generate EMG Data (Simulated EMG signals)
fprintf('Generating EMG data...\n');

for p = 1:num_participants
    for ch = 1:8
        % Create baseline EMG signal (low-amplitude noise)
        base_signal = 0.1 * randn(1, time_points);
        
        % Add muscle activity bursts to simulate object grasping and fingertip movement
        num_bursts = randi([8, 20]); % Random number of bursts
        for i = 1:num_bursts
            % Random start time and duration for each burst
            start_time = randi([1, time_points - 1000]);
            burst_duration = randi([200, 800]); % 200-800 ms muscle activity
            end_time = min(start_time + burst_duration, time_points);
            
            % Generate burst signal (high-frequency oscillations)
            burst_freq = 50 + 40*rand(); % Random frequency 50-90Hz
            burst_amp = 0.8 + 0.7*rand(); % Random amplitude
            burst_signal = burst_amp * sin(2*pi*burst_freq*(0:1/sampling_rate:(burst_duration-1)/sampling_rate));
            
            % Use hann window for smooth envelope
            envelope = hann(length(burst_signal))';
            burst_signal = burst_signal .* envelope;
            
            % Ensure signal length matches
            if length(burst_signal) > (end_time - start_time + 1)
                burst_signal = burst_signal(1:(end_time - start_time + 1));
            end
            
            % Add burst to baseline signal
            base_signal(start_time:start_time+length(burst_signal)-1) = ...
                base_signal(start_time:start_time+length(burst_signal)-1) + burst_signal;
        end
        
        emg_data(p, :, ch) = base_signal;
    end
end

%% Generate Magnetic Data (Simulated force/torque sensor signals)
fprintf('Generating magnetic data...\n');

for p = 1:num_participants
    for sensor = 1:6
        % Create baseline signal (low noise)
        base_signal = 0.05 * randn(1, time_points);
        
        % Simulate force changes during grasping
        grip_events = randi([5, 15]); % Number of grasping events
        for i = 1:grip_events
            % Random start time and duration of grasping
            start_time = randi([1, time_points - 2000]);
            grip_duration = randi([1000, 3000]); % 1-3 seconds
            end_time = min(start_time + grip_duration, time_points);
            
            % Generate grasping force signal (smooth transition)
            grip_strength = 2 + 3*rand(); % Random grip strength
            grip_signal = grip_strength * (1 - exp(-(0:grip_duration-1)/200)) .* ...
                         exp(-(0:grip_duration-1)/(grip_duration/2));
            
            % Ensure signal length matches
            if length(grip_signal) > (end_time - start_time + 1)
                grip_signal = grip_signal(1:(end_time - start_time + 1));
            end
            
            % Add grip signal to baseline
            base_signal(start_time:start_time+length(grip_signal)-1) = ...
                base_signal(start_time:start_time+length(grip_signal)-1) + grip_signal;
        end
        
        magnetic_data(p, :, sensor) = base_signal;
    end
end

%% Data Filtering (remove environmental noise)
fprintf('Applying filters to remove environmental noise...\n');

% Design Butterworth bandpass filters
nyquist_freq = sampling_rate/2;

% Ensure frequencies do not exceed Nyquist
[b_eeg, a_eeg] = butter(4, [1, min(100, nyquist_freq*0.99)]/nyquist_freq, 'bandpass');
[b_emg, a_emg] = butter(4, [20, min(500, nyquist_freq*0.99)]/nyquist_freq, 'bandpass');
[b_mag, a_mag] = butter(4, [0.1, min(10, nyquist_freq*0.99)]/nyquist_freq, 'bandpass');

% Apply filters
for p = 1:num_participants
    for ch = 1:32
        eeg_data(p, :, ch) = filtfilt(b_eeg, a_eeg, eeg_data(p, :, ch));
    end
    
    for ch = 1:8
        emg_data(p, :, ch) = filtfilt(b_emg, a_emg, emg_data(p, :, ch));
    end
    
    for sensor = 1:6
        magnetic_data(p, :, sensor) = filtfilt(b_mag, a_mag, magnetic_data(p, :, sensor));
    end
end

%% Data Visualization
fprintf('Generating visualization plots...\n');

% Create time vector
time_vector = (0:time_points-1)/sampling_rate;

for p = 1:min(3, num_participants) % Show data for first 3 participants
    figure('Position', [100, 100, 1200, 800], 'Name', sprintf('Participant %d Visualization', p), 'Color', 'w');
    
    % Plot EEG
    subplot(3, 1, 1);
    plot(time_vector, eeg_data(p, :, 1));
    title(sprintf('EEG Signal (Participant %d, Channel 1)', p));
    xlabel('Time (s)'); ylabel('Amplitude (μV)');
    xlim([0, 10]); % Show first 10s
    grid on;
    
    % Plot EMG
    subplot(3, 1, 2);
    plot(time_vector, emg_data(p, :, 1));
    title(sprintf('EMG Signal (Participant %d, Channel 1)', p));
    xlabel('Time (s)'); ylabel('Amplitude (mV)');
    xlim([0, 10]);
    grid on;
    
    % Plot Magnetic data
    subplot(3, 1, 3);
    plot(time_vector, magnetic_data(p, :, 1));
    title(sprintf('Magnetic Data (Participant %d, Sensor 1)', p));
    xlabel('Time (s)'); ylabel('Force (N)');
    xlim([0, 10]);
    grid on;
    
    % Save figures
    saveas(gcf, sprintf('participant_%d_signals.png', p));
end

%% Time-Frequency Analysis (example: EEG spectrum)
figure('Position', [100, 100, 1000, 600], 'Color', 'w');
p = 1; ch = 1; % First participant, first channel

% Compute spectrum
[pxx, f] = pwelch(eeg_data(p, :, ch), 1024, 512, 1024, sampling_rate);

subplot(2, 1, 1);
plot(f, 10*log10(pxx));
title(sprintf('EEG Spectrum (Participant %d, Channel 1)', p));
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
xlim([0, 100]);
grid on;

% Highlight EEG bands
hold on;
bands = {'δ', 'θ', 'α', 'β', 'γ'};
band_ranges = [1,4; 4,8; 8,13; 13,30; 30,100];
colors = {'r', 'g', 'b', 'm', 'c'};
for i = 1:length(bands)
    x = [band_ranges(i,1), band_ranges(i,2), band_ranges(i,2), band_ranges(i,1)];
    y = [min(ylim), min(ylim), max(ylim), max(ylim)];
    patch(x, y, colors{i}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    text(mean(band_ranges(i,:)), max(ylim)-5, bands{i}, 'Color', colors{i}, 'FontWeight', 'bold');
end
hold off;

% Spectrogram
subplot(2, 1, 2);
spectrogram(eeg_data(p, :, ch), 256, 250, 256, sampling_rate, 'yaxis');
title(sprintf('EEG Time-Frequency (Participant %d, Channel 1)', p));
ylim([0, 100]);
colorbar;

% Save time-frequency figure
saveas(gcf, sprintf('participant_%d_eeg_time_frequency.png', p));

%% Save Data to Files
fprintf('Saving data to files...\n');

for p = 1:num_participants
    % Create tables
    eeg_table = array2table(squeeze(eeg_data(p, :, :)), 'VariableNames', strcat('EEG_Channel_', string(1:32)));
    emg_table = array2table(squeeze(emg_data(p, :, :)), 'VariableNames', strcat('EMG_Channel_', string(1:8)));
    magnetic_table = array2table(squeeze(magnetic_data(p, :, :)), 'VariableNames', strcat('Magnetic_Sensor_', string(1:6)));
    
    % Add time column
    time_table = table(time_vector', 'VariableNames', {'Time_s'});
    
    % Merge
    eeg_table = [time_table, eeg_table];
    emg_table = [time_table, emg_table];
    magnetic_table = [time_table, magnetic_table];
    
    % Save as CSV
    writetable(eeg_table, sprintf('participant_%d_eeg_data.csv', p));
    writetable(emg_table, sprintf('participant_%d_emg_data.csv', p));
    writetable(magnetic_table, sprintf('participant_%d_magnetic_data.csv', p));
end

%% Generate Data Report
fprintf('Generating data report...\n');

report_file = 'data_generation_report.txt';
fid = fopen(report_file, 'w');
fprintf(fid, 'Robotic Arm Control Signal Generation Report\n');
fprintf(fid, 'Generation Time: %s\n', datetime('now'));
fprintf(fid, 'Number of Participants: %d\n', num_participants);
fprintf(fid, 'Sampling Rate: %d Hz\n', sampling_rate);
fprintf(fid, 'Recording Duration: %d seconds\n', duration);
fprintf(fid, 'EEG Channels: 32\n');
fprintf(fid, 'EMG Channels: 8\n');
fprintf(fid, 'Magnetic Sensors: 6\n');
fprintf(fid, 'Generated Data Files:\n');
fprintf(fid, '  - participant_X_eeg_data.csv\n');
fprintf(fid, '  - participant_X_emg_data.csv\n');
fprintf(fid, '  - participant_X_magnetic_data.csv\n');
fprintf(fid, 'Generated Figures:\n');
fprintf(fid, '  - participant_X_signals.png\n');
fprintf(fid, '  - participant_X_eeg_time_frequency.png\n');
fclose(fid);

fprintf('All EEG, EMG, and magnetic data for participants have been generated and saved as CSV files.\n');
fprintf('The data generation code and results can be found in the supplementary materials.\n');
