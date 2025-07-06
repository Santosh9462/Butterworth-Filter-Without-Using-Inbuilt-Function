clc;
clear;
close all;

%% PARAMETERS
fs = 4000;                   % Sampling frequency in Hz
t = 0:1/fs:1;                % Time vector (1 second)
f_signal = 50;              % Frequency of the sine wave

%% SIGNAL GENERATION
x = sin(2*pi*f_signal*t);           % Clean sine wave
noise = 0.5 * randn(size(t));       % Gaussian noise
x_noisy = x + noise;                % Noisy signal

%% FILTER DESIGN
filter_order = 4;

% Lowpass filter (cutoff at 500 Hz)
cutoff_lp = 500;
[bLow, aLow] = butter(filter_order, cutoff_lp/(fs/2), 'low');

% Highpass filter (cutoff at 300 Hz)
cutoff_hp = 300;
[bHigh, aHigh] = butter(filter_order, cutoff_hp/(fs/2), 'high');

% Bandpass filter (100–600 Hz)
cutoff_bp = [100 600];
[bBand, aBand] = butter(filter_order, cutoff_bp/(fs/2), 'bandpass');

%% BUILT-IN FILTERING
y_lp_builtin = filter(bLow, aLow, x_noisy);
y_hp_builtin = filter(bHigh, aHigh, x_noisy);
y_bp_builtin = filter(bBand, aBand, x_noisy);

%% MANUAL IIR FILTER IMPLEMENTATION

% Lowpass Manual Filtering
order_lp = length(aLow) - 1;
z_lp = zeros(order_lp,1);
y_lp_manual = zeros(size(x_noisy));
for n = 1:length(x_noisy)
    x_in = x_noisy(n);
    y_out = bLow(1)*x_in + z_lp(1);
    for i = 1:order_lp-1
        z_lp(i) = bLow(i+1)*x_in + z_lp(i+1) - aLow(i+1)*y_out;
    end
    z_lp(order_lp) = bLow(end)*x_in - aLow(end)*y_out;
    y_lp_manual(n) = y_out;
end

% Highpass Manual Filtering
order_hp = length(aHigh) - 1;
z_hp = zeros(order_hp,1);
y_hp_manual = zeros(size(x_noisy));
for n = 1:length(x_noisy)
    x_in = x_noisy(n);
    y_out = bHigh(1)*x_in + z_hp(1);
    for i = 1:order_hp-1
        z_hp(i) = bHigh(i+1)*x_in + z_hp(i+1) - aHigh(i+1)*y_out;
    end
    z_hp(order_hp) = bHigh(end)*x_in - aHigh(end)*y_out;
    y_hp_manual(n) = y_out;
end

% Bandpass Manual Filtering
order_bp = length(aBand) - 1;
z_bp = zeros(order_bp,1);
y_bp_manual = zeros(size(x_noisy));
for n = 1:length(x_noisy)
    x_in = x_noisy(n);
    y_out = bBand(1)*x_in + z_bp(1);
    for i = 1:order_bp-1
        z_bp(i) = bBand(i+1)*x_in + z_bp(i+1) - aBand(i+1)*y_out;
    end
    z_bp(order_bp) = bBand(end)*x_in - aBand(end)*y_out;
    y_bp_manual(n) = y_out;
end

%% FREQUENCY DOMAIN ANALYSIS
NFFT = length(t);
NFFT_half = floor(NFFT/2);  % Ensure integer indexing for the frequency axis
f_axis = (0:NFFT_half) * fs / NFFT;

% Compute FFTs
Y_lp_manual = fft(y_lp_manual, NFFT);
Y_hp_manual = fft(y_hp_manual, NFFT);
Y_bp_manual = fft(y_bp_manual, NFFT);
Y_input = fft(x, NFFT);
Y_noise_signal = fft(x_noisy, NFFT);

% Plot frequency spectra
figure;
subplot(2,3,1);
plot(f_axis, abs(Y_bp_manual(1:NFFT_half+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Bandpass Filtered (Manual)');

subplot(2,3,2);
plot(f_axis, abs(Y_hp_manual(1:NFFT_half+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Highpass Filtered (Manual)');

subplot(2,3,3);
plot(f_axis, abs(Y_lp_manual(1:NFFT_half+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Lowpass Filtered (Manual)');

subplot(2,3,4);
plot(f_axis, abs(Y_input(1:NFFT_half+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Input Sine Wave');

subplot(2,3,5);
plot(f_axis, abs(Y_noise_signal(1:NFFT_half+1)));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Noisy Signal');

%% TIME DOMAIN PLOTS
figure;
subplot(2,1,1);
plot(t, x);
title('Original Clean Signal'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(t, x_noisy);
title('Noisy Signal'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;

figure;
subplot(3,1,1);
plot(t, y_lp_manual);
title('Lowpass Filtered (Manual)'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(3,1,2);
plot(t, y_hp_manual);
title('Highpass Filtered (Manual)'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(3,1,3);
plot(t, y_bp_manual);
title('Bandpass Filtered (Manual)'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;

%% FREQUENCY RESPONSE OF FILTERS
[H_lp, f_resp] = freqz(bLow, aLow, 1024, fs);
[H_hp, ~] = freqz(bHigh, aHigh, 1024, fs);
[H_bp, ~] = freqz(bBand, aBand, 1024, fs);
[x_original, ~] = freqz(x, 1, 1024, fs);
[x_noise, ~] = freqz(x_noisy, 1, 1024, fs);

figure;
plot(f_resp, 20*log10(abs(H_lp)), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Lowpass Filter Frequency Response');
legend('Lowpass @ 500 Hz'); grid on; xlim([0 fs/2]);

figure;
plot(f_resp, 20*log10(abs(H_hp)), 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Highpass Filter Frequency Response');
legend('Highpass @ 300 Hz'); grid on; xlim([0 fs/2]);

figure;
plot(f_resp, 20*log10(abs(H_bp)), 'g', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Bandpass Filter Frequency Response');
legend('Bandpass 100–600 Hz'); grid on; xlim([0 fs/2]);

figure;
plot(f_resp, 20*log10(abs(x_original)), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Sine Wave Frequency Response');
legend('Sine Wave'); grid on; xlim([0 fs/2]);

figure;
plot(f_resp, 20*log10(abs(x_noise)), 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Noise Signal Frequency Response');
legend('Noise Signal'); grid on; xlim([0 fs/2]); ylim([-80 5]);

%% MANUAL VS BUILT-IN COMPARISON (Frequency Domain)
Y_lp_builtin = fft(y_lp_builtin, NFFT);
Y_hp_builtin = fft(y_hp_builtin, NFFT);
Y_bp_builtin = fft(y_bp_builtin, NFFT);

figure;
subplot(3,1,1);
plot(f_axis, 20*log10(abs(Y_lp_manual(1:NFFT_half+1))), 'b'); 
hold on;
plot(f_axis, 20*log10(abs(Y_lp_builtin(1:NFFT_half+1))), 'r--');
title('Lowpass Filtered: Manual vs Built-in'); xlabel('Freq (Hz)'); ylabel('Mag (dB)');
legend('Manual', 'Built-in'); grid on;

subplot(3,1,2);
plot(f_axis, 20*log10(abs(Y_bp_manual(1:NFFT_half+1))), 'g'); 
hold on;
plot(f_axis, 20*log10(abs(Y_bp_builtin(1:NFFT_half+1))), 'k--');
title('Bandpass Filtered: Manual vs Built-in'); xlabel('Freq (Hz)'); ylabel('Mag (dB)');
legend('Manual', 'Built-in'); grid on;

subplot(3,1,3);
plot(f_axis, 20*log10(abs(Y_hp_manual(1:NFFT_half+1))), 'c'); 
hold on;
plot(f_axis, 20*log10(abs(Y_hp_builtin(1:NFFT_half+1))), 'm--');
title('Highpass Filtered: Manual vs Built-in'); xlabel('Freq (Hz)'); ylabel('Mag (dB)');
legend('Manual', 'Built-in'); grid on;

%% DISPLAY FILTER COEFFICIENTS
disp('Lowpass Filter Coefficients:');
disp('bLow:'), disp(bLow);
disp('aLow:'), disp(aLow);

disp('Highpass Filter Coefficients:');
disp('bHigh:'), disp(bHigh);
disp('aHigh:'), disp(aHigh);

disp('Bandpass Filter Coefficients:');
disp('bBand:'), disp(bBand);
disp('aBand:'), disp(aBand);

%% DISPLAY FFT COEFFICIENTS OF SIGNALS
X_fft = fft(x, NFFT);
Noise_fft = fft(noise, NFFT);

disp('FFT Coefficients (first 10 values for inspection):');
disp('Clean Sine Wave FFT:'), disp(X_fft(1:10).');
disp('Noise FFT:'), disp(Noise_fft(1:10).');
disp('Lowpass Output FFT:'), disp(Y_lp_manual(1:10).');
disp('Highpass Output FFT:'), disp(Y_hp_manual(1:10).');
disp('Bandpass Output FFT:'), disp(Y_bp_manual(1:10).');



  




