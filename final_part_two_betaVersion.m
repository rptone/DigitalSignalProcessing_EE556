%richard chester proppe
%professor F. J. Harris
%EE556
%projectCompletionDate: 05/08/2017



%%%%%%%%%%%%%%%%%%%%%%%%%   part_a_and_b   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%%%%%%%%%%%%%%%    create kaiser window operators   %%%%%%%%%%%%%%%%%%%%%%

%lowpass halfband window
L = 2048;
Beta = 10;
window = kaiser(L,Beta)';
window = window/(sum(window));

%bandpass window low, high
L2 = 1024;
Beta2 = 10;
window2 = kaiser(L2,Beta2)';
window2 = window2/(sum(window2));

%%%%%%%%%%%%%%%%%%%%%%  create composite signal  %%%%%%%%%%%%%%%%%%%%%%%%%%

%%Dual Tone Multiple Frequencies (DTMF), in Hertz
f_r1 = 697;
f_22 = 770;
f_r3 = 852;
f_r4 = 941;
f_c1 = 1209;
f_c2 = 1336;
f_c3 = 1477;
f_c4 = 1633;

%sampling frequency
fs = 8600;
%number of samples
N = 19999;
%discrete time sampling interval
%start at 0 and go to (8600 - 0.43), counting ever 0.43 frequencies
n_DTMF = 0:N;

%declare discrete time sin functions 
%Dual Tone Multiple Frequencies (DTMF) signals
%keypad rows
r1_DTMF = cos((2*pi*f_r1/fs)*n_DTMF);
r2_DTMF = cos((2*pi*f_22/fs)*n_DTMF);
r3_DTMF = cos((2*pi*f_r3/fs)*n_DTMF);
r4_DTMF = cos((2*pi*f_r4/fs)*n_DTMF);
%keypad columns
c1_DTMF = cos((2*pi*f_c1/fs)*n_DTMF);
c2_DTMF = cos((2*pi*f_c2/fs)*n_DTMF);
c3_DTMF = cos((2*pi*f_c3/fs)*n_DTMF);
c4_DTMF = cos((2*pi*f_c4/fs)*n_DTMF);
%sum of all rows and columns = composite signal
compositeSignal = r1_DTMF + r2_DTMF + r3_DTMF + r4_DTMF + c1_DTMF + ...
    c2_DTMF + c3_DTMF + c4_DTMF;
%plot
figure(1)
subplot(2,1,1)
%plot: composite signal
plot(compositeSignal)
axis([0 length(compositeSignal) min(compositeSignal) max(compositeSignal)])
title('Composite Signal = Sum of Eight DTMF Sinusoids')
xlabel('Samples = 20000')
grid on
%plot
subplot(2,1,2)
compositeSignal_fft = fftshift(20*log10(abs(fft(compositeSignal,1024))));
plot(compositeSignal_fft)
axis([0 length(compositeSignal_fft) min(compositeSignal_fft) ...
    max(compositeSignal_fft)])
title('FFT of the Composite Signal')
xlabel('Transform Length = 1024')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%    create filters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%   lowpass halfband filter parameters  %%%%%%%%%%%%%%%%%%%
n_LPHB = 6;
Rp_LPHB = 0.02;
Rs_LPHB = 60;
%filter coefficients, numerator and denmominator
[b, a] = remez(15,[0 1650 2670 4300]/(4300), [1 1 0 0], [20 1]);
%zero vector
zero_vector = zeros(1,400);
zero_vector(1) = 1;
%load filter coefficients and zero vector into array h
h = filter(b,a,zero_vector);
%plot
figure(2)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, fftshift(20*log10(abs(fft(h,1024)/20))))
grid on
title('Frequency Response: Remez Lowpass Halfband Filter, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%stem plot of DTMF signals
x1 = [697,770,852,941,1209,1336,1477,1633];
y1 = [-60,-60,-60,-60,-60,-60,-60,-60];
stem(x1, y1, 'r');
hold off
axis([0 4300 -60 10])
grid on
%plot
figure(3)
%pole zero diagram
zplane(b,a);
title('Pole and Zero Diagram: Remez Lowpass Halfband Filter')
%plot
figure(4)
%time responses
plot(h);
axis([0 60 min(h) max(h)])
xlabel('Time in Seconds');
ylabel('Attenuation in Decibles');
title('Time Response of the Remez Half Band Low Pass Filter')

%%%%%%%%%%%%%%%   bandpass lowband filter paramters   %%%%%%%%%%%%%%%%%%%%%

Rp_BPLB = 0.02;
Rs_BPLB = 50;
%filter coefficients
%calculate filter order
n_BPLB = round((4300/250)*(50/22));
%filter coefficients
[B_BPLB,A_BPLB] = remez(n_BPLB,[0 440 690 950 1200 (2150)]/(2150),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector_BPLB = zeros(1,400);
zero_vector_BPLB(1) = 1;
%load filter coefficients and zero vector into array h_BLB
h_BPLB = filter(B_BPLB,A_BPLB,zero_vector_BPLB);

%plot
figure(5);
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*4300,fftshift(20*log10(abs(fft(h_BPLB,1024)/40))));
hold on 
%stem plot of DTMF signals
x1 = [697,770,852,941,1209,1336,1477,1633];
y1 = [-80,-80,-80,-80,-80,-80,-80,-80];
stem(x1, y1, 'r');
hold off
title('Frequency Response: Remez Bandpass Lowband Filter, fs = 4300 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
axis([0 4300 -80 10])
grid on
figure(6)
%pole zero diagram bandpass lowband
zplane(B_BPLB,A_BPLB);
title('Pole and Zero Diagram: Remez Bandpass Lowband Filter')
figure(7)
%time response
plot(h_BPLB);
axis([0 100 min(h_BPLB) max(h_BPLB)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Time Response of the Remez Bandpass Lowband Filter')

%%%%%%%%%%%%%%%   bandpass highband filter paramters   %%%%%%%%%%%%%%%%%%%%

n_BPHB = 10;
Rp_BPHB = 0.02;
Rs_BPHB = 50;
%filter order
nB_highR = round((4300/250)*(50/22)); 
%filter coefficients
[B_BPHB,A_BPHB] = remez(nB_highR,[0 950 1200 1640 1890 (2150)]/(2150),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector_BPHB = zeros(1,400);
zero_vector_BPHB(1) = 1;
%load filter coefficients and zero vector into array h_BLB
h_BPHB = filter(B_BPHB,A_BPHB,zero_vector_BPHB);

%plot
figure(8);
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*4300,fftshift(20*log10(abs(fft(h_BPHB,1024)/40))));
grid on 
axis([0 4300 -80 10])
title('Frequency Response: Remez Bandpass Highband Filter, fs = 4300Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%stem plot of DTMF signals
x2 = [697,770,852,941,1209,1336,1477,1633];
y2 = [-80,-80,-80,-80,-80,-80,-80,-80];
stem(x2, y2, 'r');
hold off
%plot
figure(9)
%pole zero diagram
zplane(B_BPHB,A_BPHB);
title('Pole and Zero Diagram: Remez Bandpass Highband Filter')
%plot
figure(10)
%time response
plot(h_BPHB);
axis([0 60 min(h_BPHB) max(h_BPHB)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Time Response of the Remez Bandpass Highband Filter')

%%%%%%%%%%%%%%%%    first filtering: lowpass halfband    %%%%%%%%%%%%%%%%%%

%output signal: lowpass halfband
filter_output_LPHB = filter(b,a,compositeSignal);

%%%%%%%%%%%%%%%%%%%%%%%   first downsampling   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%downsample filter output 1:2
downsampled_output_LPHB = filter_output_LPHB(1:2:20000);
%plot
figure(11)
subplot(2,1,1)
plot(downsampled_output_LPHB)
axis([0 length(downsampled_output_LPHB) min(downsampled_output_LPHB) ...
    max(downsampled_output_LPHB)])
title('Downsampled Output: Remez Lowpass Halfband Filter')
ylabel('Amplitude')
xlabel('Samples = 10000')
hold on
%fft of downsampled output
fft_downsampled_output_LPHB = ...
    fftshift(20*log10(abs(fft(downsampled_output_LPHB(1:2048).*window))));
%plot
subplot(2,1,2)
%windowed fft
plot((-0.5:(1/2048):0.5-(1/2048))*4300, ...
    fftshift(20*log10(abs(fft(downsampled_output_LPHB(1:2048).*window)))))

axis([0 length(fft_downsampled_output_LPHB) ...
    min(fft_downsampled_output_LPHB) max(fft_downsampled_output_LPHB)])
title('Windowed FFT of Downsampled Output: Remez Lowpass Halfband Filter')
xlabel('Transform Length = 1024')
ylabel('Magnitude in decibels')
hold off

%%%%%%%%%%%%%%%   second filtering: bandpass lowband   %%%%%%%%%%%%%%%%%%%% 

%ouput signal bandpass lowband
filter_output_BPLB = filter(B_BPLB,A_BPLB,downsampled_output_LPHB);

%%%%%%%%%%%%%%   second downsampling:bandpass lowband  %%%%%%%%%%%%%%%%%%%%

%downsample filter output 1:4
downsampled_output_BPLB = filter_output_BPLB(1:4:10000);
%plot
figure(12)
subplot(2,1,1)
plot(downsampled_output_BPLB)
axis([0 length(downsampled_output_BPLB) min(downsampled_output_BPLB) ...
    max(downsampled_output_BPLB)])
title('Downsampled Output: Remez Bandpass Lowband Filter')
ylabel('Amplitude')
xlabel('Samples = 2500')
hold on
%plot
subplot(2,1,2)
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
title('Windowed FFT of Downsampled Output: Remez Bandpass Lowband Filter')
xlabel('Transform Length = 1024')
ylabel('Magnitude in decibels')
hold off

%%%%%%%%%%%%%%%   second filtering: bandpass highband   %%%%%%%%%%%%%%%%%%% 

%output signal bandpass highband
filter_output_BPHB = filter(B_BPHB,A_BPHB,downsampled_output_LPHB);


%%%%%%%% %%%%%   second downsampling: bandpass highband   %%%%%%%%%%%%%%%%%

%downsample filter output 1:4
downsampled_output_BPHB = filter_output_BPHB(1:4:10000);
%plot
figure(13)
subplot(2,1,1)
plot(downsampled_output_BPHB)
axis([0 length(downsampled_output_BPHB) min(downsampled_output_BPHB) ...
    max(downsampled_output_BPHB)])
title('Downsampled Output: Remez Bandpass Highband Filter')
ylabel('Amplitude')
xlabel('Samples = 2500')
hold on
%plot
subplot(2,1,2)
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
title('Windowed FFT of Downsampled Output: Remez Bandpass Highband Filter')
xlabel('Transform Length = 1024')
ylabel('Magnitude in decibels')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%      %part_c      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%   bandpass lowband set   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%new sampling frequency, 8:1 of fs
fs_partC  = 8600/8;
%filter coefficients
[B_LO_f1,A_LO_f1] = remez(50,[0 100 130 140 170 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f1 = filter(B_LO_f1,A_LO_f1,zero_vector);
%plot
figure(14)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_LO_f1,1024)/1075)*350)));
grid on
title('Frequency Response: Remez Bandpass Filter, LOF1, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_LO_f1);
axis([0 60 min(h_LO_f1) max(h_LO_f1)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, LOF1')
grid on

%filter coefficients
[B_LO_f2,A_LO_f2] = remez(39,[0 188 218 228 258 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f2 = filter(B_LO_f2,A_LO_f2,zero_vector);
%plot
figure(15)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_LO_f2,1024)/1075)*1000)));
grid on
title('Frequency Response: Remez Bandpass Filter, LOF2, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_LO_f2);
axis([0 60 min(h_LO_f2) max(h_LO_f2)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, LOF2')
grid on

%filter coefficients
[B_LO_f3,A_LO_f3] = remez(39,[0 269 299 309 339 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f3 = filter(B_LO_f3,A_LO_f3,zero_vector);
%plot
figure(16)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_LO_f3,1024)/1075)*1000)));
grid on
title('Frequency Response: Remez Bandpass Filter, LOF3, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off 
%plot
%time response
subplot(2,1,2)
plot(h_LO_f3);
axis([0 60 min(h_LO_f3) max(h_LO_f3)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, LOF3')
hold off
grid on

%filter coefficients
[B_LO_f4,A_LO_f4] = remez(39,[0 342 372 382 412 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f4 = filter(B_LO_f4,A_LO_f4,zero_vector);
%plot
figure(17)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_LO_f4,1024)/1075)*1000)));
grid on
title('Frequency Response: Remez Bandpass Filter, LOF4, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_LO_f4);
axis([0 60 min(h_LO_f4) max(h_LO_f4)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, LOF4')
hold off
grid on

%%%%%%%%%%%%%%%%%%%%%    bandpass highband set    %%%%%%%%%%%%%%%%%%%%%%%%%

%filter coefficients
[B_HI_f1,A_HI_f1] = remez(39,[0 99 129 139 169 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f1 = filter(B_HI_f1,A_HI_f1,zero_vector);
%plot
figure(18)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_HI_f1,1024)/1075)*1000)));
grid on
title('Frequency Response: Remez Bandpass Filter, HIF1, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f1);
axis([0 60 min(h_HI_f1) max(h_HI_f1)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, HIF1')
hold off
grid on

%filter coefficients
[B_HI_f2,A_HI_f2] = remez(39,[0 226 256 266 296 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f2 = filter(B_HI_f2,A_HI_f2,zero_vector);
%plot
figure(19)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_HI_f2,1024)/1075)*1000)));
grid on
title('Frequency Response: Remez Bandpass Filter, HIF2, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f2);
axis([0 60 min(h_HI_f2) max(h_HI_f2)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, HIF2')
hold off
grid on

%filter coefficients
[B_HI_f3,A_HI_f3] = remez(39,[0 367 397 407 437 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f3 = filter(B_HI_f3,A_HI_f3,zero_vector);
%plot
figure(20)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_HI_f3,1024)/1075)*1000)));
grid on
title('Frequency Response: Remez Bandpass Filter, HIF3, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f3);
axis([0 60 min(h_HI_f3) max(h_HI_f3)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, HIF3')
hold off
grid on

%filter coefficients
[B_HI_f4,A_HI_f4] = remez(40,[0 470 500 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1], [20 1]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f4 = filter(B_HI_f4,A_HI_f4,zero_vector);
%plot
figure(21)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,fftshift(20*log10(abs(fft(h_HI_f4,1024)/1075)*210)));
grid on
title('Frequency Response: Remez Bandpass Filter, HIF4, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f4);
axis([0 60 min(h_HI_f4) max(h_HI_f4)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: Remez Bandpass Filter, HIF4')
hold off
grid on

%%%%%%%%%%%%%%%%%%%%%%%%   %part_d   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assign keypad array values
keypad_1 = r1_DTMF + c1_DTMF;
keypad_2 = r1_DTMF + c2_DTMF;
keypad_3 = r1_DTMF + c3_DTMF;
keypad_A = r1_DTMF + c4_DTMF;
keypad_4 = r2_DTMF + c1_DTMF;
keypad_5 = r2_DTMF + c2_DTMF;
keypad_6 = r2_DTMF + c3_DTMF;
keypad_B = r2_DTMF + c4_DTMF;
keypad_7 = r3_DTMF + c1_DTMF;
keypad_8 = r3_DTMF + c2_DTMF;
keypad_9 = r3_DTMF + c3_DTMF;
keypad_C = r3_DTMF + c4_DTMF;
keypad_STAR = r4_DTMF + c1_DTMF;
keypad_0 = r4_DTMF + c2_DTMF;
keypad_POUND = r4_DTMF + c3_DTMF;
keypad_D = r4_DTMF + c4_DTMF;

%declare data array: myInitials 
myInitials_signal = zeros(1,120000);

%loop test variable, flag
flag = 0;

%column address variables
k1 = 1:20000;
c1 = 1:20000;
for loop = 1:2
    if loop == 1
        %k>1 && k<=20000
        myInitials_signal(1,k1) = keypad_7(1,c1);
        flag = 1;
    end
end

%column address variables
k2 = 40001:60000;
c2 = 1:20000;
for loop = 1:2
    if loop == 1
        %k>40000 && k<=60000
        myInitials_signal(1,k2) = keypad_2(1,c2);
        flag = 2;
    end 
end

%column address variables
k3 = 80001:100000;
c3 = 1:20000;
for loop = 1:2
    if loop == 1
        %k>80000 && k<=100000
        myInitials_signal(1,k3) = keypad_7(1,c3);
        flag = 3;
    end
end

%plot
figure(22)
subplot(2,1,1)
plot(myInitials_signal)
title('myInitials signal')
ylabel('Magnitude')
ylabel('Samples = 120000')
grid on
%plot
subplot(2,1,2)
fft_myInitials_signal=fftshift(20*log10(abs(fft(myInitials_signal,1024))));
plot(compositeSignal_fft)
axis([0 length(fft_myInitials_signal) min(fft_myInitials_signal) ...
    max(fft_myInitials_signal)])
title('FFT of myInitials  signal')
xlabel('Transform Length = 1024')
grid on

%%%%%%%%%%%%%%%    begin filtering "myInitials_signal"     %%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%    first filtering: lowpass halfband    %%%%%%%%%%%%%%%%%%

%output signal: lowpass halfband
filter_output_LPHB2 = filter(b,a,myInitials_signal);

%%%%%%%%%%%%%%%%%%%%%%%   first downsampling   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%downsample filter output 1:2
downsampled_output_LPHB2 = filter_output_LPHB2(1:2:20000);
%plot
figure(23)
subplot(2,1,1)
plot(downsampled_output_LPHB2)
axis([0 length(downsampled_output_LPHB2) min(downsampled_output_LPHB2) ...
    max(downsampled_output_LPHB2)])
title('Downsampled Output: Lowpass Halfband Filter')
ylabel('Amplitude')
xlabel('Samples = 10000')
hold on
%fft of downsampled output
fft_downsampled_output_LPHB2 = ...
    fftshift(20*log10(abs(fft(downsampled_output_LPHB2(1:2048).*window))));
%plot
subplot(2,1,2)
%windowed fft
plot((-0.5:(1/2048):0.5-(1/2048))*4300, ...
    fftshift(20*log10(abs(fft(downsampled_output_LPHB2(1:2048).*window)))))

axis([0 length(fft_downsampled_output_LPHB2) ...
    min(fft_downsampled_output_LPHB2) max(fft_downsampled_output_LPHB2)])
title('Windowed FFT of Downsampled Output: Lowpass Halfband Filter')
xlabel('Transform Length = 1024')
ylabel('Magnitude in decibels')
hold off

%%%%%%%%%%%%%%%   second filtering: bandpass lowband   %%%%%%%%%%%%%%%%%%%% 

%ouput signal bandpass lowband
filter_output_BPLB2 = filter(B_BPLB,A_BPLB,downsampled_output_LPHB2);

%%%%%%%%%%%%%%   second downsampling:bandpass lowband  %%%%%%%%%%%%%%%%%%%%

%downsample filter output 1:4
downsampled_output_BPLB2 = filter_output_BPLB2(1:4:10000);
%plot
figure(24)
subplot(2,1,1)
plot(downsampled_output_BPLB2)
axis([0 length(downsampled_output_BPLB2) min(downsampled_output_BPLB2) ...
    max(downsampled_output_BPLB2)])
title('Downsampled Output: Bandpass Lowband Filter')
ylabel('Amplitude')
xlabel('Samples = 2500')
hold on
%plot
subplot(2,1,2)
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB2(1:1024).*window2)))))
title('Windowed FFT of Downsampled Output: Bandpass Lowband Filter')
xlabel('Transform Length = 1024')
ylabel('Magnitude in decibels')
hold off

%%%%%%%%%%%%%%%   second filtering: bandpass highband   %%%%%%%%%%%%%%%%%%% 

%output signal bandpass highband
filter_output_BPHB2 = filter(B_BPHB,A_BPHB,downsampled_output_LPHB2);


%%%%%%%% %%%%%   second downsampling: bandpass highband   %%%%%%%%%%%%%%%%%

%downsample filter output 1:4
downsampled_output_BPHB2 = filter_output_BPHB2(1:4:10000);
%plot
figure(25)
subplot(2,1,1)
plot(downsampled_output_BPHB2)
axis([0 length(downsampled_output_BPHB2) min(downsampled_output_BPHB2) ...
    max(downsampled_output_BPHB2)])
title('Downsampled Output: Bandpass Highband Filter')
ylabel('Amplitude')
xlabel('Samples = 2500')
hold on
%plot
subplot(2,1,2)
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB2(1:1024).*window2)))))
title('Windowed FFT of Downsampled Output: Bandpass Highband Filter')
xlabel('Transform Length = 1024')
ylabel('Magnitude in decibels')
hold off


%%%%%%%%%%%%%%%%%%%%%   bandpass lowband set   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%new sampling frequency, 8:1 of fs
fs_partC  = 8600/8;
%filter coefficients
[B_LO_f1,A_LO_f1] = remez(50,[0 100 130 140 170 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f1 = filter(B_LO_f1,A_LO_f1,zero_vector);
%plot
figure(26)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_LO_f1,1024)/1075)*350)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, LOF1, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_LO_f1);
axis([0 60 min(h_LO_f1) max(h_LO_f1)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, LOF1')
grid on

%filter coefficients
[B_LO_f2,A_LO_f2] = remez(39,[0 188 218 228 258 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f2 = filter(B_LO_f2,A_LO_f2,zero_vector);
%plot
figure(27)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_LO_f2,1024)/1075)*1000)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, LOF2, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_LO_f2);
axis([0 60 min(h_LO_f2) max(h_LO_f2)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, LOF2')
grid on

%filter coefficients
[B_LO_f3,A_LO_f3] = remez(39,[0 269 299 309 339 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f3 = filter(B_LO_f3,A_LO_f3,zero_vector);
%plot
figure(28)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_LO_f3,1024)/1075)*1000)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, LOF3, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off 
%plot
%time response
subplot(2,1,2)
plot(h_LO_f3);
axis([0 60 min(h_LO_f3) max(h_LO_f3)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, LOF3')
hold off
grid on

%filter coefficients
[B_LO_f4,A_LO_f4] = remez(39,[0 342 372 382 412 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f4 = filter(B_LO_f4,A_LO_f4,zero_vector);
%plot
figure(29)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_LO_f4,1024)/1075)*1000)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, LOF4, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_LO_f4);
axis([0 60 min(h_LO_f4) max(h_LO_f4)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, LOF4')
hold off
grid on

%%%%%%%%%%%%%%%%%%%%%    bandpass highband set    %%%%%%%%%%%%%%%%%%%%%%%%%

%filter coefficients
[B_HI_f1,A_HI_f1] = remez(39,[0 99 129 139 169 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f1 = filter(B_HI_f1,A_HI_f1,zero_vector);
%plot
figure(30)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_HI_f1,1024)/1075)*1000)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, HIF1, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f1);
axis([0 60 min(h_HI_f1) max(h_HI_f1)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, HIF1')
hold off
grid on

%filter coefficients
[B_HI_f2,A_HI_f2] = remez(39,[0 226 256 266 296 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f2 = filter(B_HI_f2,A_HI_f2,zero_vector);
%plot
figure(31)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_HI_f2,1024)/1075)*1000)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, HIF2, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f2);
axis([0 60 min(h_HI_f2) max(h_HI_f2)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, HIF2')
hold off
grid on

%filter coefficients
[B_HI_f3,A_HI_f3] = remez(39,[0 367 397 407 437 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1 0 0], [20 1 20]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f3 = filter(B_HI_f3,A_HI_f3,zero_vector);
%plot
figure(32)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_HI_f3,1024)/1075)*1000)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, HIF3, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f3);
axis([0 60 min(h_HI_f3) max(h_HI_f3)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, HIF3')
hold off
grid on

%filter coefficients
[B_HI_f4,A_HI_f4] = remez(40,[0 470 500 (fs_partC/2)]/(fs_partC/2),...
    [0 0 1 1], [20 1]);
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f4 = filter(B_HI_f4,A_HI_f4,zero_vector);
%plot
figure(33)
subplot(2,1,1)
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*8600/8,...
    fftshift(20*log10(abs(fft(h_HI_f4,1024)/1075)*210)));
grid on
title('Frequency Response: My Initials, Remez Bandpass Filter, HIF4, fs = 8600/8 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
plot(h_HI_f4);
axis([0 60 min(h_HI_f4) max(h_HI_f4)])
xlabel('Time in Seconds');
ylabel('Attenuation in decibels');
title('Impulse Response: My Initials, Remez Bandpass Filter, HIF4')
hold off
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end of EE556  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%thank you

