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
[b,a] = ellip(n_LPHB, Rp_LPHB, Rs_LPHB,(1650/(8600/2)), 'low');
%zero vector
zero_vector = zeros(1,200);
zero_vector(1) = 1;
%load filter coefficients and zero vector into array h
h = filter(b,a,zero_vector);
%zero, poles, and gain
[z,p,k] = ellip( n_LPHB, Rp_LPHB, Rs_LPHB, (1650/(8600/2)), 's');
%construct transfer function
[Be,Ae] = zp2tf(z,p,k);
%plot
figure(2)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h,1024)))))
grid on
title('Frequency Response: Lowpass Halfband Filter, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%stem plot of DTMF signals
x1 = [697,770,852,941,1209,1336,1477,1633];
y1 = [-80,-80,-80,-80,-80,-80,-80,-80];
stem(x1, y1, 'r');
hold off
axis([0 4300 -80 10])
grid on
%plot
figure(3)
%pole zero diagram
zplane(b,a);
title('Pole and Zero Diagram: Lowpass Halfband Filter')
%plot
figure(4)
%time responses
subplot(2,1,1)
step(Be,Ae)
title('Step Response: Lowpass Halfband Filter')
grid on 
hold on
%plot
subplot(2,1,2)
impulse(Be,Ae)
title('Impulse Response: Lowpass Halfband Filter')
grid on

%%%%%%%%%%%%%%%   bandpass lowband filter paramters   %%%%%%%%%%%%%%%%%%%%%

fs_BPLB = 4300;
n_BPLB = 10;
Rp_BPLB = 0.02;
Rs_BPLB = 50;
%filter coefficients
[B_BPLB,A_BPLB] = ellip((n_BPLB/2),Rp_BPLB,Rs_BPLB,[690 950]/(4300/2), ...
    'bandpass');
%zero vector
zero_vector_BPLB = zeros(1,400);
zero_vector_BPLB(1) = 1;
%load filter coefficients and zero vector into array h_BLB
h_BPLB = filter(B_BPLB,A_BPLB,zero_vector_BPLB);
%plot
figure(5);
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*4300, ...
    fftshift(20*log10(abs(fft(h_BPLB,1024)))));
hold on 
%stem plot of DTMF signals
x1 = [697,770,852,941,1209,1336,1477,1633];
y1 = [-80,-80,-80,-80,-80,-80,-80,-80];
stem(x1, y1, 'r');
hold off
title('Frequency Response: Bandpass Lowband Filter, fs = 4300 Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
axis([0 4300 -80 10])
grid on
figure(6)
%zeros, poles, and gain 
[Z_BPLB,P_BPLB,K_BPLB] = ellip(n_BPLB, Rp_BPLB, Rs_BPLB,(690/(4300/2)), ...
    's');
%pole zero diagram bandpass lowband
zplane(B_BPLB,A_BPLB);
title('Pole and Zero Diagram: Bandpass Lowband Filter')
figure(7)
%time response
%transfer function coefficients
[Be1,Ae1] = zp2tf(Z_BPLB,P_BPLB,K_BPLB);
%plot
subplot(2,1,1)
step(Be1,Ae1)
title('Step Response: Lowpass Halfband Filter')
grid on
%plot
subplot(2,1,2)
impulse(Be1,Ae1)
title('Impulse Response: Lowpass Halfband Filter')
grid on

%%%%%%%%%%%%%%%   bandpass highband filter paramters   %%%%%%%%%%%%%%%%%%%%

n_BPHB = 10;
Rp_BPHB = 0.02;
Rs_BPHB = 50;
%filter coefficients
[B_BPHB,A_BPHB] = ellip((n_BPHB/2),Rp_BPHB,Rs_BPHB,[1200 1640]/(4300/2),...
    'bandpass');
%zero vector
zero_vector_BPHB = zeros(1,400);
zero_vector_BPHB(1) = 1;
%load filter coefficients and zero vector into array h_BLB
h_BPHB = filter(B_BPHB,A_BPHB,zero_vector_BPHB);
%plot
figure(8);
%frequency response
plot((-0.5:1/1024:0.5-(1/1024))*4300, ...
    fftshift(20*log10(abs(fft(h_BPHB,1024)))));
grid on 
axis([0 4300 -80 10])
title('Frequency Response: Bandpass Highband Filter, fs = 4300Hz')
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
%zeros, poles, and gain
[Z_BPHB,P_BPHB,K_BPHB] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,(1200/(4300/2)),...
    's');
%construct transfer function bandpass highband
[Be2,Ae2] = zp2tf(Z_BPHB,P_BPHB,K_BPHB);
%pole zero diagram
zplane(B_BPHB,A_BPHB);
title('Pole and Zero Diagram: Bandpass Highband Filter')
%plot
figure(10)
%time response
subplot(2,1,1)
step(Be2,Ae2)
title('Step Response: Lowpass Halfband Filter')
grid on
subplot(2,1,2)
impulse(Be2,Ae2)
title('Impulse Response: Lowpass Halfband Filter')
grid on

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
title('Downsampled Output: Lowpass Halfband Filter')
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
title('Windowed FFT of Downsampled Output: Lowpass Halfband Filter')
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
title('Downsampled Output: Bandpass Lowband Filter')
ylabel('Amplitude')
xlabel('Samples = 2500')
hold on
%plot
subplot(2,1,2)
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB(1:1024).*window2)))))
title('Windowed FFT of Downsampled Output: Bandpass Lowband Filter')
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
title('Downsampled Output: Bandpass Highband Filter')
ylabel('Amplitude')
xlabel('Samples = 2500')
hold on
%plot
subplot(2,1,2)
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075, ...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB(1:1024).*window2)))))
title('Windowed FFT of Downsampled Output: Bandpass Highband Filter')
xlabel('Transform Length = 1024')
ylabel('Magnitude in decibels')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%      %part_c      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%   bandpass lowband set   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%filter coefficients
[B_LO_f1,A_LO_f1] = ellip(3,0.02,60,([128 138])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f1,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F1, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_LO_f1,P_LO_f1,K_LO_f1] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f1,Ae_LO_f1] = zp2tf(Z_LO_f1,P_LO_f1,K_LO_f1);
impulse(Be_LO_f1,Ae_LO_f1)
title('Impulse Response: Bandpass Filter, LO_F1')
grid on

%filter coefficients
[B_LO_f2,A_LO_f2] = ellip(3,0.02,60,([219 229])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f2,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F2, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_LO_f2,P_LO_f2,K_LO_f2] = ellip(n_BPHB, ...
    Rp_BPHB, Rs_BPHB,(1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f2,Ae_LO_f2] = zp2tf(Z_LO_f2,P_LO_f2,K_LO_f2);
impulse(Be_LO_f2,Ae_LO_f2)
title('Impulse Response: Bandpass Filter, LO_F2')
grid on

%filter coefficients
[B_LO_f3,A_LO_f3] = ellip(3,0.02,60,([301 311])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f3,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F3, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_LO_f3,P_LO_f3,K_LO_f3] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB, ...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f3,Ae_LO_f3] = zp2tf(Z_LO_f3,P_LO_f3,K_LO_f3);
impulse(Be_LO_f3,Ae_LO_f3)
title('Impulse Response: Bandpass Filter, LO_F3')
hold off
grid on

%filter coefficients
[B_LO_f4,A_LO_f4] = ellip(3,0.02,60,([373 383])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f4,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F4, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_LO_f4,P_LO_f4,K_LO_f4] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB, ...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f4,Ae_LO_f4] = zp2tf(Z_LO_f4,P_LO_f4,K_LO_f4);
impulse(Be_LO_f4,Ae_LO_f4)
title('Impulse Response: Bandpass Filter, LO_F4')
hold off
grid on

%%%%%%%%%%%%%%%%%%%%%    bandpass highband set    %%%%%%%%%%%%%%%%%%%%%%%%%

%filter coefficients
[B_HI_f1,A_HI_f1] = ellip(3,0.02,60,([128 138])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f1,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F1, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_HI_f1,P_HI_f1,K_HI_f1] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB, ...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f1,Ae_HI_f1] = zp2tf(Z_HI_f1,P_HI_f1,K_HI_f1);
impulse(Be_HI_f1,Ae_HI_f1)
title('Impulse Response: Bandpass Filter, HI_F1')
hold off
grid on

%filter coefficients
[B_HI_f2,A_HI_f2] = ellip(3,0.02,60,([255 265])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f2,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F2, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_HI_f2,P_HI_f2,K_HI_f2] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB, ...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f2,Ae_HI_f2] = zp2tf(Z_HI_f2,P_HI_f2,K_HI_f2);
impulse(Be_HI_f2,Ae_HI_f2)
title('Impulse Response: Bandpass Filter, HI_F2')
hold off
grid on

%filter coefficients
[B_HI_f3,A_HI_f3] = ellip(3,0.02,60,([397 407])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f3,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F3, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_HI_f3,P_HI_f3,K_HI_f3] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB, ...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f3,Ae_HI_f3] = zp2tf(Z_HI_f3,P_HI_f3,K_HI_f3);
impulse(Be_HI_f3,Ae_HI_f3)
title('Impulse Response: Bandpass Filter, HI_F3')
hold off
grid on

%filter coefficients
[B_HI_f4,A_HI_f4] = ellip(3,0.02,60,([511 521])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f4,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F4, fs = 8600Hz')
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
%zeros, poles, and gain
[Z_HI_f4,P_HI_f4,K_HI_f4] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB, ...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f4,Ae_HI_f4] = zp2tf(Z_HI_f4,P_HI_f4,K_HI_f4);
impulse(Be_HI_f4,Ae_HI_f4)
title('Impulse Response: Bandpass Filter, HI_F4')
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

%filter coefficients
[B_LO_f1,A_LO_f1] = ellip(3,0.02,60,([128 138])/(8600/2),'bandpass');
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
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f1,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F1, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB2(1:1024).*window2)))))
hold off
%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_LO_f1,P_LO_f1,K_LO_f1] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f1,Ae_LO_f1] = zp2tf(Z_LO_f1,P_LO_f1,K_LO_f1);
impulse(Be_LO_f1,Ae_LO_f1)
title('Impulse Response: Bandpass Filter, LO_F1')
grid on

%filter coefficients
[B_LO_f2,A_LO_f2] = ellip(3,0.02,60,([219 229])/(8600/2),'bandpass');
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f2 = filter(B_LO_f2,A_LO_f2,zero_vector);
%plot
figure(26)
subplot(2,1,1)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f2,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F2, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB2(1:1024).*window2)))))
hold off

%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_LO_f2,P_LO_f2,K_LO_f2] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f2,Ae_LO_f2] = zp2tf(Z_LO_f2,P_LO_f2,K_LO_f2);
impulse(Be_LO_f2,Ae_LO_f2)
title('Impulse Response: Bandpass Filter, LO_F2')
grid on

%filter coefficients
[B_LO_f3,A_LO_f3] = ellip(3,0.02,60,([301 311])/(8600/2),'bandpass');
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f3 = filter(B_LO_f3,A_LO_f3,zero_vector);
%plot
figure(27)
subplot(2,1,1)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f3,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F3, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB2(1:1024).*window2)))))
hold off 

%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_LO_f3,P_LO_f3,K_LO_f3] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f3,Ae_LO_f3] = zp2tf(Z_LO_f3,P_LO_f3,K_LO_f3);
impulse(Be_LO_f3,Ae_LO_f3)
title('Impulse Response: Bandpass Filter, LO_F3')
hold off
grid on

%filter coefficients
[B_LO_f4,A_LO_f4] = ellip(3,0.02,60,([373 383])/(8600/2),'bandpass');
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_LO_f4 = filter(B_LO_f4,A_LO_f4,zero_vector);
%plot
figure(28)
subplot(2,1,1)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_LO_f4,1024)))))
grid on
title('Frequency Response: Bandpass Filter, LO_F4, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPLB2(1:1024).*window2)))))
hold off

%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_LO_f4,P_LO_f4,K_LO_f4] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_LO_f4,Ae_LO_f4] = zp2tf(Z_LO_f4,P_LO_f4,K_LO_f4);
impulse(Be_LO_f4,Ae_LO_f4)
title('Impulse Response: Bandpass Filter, LO_F4')
hold off
grid on

%%%%%%%%%%%%%%%%%%%%%    bandpass highband set    %%%%%%%%%%%%%%%%%%%%%%%%%

%filter coefficients
[B_HI_f1,A_HI_f1] = ellip(3,0.02,60,([128 138])/(8600/2),'bandpass');
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f1 = filter(B_HI_f1,A_HI_f1,zero_vector);
%plot
figure(29)
subplot(2,1,1)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f1,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F1, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB2(1:1024).*window2)))))
hold off

%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_HI_f1,P_HI_f1,K_HI_f1] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f1,Ae_HI_f1] = zp2tf(Z_HI_f1,P_HI_f1,K_HI_f1);
impulse(Be_HI_f1,Ae_HI_f1)
title('Impulse Response: Bandpass Filter, HI_F1')
hold off
grid on

%filter coefficients
[B_HI_f2,A_HI_f2] = ellip(3,0.02,60,([255 265])/(8600/2),'bandpass');
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f2 = filter(B_HI_f2,A_HI_f2,zero_vector);
%plot
figure(30)
subplot(2,1,1)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f2,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F2, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB2(1:1024).*window2)))))
hold off

%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_HI_f2,P_HI_f2,K_HI_f2] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f2,Ae_HI_f2] = zp2tf(Z_HI_f2,P_HI_f2,K_HI_f2);
impulse(Be_HI_f2,Ae_HI_f2)
title('Impulse Response: Bandpass Filter, HI_F2')
hold off
grid on

%filter coefficients
[B_HI_f3,A_HI_f3] = ellip(3,0.02,60,([397 407])/(8600/2),'bandpass');
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f3 = filter(B_HI_f3,A_HI_f3,zero_vector);
%plot
figure(31)
subplot(2,1,1)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f3,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F3, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB2(1:1024).*window2)))))
hold off

%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_HI_f3,P_HI_f3,K_HI_f3] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f3,Ae_HI_f3] = zp2tf(Z_HI_f3,P_HI_f3,K_HI_f3);
impulse(Be_HI_f3,Ae_HI_f3)
title('Impulse Response: Bandpass Filter, HI_F3')
hold off
grid on

%filter coefficients
[B_HI_f4,A_HI_f4] = ellip(3,0.02,60,([511 521])/(8600/2),'bandpass');
%zero vector
zero_vector = zeros(1,400);
%set address 1,1 to 1
zero_vector(1) = 1;
%load the filter coefficients and the zero vector into an array, h
h_HI_f4 = filter(B_HI_f4,A_HI_f4,zero_vector);
%plot
figure(32)
subplot(2,1,1)
%frequency response
plot((-0.5:(1/1024):(0.5 - 1/1024))*8600, ...
    fftshift(20*log10(abs(fft(h_HI_f4,1024)))))
grid on
title('Frequency Response: Bandpass Filter, HI_F4, fs = 8600Hz')
xlabel('Frequency in Hertz')
ylabel('Magnitude in decibels')
hold on
%windowed fft of downsampled bandpass lowband output
plot((-0.5:(1/1024):0.5-(1/1024))*1075,...
    fftshift(20*log10(abs(fft(downsampled_output_BPHB2(1:1024).*window2)))))
hold off

%plot
%time response
subplot(2,1,2)
%zeros, poles, and gain
[Z_HI_f4,P_HI_f4,K_HI_f4] = ellip(n_BPHB, Rp_BPHB, Rs_BPHB,...
    (1200/(4300/2)),'s');
%construct transfer function bandpass highband
[Be_HI_f4,Ae_HI_f4] = zp2tf(Z_HI_f4,P_HI_f4,K_HI_f4);
impulse(Be_HI_f4,Ae_HI_f4)
title('Impulse Response: Bandpass Filter, HI_F4')
hold off
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end of EE556  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%thank you

