clear
clc
echo off;
format long;

%initial value/////////////////////////////////////
clk=6.144*10^6;         %CLK = 6.144[MHz]�@�����g��
Ts=1/clk;               %Ts = 162.76[ns]
N = 65536;             %FFT sample number

%Input Amplitude//////////////////////////////////////////////
Ampl=2^27;             % AudioData Input 

%Input_1//////////////////////////////////////////////////////
%Amp1=Ampl;              % -0dBFS
%Amp1=Ampl/2;             % -6dBFS
%Amp1=Ampl/10;           %-20dBFS
%Amp1=Ampl/100;          %-40dBFS
%Amp1=Ampl/1000;         %-60dBFS
%Amp1=Ampl/10000;        %-80dBFS 
%Amp1=Ampl/100000;       %-100dBFS
Amp1=Ampl/1000000;      %-120dBFS


osr=128;                 %over sampling rate
suso = 20;
Offset = 300;
FSC=2^27;                % Quantizer FSC   2~27=134217728
M=8;                     %number of speaker

%input frequency////////////////////////////////////
simtime = (N+Offset)*Ts;          %simulation time
fB = ceil((N)/(2*osr));             %
f0 = floor(107/256*fB);           %f=10[kHz]
finrad = 2*pi/(N)*f0*clk;                   


%delsig constant//////////////////////////////////////////
g1=0;        
a1=1;
a2=3;
a3=3;
b1=1;
c1=1;
c2=1;
c3=1;


%NSDEM constant//////////////////////////////////////////
ns_fb=g1;
SAT=127;             
A=1;
B2=8;
C=32;
c11=1;
c22=1;
c33=1;


%Error /////////////////////////////////////////////
%�X�s�[�J�̃~�X�}�b�`��z��
err = [0.972287782947746 0.967359789827757 0.961453278144533 0.982214763022800 0.933245475470565 1.03214098075439 1.01733533671518 0.999332404945249];  %    10%


%simlink model//////////////////////////////////////
%���s����simlink���f�����w��
sim('AVIC_model');


%nsdem//////////////////////////////////////////////
%simlink���f�����烏�[�N�X�y�[�X�ɃG�N�X�|�[�g�����f�[�^����������
%OUTp�AOUTn�͂��ꂼ��X�s�[�J���쓮����M��
%1�̃X�s�[�J��OUTp�AOUTn��2�̃f�W�^���M���ŋ쓮�����
for n = 1:M
 out1(:,n)=OUTp(:,n)-OUTn(:,n);
 out1(:,n)=out1(:,n)*err(n);
end


%�~�X�}�b�`���l�������e�X�s�[�J�쓮�M�������Z���A�t���X�P�[����1�ɂȂ�悤�ɒ���
out8_1=out1(:,1)+out1(:,2)+out1(:,3)+out1(:,4)+out1(:,5)+out1(:,6)+out1(:,7)+out1(:,8);
out8_1=(out8_1)/8;


window = hanning(N);
sim_n_ns = out8_1(Offset:N+Offset-1);
simfft_n_ns = fft(sim_n_ns .* hanning(N));
simfft_norm_ns = simfft_n_ns/(N/4);
simfftdb_norm_ns = 20*log10(abs(simfft_norm_ns));

fre = (0:length(simfftdb_norm_ns)-1)/length(simfftdb_norm_ns)*clk;



figure(1);
semilogx(fre,simfftdb_norm_ns,'b-');
axis([clk/100000 clk -200 20]);
ylabel('gain[dB]');
xlabel('frequency[Hz]');
grid on;
hold on;


%SNR�̎Z�o///////////////////////////////////////////////////////// AmpdB =
%20*log10(Amp);                                        %��͐U����dB�\��

hwfft = simfft_norm_ns(1:ceil(fB));

signalBins = [108];                       %hanning���ɂ��M����g����
signalBins = signalBins(signalBins>0);
signalBins = signalBins(signalBins<=length(hwfft));
s = norm(hwfft(signalBins)); %��͐M���̓d�͂��Z�o 

n1 = norm(simfft_norm_ns(5:107-suso,1));
n2 = norm(simfft_norm_ns(110+suso:213,1));

n = n1+n2;

if n==0
    snr = Inf;
else
    snr = 20*log10(s/n)                                    %SNR��\�L
end
