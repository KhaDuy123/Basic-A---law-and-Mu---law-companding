clc;
clear;
% 1. Load speech signal
Fs = 4000;
[mSpeech,Fs] = audioread("MaleSpeech-16-4-mono-20secs.wav"); %vốn dĩ từ đây đã được lượng tử hóa rồi nên ở đây ta coi như một tín hiệu analog giả
% sound(mSpeech,Fs)
% Consider the speech signal in 1.5s
t = 0:1/Fs:1.5;
plot(t,mSpeech(1:length(t)),'LineWidth',2);
hold on
% 2. Quantize the sample signal
L = 16; %the number of quantization levels
V_p = 0.5625; %the peak voltage of signal
% Determine the single quantile interval ?-wide
q = 2*V_p/(L-1); % Use the exact equation %0.0375/2 = 0.01875

s_q_2 = quan_uni(mSpeech(1:length(t)),q); % Uniform quantization
% Plot the sample signal and the quantization signal
plot(t,s_q_2,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r');

% 3. Calculate the average quantization noise power,...
% the average power of the sample signal and SNR
e_uni = mSpeech(1:length(t)) - s_q_2; % error between sample signal and quantized signal
pow_noise_uni = 0;
pow_sig = 0;
for i = 1:length(t)
    pow_noise_uni = pow_noise_uni + e_uni(i)^2;
    pow_sig = pow_sig + mSpeech(i)^2;
end

pow_noise_uni = pow_noise_uni / 12; %cai minh them vo
pow_sig = pow_sig/length(t); % cái mình thêm vô
SNR_a_uni =   pow_sig/pow_noise_uni ;



%--------compression-------------
% 5. Compress the sample signal ‘mSpeech’
mu = 255;  %or A = 87.6; use the standard value %cai minh them vo
A = 87.6;
y_max = V_p;
x_max = V_p;

s_c_5 = [];

%s_c_5 = compand(mSpeech(1:length(t),1),255,max(mSpeech(1:length(t),1)),'a/compressor');
%Doi voi mu

for i = 1:length(t)
   if (mSpeech (i,1) >= 0) 
       s_c_5(1,i) = y_max* (log(1 + mu * abs(mSpeech(i,1))/x_max )) / (log(1+mu));
   else
       s_c_5(1,i) = (-1)*y_max*(log(1 + mu * abs(mSpeech(i,1))/x_max )) / (log(1+mu));
   end
end
%}
%Doi voi A
c = 0;
%{
for i = 1:length(t)
    if ( abs(mSpeech(i,1))/x_max >= 0 &&  abs(mSpeech(i,1))/x_max <= 1/A )
        s_c_5(1,i) =  sign(mSpeech(i,1)) *y_max * (A * abs(mSpeech(i,1) / x_max)) / (1 + log(A));
    elseif ( abs(mSpeech(i,1))/x_max > 1/A &&  abs(mSpeech(i,1))/x_max <= 1)
        s_c_5(1,i) = sign(mSpeech(i,1)) *y_max* (1 + log(A * abs(mSpeech(i,1)) /x_max ) ) / (1 + log(A));
    end
end
%}
plot(t,s_c_5);
s_q_6 = quan_uni(s_c_5,q);
plot(t,s_q_6,'b^','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
%s_e_7 = compand(s_c_5,255,max(mSpeech(1:length(t),1)),'mu/expander');
s_e_7 = sign(s_q_6(1:length(t))) .* x_max.*(1/mu) .* ( (1+mu).^(abs(s_q_6(1:length(t))) / y_max) - 1 ); %for mu law

%{
for i = 1:length(t)
    if  ( abs(s_q_6(1,i)) / y_max >= 0 && abs(s_q_6(1,i)) / y_max < (1 / (1+log(A)) ))
        s_e_7(1,i) = x_max * sign(s_q_6(1,i)) * ( abs(s_q_6(1,i) / y_max * (1 + log(A)) ) ) / A;
    elseif ( abs(s_q_6(1,i)) / y_max < 1 && abs(s_q_6(1,i)) / y_max >= (1 / (1+log(A)) ))
        s_e_7(1,i) = x_max * sign(s_q_6(1,i)) * ( exp(abs(s_q_6(1,i)) / y_max * (1 + log(A)) -1) ) / A;
    end
end
%}
plot(t,s_e_7,'g*','MarkerSize',6,'MarkerEdgeColor','g','MarkerFaceColor','g');
legend('Sample signal','Uniform quantized values','Compresssignal','Compress quantized values','Nouniform quantized values');
xlim([0.53 0.58])
%ylim([-0.5625 0.5625])


e_com = mSpeech(1:length(t)) - s_e_7; % error between sample signal and quantized signal
pow_noise_com = 0;
for i = 1:length(t)
    pow_noise_com = pow_noise_com + e_com(i)^2;
end

pow_noise_com = pow_noise_com / 12;
SNR_a_com = pow_sig/pow_noise_com;







