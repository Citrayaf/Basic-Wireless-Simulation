clear all
clc

SNR=0:2:30;
numerik_snr=10.^(SNR/10);

for i=1:length(SNR)
    fprintf('SNR: %d\n', SNR(i));
    BER_bpsk_awgn(i)=0.5*erfc(sqrt(numerik_snr(i)));
end

semilogy(SNR,BER_bpsk_awgn,'-r*','linewidth',1)
