clear all;
clear all figures;
%% 7.1 Samples at the output of a detector: Gaussian white noise
%vector of noise
%--> i n p u t s a j u s t a b l e s
N_samples=10001;
rep=100;
PRF=1*10^3 ;%kHz
pot_noise=1; %1W
mean_=0;
%--> e x e c u c i ó
sample_factor_noise=randn(1,N_samples)+mean_;
Pot_n=sum(abs(sample_factor_noise).^2)/N_samples;
noise=sqrt(pot_noise).*sample_factor_noise./sqrt(Pot_n);
figure(1);
plot(abs(fft(noise)));
figure(2);
histogram(noise,'Normalization','pdf');
xlabel('Noise');
ylabel('PDF of noise');
%matrix of noise
%--> i n p u t s a j u s t a b l e s
N_samples=10001;
PRF=1*10^3 ;%kHz
pot_noise=1; %1W
rep=100;
mean_=0;
%--> e x e c u c i ó
sample_noise=randn(rep,N_samples)+mean_;
for i=1:rep
    Pot=sum(abs(sample_noise(i,:)).^2)/N_samples;
    noise_2(i,:)=sqrt(pot_noise)*sample_noise(i,:)./sqrt(Pot);
end
noise_fft=abs(fft(noise_2,[],2));
noise_mean=mean(noise_fft(:,(1:N_samples)));
figure(3);
freq=(0:N_samples-1)/N_samples*PRF;
plot(noise_mean);
figure(4);
histogram(noise_mean,'Normalization','pdf');
xlabel('Noise');
ylabel('PDF of noise');