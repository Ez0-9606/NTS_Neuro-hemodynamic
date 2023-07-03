clc;
close all;
clearvars -except hrALL

data = hrALL(1:end-1,2);

Fs = 10;
dt = 5; %[s]
window = round(dt*Fs);
% g = ones(window,1);
g = hann(window,"periodic");    
noverlap = round(0.5*dt*Fs); % 50% overlap of window
nfft = max(256,2^ceil(log2(window)));
flimit = 5;
n = round(nfft/Fs*flimit);

figure;
[sp, fp, tp, psp] = spectrogram(data, g, noverlap, nfft, Fs, 'power');
imagesc(tp, fp(1:n), 10*log10(psp(1:n,:)));hold on;
view(0, 270);
xline(20, '-',{'Stim-On'}, 'Color', 'r', 'LineWidth', 2);
xline(80, '-',{'Stim-Off'}, 'Color', 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 1]);

figure; 
[s, f, t] = stft(data, Fs, Window =g, OverlapLength = noverlap, FFTLength = nfft, FrequencyRange = "onesided");
p = abs(s).^2/(Fs*norm(g)^2);
p(2:end-1) = 2*p(2:end-1);
imagesc(t', f(1:n), 10*log10(enbw(g, Fs)*p(1:n,:))); hold on;
view(0,270);
xline(20, '-',{'Stim-On'}, 'Color', 'r', 'LineWidth', 2);
xline(80, '-',{'Stim-Off'}, 'Color', 'r', 'LineWidth', 2);
% ylim([0 1]);

%%


%%
psd = psp(1:n,:);
lim = round(nfft/Fs*0.75);
lfpw = sum(psd(2:lim,:), 1);
hfpw = sum(psd(lim+1:end,:), 1);

figure;
subplot(2,1,1);
plot(t, lfpw); hold on;
plot(t, hfpw);
subplot(2,1,2);
plot(t, lfpw./hfpw);

return;
%%
figure;
plot(2*abs(s(:,1)).^2/(Fs*norm(g)^2)./psp(:,1));
title('Spectrogram S vs P: PSD');
%%
figure;
plot(periodogram(data(1:window), g, nfft, Fs), 'LineWidth', 3); hold on;
plot(psp(:,1));
title('Periodogram vs Spectrogram: PSD');
%%
figure;
plot(abs(fft(data(1:window).*g, nfft)), 'LineWidth', 3); hold on;
plot(abs(s(:,1)));
title('FFT vs Spectrogram: STFT');
%%
dft = fft(data(1:window).*g, nfft);
% 
% figure;
% plot(abs(dft), 'LineWidth',3); hold on;
% title('FFT vs Cus-DFT');

%%
psd = dft.*conj(dft)/norm(g)^2/Fs;
psd = psd(1:nfft/2+1);
psd(2:end-1) = 2*psd(2:end-1);

figure;
plot(psd, 'LineWidth', 3); hold on;
plot(psp(:,1));
title('FFT vs Spectrogram: PSD');

%%
% cusper = zeros(nfft/2+1, 1);
% for i = 1:nfft/2+1
%     cusper(i) = 1/Fs/sum(g.^2)*abs(sum([g.*data(1:window);zeros(nfft-window,1)].*exp(-1i*2*pi/nfft*(i-1)*(0:nfft-1)')))^2;
% end
% cusper(2:end-1) = 2*cusper(2:end-1);
% 
% figure;
% % plot(psd, 'LineWidth', 3); hold on;
% plot(cusper, 'LineWidth', 3); hold on;
% plot(periodogram(data(1:window), g, nfft, Fs));
% title('Cus-Periodogram vs Periodogram: PSD');

%%
figure;
[pxx,~] = pwelch(data(1:window), g, noverlap, nfft, Fs);
plot(pxx, 'LineWidth', 3); hold on;
plot(psp(:,1));
title('Welch vs Spectrogram: PSD');