pkg load signal;
clear;

[audio_data] = audioread("audio.ogg", "ogg");

frequencia_amostragem = 44100;

fs = frequencia_amostragem;  % Frequência de amostragem
Omega_p = 200	 * (2 * pi)  ;  % Frequência de corte da banda de passagem (em rad/s)
Omega_s=  400	 * (2 * pi);  % Frequência de corte da banda de rejeição (em rad/s)

ep = 0.05; % Atenuação máxima permitida na banda de passagem
es = 0.05; % Atenuação mínima desejada na banda de rejeição

% Cálculo da ordem do filtro

ordem = ceil(((log10(((1-(1-ep)^2)/(1-ep)^2)/(((1-es^2)/(es^2)))) / log10(Omega_p/Omega_s))/ 2));
wc = (Omega_s/(((1-es^2)/es^2)^(1/(2*ordem))));
wn = wc / fs / 2;
% Projetando o filtro Butterworth
[b, a] = butter(ordem, wn, "high");

audio_filtrado = zeros(size(audio_data));  % Inicializa o sinal filtrado com zeros

% Aplicação do filtro
for n = 1:length(audio_data)
    audio_filtrado(n) = b(1) * audio_data(n);
    for k = 2:length(b)
        if n >= k
            audio_filtrado(n) = audio_filtrado(n) + b(k) * audio_data(n-k+1) - a(k) * audio_filtrado(n-k+1);
        end
    end
end

% Reproduzindo o áudio original
disp('Reproduzindo áudio original...');
%sound(audio_data, frequencia_amostragem);

% Plotando o gráfico de amplitude x tempo do áudio original
figure;
subplot(2, 1, 1);
tempo = (0:length(audio_data)-1) / frequencia_amostragem;
plot(tempo, audio_data);
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Amplitude x Tempo - Original');
grid on;

% Reproduzindo o áudio filtrado
disp('Reproduzindo áudio filtrado...');
sound(audio_filtrado, frequencia_amostragem);

% Plotando o gráfico de amplitude x tempo do áudio filtrado
subplot(2, 1, 2);
plot(tempo, audio_filtrado);
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Amplitude x Tempo - Filtrado');
grid on;

% Ajustando a escala dos eixos do gráfico do áudio filtrado
ylim([-1, 1]);  % Ajuste de acordo com os limites desejados

% Ajustando a escala dos eixos do gráfico do áudio original
subplot(2, 1, 1);
ylim([-1, 1]);  % Ajuste de acordo com os limites desejados
% Plotando o gráfico de frequência do áudio original
figure;
subplot(2, 1, 1);
N = length(audio_data);
f = (0:N-1) * (frequencia_amostragem / N);
amplitude = abs(fft(audio_data));
% Limitando a faixa de frequências
limite_inferior = 0;
limite_superior = 1000;
indice_limite_inferior = round(limite_inferior * N / frequencia_amostragem) + 1;
indice_limite_superior = round(limite_superior * N / frequencia_amostragem) + 1;
plot(f(indice_limite_inferior:indice_limite_superior), amplitude(indice_limite_inferior:indice_limite_superior));
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro de Frequência - Original');
grid on;
ylim([0, 10000]);  % Ajustando o eixo da amplitude

% Plotando o gráfico de frequência do áudio filtrado
subplot(2, 1, 2);
amplitude_filtrado = abs(fft(audio_filtrado));
% Limitando a faixa de frequências
plot(f(indice_limite_inferior:indice_limite_superior), amplitude_filtrado(indice_limite_inferior:indice_limite_superior));
xlabel('Frequência (Hz)');
ylabel('Amplitude');
title('Espectro de Frequência - Filtrado');
grid on;
ylim([0, 10000]);  % Ajustando o eixo da amplitude

