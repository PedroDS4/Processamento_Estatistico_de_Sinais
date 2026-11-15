%% Projeto: Detecção de Alvo em Clutter via Modelagem AR e Filtro Branqueador
% Aluno: Pedro Arthur Oliveira dos Santos


fontesize = 15;

%% Sistema AR

poles = [0.92*exp(-1i*0.15*pi), 0.92*exp(1i*0.15*pi), 0.85*exp(-1i*0.42*pi), 0.85*exp(1i*0.42*pi)];


D_z = poly(poles);
N_z = 1;





%% Gerando o Sinal x(n) = s(n) + c(n) + v(n)

% s(n): Senoide de fase aleatória

% A_s = 0.1;
w_true = 0.4*pi;

A_s = estimar_amplitude(w_true, poles);

% Eixo de tempo discreto
N = 16384;
n = (0:N-1);

% Variável Aleatória de Fase
phi_s = -pi + (2*pi)*rand(1);
s_n = A_s * cos(w_true*n + phi_s);





% c(n): Clutter -> Processo WSS de banda estreita
sigma_c = 1;

c_n_pre = sqrt(sigma_c)*randn(1, N);


c_n = filter(N_z, D_z, c_n_pre);



figure()
plot(c_n_pre/max(abs(c_n_pre)), 'Color', 'red', 'linewidth', 4); hold on;
plot(c_n/max(abs(c_n)), 'Color', 'blue', 'linewidth', 2);
xlabel('[n]');
ylabel('Amplitude');
legend('Ruído Branco', 'Ruído Branco de Banda Estreita');
set(gca, 'fontsize', fontesize);
grid on

% v(n): Ruído Branco
sigma_v = 0.3;

v_n = sqrt(sigma_v)*randn(1, N);



% x(n): s(n) + c(n) + v(n)

x_n = s_n + c_n + v_n;

% x_n = s_n + v_n;

%% Estimação do Periodograma bruto
Nfft = 1024;

Pxx = estimar_periodograma(x_n, Nfft);

dt = 1e-3;

df = 1/(Nfft*dt);
f = (0:Nfft-1)*df;

omega = 2*pi*f*dt;


figure()
plot(omega, abs(Pxx), 'color', 'black', 'linewidth', 3);
xlabel('w (rad/amostra)');
ylabel('Magnitude');
set(gca, 'fontsize', fontesize);
title('Pxx(w) Densidade Espectral de Potência Estimada')
grid on



%% Estimando a ordem do Filtro Branqueador



% Ordem máxima a ser procurada
p_max = 10;
sigma_c = zeros(p_max, 1);

r_x = estimar_r_x(x_n, p_max);

figure()
plot(r_x ,'color', 'blue', 'LineWidth', 3);
xlabel('k');
ylabel('Sequência de Autocorrelação de x[n]');
set(gca, 'fontsize', fontesize);
grid on;


% Resolvendo YW para cada valor de p até p_max
for p = 1:p_max
    
    r_x_p = r_x(2:p+1);

    R_x_p = toeplitz(r_x(1:p));
    
    a_p = R_x_p\(-r_x_p);
    
    b_0 = sqrt([1; a_p]' * conj(r_x(1:p+1)));

    c_n_p = filter(b_0, [1; a_p], c_n_pre);

    sigma_c(p) =  b_0^2;
end




p_vector = (1:p_max);

figure();
plot(p_vector, sigma_c, 'Color', 'red', 'LineWidth', 4);
xlabel('$p$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\sigma_v$', 'Interpreter', 'latex', 'FontSize', 12);
title('Vari\^ancia do Erro de Predi\c{c}\~ao', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 1.5*fontesize); 
grid on;


% Critério do Cotovelo (queda relativa < 5%)
limiar = 0.0165;
% limiar = 0.001;
queda_relativa = [1; -diff(sigma_c) ./ sigma_c(1:end-1)];

ordem_escolhida = find(queda_relativa < limiar, 1, 'first') - 1;
if isempty(ordem_escolhida) || ordem_escolhida < 1
    ordem_escolhida = p_max;
end
fprintf('Ordem escolhida pelo critério do cotovelo: p = %d\n', ordem_escolhida);


p = ordem_escolhida;

%% Implementação do Filtro Branqueador
% Com a ordem p escolhida, estimar o filtro branqueador para reverter o
% clutter

% [~, p] = min(abs(diff(sigma_c)./sigma_c(1:p_max-1) - 0.01)); % ordem escolhida na Tarefa anterior a partir do cotovelo da curva
% p  = 3;
[x_n_w, a_p_v] = estimar_AR_e_filtrar(r_x, x_n, p);



%% Densidade Espectral de Potência de x(n) filtrado
Pxx_w = estimar_periodograma(x_n_w, Nfft);





figure()
plot(omega, abs(Pxx), 'color', 'black', 'linewidth', 5); hold on;
plot(omega, abs(Pxx_w), 'color', 'green', 'linewidth', 3);
xline(w_true, '--r' );
xline(2*pi - w_true, '--r' );
xlabel('w (rad/amostra)');
ylabel('Magnitude');
set(gca, 'fontsize', fontesize);
legend('Sinal Original', 'Sinal Filtrado')
title('Densidade Espectral de Potência Pxx(w) ')
grid on



%%  Estimando a frequência do alvo.


[~, iw_est] = max(abs(Pxx_w));
w_est = omega(iw_est);


[~, iw_wrong] = max(abs(Pxx));
w_wrong = omega(iw_wrong);



erro = (w_est - w_true)/w_true;
fprintf('Frequência Verdadeira %.4f versus Estimada %.4f \n', w_true, w_est);
fprintf('O erro relativo entre a frequência original e a estimada após o filtro é de %.3f porcento \n', erro*100);



%% Verificação numérica da invariância do CFAR local


cfar_x_n = calcular_cfar(Pxx, omega, w_true);

cfar_x_n_w = calcular_cfar(Pxx_w, omega, w_true);


fprintf('Razão cfar antes da filtragem cfar = %.2f dB \n', cfar_x_n);
fprintf('Razão cfar após a filtragem cfar = %.2f dB \n', cfar_x_n_w);


%% Quantificando o ganho de detecção com um limiar global


limiar = 0.9;


n_before = length(Pxx(abs(Pxx/max(abs(Pxx)))>limiar));
n_after = length(Pxx_w(abs(Pxx_w/max(abs(Pxx_w)))>limiar));

perc_before = 100*n_before/Nfft;
perc_after = 100*n_after/Nfft;

fprintf('Porcentagem de frequencias acima do limiar antes da filtragem = %.2f \n', perc_before);
fprintf('Porcentagem de frequencias acima do limiar após a filtragem = %.2f \n', perc_after);




%% Simulações de Monte-Carlo para Validação


MC = 300;

disp('Iniciando as Simulações de Monte Carlo')
erro_w = zeros(MC, 1);

for m = 1:MC

    % A_s = 0.5;
    % w_true = 0.4*pi;

    % Eixo de tempo discreto
    N = 16384;
    n = (0:N-1);

    % Variável Aleatória de Fase
    phi_s = -pi + (2*pi)*rand(1);
    % phi_s
    s_n = A_s * cos(w_true*n + phi_s);

    sigma_c = 1;

    c_n_pre = sqrt(sigma_c)*randn(1, N);


    c_n = filter(N_z, D_z, c_n_pre);


    % v(n): Ruído Branco
    sigma_v = 0.3;

    v_n = sqrt(sigma_v)*randn(1, N);


    % x(n): s(n) + c(n) + v(n)

    x_n = s_n + c_n + v_n;


    Nfft = 1024;

    Pxx = estimar_periodograma(x_n, Nfft);


    % r_x estimado pela covariancia
    r_x = estimar_r_x(x_n, p);


    x_n_w = estimar_AR_e_filtrar(r_x, x_n, p);


    Pxx_w = estimar_periodograma(x_n_w, Nfft);


    [~, iw_est] = max(abs(Pxx_w));
    w_est = omega(iw_est);

    erro_w(m) = (w_est - w_true)/w_true;



    w_sorted = sort(omega);

    [~ , bin_alvo ] = min (abs ( w_sorted - w_true ) ) ; % indice do alvo na grade ordenada

    guarda = 3; janela = 40; % celulas de guarda e de referencia ( ajuste se quiser )

    idx_janela = mod (( bin_alvo - janela : bin_alvo + janela ) - 1 , Nfft ) + 1;

    idx_guarda = mod (( bin_alvo - guarda : bin_alvo + guarda ) - 1 , Nfft ) + 1;

    idx_ref = setdiff ( idx_janela , idx_guarda ) ; % referencia =  janela MENOS as celulas de guarda


    med_x_n = median(Pxx(idx_ref));

    cfar_x_n(m) = Pxx(bin_alvo)/med_x_n;

    med_x_n_w = median(Pxx_w(idx_ref));

    cfar_x_n_w(m) = Pxx_w(bin_alvo)/med_x_n_w;

end


erro_w_med = mean(erro_w);

fprintf('Erro médio de estimação de frequencia  = %.5f', erro_w_med);


figure()
hist(erro_w, 20)
xlabel('Erro de Estimação da Frequência do Alvo');
ylabel('Frequencia');
set(gca, 'fontsize', fontesize);
title('Histograma');
grid on;


figure()
plot(erro_w, 'color', 'red', 'linewidth', 3);
xlabel('Simulation');
ylabel('Frequency Estimation Error');
set(gca, 'fontsize', fontesize);
grid on;



figure()
plot(cfar_x_n, 'color', 'black', 'linewidth', 4); hold on;
plot(cfar_x_n_w, 'color', 'green', 'linewidth', 3);
xlabel('Simulação');
ylabel('CFAR (dB)');
title('Razão CFAR')
legend('Original', 'Pós Filtragem');
set(gca, 'fontsize', fontesize);
grid on;








%% Funções

% Estimação do Periodograma
function Pxx = estimar_periodograma(x_n, Nfft)

    N = length(x_n);
    
    win = hamming (Nfft)'; % qualquer janela reduz vazamento espectral
    
    winnorm = sum( win .^2) / Nfft; % normalizacao de potencia da janela
    
    nblocks = floor ( N / Nfft );
    
    Pxx = zeros (1 , Nfft);
    
    for b = 1: nblocks
        seg = x_n (( b -1) * Nfft +1 : b * Nfft);
    
        Xb = fft ( seg .* win ) ;
    
        Pxx = Pxx + abs( Xb ) .^2/( Nfft * winnorm ) ;
    end
    
    Pxx = Pxx / nblocks ; % periodograma de Welch ( media entre blocos )

end



% CFAR
function cfar = calcular_cfar(Pxx, omega, omega_true)

    Nfft = length(omega);
    
    omega_sorted = sort(omega);
    
    [~ , bin_alvo ] = min (abs ( omega_sorted - omega_true ) ) ; % indice do alvo na grade ordenada
    
    guarda = 3; janela = 40; % celulas de guarda e de referencia ( ajuste se quiser )
    
    idx_janela = mod (( bin_alvo - janela : bin_alvo + janela ) - 1 , Nfft ) + 1;
    
    idx_guarda = mod (( bin_alvo - guarda : bin_alvo + guarda ) - 1 , Nfft ) + 1;
    
    idx_ref = setdiff ( idx_janela , idx_guarda ) ; % referencia =  janela MENOS as celulas de guarda
    
    
    med_x_n = median(Pxx(idx_ref));
    
    cfar = 20*log10(Pxx(bin_alvo)/med_x_n);
end



function [x_n_w, a_hat_full] = estimar_AR_e_filtrar(r_x, x_n, p)

    
    kidx = (p+1):(2*p); % defasagens k = p+1 ,... ,2p -- nunca k = 0;
    
    Mmod = zeros(p);
    rhsmod = zeros(p, 1);
    
    
    for i = 1:p
        kk = kidx(i);
        for l = 1:p
            Mmod (i , l) = r_x(kk - l + 1); % <-- r_x na defasagem correspondente (ver equacao acima )
        end
        rhsmod (i) = -r_x(kk + 1) ;
    end
    
    
    a_hat = ( Mmod \ rhsmod )';
    
    a_hat_full = [1 , a_hat ] % A_p (z) = 1 + a_hat (1) z^ -1 + ... + a_hat (p)z^ -p;
    
    b_0 = sqrt(a_hat_full(1:p) * conj(rhsmod));

    x_n_w = filter ( a_hat_full , 1 , x_n); % sinal branqueado

end



function r_x = estimar_r_x(x_n, p)

    % r_x estimado pela covariancia
    r_x = zeros(2*p + 1, 1);
    mu_x = mean(x_n);
    N = length(x_n);
    for k = 0:2*p
        x_k = circshift(x_n, k);
        c_0 = cov(x_n);
        cov_x = (sum((x_n-mu_x) .* (x_k-mu_x)))/N;
        r_x(k+1) = cov_x;
    end

end



function A = estimar_amplitude(w_alvo, poles)
    
    z = exp(1i*w_alvo);
    p1 = (1 - z*(poles(1) + poles(2)) + abs(poles(1)*z)^2 );
    p2 = (1 - z*(poles(3) + poles(4)) + abs(poles(4)*z)^2 );
    
    D_w = p1*p2;

    H_w_2 = abs(1/D_w)^2;

    A = sqrt(2 * H_w_2/(150*pi));

end



