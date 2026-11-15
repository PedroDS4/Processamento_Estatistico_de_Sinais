<div align="center">

# Relatório Técnico
### Projeto: Detecção de Alvo em Clutter via Modelagem AR e Filtro Branqueador

**Autor:** Pedro Arthur Oliveira dos Santos  
**Orientador:** Prof. Luiz Felipe de Queiroz Silveira  
**Instituição:** Universidade Federal do Rio Grande do Norte (UFRN) / Pós-Graduação em Engenharia Elétrica e Computação  

</div>

---

## Introdução

O objetivo deste relatório é descrever as etapas para a realização do projeto de detecção de alvo em um sinal de radar contaminado, utilizando os conceitos de filtragem de processos estocásticos.

Nas próximas seções serão descritas detalhadamente cada uma das etapas até o objetivo final do projeto.

---

## Modelagem AR e Seleção da Ordem

### Modelagem do Processo

O sinal recebido no radar é dado por:

$$x[n] = s[n] + c[n] + v[n]$$

onde o sinal $s[n] = A_S \cos(\omega_0 n + \phi)$ é um processo harmônico, $\phi$ é uma variável aleatória uniforme no intervalo $[-\pi, \pi]$, $v[n]$ é um ruído branco de média zero, e $c[n]$ é um processo de ruído branco de média zero filtrado por um sistema todo polo de ordem $4$, ou $\text{AR}(4)$.

Dessa forma, o clutter no domínio $Z$ é dado por:

$$C(z) = H(z) V_c(z)$$

Onde $H(z)$ é um sistema todo-polo. A densidade espectral de potência do clutter pode ser escrita como:

$$P_c(e^{j\omega}) = |H(e^{j\omega})|^2 \sigma_c^2 = |H(e^{j\omega})|^2$$

Como $\sigma_c^2 = 1$, podemos então estimar a potência do sinal conhecendo apenas os polos do sistema. Assim, um sistema de quarta ordem com os polos:

$$\begin{cases}
p_1 = 0.85e^{j0.42\pi} \\
p_2 = 0.85e^{-j0.42\pi} \\
p_3 = 0.92e^{j0.15\pi} \\
p_4 = 0.92e^{-j0.15\pi}
\end{cases}$$

A potência teórica do processo harmônico é dada por:

$$P_s(\omega)|_{\omega = \omega_0} = \frac{\pi |A|^2}{2} (\delta(\omega_0 - \omega_0) + \delta(\omega_0 + \omega_0)) = \frac{\pi |A|^2}{2}$$

A potência do processo harmônico é determinada a partir da potência do clutter na frequência dominante como:

$$P_s(\omega_0) = \frac{P_c(\omega_0)}{150}$$

Assim temos a seguinte equação:

$$\frac{\pi |A|^2}{2} = \frac{|H(e^{j\omega})|^2}{150}$$

sendo a resposta ao impulso no domínio $Z$ do filtro $\text{AR}(4)$ dada por:

$$H(z) = \frac{1}{(1 - p_1z^{-1})(1 - p_1^*z^{-1})(1 - p_2z^{-1})(1 - p_2^*z^{-1})}$$

Após substituir os valores dos polos, a amplitude de sinal procurada é dada por:

$$A = \sqrt{\frac{2|H(e^{j\omega})|^2}{150\pi}} = 0.0193$$

Então, a amplitude total do sinal do alvo é menor que $2\%$ da variância do clutter, sendo o alvo completamente mascarado no domínio do tempo.

### Determinação da Ordem

Para determinar a ordem do modelo AR, assumindo um cenário real onde não a conhecemos, as equações de Yule-Walker serão resolvidas para vários valores de $p$ no intervalo $[1, p_{max}]$, e o valor da ordem escolhido ficará no cotovelo da curva do erro de predição em função de $p$.

Resolveu-se então a equação de Yule-Walker:

$$r_x[k] + \sum_{l = 1}^{p} a_p(l)r_x[k-l] = 0$$

Para os valores de $p$ no intervalo $[0, p+1]$.

```matlab
% Resolvendo YW para cada valor de p ate p_max
for p = 1:p_max
    r_x_p = r_x(2:p+1);
    R_x_p = toeplitz(r_x(1:p));
    a_p = R_x_p \ (-r_x_p);
    b_0 = sqrt([1; a_p]' * conj(r_x(1:p+1)));
    c_n_p = filter(b_0, [1; a_p], c_n_pre);
    sigma_c(p) = b_0^2;
end
