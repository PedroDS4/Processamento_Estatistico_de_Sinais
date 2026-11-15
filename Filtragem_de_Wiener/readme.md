##**Filtragem com Filtro de Wiener**
Considere o problema aonde temos um sinal contaminado com ruído aditivo da forma

$$
x[n] = d[n] + v[n]
$$

Queremos encontrar uma estimativa para o sinal $d[n]$, dessa forma é preciso conhecer algumas estatísticas desses sinais.

As equações de Wiener-Hopf são dadas por

$$
R_x w = r_{xd}
$$

também existe o formato onde se conhecem apenas as correlações, assumindo o ruído descorrelacionado com o sinal $d[n]$. Nesse caso dada por

$$
(R_d + R_v)w = r_{d}
$$




**Exemplo 1**:
Seja $d[n]$ um processo AR de orden 1, com $r_d[k] = \alpha^{|k|}$, suponha que $d[n]$ é observado na presença de ruído branco $v[n]$ com $\sigma_v^2$.

Vamos descobrir qual filtro $H(z)$ gera esse processo AR, utilizando as equações com as densidades espectrais de potência, temos

$$
P_d(z) = H(z) H^*(1/z^*) \sigma_n^2
$$

calculando a transformada de Fourier da autocorrelação e reescrevendo a função de autocorrelação de $d$, temos

$$
P_d(z) = TZ\{\ \alpha^{k}u[k] + \alpha^{-k}u[-k-1] \}\
$$

Dessa forma, temos

$$
P_d(z) = \frac{1}{1 - \alpha z^{-1}} - \frac{1}{1 - \frac{1}{\alpha} z^{-1}}
$$

Multiplicando o segundo termo da soma em cima e em baixo por $-\alpha z$:

$$
P_d(z) = \frac{1}{1 - \alpha z^{-1}} + \frac{\alpha z}{1 - \alpha z}
$$

distribuindo e simplificando, temos

$$
P_d(z) = \frac{1 - \alpha z + \alpha z - \alpha^2}{(1 - \alpha z^{-1})(1 - \alpha z)} = \frac{1 - \alpha^2}{(1 - \alpha z^{-1})(1 - \alpha z)}
$$

Dessa forma temos que

$$
\begin{cases}
H(z) = \frac{1}{1 - \alpha z^{-1}} \\
\sigma_n^2 = 1-\alpha^2
\end{cases}
$$

A autocorrelação do ruído será dada por

$$
r_v[k] = \sigma_v^2 \delta[k]
$$

Podemos então calcular agora as matrizes
$$
R_d = toeplitz(r_d[0:1])
$$

Dessa forma, temos as seguintes equações

$$
\begin{bmatrix}
1 + \sigma_v^2     & \alpha  \\
\alpha & 1 + \sigma_v^2
\end{bmatrix}
\begin{bmatrix}
w_0 \\
w_1
\end{bmatrix}
=
\begin{bmatrix}
1 \\
\alpha
\end{bmatrix}
$$

Solucionando isso, obtemos

$$
\begin{bmatrix}
w_0 \\
w_1
\end{bmatrix}
=
\begin{bmatrix}
\frac{1+\sigma_v^2 - \alpha^2}{(1+\sigma_v^2)^2 - \alpha^2} \\
\frac{\alpha \sigma_v^2}{(1+\sigma_v^2)^2 - \alpha^2}
\end{bmatrix}
$$

Abaixo essas equações serão resolvidas computacionalmente com essas fórmulas, e também depois utilizando estimativas da matriz de autocorrelação, e também pela matriz de toeplitz
