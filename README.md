# levy-option-pricing
Matlab implementation of option pricing and calibration methods based on Levy stochastic processes.

The following is a chapter of my master thesis on bid-ask calibration of Levy models. It makes use of some of the provided code.

## Comparison of Pricing Algorithms

The pricing algorithms developed in this chapter all aim to allow for an
efficient calculation of multiple European call options written on the
same underlying. While the Fourier transform based algorithms gradually
improve the theoretical computational effort by advancing to the FFT and
FRFT algorithms respectively, the COS method utilises the rather fast
convergence property of the cosine series expansion. In this section we
examine the actual performance capabilities of our own MATLAB
implementation of the following four pricing algorithms:

  - pFT: naive Fourier transform pricing as in Section
    [\[subsec:numerical\_algorithms\]](#subsec:numerical_algorithms)

  - pFFT: Fourier pricing based on the FFT algorithm, see Section
    [\[subsec:numerical\_algorithms\]](#subsec:numerical_algorithms)

  - pFRFT: Fourier pricing based on the FRFT algorithm, see Section
    [\[subsec:numerical\_algorithms\]](#subsec:numerical_algorithms)

  - pCOS: COS pricing method as in Section
    [\[sec:cos\_method\]](#sec:cos_method)

For details on the general MATLAB implementation architecture refer to
Appendix [\[app:matlab\_implementation\]](#app:matlab_implementation).
In order to make a technically rapid calibration procedure possible in
the following
Chapter [\[chp:model\_calibration\]](#chp:model_calibration), we seek
the pricing algorithm yielding the lowest pricing error given a fixed
run-time. For this purpose we consider a Black-Scholes model with stock
price \(S_0 = 2\) and volatility parameter \(\sigma = 0.2\). Assume an
interest rate of \(r = 0.01\). Consider 100 European call options with
maturity \(T = 3\) and strike prices evenly spreading from \(1.8\) to
\(2.2\) written on this stock. For each pricing run with one of the
algorithms above and a certain parametrisation we measure the pricing
error over all 100 options as follows: Since we assumed a Black-Scholes
model, the analytical option prices may be obtained via the
Black-Scholes formula. The pricing error is then defined as the sum of
all absolute errors between the Black-Scholes price for an option and
the price calculated with one of the pricing algorithms. For each
pricing run of 100 options we measure the CPU time spent on the
calculation. To avoid influences of machine background tasks we average
all pricing errors and run-times over ten independent runs.
Figure [\[fig:pricing\_errors\_vs\_cputime\]](#fig:pricing_errors_vs_cputime)
depicts the results received by varying the number of summands
\(N = 2^n\) for \(2 \leq n \leq 15\) for each pricing algorithm and
storing the corresponding cpu time as well as the pricing error. The
remaining parameters for the Fourier transform based pricing algorithms
are \(\alpha = -4\) and \(\eta = 0.25\). The parameter \(\lambda\) was
chosen adaptively by the respective pricing algorithms.

![image](images/pricing_error_vs_cputime.pdf)
<span id="fig:pricing_errors_vs_cputime" label="fig:pricing_errors_vs_cputime">\[fig:pricing\_errors\_vs\_cputime\]</span>

First we notice that for both pCOS as well as pFT the pricing error
quickly approaches a very small constant level as \(N\) and thus the CPU
time increase. The methods pFFT and pFRFT however show a linear decrease
in the pricing error for increasing CPU time spent on the log-log scale.
This is especially due to pFT and pCOS (approximately) evaluating the
call option prices at the exact specified exercise prices whereas pFFt
and pFRFT evaluate prices along an equidistant log-strike grid. Thus the
resulting option prices of pFFT and pFRFT lie on a grid of exercise
prices with increasing mesh-size. The hence required final interpolation
of option prices at the specified strike grid adds another source of
pricing error besides the error occurring due to a finite \(N\). We
found that log-linear interpolation as in equation
[\[eq:log\_linear\_interpolation\]](#eq:log_linear_interpolation) yields
more stable and accurate prices than other interpolation techniques.

Secondly, to check whether the CPU time actually is reduced by using FFT
and FRFT based pricing algorithms instead of the ordinary pFT we
consider Figure [\[fig:cputime\_vs\_N\]](#fig:cputime_vs_N).

![image](images/cputime_vs_n.pdf)
<span id="fig:cputime_vs_N" label="fig:cputime_vs_N">\[fig:cputime\_vs\_N\]</span>

Indeed the FFT based pricing implementation requires less CPU time for a
fixed \(N\) than the FT based implementation. The reason why this does
not hold for the first three \(N\) is that in these cases the cost of
calling the recursion of the FFT algorithm on the machine level
predominates the effect of saving basic numerical operations. This is a
general issue when dealing with the implementation of numerical
procedures. Obviously for a fixed number of summands \(N\) the pFRFT
method is slower than the pFFT method since the FRFT algorithm applies
two inverse and one ordinary FFT, c.f. Section
[\[subsec:numerical\_algorithms\]](#subsec:numerical_algorithms).

Finally, we choose the pCOS pricing method as the most appropriate one
for our purposes. The pFT method is even more accurate for a fixed level
of spent CPU time as we deduce from
Figure [\[fig:pricing\_errors\_vs\_cputime\]](#fig:pricing_errors_vs_cputime),
but this takes place on an already extremely small level of overall
pricing error. The pCOS implementation however accomplishes this
precision with a much lower CPU time to be spent on the computations.

## Insights into the COS Method

Since we identified the COS method as the best-suited out of the four
presented numerical algorithms to solve the pricing problem for
arbitrary Lévy models, we now address some implementational issues. We
saw above that the COS method requires a lot less computational effort
than the FFT or FRFT algorithms due to reaching a comparable accuracy
with less summands. A natural question now is: How far can the number
\(N\) of summands be decreased to still get a sufficiently accurate
option value? In practise one could require the numerical result to be
correct up to two digits after the decimal point since this corresponds
to the smallest payable Euro amount. To analyse the effect of
diminishing \(N\) we perform a simulation study and apply the COS method
for both the Black-Scholes and the CGMY model to an initially fixed
randomly generated parameter set of size \(100\). For each of the
\(100\) parameter sets each single parameter is independently uniformly
sampled from a predetermined bounded interval, c.f.
Table [1](#tab:parameters).

<div id="tab:parameters">

|               |                  |
| :-----------: | :--------------: |
| **Parameter** |    **Domain**    |
|     \(r\)     |     \(0.01\)     |
|    \(S_0\)    |      \(5\)       |
|  \(\sigma\)   | \([0.05, 0.95]\) |
|     \(C\)     |  \([0.01, 2]\)   |
|     \(G\)     |    \([1, 5]\)    |
|     \(M\)     |    \([1, 5]\)    |
|     \(Y\)     |  \([0.1, 0.9]\)  |
|     \(T\)     |    \([1, 3]\)    |
|     \(K\)     |    \([4, 6]\)    |

Parameters domains for COS pricing in both the Black-Scholes and the
CGMY model.

</div>

The market parameters \(r\) and \(S_0\) as well as the option parameter
intervals for \(T\) and \(K\) are also depicted in
Table [\[tab:parameters\]](#tab:parameters). The strike range from
\(4\) to \(6\) corresponds to a \(20 \%\) spreading around the current
stock price, which covers relevant options for applications. As
reference prices we use the results of direct integration of
[\[eq:fourier\_pricing\_formula\]](#eq:fourier_pricing_formula) with
\(\alpha = -1.1\) via the MATLAB function *quadgk* in the CGMY case and
formula [\[eq:black\_scholes\_formula\]](#eq:black_scholes_formula) in
the Black-Scholes case. The truncation bounds for both the Black-Scholes
and the CGMY model are tabulated in Appendix
[\[app:cumulant\_derivation\]](#app:cumulant_derivation). For various
numbers of summands, ranging from \(10\) to \(20,000\), we apply the COS
method for both Lévy models to the random parameter sets, determine the
resulting \(100\) absolute pricing errors and plot the maximal, average
and minimal of these errors over the corresponding \(N\). The outcome is
depicted in
Figure [\[fig:pricing\_error\_vs\_cos\_n\]](#fig:pricing_error_vs_cos_n).

![image](images/pricing_error_vs_cos_n.pdf)
<span id="fig:pricing_error_vs_cos_n" label="fig:pricing_error_vs_cos_n">\[fig:pricing\_error\_vs\_cos\_n\]</span>
![image](images/pricing_error_for_parameter_set.pdf)
<span id="fig:pricing_error_for_parameter_set" label="fig:pricing_error_for_parameter_set">\[fig:pricing\_error\_for\_parameter\_set\]</span>

First note that the pricing error decay is much quicker for the
Black-Scholes model than for the CGMY model. This is due to the higher
degree of complexity of the four-parametric CGMY model compared to the
1-parametric Black-Scholes model. From the cumulant derivation in
Appendix [\[app:cumulant\_derivation\]](#app:cumulant_derivation) we get
that the Black-Scholes cumulants are zero from order three on. Thus the
truncation bounds do not change as higher order cumulants are taken into
account. For the CGMY model these higher order cumulants are especially
relevant, c.f. .
Figure [\[fig:pricing\_error\_for\_parameter\_set\]](#fig:pricing_error_for_parameter_set)
shows the pricing error for each of the parameter sets for \(20,000\)
summands in the COS method. It would be difficult to reduce these even
further as in contrast to the Black-Scholes case the reference value
does not stem from an analytical formula. It rather is itself obtain
from a numerical integration scheme. As such it is exposed to the
integration parameters, i.e. the level of the dampening coefficient
\(\alpha\). In Figure [1](#fig:pricing_error_vs_b-a) for a COS parameter
of \(N = 20,000\) we plot the resulting pricing error versus the width
\(b - a\) of the underlying integration domain. The plot indicates that
the maximal pricing error of \(0.0095\) in this case is not affected by
the size of the truncation bounds \(a\) and \(b\). Instead it reflects
some amount of inaccuracy of both COS and *quadgk* based pricing.

![image](Matlab_Images/pricing_error_vs_b-a.pdf)
<span id="fig:pricing_error_vs_b-a" label="fig:pricing_error_vs_b-a">\[fig:pricing\_error\_vs\_b-a\]</span>
