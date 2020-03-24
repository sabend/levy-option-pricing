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

![image](Matlab_Images/Pricing_Error_vs_cputime.pdf)
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

![image](Matlab_Images/cputime_vs_N.pdf)
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
