# Runtime analysis of differently scaled problems

## Setup

We run the [Calliope Euro Model](https://github.com/timtroendle/euro-calliope) without scaling (*unscaled*), with manual scaling (*scaled*) and with Gurobi's automatic scaling (*gurobi*), i.e. Gurobi ScaleFlag = 2.

## Results

The log outputs are separately linked for [*scaled*](./analysis/unscaled.log), [*unscaled*](./analysis/scaled.log) and [*gurobi*](./analysis/gurobiscaled.log). What follows is essentially an analysis of these files.

## Quick summary

It follows a summary of how much time was spent in which part of the optimization of the different instances of the model. This only considers time spent in the Gurobi backend. We run the model in each case on a time subset of 01.01. to 30.06.

|sec\\instance|	*unscaled*|	*scaled*|	*gurobi*|	
|---|---|---|---|
|Preprocess |	61 | 69 |	61 |
|Barrier | 5609 | 4517 | 2782 |
|DPushes | 5241 | 25 | 1515 |
|PPushes | 971 | 46 | 1388 |
|Simplex | 2020 | 870 | 1552 |
|||||
|total| 13847 | 5462 | 7241 |
|lost in restarts| 4650 | 700 | 0 |

All times are given in seconds. 
- Preprocess entails everything Gurobi does before starting barrier, in particular removing redundant rows and cols in the constraint matrix (c.f. [seminal paper about lp preprocessing](https://link.springer.com/article/10.1007/BF01580428).
- Barrier is the application of the interior point method that finds an optimal point x. This point is not guaranteed to be a basic solution, on the contrary it is very likely to be a midface solution. To recover a basic solution with the optimal objective value next the crossover phase is started.
- First we want to find a square submatrix AB of the matrix A with full rank (a basis for A) such that the induced basic feasible solution xB = AB.inv (b - AN xN) is close to x (here AN is the complement of AB and xN the complement of xB. xN has all variables pushed to the bounds). For this we need to push some of the variables of x to their bounds. This happens in DPushes for the dual variables and in PPushes for the primal variables (why necessary for both? do we use primal-dual simplex for cleanup?).
- Lastly, we run a simplex starting from the "advanced basis" we just found to reoptimize xB.
In DPushes the dual variables are pushed to their bounds, i.e. 

## Observations

Caveat: what follows is derived from the log files and not directly visible in the above results
- *scaled* is considerably faster at finding a basis than the other two models!
- There is no model that is fastest in all parts of the algorithm.
- *unscaled* loses a lot of time because it needs to restart the crossover after about 2500s
- Both *scaled* and *unscaled* need to restart barrier method!
- both *scaled* and *gurobi* need to tighten their Markowitz tolerance, *gurobi* to 0.5, *scaled* to 0.25

My best guess is that the runtime of the Barrier and the Simplex method cannot be consistently improved by scaling the problem (c.f. also [this paper](https://link.springer.com/article/10.1007/s10589-011-9420-4)) although the literature does not agree on this. UPDATE: considering that barrier also needs to be restarted for the *unscaled* model we might need to reconsider this. It now seems that the largest part of the runtime difference is due to restarts. What are the exact conditions that barrier is restarted? 
Moreover, it seems that the basis recovery phase of crossover (DPushes and PPushes) depend on how well scaled the problem is.

Runtime plots of the DPushes phase of the algo.
For *unscaled*, *scaled* and *gurobi*. We plot seconds on the x-axis (not amount of seconds but the n-th second in the whole model run, i.e. to get seconds since start of DPushes subtract the smallest x-value) and amount of variables that need to be pushed on the y-axis. For readability we only show a subset of datapoints.


![unscaled](./analysis/scaled.png)
![scaled](./analysis/unscaled.png)
![gurobi](./analysis/gurobi.png)

Most interesting points:
- *unscaled* needs to restart the DPushes phase, probably due to numerical issues. Hypothesis: in order to push variables to their bounds the algo described by [Bixby](https://scholarship.rice.edu/bitstream/handle/1911/101733/TR91-32.pdf?sequence=1&isAllowed=y). This proceeds in phases (as visible in the log files). In each phase a set of desirable variables are chosen (desirable w.r.t. how close to their bounds they already are) and pushed if possible. Whether or not a variable can be pushed depends on condition number of the to-be-constructed basis (this is where Markowitz tolerance comes into play). If no basis can be constructed because too many variables are discarded, the process needs to be retried with a tightened (higher) tolerance.
- At least in the case of *unscaled* (but I conjecture this is the case in general) the first Dpushes are much faster than the last ones. This might be due to the fact that in the beginning it is simple to choose a large set of desirable variables to push but it gets continuously harder when there are already many variables in the basis (new variables / columns need to be sufficiently independent).

## Unscaled Log with commentaries (to be continued)

> [2020-06-11 09:57:17] DEBUG    Parameter OutputFlag unchanged  
> [2020-06-11 09:57:17] DEBUG    Value: 1  Min: 0  Max: 1  Default: 1  
> [2020-06-11 09:57:17] DEBUG    Changed value of parameter LogFile to /scratch/124306203.tmpdir/tmpok53brrv.log  
> [2020-06-11 09:57:17] DEBUG    Prev: gurobi.log  Default:  
> [2020-06-11 09:57:17] DEBUG    Parameter Crossover unchanged  
> [2020-06-11 09:57:17] DEBUG    Value: -1  Min: -1  Max: 5  Default: -1  
> [2020-06-11 09:57:17] DEBUG    Changed value of parameter FeasibilityTol to 1e-05  
> [2020-06-11 09:57:17] DEBUG    Prev: 1e-06  Min: 1e-09  Max: 0.01  Default: 1e-06  
> [2020-06-11 09:57:17] DEBUG    Changed value of parameter Method to 2  
> [2020-06-11 09:57:17] DEBUG    Prev: -1  Min: -1  Max: 5  Default: -1  
> [2020-06-11 09:57:17] DEBUG    Changed value of parameter OptimalityTol to 1e-05  
> [2020-06-11 09:57:17] DEBUG    Prev: 1e-06  Min: 1e-09  Max: 0.01  Default: 1e-06  
> [2020-06-11 09:57:17] DEBUG    Changed value of parameter Threads to 4  
> [2020-06-11 09:57:17] DEBUG    Prev: 0  Min: 0  Max: 1024  Default: 0  
> [2020-06-11 09:57:17] DEBUG    Optimize a model with 10845632 rows, 8413076 columns and 28359019 nonzeros  
  
> [2020-06-11 09:57:17] DEBUG    Coefficient statistics:  
> [2020-06-11 09:57:17] DEBUG    Matrix range     [1e-08, 2e+05]  
> [2020-06-11 09:57:17] DEBUG    Objective range  [1e+00, 1e+00]  
> [2020-06-11 09:57:17] DEBUG    Bounds range     [5e+01, 2e+06]  
> [2020-06-11 09:57:17] DEBUG    RHS range        [2e-01, 4e+07]  
> [2020-06-11 09:57:17] DEBUG    Warning: Model contains large matrix coefficient range  
> [2020-06-11 09:57:17] DEBUG    Consider reformulating model or setting NumericFocus parameter  
            	      	       to avoid numerical issues.			        
  
Gurobi warns us about bad conditioning of matrix  
  
> [2020-06-11 09:57:22] DEBUG    Presolve removed 7438998 rows and 4590000 columns (presolve time = 6s) ...  
> [2020-06-11 09:57:27] DEBUG    Presolve removed 8431809 rows and 5557257 columns (presolve time = 10s) ...  
> [2020-06-11 09:57:29] DEBUG    Presolve removed 8431809 rows and 5582505 columns  
> [2020-06-11 09:57:29] DEBUG    Presolve time: 12.87s  
> [2020-06-11 09:57:29] DEBUG    Presolved: 2413823 rows, 2830571 columns, 9098649 nonzeros  
> [2020-06-11 09:58:16] DEBUG    Ordering time: 0.00s  
  
In the presolve step redundant constraints and variables are removed. Moreover, variables get fixed at their bounds if possible  
  
> [2020-06-11 09:58:26] DEBUG    Barrier statistics:  
> [2020-06-11 09:58:26] DEBUG    Dense cols : 316  
> [2020-06-11 09:58:26] DEBUG    AA' NZ     : 7.114e+06  
> [2020-06-11 09:58:26] DEBUG    Factor NZ  : 5.573e+07 (roughly 2.6 GBytes of memory)  
> [2020-06-11 09:58:26] DEBUG    Factor Ops : 1.448e+10 (roughly 5 seconds per iteration)  
> [2020-06-11 09:58:26] DEBUG    Threads    : 4  
  
In each step of the barrier method AA' needs to be factored (how exactly?)  
How well this can be done depends on its sparsity and size.  
  
> [2020-06-11 09:58:35] DEBUG    Objective			   Residual  
> [2020-06-11 09:58:35] DEBUG   Iter Primal         Dual             Primal   Dual      Compl      Time  
> [2020-06-11 09:58:36] DEBUG    0   2.53024472e+14 -2.27949438e+15  1.16e+10 6.03e+02  4.73e+10    79s  
> [2020-06-11 09:58:46] DEBUG    1   2.07660823e+14 -1.18328109e+15  7.82e+09 8.79e+03  2.36e+10    90s  
> [2020-06-11 09:58:58] DEBUG    2   1.99818127e+14 -9.75190457e+14  7.16e+09 4.53e+03  2.00e+10   101s  
> [2020-06-11 09:59:12] DEBUG    3   1.63398081e+14 -1.00869916e+15  5.73e+09 1.32e+03  1.61e+10   115s  
...  
> [2020-06-11 10:33:33] DEBUG    179   1.03889559e+11  1.03878701e+11  5.84e-04 5.36e-07  1.75e+00  2177s  
> [2020-06-11 10:33:50] DEBUG    180   1.03889122e+11  1.03878847e+11  5.50e-04 5.39e-07  1.66e+00  2193s  
> [2020-06-11 10:34:01] DEBUG    181   1.03889045e+11  1.03879004e+11  5.44e-04 4.95e-07  1.62e+00  2204s  
> [2020-06-11 10:34:11] DEBUG    182   1.03888860e+11  1.03879102e+11  5.30e-04 4.57e-07  1.58e+00  2214s  
  
note that here the progress of the barrier method is set back: complementary slackness is back at initial 4.73e+10  
  
> [2020-06-11 10:34:30] DEBUG    183   2.53024472e+14 -2.27949438e+15  1.16e+10 6.03e+02  4.73e+10  2233s  
> [2020-06-11 10:34:42] DEBUG    184   2.16453812e+14 -1.82787438e+15  8.15e+09 1.03e+04  3.29e+10  2245s  
> [2020-06-11 10:34:55] DEBUG    185   1.97377174e+14 -1.55622510e+15  6.90e+09 2.42e+03  2.55e+10  2259s  
> [2020-06-11 10:35:07] DEBUG    186   1.65792608e+14 -1.44971309e+15  5.65e+09 1.32e+03  2.07e+10  2270s  
...  
> [2020-06-11 11:29:51] DEBUG    441   1.03882354e+11  1.03882107e+11  1.62e-05 9.30e-07  3.98e-02  5554s  
> [2020-06-11 11:30:08] DEBUG    442   1.03882342e+11  1.03882109e+11  1.52e-05 2.57e-06  3.76e-02  5572s  
> [2020-06-11 11:30:25] DEBUG    443   1.03882336e+11  1.03882114e+11  1.46e-05 2.39e-06  3.58e-02  5588s  
> [2020-06-11 11:30:46] DEBUG    444   1.03882329e+11  1.03882117e+11  7.25e-05 7.49e-06  3.42e-02  5610s  
> [2020-06-11 11:30:47] DEBUG    Barrier solved model in 444 iterations and 5609.83 seconds  
> [2020-06-11 11:30:47] DEBUG    Optimal objective 1.03882329e+11  
  
Barrier method has finished and yields an optimal non-basic x  
  
> [2020-06-11 11:30:48] DEBUG    Crossover log...  
> [2020-06-11 11:30:49] DEBUG    1032503 DPushes remaining with DInf 0.0000000e+00              5613s  
> [2020-06-11 11:30:52] DEBUG    875018 DPushes remaining with DInf 0.0000000e+00              5616s  
> [2020-06-11 11:30:57] DEBUG    872035 DPushes remaining with DInf 0.0000000e+00              5621s  
> [2020-06-11 11:31:02] DEBUG    869797 DPushes remaining with DInf 0.0000000e+00              5625s  
...  
> [2020-06-11 12:12:21] DEBUG    3426 DPushes remaining with DInf 0.0000000e+00              8104s  
> [2020-06-11 12:12:27] DEBUG    2680 DPushes remaining with DInf 0.0000000e+00              8111s  
> [2020-06-11 12:12:33] DEBUG    1923 DPushes remaining with DInf 0.0000000e+00              8116s  
> [2020-06-11 12:12:37] DEBUG    0 DPushes remaining with DInf 0.0000000e+00              8120s  
> [2020-06-11 12:12:37] DEBUG    Restart crossover...  
> [2020-06-11 12:12:39] DEBUG    1032503 DPushes remaining with DInf 0.0000000e+00              8123s  
> [2020-06-11 12:12:43] DEBUG    881068 DPushes remaining with DInf 0.0000000e+00              8127s  
> [2020-06-11 12:12:48] DEBUG    873666 DPushes remaining with DInf 0.0000000e+00              8132s  
> [2020-06-11 12:12:52] DEBUG    872174 DPushes remaining with DInf 0.0000000e+00              8135s  
...  
> [2020-06-11 12:58:00] DEBUG    2720 DPushes remaining with DInf 0.0000000e+00             10843s  
> [2020-06-11 12:58:06] DEBUG    1967 DPushes remaining with DInf 0.0000000e+00             10850s  
> [2020-06-11 12:58:08] DEBUG    1171 DPushes remaining with DInf 0.0000000e+00             10852s  
> [2020-06-11 12:58:10] DEBUG    0 DPushes remaining with DInf 2.7334562e-07             10854s  
  
In order to get a basic xB close to x returned by barrier, some of its variables are pushed to their bounds (these become the non-basic variables). The remaining variables that aren't at their bounds induce the indices of the basis B we use for starting the simplex.  
Here Pushing dual variables needs to be restarted. what is changed between the two runs??? it doesn't say anything about Markowitz tolerance here, only in the other examples where DPushes isn't restarted. But Bixby paper implies that the second run proceeds with a different singularity tolerance \tau.  
  
> [2020-06-11 12:58:11] DEBUG    2334769 PPushes remaining with PInf 1.0487462e+03             10854s  
> [2020-06-11 12:58:52] DEBUG    1680983 PPushes remaining with PInf 1.1602616e+04             10896s  
> [2020-06-11 12:59:15] DEBUG    1501353 PPushes remaining with PInf 1.0264356e+04             10919s  
> [2020-06-11 12:59:51] DEBUG    1469589 PPushes remaining with PInf 9.8895391e+03             10955s  
...  
> [2020-06-11 13:14:09] DEBUG    13896 PPushes remaining with PInf 3.3107697e+00             11813s  
> [2020-06-11 13:14:13] DEBUG    7248 PPushes remaining with PInf 8.7386683e-01             11816s  
> [2020-06-11 13:14:17] DEBUG    372 PPushes remaining with PInf 1.9073876e+04             11821s  
> [2020-06-11 13:14:22] DEBUG    0 PPushes remaining with PInf 1.8688597e+04             11825s  
> [2020-06-11 13:14:22] DEBUG    Push phase complete: Pinf 1.8688597e+04, Dinf 6.9700614e+09  11825s  
  
Same for primal variables, here no problems.  
  
> [2020-06-11 13:14:23] DEBUG    Iteration    Objective       Primal Inf.    Dual Inf.      Time  
> [2020-06-11 13:14:23] DEBUG    3365177    1.0388224e+11   0.000000e+00   6.970061e+09  11827s  
> [2020-06-11 13:14:56] DEBUG    3365647    1.0388224e+11   0.000000e+00   5.088006e+10  11859s  
> [2020-06-11 13:15:28] DEBUG    3366377    1.0388224e+11   0.000000e+00   2.621673e+09  11891s  
> [2020-06-11 13:16:12] DEBUG    3367047    1.0388223e+11   0.000000e+00   5.671475e+09  11935s  
...  
> [2020-06-11 13:47:24] DEBUG    3414608    1.0388212e+11   0.000000e+00   9.387971e+01  13808s  
> [2020-06-11 13:47:47] DEBUG    3414930    1.0388218e+11   1.837414e+02   0.000000e+00  13830s  
> [2020-06-11 13:47:55] DEBUG    3415040    1.0388218e+11   0.000000e+00   0.000000e+00  13838s  
> [2020-06-11 13:48:04] DEBUG    3415040    1.0388218e+11   0.000000e+00   0.000000e+00  13847s  
  
> [2020-06-11 13:48:04] DEBUG    Solved in 3415040 iterations and 13847.15 seconds  
> [2020-06-11 13:48:04] DEBUG    Optimal objective  1.038821755e+11  
> [2020-06-11 13:48:04] DEBUG    Warning: unscaled primal violation = 0.000488307 and residual = 0.000488307  
