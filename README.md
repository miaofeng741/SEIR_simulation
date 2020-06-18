Code instructions
====
SEIR_simulation_model
--------
Wuhan city, the epicentre of COVID-19 in China, has experienced an epidemic outbreak since Dec 2019.

We developed an extended SEIR model using epidemiological and hospital admission data. 

A time-varying transmission rate was estimated, and the epidemic in Wuhan was reconstructed.Here, our procedure is divided into two parts. The first part is used to estimate the transmission rate, and the second part is used to simulate the situation of the epidemic situation when the hospital construction plan is changed.

Refer to the paper for more details on the model.

Transmission Rate 
-----
We collected public epidemiological data and developed an extended susceptible-exposed-infectious-removed (SEIR) model (see Methods and Supplementary Information in article). The time-varying transmission rate was estimated using Monte Carlo simulation and a time series moving average method (see Supplementary Information), and the epidemic outbreak in Wuhan was reconstructed.

In this model, we took into account factors such as population movement, person-to-person transmission, hospital admission, patient recovery and hospital treatment to determine the transmission rate.

See how to implement this part of the model in solve_beta_MTKL_simulation.py.

Sensitivity Analysis
-------
At the end of the model, we analyze the sensitivity of some necessary parameters.

The experiment of sensitivity analysis can be roughly divided into two parts.

The first part is about some sensitivity analysis in the process of solving beta.We studied the influence degree of initial value, incubation period and self-cure preiod on the transmission rate. The specific scheme is as follows:
* Select different program initial values:set1~set10.
* Changed the incubation period parameter from lognormal(1.443,0.64) to gamma(1.975,2.633)、lognormal(1.335,0.683) or lognormal(1.712,0.537).
* Changed the self-cure preiod parameter from gamma(22.563,0.842) to gamma(1.975,2.633)、lognormal(1.335,0.683) or lognormal(1.712,0.537).

This part of simulation experiment is realized by changing the relevant parameters in solve_beta_MTkl_simulation.py

The second part is about some sensitivity analysis when changing hospital protocols.In this part, we conducted sensitivity analysis on the transmission rate, the proportion of mild and severe patients, and time-related parameters.
* Three different transmission rates were adopted(Please refer to the paper for the specific way).
* Change the proportion of patients with mild or severe diseases from 0.81/0.19 to 0.86/0.14 or 0.76/0.24.
* The incubation period and the self-healing period were extended or shortened by three days, respectively, or their distribution patterns were changed.





