Code instructions
====


# SEIR_simulation_model
--------
Wuhan city, the epicentre of COVID-19 in China, has experienced an epidemic outbreak since Dec 2019.

We developed an extended SEIR model using epidemiological and hospital admission data. 

A time-varying transmission rate was estimated, and the epidemic in Wuhan was reconstructed.Here, our procedure is divided into two parts. The first part is used to estimate the transmission rate, and the second part is used to simulate the situation of the epidemic situation when the hospital construction plan is changed.

Refer to the paper for more details on the model.

# Transmission Rate 
-----
We collected public epidemiological data and developed an extended susceptible-exposed-infectious-removed (SEIR) model (see Methods and Supplementary Information in article). The time-varying transmission rate was estimated using Monte Carlo simulation and a time series moving average method (see Supplementary Information), and the epidemic outbreak in Wuhan was reconstructed.

In this model, we took into account factors such as population movement, person-to-person transmission, hospital admission, patient recovery and hospital treatment to determine the transmission rate.

See how to implement this part of the model in solve_beta_MTKL_simulation.py.

# Sensitivity Analysis
-------
At the end of the model, we analyze the sensitivity of some necessary parameters.

The experiment of sensitivity analysis can be roughly divided into two parts.

The first part is about some sensitivity analysis in the process of solving beta.

The second part is about some sensitivity analysis when changing hospital protocols.





