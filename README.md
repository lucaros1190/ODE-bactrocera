# ODE-bactrocera
ODE model to simulate Bactrocera oleae life cycle associated with the publication of Rossini et al. (2022). Please if you use the data anche the code in this page cite the following papers:

[L. Rossini, N. Bono-Rossello, S. Speranza, E. Garone (2021). A general ODE-based model to describe the physiological age structure of ectotherms: Description and application to _Drosophila suzukii_ ](https://doi.org/10.1016/j.ecolmodel.2021.109673)

[L. Rossini, O. A. Bruzzone, M. Contarini, L. Bufacchi, S. Speranza (2022). A physiologically based ODE model for an old pest: Modeling life cycle and population dynamics of _Bactrocera oleae_ (Rossi).](https://doi.org/10.3390/agronomy12102298)

Insert the parameters (rate functions parameters and initial conditions of the model) into 'Parameters.py'
Fill 'Experimental_data.csv' with temperature values and monitoring data following the example.

Run 'python3 ODE-Solver.py' to run the script.

If you want to change the ODE model or the functions involved, just change the 'ODEbactro.py' file.
