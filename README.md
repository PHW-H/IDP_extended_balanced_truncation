## Buildup of this web site

This site contains the MATLAB code from the Integrated design project of P.H.W. Hogendoorn. As well as the Simulink model used to simulate the reduced-order models and all data sets generated and used for validating the theoretical framework by Borja, Scherpen, and Fujimoto for extended balanced truncation.

### Matlab Models

Here all MATLAB models used for the intergration project can be found
Main matlab code:
- [Main Matlab Model](RLC_system_Pancras_version.m)

Random model generator
- [Random Generator](Random_model_generator.m)

Port Hamiltonian model constructors:
- [Model type 1](Modeltype41.m)
- [Model type 2](Modeltype42.m)
- [Model type 3](Modeltype43.m)


### Simulink Model 

- [Simulink Model](balanced_modelreduction_rlc.slx)

### Used Models and results

| Model download | Model Type | input | Reduction |  Model dimentions | Deviation from orgional General | Deviation from orgional extended | Simulink result |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| [Model 1](Modeltype41.m)| 1 | step 50 | 25% | C:20 L:20 | 1.6685e-7 | 0 | git status |
| [Model 2](Modeltype41.m)| 1 | step 50 | 25% | C:40 L:40 | 2.1235e-5 | 0 | git status |
| [Model 3](Modeltype41.m)| 1 | step 50 | 25% | C:60 L:60 | 2.1920e-5 | 5.4210e-20 | git status |
| [Model 4](Modeltype41.m)| 1 | sinusoidal | 25% | C:40 L:40 | 3.7724e-4 | 2.1684e-18 | git status |

### The Paper can be downloaded here

[IDP by Pancras Hogendoorn](Reductionoflarge_scaleelectricalmodels (1).pdf)


### The Paper can be downloaded here

[IDP by Pancras Hogendoorn](Reductionoflarge_scaleelectricalmodels (1).pdf)

