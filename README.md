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
| :---: | :---: | :---: | :---: | :--------------: | :--------------: | :-----: |
| [Model 1](Model_1.m)| 1 | step 50 | 25% | C:20 L:20 | 1.6685e-7 | 0 | [Image](Model_1.png) |
| [Model 2](Model_2.m)| 1 | step 50 | 25% | C:40 L:40 | 2.1235e-5 | 0 | [Image](Model_2.png) |
| [Model 3](Model_3.m)| 1 | step 50 | 25% | C:60 L:60 | 2.1920e-5 | 5.4210e-20 | [Image](Model_3.png) |
| [Model 4](Model_4.m)| 1 | sinusoidal | 25% | C:40 L:40 | 3.7724e-4 | 2.1684e-18 | [Image](Model_4.png) |
| [Model 5](Model_5.m)| 1 | step 50 | 50% | C:20 L:20 | 4.3211e-4 | 5.7762e-13 | [Image](Model_5.png) |
| [Model 6](Model_6.m)| 1 | step 50 | 50% | C:40 L:40 | 7.2875e-4 | 6.7763e-21 | [Image](Model_6.png) |
| [Model 7](Model_7.m)| 1 | sinusoidal | 50% | C:40 L:40 | 2.1920e-5 | 1.0842e-19 | [Image](Model_7.png) |
| [Model 8](Model_8.m)| 1 | sinusoidal | 25% | C:100 L:100 | 5.7047e-4 | 4.3368e-19 | [Image](Model_8.png) |

| Model download | Model Type | input | Reduction |  Model dimentions | Deviation from orgional General | Deviation from orgional extended | Simulink result |
| :---: | :---: | :---: | :---: | :--------------: | :--------------: | :-----: |
| [Model 9](Model_9.m)| 2 | step 50 | 25% | C:20 L:20 | 1.6685e-7 | 0 | [Image](Model_9.png) |
| [Model 2](Model_2.m)| 2 | step 50 | 25% | C:40 L:40 | 2.1235e-5 | 0 | [Image](Model_2.png) |
| [Model 3](Model_3.m)| 2 | step 50 | 25% | C:60 L:60 | 2.1920e-5 | 5.4210e-20 | [Image](Model_3.png) |
| [Model 4](Model_4.m)| 2 | sinusoidal | 25% | C:40 L:40 | 3.7724e-4 | 2.1684e-18 | [Image](Model_4.png) |
| [Model 9](Model_9.m)| 2 | step 50 | 50% | C:20 L:20 | 4.3211e-4 | 5.7762e-13 | [Image](Model_9.png) |
| [Model 6](Model_6.m)| 2 | step 50 | 50% | C:40 L:40 | 7.2875e-4 | 6.7763e-21 | [Image](Model_6.png) |
| [Model 7](Model_7.m)| 2 | sinusoidal | 50% | C:40 L:40 | 2.1920e-5 | 1.0842e-19 | [Image](Model_7.png) |
| [Model 8](Model_8.m)| 2 | sinusoidal | 25% | C:100 L:100 | 5.7047e-4 | 4.3368e-19 | [Image](Model_8.png) |

### error Grapes

### error of Model 1

A model of 100 Capasitors and 100 Inductors was used.

The used model can be downloaded here:

| [Model 10](Model_10.m) | [Model 12](Model_12.m) | [Model 13](Model_13.m) | [Model 14](Model_14.m) | [Model 15](Model_15.m) |


error generalised balanced:

<img src="error_type1_general.png">

error extended balanced:

<img src="error_type1_extanded.png">

####Model 2

[Model 9](Model_9.m)

####Model 3

[Model 11](Model_11.m)

### The Paper can be downloaded here

[IDP by Pancras Hogendoorn](Reductionoflarge_scaleelectricalmodels.pdf)
