# PFC-CA1-CA3Model_Taxidisetal
Cortical-CA3-CA1 network model by [Taxidis et al]

## Introduction
During an Erasmus internship on which I based my Master thesis, I developed this model to fit the experimental data obtained by the University of Genoa. 
I have studied the behaviour of the prefrontal cortex and the hippocampal models separately before connecting them. Here is the code for the complete model.

### File structure
- `PFC-HPC.py` is the main file.
- `Parameters.py` contains all the parameters of the model.
- `Equations.py` contains all the equations of the model.
- Files that generate the matrices which define the connections in the network. There are two types of connectivity rules: random and gaussian. Following are the files that create random connectivity matrices, in which each cell has the same probability to be connected to any other cell in the network:
    - `Connectivity_HP_Random.py` creates the connection matrices for the CA1 and CA3.
    - `Connectivity_PFC_Random.py` creates the connection matrices for the prefrontal cortex.
    - `Connectivity_HPC_PFC_Random.py` creates the connection matrices from the CA1 and CA3 to the prefrontal cortex.
    - `Connectivity_PFC_HPC_Random.py` creates the connection matrices from the prefrontal cortex to the CA1 and CA3.
- Following are the files that create gaussian connectivity matrices, in which each cell is more likely to connect to its neighbours than to farther away cells:
    - `Connectivity_HP_Gauss.py` creates the connection matrices for the CA1 and CA3.
    - `Connectivity_PFC_Gauss.py` creates the connection matrices for the prefrontal cortex.
    - `Connectivity_HPC_PFC_Gauss.py` creates the connection matrices from the CA1 and CA3 to the prefrontal cortex.
    - `Connectivity_PFC_HPC_Gauss.py` creates the connection matrices from the prefrontal cortex to the CA1 and CA3.

## Usage

### Preparing the simulation environment

The model is run using [Brian2] simulator, which is based on Python language. 

1. Install brian2 package:
```
$ conda install -c conda-forge brian2
```

2. Install other useful python packages
```
$ conda install matplotlib pytest ipython notebook
```

3. Install brian2tools package, which contains several useful functions to visualize Brian 2 simulations and recordings
```
$ conda install -c brian-team brian2tools
```

### Execution

1. Create the connectivity matrices. There are two types of connectivity rules: random and gaussian. Choose the preferred one and run the files accordingly.
- To create a randomly connected network run:
```
$ python Connectivity_HP_Random.py
$ python Connectivity_PFC_Random.py
$ python Connectivity_HPC_PFC_Random.py
$ python Connectivity_PFC_HPC_Random.py
```
- To create a network in which close cells are more likely to be connected to each other that farther away cells run:
```
$ python Connectivity_HP_Gauss.py
$ python Connectivity_PFC_Gauss.py
$ python Connectivity_HPC_PFC_Gauss.py
$ python Connectivity_PFC_HPC_Gauss.py
```

2. Run the simulation:
```
$ python PFC-HPC.py
```
### Expected output
A successful simulation outputs:
- `V_pc.mat` a matrix containing the voltage values of 6 cortical pyramidal cells 
- `V_ic.mat` a matrix containing the voltage values of 6 cortical interneurons 
- `V_p1.mat` a matrix containing the voltage values of 6 CA1 pyramidal cells 
- `V_p3.mat` a matrix containing the voltage values of 6 CA3 pyramidal cells 
- `V_ih.mat` a matrix containing the voltage values of 6 hippocampal interneurons (they are the same for CA1 and CA3)
- `Spikes_pc.mat` a matrix containing the time and the cell index of cortical pyramidal neurons' spikes
- `Spikes_ic.mat` a matrix containing the time and the cell index of cortical  interneurons' spikes
- `Spikes_p1.mat` a matrix containing the time and the cell index of CA1 pyramidal neurons' spikes
- `Spikes_p3.mat` a matrix containing the time and the cell index of CA3 pyramidal neurons' spikes
- `Spikes_ih.mat` a matrix containing the time and the cell index of hippocampal interneurons' spikes

## Contact information

Eleonora Bernasconi ely97ber [at] gmail.com

<!-- simulators -->
[BRIAN2]: https://brian2.readthedocs.io/en/stable/index.html
<!-- references -->
[Taxidis et al]: https://www.frontiersin.org/articles/10.3389/fncom.2013.00003/full
