# PFC-CA1-CA3Model_Taxidisetal
Cortical-CA3-CA1 network model by [Taxidis et al]

## Introduction
During an Erasmus internship on which I based my Master thesis, I developed this model to fit the experimental data obtained by the University of Genoa. 
I have studied the behaviour of the prefrontal cortex and the hippocampal models separately before connecting them. Here is the code for the model of the prefrontal cortex only.

### File structure
- `OnlyCortexconnected.py` is the main file.
- `Parameters.py` contains all the parameters of the model.
- `Equations.py` contains all the equations of the model.
- Files that generate the matrices which define the connections in the network. There are two types of connectivity rules: random and gaussian. Following is the file that creates random connectivity matrices, in which each cell has the same probability to be connected to any other cell in the network:
    - `Connectivity_PFC_Random.py` creates the connection matrices for the prefrontal cortex.
- Following is the file that creates gaussian connectivity matrices, in which each cell is more likely to connect to its neighbours than to farther away cells:
    - `Connectivity_PFC_Gauss.py` creates the connection matrices for the prefrontal cortex.
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
$ python Connectivity_PFC_Random.py
```
- To create a network in which close cells are more likely to be connected to each other that farther away cells run:
```
$ python Connectivity_PFC_Gauss.py
```

2. Run the simulation:
```
$ python OnlyCortexconnected.py
```
### Expected output
A successful simulation outputs:
- `V_pc.mat` a matrix containing the voltage values of 6 cortical pyramidal cells 
- `V_ic.mat` a matrix containing the voltage values of 6 cortical interneurons 
- `Spikes_pc.mat` a matrix containing the time and the cell index of cortical pyramidal neurons' spikes
- `Spikes_ic.mat` a matrix containing the time and the cell index of cortical  interneurons' spikes

## Contact information

Eleonora Bernasconi ely97ber [at] gmail.com

<!-- simulators -->
[BRIAN2]: https://brian2.readthedocs.io/en/stable/index.html
<!-- references -->
[Taxidis et al]: https://www.frontiersin.org/articles/10.3389/fncom.2013.00003/full
