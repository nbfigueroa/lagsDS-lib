# lagsDS-lib
This package provides a C++ library for evaluation of the Locally Active Globally Stable (LAGS) Dynamical System, proposed in [1].

Descriptino of LAGS...

## Installation
Do the following steps:
* In your catkin src directory clone the repository
```
$ git clone https://github.com/nbfigueroa/lagsDS-lib.git
```
* wstool gets all other git repository dependencies, after the following steps you should see extra catkin 
  packages in your src directory.
```
$  wstool init
$  wstool merge lagsDS-lib/dependencies.rosinstall 
$  wstool up 
```
* Query and installs all libraries and packages 
```
$ rosdep install --from-paths . --ignore-src --rosdistro indigo 
```

### Usage
...

**References**     
> [1] Figueroa, N. and Billard, A. (2018) Locally Active Globally Stable Dynamical Systems: Theory, Learning and Experiments. In preparation.      

**Contact**: [Nadia Figueroa](http://lasa.epfl.ch/people/member.php?SCIPER=238387) (nadia.figueroafernandez AT epfl dot ch)

**Acknowledgments**
This work was supported by the EU project [Cogimon](https://cogimon.eu/cognitive-interaction-motion-cogimon) H2020-ICT-23-2014.
