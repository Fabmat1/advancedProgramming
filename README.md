# Exercise submissions for advaced programming
### Prerequisites

Please install ```cmake``` by running

```
sudo apt-get update
sudo apt-get install cmake
```

### Usage
This project uses ```cmake``` to function, to build any particular executeable first navigate to the root directory of
the project and run 

```
cmake ./
```

After the program has finished, build any target by running

```
make <target_name>
```

and run it by navigating to the output folder and executing the script

```
cd output/
./<target_name>
```

The target names and their functions for each project are found below.

### Targets
#### Exercise 1

* ```quickpi``` - Script for computing pi with AVX-2 instructions and multiprocessing
* ```quickpi_512``` - Script for computing pi with AVX-512 instructions and multiprocessing **only works on compatible CPUs!**

#### Dark Matter Bonus Exercise

This consists of three files: Two python files and one C++ file. Instructions on reconstructing the exercise are:

* Move the ```GalaxyFromIllustrisTNG50Dark_DM_Subhalo852966.txt``` file into the ```Bonus_DM``` directory 
* Run ```task_1_2_8.py``` to make the plot for tasks 1,2 and 8. The result for task 8 is calculated later, but the result is hardcoded into ```task_1_2_8.py``` to avoid repetition.
* Make the target ```bonus_DM``` and run it
* Run ```task3ff.py``` to make neccessary calculations for task 3 through 7

#### Stars Bonus Exercise

This exercise requires that you install the Eigen library for c++!<br><br>
https://eigen.tuxfamily.org/
<br>
How to use the scripts for this exercise:
* Move the ```GalaxyFromIllustrisTNG50_Stars_Subhalo521803.txt``` file into the ```Bonus_Stars``` directory
* Make the target ```bonus_stars``` and run it
* Run ```bonus_stars_plotting.py``` to create all plots, plots are saved into the Bonus_Stars directory by default

#### Exercise 3: N-Body Simulations

* Make the target ```quickNbody``` and run it
* Run ```Nbody_sims.py``` to create all plots
