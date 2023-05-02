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