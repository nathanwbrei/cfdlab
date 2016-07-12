# Free Surfaces with LBM D3Q19
Project for Computational Fluid Dynamics Lab course. Group 6.

## Compilation.

Requirements: UNIX based OS.

Use ```make``` command to compile:
```
>> make
```

There is debug mode available. This mode includes additional checks for velocity and density, computes total mass in the domain and checks that there is no FLUID cell that is neighboring with a GAS cell directly. Also it calculates a time spent on the following parts: streaming, collision, flag updating and boundary treatment. 
```
>> make debug
```
## Usage
To run the program use
```
>> ./lbsim <exampleName> <number of processes>
```
For each example it is required to have ```<exampleName>.dat``` and ```<exampleName>.pgm``` files. They have to be placed in the ```examples``` directory.

Output is vtk-files and they can be found in the ```vtk-output``` directory (if thare is no such directory it will be created).

### Description of input files
We added a few new fields to the ```input.dat``` file:

| Parameter | Description |
| --- | --- |
|radius| Radius of a droplet|
|exchange_factor| Factor to increase speed of a mass exchnge (mainly to get bore beautiful videos) |
|forces_x| Gravity force in the ```x``` direction |
|forces_y| Gravity force in the ```y``` direction |
|forces_z| Gravity force in the ```z``` direction |

## What we have done
-1. Add test function
0. Gravity added
1. new flags
2. mass and fractions
3. DF from gas cells
4. Flag update
5. Parallelization with OpenMP
6. Other optimizations
