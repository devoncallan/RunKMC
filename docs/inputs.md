# RunKMC Input File Specification

The RunKMC input file is a plain text file (`.txt`) consisting of **four required sections**, each terminated with an `end` keyword. 

1. [Parameters](#1-parameters-section)
2. [Species](#2-species-section)
3. [Rate Constants](#3-rate-constants-section)
4. [Reactions](#4-reactions-section)

## 1. Parameters Section
Defines the simulation parameters for the Kinetic Monte Carlo simulation.
### **Example:**
```
parameters
    num_units = 1e10
    termination_time = 50,000
    analysis_time = 500.0
end
```

### **Required parameters:**
- `num_units`: `integer`
    - Total number of particles/molecules
- `termination_time`: `float`
    - When to stop the simulation
- `analysis_time`: `float`
    - How frequently to output/analyze data

Note: the units for `termination_time` and `analysis_time` are arbitrary but should be consistent.

## 2. Species Section
Defines all chemical species in the system with 
### **Example:**
```
species
    # Syntax
    # {type} {name} [C0]={concentration} [FW]={molecular_weight} {optional_params}
    
    I   AIBN    [C0]=0.01   FW=164.2   f=0.50
    U   R       [C0]=0.0    FW=82.1
    M   A       [C0]=2.1    FW=104.2
    M   B       [C0]=0.9    FW=100.1
    P   PA
    P   PB
    P   PA|PB
end
```
### Species Types
- `I` = Initiator species
- `M` = Monomer species
- `U` = Generic unit species
- `P` = Polymer species


### Species Options:
- `[C0]` = Initial concentration (mol/L). If not set, defaults to zero.
- `[FW]` = Formula weight (g/mol). Mostly used for monomers (`M`) for molecular weight calculations.
- `f` = Initiator efficiency (0, 1]. Required for all initiators (`I`).

### Parameter Groups:
Another way to define a polymer species is by defining it as a 
```
P PA
P PB
P P PA|PB
```

What's the difference between defining `PA` and `P[A]`?

- `P[X]` = Polymer ending with X (`~~X`)
- `P[X.Y]` = Polymer ending with XY (`~~XY`)
- `P[X.Y.Z]` = Polymer ending with XYZ (`~~XYZ`)
  
This notation is used for specifying sequence, which is necessary for considering depropagation reactions. Currently, it is not strictly necessary for systems without depropagation.

**Wildcard notation**
  
- `P[-.A] P[R.A]|P[A.A]|P[B.A]`
  
This statement says that `P[R.A]`, `P[A.A]`, `P[B.A]` are all considered `P[-.A]`. So if a reaction involves `P[-.A]`, any one of the above species can be used.

### **Examples:**
```
I   AIBN    [C0]=0.01   FW=X    f=0.5
M   A       [C0]=8.7    FW=X
M   B       [C0]=1.3    FW=X
U   R       
P   PA
P   PB
P   D
```

## 3. Rate Constants Section
Defines all kinetic parameters:
```
rateconstants
    kd      = {float}   # Initiator decomposition
    kpAA    = {float}   # Propagation PA + A
    ktcAA   = {float}   # Termination (comb.) of PA + PA
    ...
end
```

## 4. Reactions Section
Defines the reaction network. Suppored reactions can be found [here](./reactions.md).

```
reactions
    ID AIBN -kd-> R + R
    IN R + A -kpAA-> P
    ...
end
```