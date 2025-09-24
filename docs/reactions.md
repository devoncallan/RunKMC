## Supported Reactions


How to specify a reaction:
```
{RXN_ID} {REACTANTS} -{RATE COEFFICIENT}-> {PRODUCTS}
```

Initiator Decomposition
```
ID I -kd-> R + R 
```

Initiation
```
IN R + A -kpA-> PA
```

Thermal Initiation
```
TH 3A -kthA-> 2PA
```

Propagation
```
PR PA + B -kpAB-> PB
```

Termination by combination
```
TC PA + PB -ktcAB-> D
```

Termination by dispraportionation
```
TD PA + PA -ktdAA-> D + D
```

Chain transfer to monomer
```
TM PA + B -ktrAB-> D + PB
```

