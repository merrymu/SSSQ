This repository implements the algorithms in 'Secure Spatial Skyline Queries On Encrypted Data'.

We adopted 'Qhull' at http://www.qhull.org/ to generate rc.txt for storing the information of Voronoi diagram for the data points. 

We referred to the code at https://github.com/ruizhang2015/homework/SNN\_code/MinMax for reproducing MinMax to complete our project. 

We employed two distinct schemes: the practical ORE scheme, pORE at https://github.com/kevinlewi/fastore, and the multi-client ORE scheme, m-ORE at https://github.com/collisionl/m-ORE, designed specifically for multi-source scenarios.

The program can be run with the following command:
```./sssq db50k.txt rc50k.txt```

Execute the 'make' command before running.
