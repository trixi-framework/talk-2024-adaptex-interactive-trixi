# Interactive Trixi.jl demo

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

ADAPTEX Annual meeting, 14th of May, 2024.

To reproduce the live demo of [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
shown in the talk, perform the following steps:

### Install Julia
Go to https://julialang.org/downloads and download the latest stable version of Julia (this
repository was created with Julia v1.10.3).

### Get reproducibility repository
Clone this reproducibility repository by executing
```shell
git clone https://github.com/trixi-framework/talk-2024-adaptex-interactive-trixi.git
```

### Install required packages
Go to the cloned repository and execute
```shell
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```
**This will take some time!**

### Start the demo

#### Using Jupyter

Go to the cloned repository and execute
```shell
JULIA_NUM_THREADS=4 jupyter-lab demo.ipynb
```
(multiple threads are not required, but will be helpful later on)

#### Using Julia directly

To generate the mesh:
```shell
julia --project=. scripts/mountain_dome.jl
```

To run the simulation:
```shell
JULIA_NUM_THREADS=4 julia --project=. -e scripts/mountain_gravity_wave.jl
```

## Authors
This repository was initiated by
[Benedict Geihe](https://www.mi.uni-koeln.de/NumSim/dr-benedict-geihe/).

The HOHQMesh Julia scripts were created with great help from 
[Andrew R. Winters](https://liu.se/en/employee/andwi94)!


## License
The contents of this repository are licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
