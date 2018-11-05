# Rational Agent Access Model

The Rational Agent Access Model (RAAM) is a greedy optimization framework for calculating the accessibility of public resources.  The model is based on the intuition that agents minimize a cost function that is composed of a congestion and travel time pieces:

* min<sub>ℓ</sub> [(d<sub>ℓ</sub>/s<sub>ℓ</sub>)/ρ + t<sub>rℓ</sub>/τ]

The congestion is simply the demand for the resource over the supply, converted to a cost by a (fixed) factor ρ.  In our case -- primary healthcare in the United States, ρ is the national patient to physician ratio.  The travel cost from a residence to a resource location, t<sub>rℓ</sub>, is normalized by a configurable parameter, τ.

RAAM treats each demand location as a single agent, and on subsequent iterations shifts demand from the most expensive (used) to the least expensive available location.  The algorithm terminates with a single cost at each location (no cheaper cost anywhere).  RAAM can also be configured to allow agents to shift their demand to another agent, through "tunnels."

RAAM is implemented as a series of c++ classes, the highest level of which is exposed with cython, to python.

If c++ compiles are installed, you can build the `so` via

```
python setup.py build_ext --inplace
```

An example script and data for Chicago (2010 Census Tracts), `chicago.py`, suggests the basic functionality and file formats.

This is part of a project on accesssibility by James Saxon and Dan Snow, at the Center for Spatial Data Science and the Harris School, of the University of Chicago.  Contact jsaxon@uc for more information, or if you're interested in using the tool.

A Docker container for creating the origin-destination matrices used in RAAM (t<sub>rℓ</sub>) can be found at [routing-container](https://github.com/JamesSaxon/routing-container).
