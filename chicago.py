#!/usr/bin/env python 

import pyraam 
import math

tau, tol = 60, 0.01

g = pyraam.graph(tau, tol)

g.load_agents   ("data/cook_county_pop.csv")
g.load_resources("data/cook_county_doc.csv")
g.load_edges    ("data/cook_county_times.csv")

g.allocate_min_fixed()
g.equalize_use(2500, 250, 50, 0, True)

g.write("cook_county_test.csv")


