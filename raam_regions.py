#!/usr/bin/env python 

import pyraam  

tau, tol = 15, 0.01

default_region = -80

g = pyraam.graph(tau, tol, default_region)

g.load_regions  ("data/pcsa_frac.csv", False)

g.load_agents   ("data/pop.csv", True)
g.load_resources("data/doc.csv", True)
g.load_edges    ("data/self_times.csv")
g.load_edges    ("data/times.csv")

g.allocate_min_fixed()
g.equalize_use(1000, 250, 100, verbose = True)

for t in range(180, 301, 15):
  g.set_tau(t)
  g.equalize_use(1000, 250, 100, verbose = True)
  g.write_frac_in_region("f_in_region/{:03d}.csv".format(t))


