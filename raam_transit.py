#!/usr/bin/env python 

import pyraam 

tau, tol = 60, 0.01

g = pyraam.graph(tau, tol)

g.load_agents   ("data/us_car_pop.csv")
g.load_resources("data/doc.csv")

# Driving times
g.load_edges    ("data/self_times.csv", 10, 0)
g.load_edges    ("data/times.csv",      10, 0)

# Public transportation.
g.load_edges    ("data/otp.csv",        10, 1)

g.allocate_min_fixed()
g.equalize_use(1000, 250, 100, verbose = True)
g.write("full_us_transit.csv")

