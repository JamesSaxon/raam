#!/usr/bin/env python 

import pyraam 
import math

tau = 60

g = pyraam.graph(tau, 0.01)

g.load_agents   ("data/il_pop.csv")
g.load_resources("data/il_doc.csv")
g.load_edges    ("data/il_times.csv")

g.allocate_min_fixed()
g.equalize_use(2500, 250, 50, 0, True)

g.write_frac_in_region("il_tun_ex/frac_in_region.csv")

# Must remove from tunnels 99'ers --
# non-zero tunnel cap but no population!
# g.load_tunnels  ("data/il_tunnels.csv")
g.load_tunnels  ("data/il_test_tunnel.csv")

steps, decay = 1, 500
init_loc, init_tun = 250, 10
for x in range(100):

  g.write("il_tun_ex/{:02d}.csv".format(x))

  redux = 2.**(x * steps / decay)

  g.equalize_use(steps, math.ceil(init_loc / redux), decay,
                        math.ceil(init_tun / redux), True)


g.write_frac_in_region("il_tun_ex/frac_in_region_tunnels.csv")

