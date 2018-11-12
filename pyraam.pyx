# distutils: language=c++
# distutils: include_dirs = /usr/local/include

# % # cython: profile=True
# % # cython: linetrace=True

from libcpp cimport bool
from libcpp.string cimport string
# from libcpp.vector cimport vector
# from libcpp.unordered_map cimport unordered_map
# from libcpp.unordered_set cimport unordered_set

cdef extern from "raam.h" namespace "RAAM" :

    cdef cppclass Graph:

        Graph() except + 
        Graph(Graph) except + 
        Graph(double, double, int) except +

        int load_regions  (string, bool)
        int load_agents   (string, bool)
        int load_resources(string, bool)
        int load_edges    (string, int, int)
        int load_tunnels  (string)

        int increment_resource_supply(long long)
        int decrement_resource_supply(long long)
        int change_resource_supply(int, int)

        void set_transform(float, float, float, bool, float)

        void allocate_min_fixed(bool)
        void set_tau(double)
        void set_tol(double)

        void equalize_use(int, int, double, int, bool)

        double get_cost()

        void write(string)
        void write_frac_in_region(string)

                           
cdef class graph:
    """Docstring for class pyraam graph class."""

    cdef Graph c_graph #: This is an instance of the C++ object that we're wrapping.

    def __cinit__(self, double tau, double tol, int default_region = -1):
        self.c_graph = Graph(tau, tol, default_region)

    def load_regions(self, f, bool load_tau):
        """
        Blah blah blah.

        Parameters
        ---------
        name
            A string to assign to the `name` instance attribute.

        """

        cdef string cf = str.encode(f)
        return self.c_graph.load_regions(cf, load_tau)

    def load_agents(self, f, bool regions = False):
        cdef string cf = str.encode(f)
        return self.c_graph.load_agents(cf, regions)

    def load_resources(self, f, bool regions = False):
        cdef string cf = str.encode(f)
        return self.c_graph.load_resources(cf, regions)

    def load_edges(self, f, int mult = 1, int add = 0):
        cdef string cf = str.encode(f)
        return self.c_graph.load_edges(cf, mult, add)

    def load_tunnels(self, f):
        cdef string cf = str.encode(f)
        return self.c_graph.load_tunnels(cf)

    def decrement_resource_supply(self, long long resource):
        return self.c_graph.decrement_resource_supply(resource)

    def increment_resource_supply(self, long long resource):
        return self.c_graph.increment_resource_supply(resource)

    def change_resource_supply(self, int resource, int delta):
        return self.c_graph.change_resource_supply(resource, delta)

    def set_tol(self, double tol):
        self.c_graph.set_tol(tol)

    def set_tau(self, double tau):
        self.c_graph.set_tau(tau)

    def set_transform(self, float scale = 1, float offset = 0, float power = 1, bool log = False, float log_base = 2):
        self.c_graph.set_transform(scale, offset, power, log, log_base)

    def allocate_min_fixed(self, bool verbose = False):

        self.c_graph.allocate_min_fixed(verbose)
        
    def equalize_use(self, int cycles = 1, int max_moves = 0, double decay = 0, int max_tunnel_moves = 0, bool verbose = False):
        self.c_graph.equalize_use(cycles, max_moves, decay, max_tunnel_moves, verbose)
        
    def get_cost(self):
        return self.c_graph.get_cost()

    def write(self, f):
        cdef string cf = str.encode(f)
        self.c_graph.write(cf)
        
    def write_frac_in_region(self, f):
        cdef string cf = str.encode(f)
        self.c_graph.write_frac_in_region(cf)
        
