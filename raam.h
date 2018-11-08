#include <assert.h>

#include <iterator>
#include <map>

#include <math.h>
#include <limits>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream> 
#include <fstream>
#include "stdio.h"
#include "string"
#include "math.h"
#include <vector>
#include <iomanip>
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm> 

#include "csv.h"


#define DEFAULT_REGION -1
#define RHO 1.315e3

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::vector;
using std::string;
using std::stringstream;
using std::ofstream;

namespace RAAM {



const float FMAX = std::numeric_limits<float>::max();
const float FINF = std::numeric_limits<double>::infinity();

class Graph;
class Region;
class Agent;
class Resource;
class Edge;
class Tunnel;

class Graph {

  public:
    
    Graph(float tau = 15, float tol = 0.01, long long default_region = DEFAULT_REGION);

    unsigned int load_regions  (string fname, bool load_tau = false);
    unsigned int load_resources(string fname, bool load_regions = false);
    unsigned int load_agents   (string fname, bool load_regions = false);
    unsigned int load_edges    (string fname, int mult = 1, int add = 0);
    unsigned int load_tunnels  (string fname);

    void new_region  (long long _id, float _frac = 0.5, float _tau = 0);
    void new_agent   (long long _id, int _pop, long long _region = DEFAULT_REGION);
    void new_resource(long long _id, int _doc, long long _region = DEFAULT_REGION);
    void new_edge    (long long _agent, long long _resouce, float cost);
    void new_tunnel  (long long origin, long long destination, float capacity);

    unsigned int change_resource_supply(long long id, int delta);
    unsigned int increment_resource_supply(long long _id) { return change_resource_supply(_id, +1); }
    unsigned int decrement_resource_supply(long long _id) { return change_resource_supply(_id, -1); }

    float transform_cost(float cost);
    void set_transform(float scale = 1, float offset = 0, float power = 1, bool log = false, float log_base = 2);

    void allocate_min_fixed(bool verbose = false);
    void set_tol(float t);
    void set_tau(float a);
    void tune_tau(float delta = 1, float min = 0.1, float max = 120);
    bool equalize_use(unsigned int cycles = 1, unsigned int max_moves = 0, float decay = 0,
                      unsigned int max_tunnel_moves = 0, bool verbose = true);

    double get_cost();
    float get_error();

    void write(string filename);

    void write_frac_in_region(string filename);

  private:

    float _tau, _tol;

    bool _transform;
    float _transform_scale, _transform_offset, _transform_power;
    bool _transform_log;
    float _transform_log_base;

    map<long long, unsigned int> _addrA;
    map<unsigned int, long long> _idenA;

    map<long long, unsigned int> _addrR;
    map<unsigned int, long long> _idenR;

    map<long long, unsigned int> _addrReg;

    vector<Region*>   _regions;
    vector<Agent*>    _agents;
    vector<Resource*> _resources;
    vector<Edge*>     _edges;
    vector<Tunnel*>   _tunnels;

};

class Region {

  public:

    Region();
    Region(long long, float, float);

    void add_agent(Agent* a) { _agents.push_back(a); }

    long long get_id() { return _id; }
    float get_frac()   { return _frac; }
    float get_tau()  { return _tau; }

    unsigned int get_demand_in_region();
    float get_frac_in_region();
    void  set_tau(float a);
    void  tune_tau(float delta = 1, float min = 0, float max = 120);

  private:

    long long _id;
    float     _frac;
    float     _tau;

    vector<Agent*> _agents;

};

class Agent {

  public:
    
    Agent();
    Agent(long long id, int demand, long long region, float tau, float tol);

    void add_edge(Edge* e) { _edges.push_back(e); }
    void add_tunnel(Tunnel* t) { _tunnels.push_back(t); }

    string get_string();

    long long get_id()        { return _id; }
    int       get_demand()    { return _demand; }
    int       get_n_choices() { return _edges.size(); }

    float     get_tau()        { return _tau; }
    void      set_tau(float a) { _tau = a; }

    float     get_tol()        { return _tol; }
    void      set_tol(float a) { _tol = a; }

    unsigned int get_region() { return _region; }

    float     get_avg_cost();
    float     get_avg_fixed_cost();
    float     get_avg_supply_cost();

    float     get_avg_local_cost();
    float     get_avg_local_fixed_cost();
    float     get_avg_local_supply_cost();

    int       get_total_in_region();
    float     get_frac_in_region();

    Edge* get_min_fixed_edge();
    Edge* get_max_cost_edge();
    void allocate_min_fixed(bool verbose = false);

    bool equalize_use(unsigned int max_moves = 0);
    bool equalize_tunnels(unsigned int max_moves = 0);

    void return_demand(unsigned int d);
    void remove_demand(unsigned int d);

    void add_visitors(unsigned int nv);
    void remove_visitors(unsigned int nv);

  private:

    long long _id;
    long long _region;

    unsigned int _demand;
    unsigned int _absent;
    unsigned int _visitors;

    float _tau, _tol;
    bool _disconnected;

    vector<Edge*> _edges;
    vector<Tunnel*> _tunnels;

};

class Resource { 

  public:

    Resource();
    Resource(long long id, unsigned int supply, long long region);

    long long get_id()        { return _id; }
    int       get_supply()    { return _supply; }
    int       get_demand()    { return _demand; }
    int       get_n_agents()  { return _edges.size(); }
    int       get_rDS()       { return _demand / _supply; }

    unsigned int get_region() { return _region; }

    void add_edge(Edge* e) { _edges.push_back(e); }

    void change_demand(int d_use) { _demand += d_use; }

    unsigned int change_supply(int delta);

  private:

    long long _id;
    long long _region;
    unsigned int _supply;
    unsigned int _demand;
    vector<Edge*> _edges;

};

class Edge {

  public:
    
    Edge();
    Edge(Agent* agent, Resource* resource, float cost, bool in_region);

    long long get_agent_id()    { return _agent   ->get_id(); }
    long long get_resource_id() { return _resource->get_id(); }

    Agent*    agent   () { return _agent; }
    Resource* resource() { return _resource; }

    unsigned int  get_use() { return _use; }
    void set_use(unsigned int use);
    void change_use(int d_use);

    bool is_used() { return _use > 0; }
    bool is_in_region() { return _in_region; }

    float fixed_cost()  { return _cost; }
    float supply_cost() { return _resource->get_rDS() / RHO; }

    float total_cost(float tau);

  private:

    Agent* _agent;
    Resource* _resource;
    float _cost;
    unsigned int _use;

    bool  _in_region;

};


class Tunnel {

  public:
    
    Tunnel();
    Tunnel(Agent* origin, Agent* destination, unsigned int capacity);

    long long get_destination_id() { return _destination->get_id(); }
    Agent*    origin() { return _origin; }
    Agent*    destination() { return _destination; }

    unsigned int get_capacity() { return _capacity; }

    unsigned int get_use() { return _use; }
    unsigned int get_free_capacity() { return _capacity - _use; }

    void set_use(unsigned int use);
    void change_use(int d_use);

    float get_local_cost_at_destination() { return _destination->get_avg_local_cost(); }
    float get_local_fixed_cost_at_destination() { return _destination->get_avg_local_fixed_cost(); }
    float get_local_supply_cost_at_destination() { return _destination->get_avg_local_supply_cost(); }

  private:

    Agent* _origin;
    Agent* _destination;

    unsigned int _capacity;
    unsigned int _use;

};


Graph::Graph(float tau, float tol, long long default_region) :
  _tau(tau), _tol(tol),
  _transform(false), _transform_scale(1), _transform_offset(0), _transform_power(1),
  _transform_log(false), _transform_log_base(2) {
  new_region(default_region, 0.5); 
}

unsigned int Graph::load_regions(string fname, bool load_tau) {

  unsigned int nloaded(0);

  long long region; float frac_in_region, tau;

  cout << "Loading regions from " << fname << " ... ";
  if (load_tau) {
    io::CSVReader<3> region_csv(fname);
    region_csv.read_header(io::ignore_extra_column, "region", "frac_in_region", "tau");
    while(region_csv.read_row(region, frac_in_region, tau)){
      new_region(region, frac_in_region, tau);
      nloaded++;
    }
  } else {
    io::CSVReader<2> region_csv(fname);
    region_csv.read_header(io::ignore_extra_column, "region", "frac_in_region");
    while(region_csv.read_row(region, frac_in_region)){
      new_region(region, frac_in_region);
      nloaded++;
    }
  }
  cout << "done." << endl;

  return nloaded;

}

unsigned int Graph::load_edges(string fname, int mult, int add) {

  unsigned int nloaded(0);

  cout << "Loading edges from " << fname << " ... ";
  io::CSVReader<3> times_csv(fname);
  times_csv.read_header(io::ignore_extra_column, "origin", "destination", "cost");
  long long origin; long long destination; float cost;
  while (times_csv.read_row(origin, destination, cost)) {

    new_edge(origin * mult + add, destination, transform_cost(cost));
    nloaded++;

    if (!(nloaded % 1000000)) {
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
      cout << "Loaded " << nloaded / 1000000 << "M edges from " << fname << " ... ";
    }
  }
  cout << "done." << endl;

  return nloaded;

}


unsigned int Graph::load_tunnels(string fname) {

  unsigned int nloaded(0);

  cout << "Loading tunnels from " << fname << " ... ";
  io::CSVReader<3> tunnels_csv(fname);
  tunnels_csv.read_header(io::ignore_extra_column, "origin", "destination", "capacity");
  long long origin; long long destination; float capacity;
  while (tunnels_csv.read_row(origin, destination, capacity)) {
    new_tunnel(origin, destination, capacity);

    if (origin == destination) continue;

    nloaded++;

    if (!(nloaded % 100000)) {
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
      cout << "Loaded " << nloaded / 1000000 << "M tunnels from " << fname << " ... ";
    }
  }
  cout << "done." << endl;

  return nloaded;

}



unsigned int Graph::load_agents(string fname, bool load_regions) {

  cout << "Loading agents from " << fname << " ... ";

  long long geoid, region; int pop;

  unsigned int nloaded(0);

  if (load_regions) {
    io::CSVReader<3> region_csv(fname);
    region_csv.read_header(io::ignore_extra_column, "geoid", "pop", "region");
    while(region_csv.read_row(geoid, pop, region)){
      new_agent(geoid, pop, region);
      nloaded++;
    }
  } else {
    io::CSVReader<2> region_csv(fname);
    region_csv.read_header(io::ignore_extra_column, "geoid", "pop");
    while(region_csv.read_row(geoid, pop)){
      new_agent(geoid, pop);
      nloaded++;
    }
  }

  cout << "done." << endl;
  return nloaded;

}

unsigned int Graph::load_resources(string fname, bool load_regions) {

  cout << "Loading resources from " << fname << " ... ";

  long long geoid, region; int supply;

  unsigned int nloaded(0);

  if (load_regions) {
    io::CSVReader<3> region_csv(fname);
    region_csv.read_header(io::ignore_extra_column, "geoid", "supply", "region");
    while(region_csv.read_row(geoid, supply, region)){
      new_resource(geoid, supply, region);
      nloaded++;
    }
  } else {
    io::CSVReader<2> region_csv(fname);
    region_csv.read_header(io::ignore_extra_column, "geoid", "supply");
    while(region_csv.read_row(geoid, supply)){
      new_resource(geoid, supply);
      nloaded++;
    }
  }

  cout << "done." << endl;
  return nloaded;

}


void Graph::new_region(long long id, float frac, float tau) {

  assert(_addrReg.find(id) == _addrReg.end());

  _addrReg[id] = _regions.size();

  if (tau <= 0) tau = _tau;

  _regions.push_back(new Region(id, frac, tau));

}

void Graph::new_agent(long long id, int demand, long long region) {

  assert(_addrA.find(id) == _addrA.end());
  assert(_addrReg.find(region) != _addrReg.end());

  _addrA[id] = _agents.size();
  _idenA[_agents.size()] = id;
  
  float tau = _regions[_addrReg.find(region)->second]->get_tau();

  _agents.push_back(new Agent(id, demand, _addrReg[region], tau, _tol));

  _regions[_addrReg.find(region)->second]->add_agent(_agents.back());

}

void Graph::new_resource(long long id, int supply, long long region) {

  assert(_addrR.find(id) == _addrR.end());
  if (_addrReg.find(region) == _addrReg.end()) cout << region << endl;

  assert(_addrReg.find(region) != _addrReg.end());

  _addrR[id] = _resources.size();
  _idenR[_resources.size()] = id;

  _resources.push_back(new Resource(id, supply, _addrReg[region]));

}


void Graph::new_edge(long long agent_id, long long resource_id, float cost) {

  if (_addrA.find(agent_id) == _addrA.end()) cerr << agent_id << endl;

  if (_addrA.find(agent_id) == _addrA.end()) cout << agent_id << " " << resource_id << " " << cost << endl;
  assert(_addrA.find(agent_id)    != _addrA.end());

  if (_addrR.find(resource_id) == _addrR.end()) cout << agent_id << " " << resource_id << " " << cost << endl;
  assert(_addrR.find(resource_id) != _addrR.end());

  if (!_resources[_addrR[resource_id]]->get_supply()) return;

  bool in_region = (_agents[_addrA[agent_id]]->get_region() == _resources[_addrR[resource_id]]->get_region());

  _edges.push_back(new Edge(_agents[_addrA[agent_id]], _resources[_addrR[resource_id]], cost, in_region));

  _agents   [_addrA[agent_id]]   ->add_edge(_edges.back());
  _resources[_addrR[resource_id]]->add_edge(_edges.back());

}

void Graph::new_tunnel(long long origin, long long destination, float capacity) {

  if (!capacity) return;

  if (_addrA.find(origin)      == _addrA.end()) { cerr << origin      << " not yet set for tunnel!" << endl; return; }
  if (_addrA.find(destination) == _addrA.end()) { cerr << destination << " not yet set for tunnel!" << endl; return; }

  assert(_addrA.find(origin)      != _addrA.end());
  assert(_addrA.find(destination) != _addrA.end());

  _tunnels.push_back(new Tunnel(_agents[_addrA[origin]], _agents[_addrA[destination]], capacity));

  _agents[_addrA[origin]]->add_tunnel(_tunnels.back());

}

unsigned int Graph::change_resource_supply(long long resource, int delta) {

  if (_addrR.find(resource) == _addrR.end()) {
    cout << "Did not find resource " << resource << "; did not alter supply." << endl;
    return 0;
  } 

  return _resources[_addrR.find(resource)->second]->change_supply(delta);

}

float Graph::transform_cost(float cost) { 
  
  if (!_transform) return cost; 

  cost = cost * _transform_scale + _transform_offset;

  if (_transform_log) return log(cost) / log(_transform_log_base);
  
  return pow(cost, _transform_power);

}


void Graph::set_transform(float scale, float offset, float power, bool log, float log_base) {
  
  if (scale == 1 && offset == 0 && power == 1) return; 
  if (log && power != 1) {
    cout << "Warning: log and power both set.  Turning off power." << endl;
    power = 1;
  }

  if (log && offset <= 0) {
    cout << "Warning: log set with non-positive offset -- undefined for 0 cost.  Setting to offset to 1." << endl;
    offset = 1;
  }
  
  _transform = true; 
  _transform_scale = scale; 
  _transform_offset = offset; 
  _transform_power = power;

  _transform_log = log;
  _transform_log_base = log_base;

}


void Graph::allocate_min_fixed(bool verbose) {

  for (auto a : _agents) {

    a->allocate_min_fixed(verbose);

  }

}

bool Graph::equalize_use(unsigned int cycles, unsigned int max_moves, float decay,
                         unsigned int max_tunnel_moves, bool verbose) {


  bool converged;
  int ni;
  for (unsigned int ci = 0; ci < cycles; ci++) {

    unsigned int iter_max_moves(max_moves);
    unsigned int iter_max_tunnel_moves(max_tunnel_moves);
    if (decay > 0) {
      iter_max_moves = max_moves * pow(0.5, ci / decay);
      if (max_moves && iter_max_moves == 0) iter_max_moves = 1;

      iter_max_tunnel_moves = max_moves * pow(0.5, ci / decay);
      if (max_tunnel_moves && !iter_max_tunnel_moves) iter_max_tunnel_moves = 1;
    }

    converged = true;
    ni = 0;

    for (auto agent : _agents) {

      if (!agent->equalize_use(iter_max_moves) || 
          (max_tunnel_moves && !agent->equalize_tunnels(iter_max_tunnel_moves))) {
        converged = false;
      } else ni++;
    }
    
    if (verbose) {
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
      cout << "cycle " << ci << " :: " << std::setprecision(3) << std::fixed
                             << 1.*ni/_agents.size() << "   ";
      cout << std::flush;
    }

    if (converged) {
      if (verbose) cout << endl;
      return true;
    }
  }

  if (verbose) cout << endl;
  return false;


}


void Graph::set_tau(float a) { 
  
  for (auto r : _regions) r->set_tau(a);

}

void Graph::set_tol(float new_tol) {

  _tol = new_tol;

  for (auto a : _agents) a->set_tol(_tol);

}


void Graph::tune_tau(float delta, float min, float max) {

  for (auto r : _regions) r->tune_tau(delta, min, max);

}

void Graph::write(string filename) {

  ofstream ofile;
  ofile.open(filename);

  for (auto a : _agents) {
    ofile << a->get_string();
  }

  ofile.close();

}

void Graph::write_frac_in_region(string filename) {

  ofstream ofile;
  ofile.open(filename);

  for (auto r : _regions) {
    float F = r->get_frac_in_region();
    if (!std::isnan(F)) {
      ofile << r->get_id() << "," 
            << r->get_frac() << ","
            << F << endl;
    }
  }

  ofile.close();
}

double Graph::get_cost() {

  double total_cost(0);
  long long total_demand(0);
  for (auto a : _agents) {

    unsigned int N = a->get_demand();
    if (!N) continue;

    total_demand += N;
    float cost = a->get_avg_cost();
    if (std::isnan(cost)) continue;
    if (!(cost < FMAX)) continue;

    total_cost += N * a->get_avg_cost();
  }

  return total_cost / total_demand;

}

float Graph::get_error() {

  float error(0);
  for (auto r : _regions) {
    unsigned int N = r->get_demand_in_region();
    if (N) error += N * pow(r->get_frac_in_region() - r->get_frac(), 2);
  }

  return error;

}

Region::Region(long long id, float frac, float tau) :
  _id(id), _frac(frac), _tau(tau) {};

unsigned int Region::get_demand_in_region() {

  unsigned int den(0);
  for (auto a : _agents) den += a->get_demand();
  return den;

}

float Region::get_frac_in_region() {

  unsigned int num(0), den(0);
  for (auto a : _agents) {
    num += a->get_total_in_region();
    den += a->get_demand();
  }

  return 1. * num / den;

}

void Region::set_tau(float new_tau) {

  _tau = new_tau;

  for (auto a : _agents) a->set_tau(_tau);

}

void Region::tune_tau(float delta, float min, float max) {

  float fcalc = get_frac_in_region();
  if (fcalc > _frac) {
    if (_tau + delta < max) _tau += delta;
    else _tau = max;

  } else if (fcalc < _frac) {
    if (_tau - delta > min) _tau -= delta;
    else _tau = min;
  }

  set_tau(_tau);

}

Agent::Agent(long long id, int demand, long long region, float tau, float tol) : 
  _id(id), _region(region), _demand(demand), _absent(0), _visitors(0), _tau(tau), _tol(tol), _disconnected(false) { };


string Agent::get_string() {

  float avg_cost = get_avg_cost();
  float fixed_cost = get_avg_fixed_cost();
  float supply_cost = get_avg_supply_cost();

  float cost = supply_cost + fixed_cost;

  if (!(cost < FMAX) && fabs(avg_cost - supply_cost - fixed_cost) > 0.02) {
    cout << "cost=" << avg_cost << "  supply=" << supply_cost << "  fixed=" << fixed_cost << endl;
  }

  assert(!(cost < FMAX) || fabs(avg_cost - supply_cost - fixed_cost) < 0.02);

  if (std::isnan(cost)) return "";

  stringstream ss;
  ss << _id << "," << _demand << ",";
  if (cost < FMAX) {
    ss << cost << "," << fixed_cost << "," << supply_cost;
  } else ss << ",,";
  ss << endl;

  return ss.str();

}

float Agent::get_avg_supply_cost() {

  float total_cost = get_avg_local_supply_cost() * (_demand - _absent);
  if (_demand == _absent) total_cost = 0;

  unsigned int tunnel_use(0);
  for (auto t : _tunnels) {

    if (!t->get_use()) continue;

    total_cost += t->get_use() * t->get_local_supply_cost_at_destination();
    tunnel_use += t->get_use();
  }

  assert(_absent == tunnel_use);

  return total_cost / _demand;

}

float Agent::get_avg_fixed_cost() {

  float total_cost = get_avg_local_fixed_cost() * (_demand - _absent);
  if (_demand == _absent) total_cost = 0;

  unsigned int tunnel_use(0);
  for (auto t : _tunnels) {

    if (!t->get_use()) continue;

    total_cost += t->get_use() * t->get_local_fixed_cost_at_destination();
    tunnel_use += t->get_use();
  }

  assert(_absent == tunnel_use);

  return total_cost / _demand;

}

float Agent::get_avg_cost() {

  float total_cost = get_avg_local_cost() * (_demand - _absent);
  if (_demand == _absent) total_cost = 0;

  unsigned int tunnel_use(0);
  for (auto t : _tunnels) {

    if (!t->get_use()) continue;

    total_cost += t->get_use() * t->get_local_cost_at_destination();
    tunnel_use += t->get_use();
  }

  assert(_absent == tunnel_use);

  return total_cost / _demand;

}

float Agent::get_avg_local_fixed_cost() {

  float cost(0);

  if (!_edges.size()) return FINF;

  unsigned int local_demand(0);
  for (auto e : _edges) {
    if (e->is_used()) {
      local_demand += e->get_use();
      cost += e->get_use() * e->fixed_cost();
    }
  }

  assert(local_demand == _demand + _visitors - _absent);

  return cost / _tau / local_demand;

}


float Agent::get_avg_local_supply_cost() {

  float cost(0);

  if (!_edges.size()) return FINF;

  unsigned int local_demand(0);
  for (auto e : _edges) {
    if (e->is_used()) {
      local_demand += e->get_use();
      cost += e->get_use() * e->supply_cost();
    }
  }

  assert(local_demand == _demand + _visitors - _absent);

  return cost / local_demand;

}


float Agent::get_avg_local_cost() {

  float cost(0);

  if (!_edges.size()) return FINF;

  unsigned int local_demand(0);
  for (auto e : _edges) {
    if (e->is_used()) {
      local_demand += e->get_use();
      cost += e->get_use() * e->total_cost(_tau);
    }
  }

  if (local_demand != _demand + _visitors - _absent) cout << local_demand << " == " <<  _demand << " + " <<  _visitors << " - " << _absent << " == " << _demand + _visitors - _absent << endl;
  assert(local_demand == _demand + _visitors - _absent);

  return cost / local_demand;

}

int Agent::get_total_in_region() {

  int total_in_region(0);
  for (auto e : _edges) if (e->is_in_region()) {
    total_in_region += e->get_use();
  }

  return total_in_region;

}

float Agent::get_frac_in_region() {

  return 1. * get_total_in_region() / _demand;

}

Edge* Agent::get_min_fixed_edge() {

  Edge* min_edge(0);
  float min_cost(1e10);

  for (auto e : _edges) {

    float mcost = e->fixed_cost();

    if (mcost < min_cost) {
      min_cost = mcost;
      min_edge = e;
    }
  }

  return min_edge;

}

Edge* Agent::get_max_cost_edge() {

  Edge* max_edge(0);
  float max_cost(0);

  for (auto e : _edges) {

    if (!e->get_use()) continue;

    float mcost = e->total_cost(_tau);

    if (mcost > max_cost) {
      max_cost = mcost;
      max_edge = e;
    }
  }

  return max_edge;

}


void Agent::allocate_min_fixed(bool verbose) {

  for (auto e : _edges)   e->set_use(0);
  for (auto t : _tunnels) t->set_use(0);

  Edge* min_edge = get_min_fixed_edge();
  if (min_edge) min_edge->set_use(_demand);
  else _disconnected = true;

  if (verbose && _disconnected && _demand) 
    cout << "Warning: no minimum location found for " << _id << ", which is populated." << endl;

}


bool Agent::equalize_use(unsigned int max_moves) {

  Edge* min_edge(0); float min_cost(1e11);
  Edge* max_edge(0); float max_cost(0);


  for (auto e : _edges) {

    float mcost = e->total_cost(_tau);

    if (mcost < min_cost) {
      min_cost = mcost;
      min_edge = e;
    }

    if (!e->is_used()) continue;
    
    if (mcost > max_cost) {
      max_cost = mcost;
      max_edge = e;
    }
  }

  if (!min_edge) return true;
  if (!max_edge) return true;
  if (min_edge == max_edge) return true;
  if (max_cost - min_cost < _tol) return true;


  float dtot = min_edge->get_use() + max_edge->get_use();

  float dmin_other = min_edge->resource()->get_demand() - min_edge->get_use();
  float dmax_other = max_edge->resource()->get_demand() - max_edge->get_use();

  float smin = min_edge->resource()->get_supply();
  float smax = max_edge->resource()->get_supply();

  assert(smin > 0);
  assert(smax > 0);

  float tmin = min_edge->fixed_cost();
  float tmax = max_edge->fixed_cost();

  int Nc = RHO;

  int dmin = floor(
      (smin * smax / (smin + smax)) * (
        (Nc / _tau) * (tmax - tmin) +
        (dtot + dmax_other) / smax - dmin_other / smin
        )
      );
  int dmax = dtot - dmin;
  if (dmax < 0) dmax = 0; // Hard constraint...

  // int delta_min = dmin - min_edge->resource()->get_demand();
  int delta_max = max_edge->get_use() - dmax;

  if (delta_max < 0) {
    cout << "delta max = " << delta_max << "  max_edge agent=" << max_edge->get_agent_id() << "  resource=" << max_edge->get_resource_id() << endl;
    return false; // cheap fix...
  }
  assert(delta_max >= 0);

  unsigned int delta = delta_max;
  if (delta > max_edge->get_use())    delta = max_edge->get_use();
  if (max_moves && delta > max_moves) delta = max_moves;

  assert(delta >= 0);

  min_edge->change_use(+delta);
  max_edge->change_use(-delta);

  // Disconnected elements should not reach here.
  // unsigned int local_demand = 0;
  // for (auto e : _edges) local_demand += e->get_use();
  // assert(local_demand == _demand + _visitors - _absent);

  if (delta_max == 1) return true;

  // cout << "tmin=" << tmin << " tmax=" << tmax << endl;
  // cout << "change dmin " << delta << " from=" << min_edge->get_use() << " to=" << dtot-dmax << endl;
  // cout << "change dmax " << delta << " from=" << max_edge->get_use() << " to=" << dmax << endl;
  // cout << "delta_max=" << delta_max << endl;
  // cout << "smin=" << smin << "  smax=" << smax << endl;
  // cout << "max_cost=" << max_cost << "  min_cost=" << min_cost << endl;
  // cout << "Not converged=" << _id << " shift from=" << max_edge->get_resource_id() << "  to=" << min_edge->get_resource_id() << "  delta=" << delta << "  delta_max=" << delta_max << endl;

  return false;

}


bool Agent::equalize_tunnels(unsigned int max_moves) {

  if (!_tunnels.size()) return true;

  float local_cost = get_avg_local_cost();

  Tunnel* max_tunnel(0);
  float max_cost = local_cost + _tol;
  for (auto t : _tunnels) {

    if (!t->get_use()) continue;

    float tunnel_cost = t->get_local_cost_at_destination();
    if (tunnel_cost >= max_cost) {

      max_cost = tunnel_cost;
      max_tunnel = t;
    }
  }
  
  if (max_tunnel) {

    int delta = std::min(max_tunnel->get_use(), max_moves);

    max_tunnel->change_use(-delta);
    return false;
  }

  Tunnel* min_tunnel(0);
  float min_cost = local_cost - _tol;
  for (auto t : _tunnels) {

    if (!t->get_free_capacity()) continue;

    float tunnel_cost = t->get_local_cost_at_destination();

    if (tunnel_cost > FMAX) continue;

    if (tunnel_cost < min_cost) {

      min_cost = tunnel_cost;
      min_tunnel = t;
    }
  }
  
  if (min_tunnel) {

    int delta = std::min(min_tunnel->get_free_capacity(), max_moves);

    // cout << _id << "->" << min_tunnel->get_destination_id()
    //      << "   local=" << local_cost
    //      << "   cost=" << min_cost
    //      << "   >fmax=" << (min_cost > FMAX)
    //      << "   max_moves=" << max_moves
    //      << "   capacity=" << min_tunnel->get_free_capacity()
    //      << "   delta=" << delta
    //      << endl;

    min_tunnel->change_use(delta);
    return false;
  }

  return true;

}


void Agent::remove_demand(unsigned int d) {

  if (d > _demand - _absent) { 
    
    int tunnel_cap(0);
    for (auto t : _tunnels) tunnel_cap += t->get_capacity();

    cerr << _id << "  demand=" << _demand << "  tunnel_cap=" << tunnel_cap << endl;
  }

  assert(d <= _demand - _absent);

  _absent += d;

  // The location is disconnected and there are no docs.
  // The demand is not assigned to an edge, so it cannot be removed.
  // Demand removed through a tunnel will never come back,
  // and there are no other choices as to its assignment.
  // It will necessarily leave.
  if (_disconnected) d = 0;

  Edge* max_edge(0); 
  while (d > 0) {

    max_edge = get_max_cost_edge();

    assert(max_edge);
    assert(max_edge->get_use() > 0);

    unsigned int delta = std::min(max_edge->get_use(), d);
    max_edge->change_use(-delta);

    assert(delta > 0);
    d -= delta;

  }

}


void Agent::return_demand(unsigned int d) {

  assert(d <= _absent);

  Edge* min_edge = get_min_fixed_edge();
  if (min_edge) {

    _absent -= d;
    min_edge->change_use(d);

  }

}



void Agent::add_visitors(unsigned int nv) {

  Edge* min_edge = get_min_fixed_edge();
  if (min_edge) {

    _visitors += nv;
    min_edge->change_use(nv);
  }

}

void Agent::remove_visitors(unsigned int nv) {

  assert(nv <= _visitors);

  _visitors -= nv;

  Edge* max_edge(0); 
  while (nv > 0) {

    max_edge = get_max_cost_edge();

    assert(max_edge);
    assert(max_edge->get_use() > 0);

    unsigned int delta = std::min(max_edge->get_use(), nv);
    max_edge->change_use(-delta);

    assert(delta > 0);
    nv -= delta;

  }

}


Resource::Resource(long long id, unsigned int supply, long long region) : 
  _id(id), _region(region), _supply(supply), _demand(0) { };

unsigned int Resource::change_supply(int delta) {
  
  if (_supply + delta <= 0) {
    cout << "Attempting to set supply to 0 or negative.   Aborting and returning." << endl;
    return 0;
  }

  _supply += delta;

  return _supply;

}

Edge::Edge(Agent* agent, Resource* resource, float cost, bool in_region) :
  _agent(agent), _resource(resource), _cost(cost), _use(0), _in_region(in_region) { }


void Edge::change_use(int d_use) {

  assert(_use + d_use >= 0);

  _use += d_use;

  _resource->change_demand(d_use);

}

void Edge::set_use(unsigned int use) {

  int d_use = use - _use;
  _use = use;

  _resource->change_demand(d_use);

}


float Edge::total_cost(float tau) {

  return supply_cost() + fixed_cost() / tau;

}

Tunnel::Tunnel(Agent* origin, Agent* destination, unsigned int capacity) : 
  _origin(origin), _destination(destination),
  _capacity(capacity), _use(0) { }

void Tunnel::set_use(unsigned int use) {

  int d_use = use - _use;
  _use = use;

  change_use(d_use);

}

void Tunnel::change_use(int d_use) {

  assert(_use + d_use >= 0);

  if (_use + d_use > _capacity) {
    cerr << _origin->get_id() << "->" << _destination->get_id()
          << " attempting use=(" << _use << " + " << d_use << ") "
          << " above capacity=" << _capacity << endl;
  }

  assert(_use + d_use <= _capacity);

  _use += d_use;

  if (d_use > 0) {
    _origin->remove_demand(d_use);
    _destination->add_visitors(d_use);
  }

  if (d_use < 0) {
    _origin->return_demand(abs(d_use));
    _destination->remove_visitors(abs(d_use));
  }

}


}



