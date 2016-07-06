#include<random>
#include<vector>
#include<iostream>
#include<assert.h>
#include<fstream>
#include<algorithm>

using std::cout; using std::endl;
using std::vector;

int POP_X = 100;
int POP_Y = 100;

struct Symbiont {
  float donation;
  float points;

  Symbiont() : donation(-1), points(-1) {};
  Symbiont(float d, float p) : donation(d), points(p) {};
  Symbiont(float d) : donation(d), points(0.0) {};
  Symbiont(Symbiont &parent) : donation(parent.donation), points(0.0) {};
  Symbiont(const Symbiont& orig) : donation(orig.donation), points(orig.points) {};

  float update(float, int);
  void mutate(std::default_random_engine& r, float);


};

float Symbiont::update(float res, int rate) {
  //Sym receives some resources from host and can return some back
  //cout << "updating a sym, yey!" << endl;
  float returned = (res*donation)*rate;
  points += res - returned;
  return returned;
}

void Symbiont::mutate(std::default_random_engine& r, float rate){
  std::normal_distribution<double> dist(0.5, rate);
  double mutation = dist(r); // pull from dist
  donation = (2*donation)/(1+2*mutation);
  if(donation > 1) donation = 1;
  if(donation < 0) donation = 0;
  points = 0;
}




//Host struct
struct Host {
  float donation;
  float points;

  Symbiont sym;
  int cell_id;
  

  Host() =delete;
  Host(float d, float p, Symbiont s, int id) : donation(d), points(p), sym(s), cell_id(id) {};
  Host(float d, Symbiont s, int id) : donation(d), points(0.0), sym(s), cell_id(id) {};
  Host(const Host &orig) : donation(orig.donation), points(orig.points), cell_id(orig.cell_id), sym(orig.sym) {};

  void update(int);
  int chooseNeighbor(std::default_random_engine&);
  void mutate(std::default_random_engine&, double);
  void birth(Host);
  
};

void Host::birth(Host parent) {
  donation = parent.donation;
  points = 0.0;
  cell_id = -1;
  sym = Symbiont();
}

void Host::update(int sym_mult) {
  if(sym.donation != -1) {
    //Host donates resources to sym based on donation amount,
    // sym returns some of the resources if it chooses to
    //Sym always gets 5 units
    //cout << "updating a host yey!" << endl;
    //Resources not being donated to sym
    points += (5-donation*5);
    //Resources being donated to sym and possibly some back
    points += sym.update(donation*5 + 5, sym_mult);
  }
  else points++;
}

int Host::chooseNeighbor(std::default_random_engine &r){
  //TODO: Write this function yey
  int radius = 1;
  int cell_x = cell_id % POP_X;
  int cell_y = (cell_id - cell_x)/POP_X;
  int x, y;
  vector<int> neighbors;
  for(int i=(cell_x - radius); i<=(cell_x + radius); i++) {
    for(int j=(cell_y-radius); j<=(cell_y+radius); j++){
      if(i<0) x = POP_X + i;
      else if(i>POP_X) x = i-POP_X;
      else x = i;
      
      if(j<0) y = POP_Y +j;
      else if(j>=POP_Y) y = j-POP_Y;
      else y = j;
      
      neighbors.push_back(y*POP_X+x);
    }
  }
  
  std::random_shuffle(neighbors.begin(), neighbors.end());
  
  return neighbors[0];
}

void Host::mutate(std::default_random_engine& r, double rate){

  //Mutate function that moves toward .5 courtesty of Bob
  std::normal_distribution<double> dist(0.5, rate);
  double mutation = dist(r); // pull from dist
  donation = (2*donation)/(1+2*mutation);
  if (donation > 1) donation = 1;
  if (donation < 0) donation = 0;
  

  points = 0;
}


struct Population{
  vector<Host> pop;
  int seed;
  int final_update;
  int cur_update;
  float vert_rate; // Vertical transmission rate
  float mut_rate; // Mutation Rate
  int sym_mult; // Multiplication factor for symbionts
  float start_rate; // Donation rate the initial orgs start at
  std::default_random_engine engine;
  std::ofstream data_file;
  

  Population() : final_update(1) {};
  Population(int pop_count, int f, int seed);

  void init_pop(int pop_count);
  void evolve();
  void print_stats();


};


void Population::print_stats() {
  cout << cur_update << endl;
  double host_sum = 0.0;
  double sym_sum = 0.0;
  int host_count = 0;
  int sym_count = 0;
  //figure out avg host and sym donations
  for(auto org : pop){
    host_count++;
    host_sum += org.donation;
    if(org.sym.donation > -1){
      sym_count++;
      sym_sum += org.sym.donation;
    }
  }
  
  
  data_file << cur_update <<", "<< host_sum/host_count << ", " << sym_sum/sym_count << ", " << host_count << ", " << sym_count << endl <<std::flush;

  
}

Population::Population(int pop_count, int f, int seed_i) {
  final_update = f;
  cur_update = 0;
  seed = seed_i;
  data_file.open("avg_donation_"+std::to_string(seed)+"_mut"+std::to_string(mut_rate)+"_mult"+std::to_string(sym_mult)+"_vert"+std::to_string(vert_rate)+".csv", std::ofstream::ate);
  data_file << "Update, Host_Donation, Sym_Donation, Host_Count, Sym_Count" << endl << std::flush;

  std::default_random_engine engine(seed);
  init_pop(pop_count);
}

void Population::init_pop(int pop_count) {

  assert(pop_count == POP_X * POP_Y);
  for(int i=0; i<pop_count; ++i){
    Symbiont new_sym(start_rate);
    Host new_org(start_rate, new_sym, i);
    //cout << "new org!" << new_org.sym.donation << endl;
    pop.push_back(new_org);
    //cout << "inserted new org!" << pop[-1].sym.donation << endl;
    
  }
}

void Population::evolve(){
  
  std::uniform_real_distribution<double> dist(0, 1);

  for(cur_update = 0; cur_update < final_update; ++cur_update){
    //cout << "updating orgs" << endl;
    //Give everyone their points
    for(auto &org : pop) {
      //Run host and sym updates
      
      org.update(sym_mult);
      //cout << "updated an org " << org.donation << " " << org.sym.donation << endl;
      }
    
    //Go around again for reproduction
    for(auto &org : pop) {
      //See if sym reproduces
      if ((org.sym.donation > -1) && (org.sym.points>=50)){
        //cout << "Making a baby!" << endl;
        Symbiont &parent = org.sym;
        //Baby sym!
        Symbiont baby(parent);
        //Mutate both and reset and place
        baby.mutate(engine, mut_rate);
        parent.mutate(engine, mut_rate);
        int infected = org.chooseNeighbor(engine);
        pop[infected].sym = baby;
	
      }
      //See if host reproduces
      if (org.points >=50){
        //cout << "Making a host baby!" << endl;
        Host baby(org.donation, 0.0, Symbiont(), -1);
        if(dist(engine) < vert_rate){
          //Vertical transmission, TODO: put sym repro in a function
          Symbiont s_baby(org.sym);
          s_baby.mutate(engine, mut_rate);
          org.sym.mutate(engine, mut_rate);
          baby.sym = s_baby;
        }

        //Mutate both and reset and place
        baby.mutate(engine, mut_rate);
        org.mutate(engine, mut_rate);
        int squashed = org.chooseNeighbor(engine);
        baby.cell_id = squashed;
        pop[squashed] = baby;
      }
      
    }
    print_stats();
  
  }
  data_file.close();
}



int main(int argc, char *argv[]) {
  if (argc < 6) cout << "Usage: seed mut mult vert start" << endl;
  else{
  int seed = atoi(argv[1]);
  Population pop(10000, 1000, seed);
  pop.mut_rate = atof(argv[2]);
  pop.sym_mult = atoi(argv[3]);
  pop.vert_rate = atof(argv[4]);
  pop.start_rate = atof(argv[5]);
  pop.evolve();}


}
