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

  Symbiont() : donation(-2), points(-2) {};
  Symbiont(float d, float p) : donation(d), points(p) {};
  Symbiont(float d) : donation(d), points(0.0) {};
  Symbiont(Symbiont &parent) : donation(parent.donation), points(0.0) {};
  Symbiont(const Symbiont& orig) : donation(orig.donation), points(orig.points) {};

  float update(float, int);
  void mutate(std::mt19937& r, float);


};

float Symbiont::update(float res, int rate = 5) {
  //Sym receives some resources from host and can return some back
  //cout << "updating a sym, yey!" << endl;
  points += res;
  return 0.0;
}

void Symbiont::mutate(std::mt19937& r, float rate){
  //std::normal_distribution<double> dist(0.5, rate);
  //double mutation = dist(r); // pull from dist
  //donation = (2*donation)/(1+2*mutation);
  std::normal_distribution<double> dist(0, rate);
  donation += dist(r);
  
  if(donation > 1) donation = 1;
  if(donation < -1) donation = -1;
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
  int chooseNeighbor(std::mt19937&);
  void mutate(std::mt19937&, double);
  void birth(Host);
  
};

void Host::birth(Host parent) {
  donation = parent.donation;
  points = 0.0;
  cell_id = -1;
  sym = Symbiont();
}

void Host::update(int sym_mult) {
  float pool = 100.0;
  if (donation >= 0){
    //Donate resources, keep some for reproduction
    //If donation is 0, nothing is donated to sym or put into defense
    float donation_res = pool * donation;
    pool -= donation_res;
    
    
    if(sym.donation >= 0 ) {
        //Sym is nice so it won't steal
        //If donation is more than zero, some returned to host
        //If donation is zero, nothing returned and sym keeps all donated but doesn't steal
      float returned = donation_res * sym.donation;
      //Sym gets its resources
      sym.update(donation_res - returned);
      //Host gets the boosted res from symbiont
      points += returned * sym_mult;
    }
    else if (sym.donation < 0) {
      //Mean sym is going to keep all donated and steal some more!
      
      //Calculate how much the sym is stealing from the pool remaining and switch it to positive
      float stolen = sym.donation * pool * -1;
      pool -= stolen;
      //Give sym its donated and stolen resources, the meanie!
      sym.update(donation_res + stolen);
    }
    

  }
  else if (donation < 0) {
    //resources invested into defense and removed from pool, since donation is negative and pool is positive, addition is actually removing
    pool += pool*donation;

  
    if (sym.donation != -2 && sym.donation > 0){
    //Need to at some point fix how I represent an absent sym, but for now -2 means there isn't one
    //nice sym in a defensive host, so it doesn't get anything I guess?
      
    
    }
    else if (sym.donation != -2 && sym.donation < 0) {
    //mean sym in a defensive host, battle!
      if(sym.donation < donation){
        //Sym is able to steal from the remaining pool proportional to how mean it is
        //Example: sym is -0.5 and host is -0.4, difference is -0.1
        //Host had invested 0.4 into defense, so pool has 6 resources left and sym gets 0.1 of that
        float stolen = ((sym.donation - donation)* pool * -1);
        //Remove the stolen resources from the pool, stolen is a positive value now
        pool -= stolen;
        sym.update(stolen);
      
      }
      //If host puts more into defense than sym does into attack, host gets to keep the resources
    }
  }
  //Host gets to keep whatever is left in pool
  points += pool;
}

int Host::chooseNeighbor(std::mt19937 &r){

  std::uniform_int_distribution<int> dist(0, (POP_X*POP_Y) -1);
  return dist(r);
  
}

void Host::mutate(std::mt19937& r, double rate){

  //Mutate function that moves toward .5 courtesty of Bob
  //std::normal_distribution<double> dist(0.5, rate);
  //double mutation = dist(r); // pull from dist
  //donation = (2*donation)/(1+2*mutation);
  
  std::normal_distribution<double> dist(0, rate);
  donation += dist(r);
  if (donation > 1) donation = 1;
  if (donation < -1) donation = -1;
  

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
  std::mt19937 engine;
  std::ofstream data_file;
  std::ofstream sym_file;
  std::ofstream host_file;
  int rec_res = 100;

  Population() : final_update(1) {};
  Population(int pop_count, int f, int seed, float start_rate);

  void init_pop(int pop_count);
  void evolve();
  void print_stats();


};


void Population::print_stats() {
  //  cout << cur_update << endl;
  double host_sum = 0.0;
  double sym_sum = 0.0;
  int host_count = 0;
  int sym_count = 0;
  //figure out avg host and sym donations

  std::vector<int> sym_dists(21, 0);
  std::vector<int> host_dists(21, 0);


  for(auto org : pop){
    host_count++;
    host_sum += org.donation;
    int host_donate = (org.donation) * 10 + 10;
    host_dists[host_donate] += 1;
    if(org.sym.donation > -2){
      int sym_donate = (org.sym.donation * 10) + 10;
      sym_count++;
      sym_sum += org.sym.donation;

    }
  }
  
  //cout << cur_update <<", "<< host_sum/host_count << ", " << sym_sum/sym_count << ", " << host_count << ", " << sym_count << endl <<std::flush;
  data_file << cur_update <<", "<< host_sum/host_count << ", " << sym_sum/sym_count << ", " << host_count << ", " << sym_count << endl << std::flush;
  

  sym_file << cur_update << ", ";
  std::copy(sym_dists.begin(), sym_dists.end(), std::ostream_iterator<int>(sym_file, ", "));
  sym_file << endl << std::flush;
  host_file << cur_update << ", ";
  std::copy(host_dists.begin(), host_dists.end(), std::ostream_iterator<int>(host_file, ", "));
  host_file << endl << std::flush;

  
}

Population::Population(int pop_count, int f, int seed_i, float start) {
  final_update = f;
  cur_update = 0;
  seed = seed_i;
  start_rate = start;




  init_pop(pop_count);
}

void Population::init_pop(int pop_count) {

  assert(pop_count == POP_X * POP_Y);
  //cout << "start_rate" << start_rate << endl;
  for(int i=0; i<pop_count; ++i){
    Symbiont new_sym(start_rate);
    Host new_org(start_rate, new_sym, i);
    //cout << "new org!" << new_org.sym.donation << endl;
    pop.push_back(new_org);
    //cout << "inserted new org!" << pop[-1].sym.donation << endl;
    
  }
}

void Population::evolve(){
  std::string str_seed = std::to_string(seed);
  std::string str_mut = std::to_string(mut_rate);
  str_mut.erase ( str_mut.find_last_not_of('0') + 1, std::string::npos );
  std::string str_mult = std::to_string(sym_mult);
  //str_mult.erase ( str_mult.find_last_not_of('0') + 1, std::string::npos );
  std::string str_vert = std::to_string(vert_rate);
  str_vert.erase ( str_vert.find_last_not_of('0') + 1, std::string::npos );
  std::string str_start = std::to_string(start_rate);
  str_start.erase ( str_start.find_last_not_of('0') + 1, std::string::npos );
  
  
  data_file.open("avg_donation_"+str_seed+"_mut"+str_mut+"_mult"+str_mult+"_vert"+str_vert+"_start"+str_start+".csv", std::ofstream::ate);
  data_file << "Update, Host_Donation, Sym_Donation, Host_Count, Sym_Count" << endl << std::flush;
  sym_file.open("syms_donation_"+str_seed+"_mut"+str_mut+"_mult"+str_mult+"_vert"+str_vert+"_start"+str_start+".csv", std::ofstream::ate);
  sym_file << "Update, -1_-.9 -.9_-.8 -.8_-.7 -.7_-.6 -.6_-.5 -.5_-.4 -.4_-.3 -.3_-.2 -.2_-.1 -.1_0 0_.1 .1_.2 .2_.3 .3_.4 .4_.5 .5_.6 .6_.7 .7_.8 .8_.9 .9_1 1" << endl << std::flush;
  host_file.open("hosts_donation_"+str_seed+"_mut"+str_mut+"_mult"+str_mult+"_vert"+str_vert+"_start"+str_start+".csv", std::ofstream::ate);
  host_file << "Update, -1_-.9 -.9_-.8 -.8_-.7 -.7_-.6 -.6_-.5 -.5_-.4 -.4_-.3 -.3_-.2 -.2_-.1 -.1_0 0_.1 .1_.2 .2_.3 .3_.4 .4_.5 .5_.6 .6_.7 .7_.8 .8_.9 .9_1 1" << endl << std::flush;
  
  std::mt19937 engine(seed);
  std::uniform_real_distribution<double> dist(0, 1);


  for(cur_update = 0; cur_update < final_update; ++cur_update){
    //cout << cur_update << endl;
    //Give everyone their points
    for(auto &org : pop) {
      //Run host and sym updates
      
      org.update(sym_mult);
      //cout << "updated an org " << org.donation << " " << org.sym.donation << endl;
      }
    
    //Go around again for reproduction
    for(auto &org : pop) {
      //See if sym reproduces
      //cout << "sym here: " << org.sym.donation << " " << org.sym.points << endl;
      if ((org.sym.donation > -2) && (org.sym.points>=rec_res)){
        //cout << "Making a sym baby!" << endl;
        Symbiont &parent = org.sym;
        //Baby sym!
        Symbiont baby(parent);
        //Mutate both and reset and place
        baby.mutate(engine, mut_rate);
        parent.mutate(engine, mut_rate);
        int infected = org.chooseNeighbor(engine);
	if(pop[infected].sym.donation == -2)
	  //cout << "baby sym survives!" << endl;
	  pop[infected].sym = baby;
	
      }
      //See if host reproduces
      if (org.points >=(10*rec_res)){
        //cout << "Making a host baby!" << endl;
        Host baby(org.donation, 0.0, Symbiont(), -2);
        if(org.sym.points > 0 && dist(engine) < vert_rate){
          //Vertical transmission, TODO: put sym repro in a function
          Symbiont s_baby(org.sym);
          s_baby.mutate(engine, mut_rate);
          org.sym.mutate(engine, mut_rate);
          baby.sym = s_baby;
        }

	//cout << "host baby has sym: " << baby.sym.donation << endl;
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
  sym_file.close();
  host_file.close();
}



int main(int argc, char *argv[]) {
  if (argc < 6) cout << "Usage: seed mut mult vert start" << endl;
  else{
  int seed = atoi(argv[1]);

  Population pop(10000, 100000, seed, atof(argv[5]));
  pop.mut_rate = atof(argv[2]);
  pop.sym_mult = atoi(argv[3]);
  pop.vert_rate = atof(argv[4]);

  pop.evolve();}


}
