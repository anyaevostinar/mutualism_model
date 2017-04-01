#include<random>
#include<vector>
#include<iostream>
#include<assert.h>
#include<fstream>
#include<algorithm>
#include<iterator>
#include<set>
#include "../../Empirical/tools/vector.h"
#include "../../Empirical/tools/math.h"
using std::cout; using std::endl;
using std::vector;
using std::set;

int POP_X = 100;
int POP_Y = 100;
vector<int> resources = {0, 1};
int host_tasks_lim = 1;
int sym_tasks_lim = 1;

//Code from http://stackoverflow.com/questions/6942273/get-random-element-from-container
template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
  std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
  std::advance(start, dis(g));
  return start;
}


struct Symbiont {
  float donation;
  float points;
  vector<int> tasks;

  Symbiont() : donation(-2), points(0) {};
  Symbiont(float d, float p, vector<int> s) : donation(d), points(p), tasks(s) {};
  Symbiont(float d, vector<int> s) : donation(d), points(0.0), tasks(s) {};
  Symbiont(Symbiont &parent) : donation(parent.donation), points(0.0), tasks(parent.tasks) {};
  Symbiont(const Symbiont& orig) : donation(orig.donation), points(orig.points), tasks(orig.tasks) {};

  float update(float, int);
  void mutate(std::mt19937& r, float);


};

float Symbiont::update(float res, int rate = 5) {
  //Sym receives some resources from host and can return some back
  //cout << "updating a sym, yey!" << endl;
  //  cout << "Symbiont res: " << res << endl;
  points += res;
  return 0.0;
}

void Symbiont::mutate(std::mt19937& r, float rate){
  //  std::normal_distribution<double> dist(0.5, rate);
  //double mutation = dist(r); // pull from dist
  //donation = (2*donation)/(1+2*mutation);
  std::normal_distribution<double> dist(0, rate);
  double mutation = dist(r);
  donation += mutation;
  
  if(donation > 1) donation = 1;
  if(donation < -1) donation = -1;

  std::uniform_real_distribution<double> dist2(0, 1);
  double mutation2 = dist2(r);
  if (mutation2<rate){
    tasks.push_back(*select_randomly(resources.begin(), resources.end(), r));
    if (tasks.size() >= sym_tasks_lim)
      tasks.erase(std::remove(tasks.begin(), tasks.end(), *select_randomly(tasks.begin(), tasks.end(), r)), tasks.end());
    set<int> temp(tasks.begin(), tasks.end());
    tasks = vector<int>(temp.begin(), temp.end());
  }

  points = 0;
}




//Host struct
struct Host {
  float donation;
  float points;
  vector<int> tasks;

  Symbiont sym;
  int cell_id;
  

  Host() =delete;
  Host(float d, float p, Symbiont s, int id, vector<int> t) : donation(d), points(p), sym(s), cell_id(id), tasks(t) {};
  Host(float d, Symbiont s, int id, vector<int> t) : donation(d), points(0.0), sym(s), cell_id(id), tasks(t) {};
  Host(const Host &orig) : donation(orig.donation), points(orig.points), cell_id(orig.cell_id), sym(orig.sym), tasks(orig.tasks) {};

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
  tasks = parent.tasks;
}

void Host::update(int sym_mult) {
  ///  std::cout << "---------------------------" << endl;
  ///  std::cout << "start host: " << points << " behave: " << donation << endl;
  ///  std::cout << "start sym: " << sym.points << " behave: " << sym.donation << endl;
  float host_start = points;
  float sym_start = sym.points;
  //we need to have a list of resource pools of 25 each
  std::vector<float> pools;
  for (auto resource : resources){
    pools.push_back(25);
  }
  
  if (donation >= 0){
    if (sym.donation == -2) {
      for(int p=0; p<pools.size(); p++){
	if (std::find(tasks.begin(), tasks.end(), p) != tasks.end()){
	  points += pools[p];
	  pools[p] = 0;
	} else{
	  points += (1-donation) * pools[p] * 0.5; //host wastes some
	  pools[p] -= pools[p] * (1-donation);
	}
      }
    } else if (sym.donation>=0){
      for(int p=0; p<pools.size(); p++){
	if (std::find(tasks.begin(), tasks.end(), p) != tasks.end()){
	  points += pools[p];
	  pools[p] = 0;
	} else{
	  points += (1-donation) * pools[p] * 0.5; //host wastes some
	  pools[p] -= pools[p] * (1-donation);
	}
      }
      for(int p=0; p<pools.size(); p++){
	if(std::find(sym.tasks.begin(), sym.tasks.end(), p) != sym.tasks.end()){
	  points += pools[p] * donation * sym.donation;
	  pools[p] -= pools[p] * donation * sym.donation;
	  sym.update(pools[p]);
	  pools[p] = 0;
	}
      }
    } else if (sym.donation<0){
      for(int p=0; p<pools.size(); p++){
	if(std::find(sym.tasks.begin(), sym.tasks.end(), p) != sym.tasks.end()){
	  sym.update(pools[p]);
	  pools[p] = 0;
	}
	if (std::find(tasks.begin(), tasks.end(), p) != tasks.end()){
	  sym.update(pools[p]*(-sym.donation));
	  pools[p] += pools[p]*sym.donation;
	  points += pools[p]*(1-donation);
	  pools[p] = 0;
	}
      }
    } 
  } else if (donation < 0){
    if (sym.donation == -2){
      for(int p=0; p<pools.size(); p++){
	if (!(std::find(tasks.begin(), tasks.end(), p) != tasks.end())){
	  pools[p] += donation*pools[p]; //general resources wasted for speciality
	}else{
	  points += pools[p] * 0.5;
	}
      }
      
    } else if(sym.donation < 0){
      for(int p=0; p<pools.size(); p++){
	if (!(std::find(tasks.begin(), tasks.end(), p) != tasks.end())){
	  pools[p] += donation*pools[p]; //general resources wasted for speciality
	}
      }
      if (sym.donation < donation){
	for(int p=0; p<pools.size(); p++){
	  if (std::find(tasks.begin(), tasks.end(), p) != tasks.end()){
	    pools[p] += pools[p] * (sym.donation-donation); //sym wastes some
	        
	  }if(std::find(sym.tasks.begin(), sym.tasks.end(), p) != sym.tasks.end()){
	    sym.update(pools[p]*(-sym.donation));
	    pools[p] += pools[p] * sym.donation;
	  }
	      
	}
      }else{
	for(int p=0; p<pools.size(); p++){
	  if(std::find(sym.tasks.begin(), sym.tasks.end(), p) != sym.tasks.end()){
	    sym.update(pools[p]*0.5);
	    points += pools[p] * 0.5;
	    pools[p] = 0;
	  }
	}
      }
      for(int p=0; p<pools.size(); p++){
	if (std::find(tasks.begin(), tasks.end(), p) != tasks.end()){
	  points += pools[p];
	  pools[p] = 0;
	}
	if (std::find(sym.tasks.begin(), sym.tasks.end(), p) != sym.tasks.end()){
	  sym.update(pools[p]);
	  pools[p] = 0;
	}
      }
    } else if(sym.donation >=0){
      for(int p=0; p<pools.size(); p++){
	if (std::find(tasks.begin(), tasks.end(), p) != tasks.end()){
	  points += pools[p] * (1+donation);
	  pools[p] = 0;
	}else{
	  pools[p] += pools[p]*donation; //general resources spent on defense
	}
	if (std::find(sym.tasks.begin(), sym.tasks.end(), p) != sym.tasks.end()){
	  pools[p] = pools[p] * 0.5; //sym fighting to get to resources loses some
	  points += pools[p] * sym.donation;
	  pools[p] -= pools[p] * sym.donation;
	  sym.update(pools[p]);
	  pools[p]=0;
	}
      }
    } else{
      cout << "Error in host update, no matching condition for symbiont." << endl;
      
    }
  } else
    cout << "Error in host update, no matching condition for host." << endl;

  ///  std::cout << "diff host: " << points - host_start << endl;                  
  ///  std::cout << "diff sym: " << sym.points - sym_start << endl;
  ///  std::cout << "--------------------------" << endl;  

}

int Host::chooseNeighbor(std::mt19937 &r){

  constexpr int radius = 1;
  emp_assert(radius <= POP_X && radius <= POP_Y);
  emp_assert(cell_id >=0 && cell_id < POP_X * POP_Y, cell_id);
  int cell_x = cell_id % POP_X;
  int cell_y = (cell_id - cell_x)/POP_Y;
  emp::vector<int> neighbor_ids;

  emp_assert(POP_X * POP_Y > 1, POP_X, POP_Y);

  for(int i=cell_x-radius; i<=cell_x+radius; i++){
    const int x = emp::Mod(i, POP_X);

    for(int j=cell_y-radius; j<=cell_y+radius; j++){
      if(i ==0 && j == 0) continue;
      const int y = emp::Mod(j, POP_Y);

      neighbor_ids.push_back((y*POP_X+x));
    }
  }

  emp_assert(neighbor_ids.size(), neighbor_ids.size());
  return(*select_randomly(neighbor_ids.begin(), neighbor_ids.end(), r));
}

void Host::mutate(std::mt19937& r, double rate){

  //Mutate function that moves toward .5 courtesty of Bob
  //std::normal_distribution<double> dist(0.5, rate);
  //double mutation = dist(r); // pull from dist
  //donation = (2*donation)/(1+2*mutation);
  
  std::normal_distribution<double> dist(0, rate);
  double mutation = dist(r);
  donation += mutation;
  if (donation > 1) donation = 1;
  if (donation < -1) donation = -1;

  std::uniform_real_distribution<double> dist2(0, 1);
  double mutation2 = dist2(r);
  if (mutation2<rate){
    tasks.push_back(*select_randomly(resources.begin(), resources.end(), r));
    if(tasks.size() >= host_tasks_lim)
      tasks.erase(std::remove(tasks.begin(), tasks.end(), *select_randomly(tasks.begin(), tasks.end(), r)), tasks.end());
    set<int> temp(tasks.begin(), tasks.end());
    tasks = vector<int>(temp.begin(), temp.end());
  }

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
  bool movie = false;
  std::mt19937 engine;
  std::ofstream sym_map;
  std::ofstream host_map;
  std::ofstream data_file;
  std::ofstream sym_file;
  std::ofstream host_file;
  std::ofstream tasks_file;
  int rec_res = 100;


  Population() : final_update(1) {};
  Population(int pop_count, int f, int seed);

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
  int sym_donate = 10000;
  int host_donate = 100000;

  std::vector<int> sym_dists(21, 0);
  std::vector<int> host_dists(21, 0);
  int res_size = resources.size();
  std::vector<int> tasks(2*res_size, 0);

  if (movie and (cur_update % 1000==0)) {
    host_map << cur_update << ", {";
    sym_map << cur_update << ", {";
  }

  for(auto org : pop){
    host_count++;
    host_sum += org.donation;
    host_donate = (org.donation) * 10 + 10;
    host_dists[host_donate] += 1;
    for(auto task : org.tasks){
      tasks[task] += 1;
    }
    if(org.sym.donation > -2){
      sym_donate = (org.sym.donation * 10) + 10;
      sym_count++;
      sym_sum += org.sym.donation;
      sym_dists[sym_donate] +=1;
      for(auto task : org.sym.tasks) {
	tasks[task+res_size] += 1;
      }
    }
    if (movie and (cur_update % 1000 == 0)){
      if (org.sym.donation != -2) sym_donate = (org.sym.donation * 10) + 10;
      else sym_donate = -100;

      //cout << (org.cell_id) << " " << (org.cell_id % POP_X) << " " << ((org.cell_id % POP_X) == (POP_X -1)) << endl;
      if (org.cell_id % POP_X == 0) {
	host_map << "[";
	sym_map << "[";
      }
      host_map << host_donate << ", ";
      sym_map << sym_donate << ", ";
      if ((org.cell_id % POP_X) == (POP_X-1)){
	host_map <<"],";
	sym_map << "],";
      }
     
    }
  }
  if (movie and (cur_update % 1000==0)) {
    host_map << "] }\n";
    sym_map << "] }\n";
  }

  //cout << cur_update <<", "<< host_sum/host_count << ", " << sym_sum/sym_count << ", " << host_count << ", " << sym_count << endl <<std::flush;
  data_file << cur_update <<", "<< host_sum/host_count << ", " << sym_sum/sym_count << ", " << host_count << ", " << sym_count << endl << std::flush;
  

  sym_file << cur_update << ", ";
  std::copy(sym_dists.begin(), sym_dists.end(), std::ostream_iterator<int>(sym_file, ", "));
  sym_file << endl << std::flush;
  host_file << cur_update << ", ";
  std::copy(host_dists.begin(), host_dists.end(), std::ostream_iterator<int>(host_file, ", "));
  host_file << endl << std::flush;


  tasks_file << cur_update << ", ";
  std::copy(tasks.begin(), tasks.end(), std::ostream_iterator<int>(tasks_file, ", "));
  tasks_file << endl << std::flush;
  
}

Population::Population(int pop_count, int f, int seed_i) {
  final_update = f;
  cur_update = 0;
  seed = seed_i;

  init_pop(pop_count);
}

void Population::init_pop(int pop_count) {

  assert(pop_count == POP_X * POP_Y);
  //This means we have two random number distributions but that's because the state was getting messed up
  //when I had one for the Population object.... so need to fix that someday

  std::mt19937 engine(seed);


  for(int i=0; i<pop_count; ++i){
    vector<int> sym_tasks;
    for(int t= 0; t<sym_tasks_lim; ++t){
      sym_tasks.push_back(*select_randomly(resources.begin(), resources.end(), engine)); 
      //sym_tasks.push_back(1);
    }
    Symbiont new_sym(0, sym_tasks);
    vector<int> host_tasks;
    for(int t=0; t<host_tasks_lim; ++t){
      host_tasks.push_back(*select_randomly(resources.begin(), resources.end(), engine));
      //host_tasks.push_back(0);
    }
    
    Host new_org(0, new_sym, i, host_tasks);

    pop.push_back(new_org);

    
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
  if (movie){
    sym_map.open("sym_map_"+str_seed+"_vert"+str_vert+".map");
    host_map.open("host_map_"+str_seed+"_vert"+str_vert+".map");
    sym_map << "update, state" << endl << std::flush;
    host_map << "update, state" << endl << std::flush;
  }
  tasks_file.open("tasks_"+str_seed+"_mut"+str_mut+"_mult"+str_mult+"_vert"+str_vert+"_start"+str_start+".csv", std::ofstream::ate);
  tasks_file << "Update";
  for(auto resource : resources){
    tasks_file << ", Host_" << resource;
  }
  for (auto resource : resources){
    tasks_file << ", Sym_" << resource;
  }
  tasks_file << endl << std::flush;

  
  std::mt19937 engine(seed);
  std::uniform_real_distribution<double> dist(0, 1);

  //This is to disable mutation after a while to let things settle out
  bool ecological = false;
  for(cur_update = 0; cur_update < final_update+10000; ++cur_update){
    if (!ecological && cur_update >= final_update) ecological = true;
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
	if (!ecological){
        baby.mutate(engine, mut_rate);
        parent.mutate(engine, mut_rate);
	} else {baby.points = 0; parent.points=0;}

        int infected = org.chooseNeighbor(engine);
	if(pop[infected].sym.donation == -2)
	  //cout << "baby sym survives!" << endl;
	  pop[infected].sym = baby;
	
      }
      //See if host reproduces
      if (org.points >=(10*rec_res)){
        //cout << "Making a host baby!" << endl;
        Host baby(org.donation, 0.0, Symbiont(), -2, org.tasks);

        if((org.sym.points >= 0) && org.sym.donation != -2  && (dist(engine) <= vert_rate)){
          //Vertical transmission, TODO: put sym repro in a function
          Symbiont s_baby(org.sym);
	  if (!ecological) {
          s_baby.mutate(engine, mut_rate);
          org.sym.mutate(engine, mut_rate);
	  } else {
	    s_baby.points=0; org.sym.points=0;
	  }
          baby.sym = s_baby;
        }

	//cout << "host baby has sym: " << baby.sym.donation << endl;
        //Mutate both and reset and place
	if (!ecological) {
        baby.mutate(engine, mut_rate);
        org.mutate(engine, mut_rate);
	} else {baby.points=0; org.points=0;}
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
  sym_map.close();
  host_map.close();
}



int main(int argc, char *argv[]) {
  if (argc < 5) cout << "Usage: seed mut vert movie(0/1)" << endl;
  else{
  int seed = atoi(argv[1]);

  Population pop(10000, 500000, seed);
  pop.mut_rate = atof(argv[2]);
  //with division of labor, synergy no longer needs to be artificially enforced
  pop.sym_mult = 1;
  pop.vert_rate = atof(argv[3]);
  if (atoi(argv[4]) == 1){
    pop.movie = true;
    cout << "Making movie files" << endl;
  }else cout << "Not making movie files" << endl;

  pop.evolve();}


}
