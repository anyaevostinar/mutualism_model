#!/usr/bin/env python

import random
import numpy
import sys
import argparse

class Organism:
  '''A class to contain an organism'''
  def __init__(self, cellID, genome=[], parent=False, empty=False):
    self.age = 0
    self.empty = empty
    self.ID = cellID
    self.fitness = 0
    self.genome = genome
    self.symbionts = []
    self.donation = donation
    self.vert_trans = vertical_transmission
    if not self.empty:
      if len(genome):
        self.genome = genome
      elif parent:
        newGenome =[]
        for i in range(len(parent.genome)):
          newGenome.append(parent.genome[i])
        self.genome = newGenome
        self.vert_trans = parent.vert_trans
        self.donation = parent.donation
        if args.host_evolve == 1:
          self.mutate()
          parent.mutate()
        parent.fitness = 0
        parent.age = 0
        for symbiont in parent.symbionts:
          symbiont.reset()
      else:
        print "fail"

  def __repr__(self):
    info = "empty: " + str(self.empty) + ", ID: " + str(self.ID) + ", genome: " + str(self.genome) + ", fitness: " + str(self.fitness) + "\n\tSymbionts:"
    for symbiont in self.symbionts:
      info+="\n\t"+symbiont.__repr__()+"\n"
    return info

  def update(self):
    '''Updates the organism's fitness based on its age'''
    if not self.empty:
      if len(self.symbionts) > 0 and (random.random() < self.donation):
        #symbiont gets this update
        lucky_sym = random.choice(self.symbionts)

        self.fitness += lucky_sym.update()
#        print "New Host",
#        print self.fitness
        assert(self.fitness >= lucky_sym.symbiont_bonus)
        
        return True
      else:
        self.age += 1
        cur_gene = self.genome[self.age%len(self.genome)]
        self.fitness += cur_gene
        return True
    else:
      return False
      

  def mutate(self):
    if random.random() < .02:
      newGenome = self.genome
      flipBit = random.randint(0, len(newGenome)-1)
      newGenome[flipBit] = random.randint(0,1)
      self.genome = newGenome
    if vert_sd:
      self.vert_trans = numpy.random.normal(self.vert_trans, vert_sd)
    if donation_sd:
      self.donation = numpy.random.normal(self.donation, donation_sd)
    
      
  def findNeighbors(self):
    cellID = self.ID
    radius = 1
    world_x = pop_x
    world_y = pop_y
    cell_x = cellID % world_x
    cell_y = (cellID - cell_x)/world_x
    neighbor_ids = []
    for i in range(cell_x-radius, cell_x+radius+1):
      for j in range(cell_y-radius, cell_y+radius+1):
        if i<0:
          x = world_x + i
        elif i>=world_x:
          x = i-world_x
        else:
          x = i

        if j<0:
          y = world_y + j
        elif j>=world_y:
          y = j-world_y
        else:
          y = j

        neighbor_ids.append(y*world_x+x)

    return neighbor_ids

class Symbiont:
  '''A class to contain a symbiont.
  TODO: make super class with symbiont and host children'''
  def __init__(self, cellID, genome=[], parent=False, empty=False):
    self.age = 0
    self.ID = cellID
    self.empty = empty
    self.fitness = 0
    self.symbiont_bonus = 0
    self.genome = genome
    if not self.empty:
      if len(genome):
        self.genome = genome
      elif parent:
        self.genome = list(parent.genome)
        if args.sym_evolve == 1:
          self.mutate()
          parent.mutate()
        parent.fitness = 0
        parent.age = 0
        parent.symbiont_bonus = 0
      else:
        print "fail"

  def reset(self):
    '''Resets symbiont after host reproduction.'''
    self.age = 0
    self.fitness = 0
    self.symbiont_bonus = 0

  def __repr__(self):
    info = "empty: " + str(self.empty) + ", ID: " + str(self.ID) + ", genome: " + str(self.genome) + ", fitness: " + str(self.fitness)
    return info

  def update(self):
    '''Updates the organism's fitness based on its age'''
    if not self.empty:
      self.age += 1
      cur_gene = self.genome[self.age%len(self.genome)]
      self.fitness += cur_gene
      ones = cur_gene
      zeroes = 0
      if cur_gene == 0:
        for g in range(self.age%len(self.genome), 0, -1):
          if self.genome[g] == 0:
            zeroes +=1
          else:
            break
        self.fitness += zeroes
      elif cur_gene == 1:
        for h in range(self.age%len(self.genome), 0, -1):
          if self.genome[h] == 1:
            ones +=1
          else:
            break



      self.symbiont_bonus += ones
      return ones

      

  def mutate(self):
    if random.random() < .02:
      newGenome = self.genome
      flipBit = random.randint(0, len(newGenome)-1)
      newGenome[flipBit] = random.randint(0,1)
      self.genome = newGenome
    


        
class Population:
  '''A class to contain the population and do stuff'''
  def __init__(self, popsize):
    self.currentUpdate = 0
    self.orgs = []
    self.pop_size = popsize
    self.vert_trans_count = 0
    self.hor_trans_count = 0
    

    for i in range(popsize):
      self.orgs.append(self.makeOrg())

  def makeOrg(self):
    '''A function to make a new organism randomly'''
    if args.host_evolve == 1:
      randomBitArray = list(numpy.random.randint(2, size=(100,)))
    else:
      randomBitArray = [0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1]
    newOrg = Organism(len(self.orgs), genome=randomBitArray)
    if random.random() < starting_symbiont_proportion:
      #Add a symbiont! 
      if args.sym_evolve == 0:
        #Can't be random since we want every one the same
        randomBitArray2 = [1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1]
      else:
        randomBitArray2 = list(numpy.random.randint(2, size=(50,)))
      #Making symbiont more mutualistic or parasitic based on setting
      if symbiont_mutualism > 0:
        randomBitArray[0:symbiont_mutualism]=1
      elif symbiont_mutualism < 0:
        reverse = -symbiont_mutualism
        randomBitArray[0:reverse]=0
      newSymbiont = Symbiont(newOrg.ID, genome=randomBitArray2)
      newOrg.symbionts.append(newSymbiont)
    return newOrg

  def checkSymbionts(self, org):
    '''A helper function to check if the symbionts need to be reproduced.'''
    ## do the orgs symbionts get to reproduce?
    for symbiont in org.symbionts:
      if symbiont.fitness >= len(symbiont.genome):
        #Symbiont reproduction! Horizontal Transmission.
        newHostID = random.choice(org.findNeighbors())
        newHost = self.orgs[newHostID]
        if not newHost.empty:
          newSymbiont = Symbiont(newHost.ID, parent = symbiont)
          ##currently only letting there be one symbiont per host, but will change this
          if len(newHost.symbionts) >0:
            assert len(newHost.symbionts) == 1, "There should be only one symbiont"
            if random.random() < 0.5:
              #Symbiont has a 50/50 shot of repelling new symbiont
              self.hor_trans_count += 1
              newHost.symbionts[0] = newSymbiont
          else:
            newHost.symbionts.append(newSymbiont)
            self.hor_trans_count +=1
    return len(org.symbionts)
    
  def reproduceOrg(self, org):
    '''A helper function to reproduce a host organism.'''
    ## Time to reproduce
    ## Returns list of neighbor indices
    dead_neighbor = False
    neighbors = org.findNeighbors()
    for ID in neighbors:
      if self.orgs[ID].empty:
        dead_neighbor = ID
        break
    if dead_neighbor:
      position = dead_neighbor
    else:
      position = random.choice(neighbors)
    newOrg = Organism(position, parent = org)
    ##Vertical transmission of the symbiont time
    if len(org.symbionts) > 0 and (random.random() < org.vert_trans):
      self.vert_trans_count += 1
      parentSymbiont = random.choice(org.symbionts)
      newSymbiont = Symbiont(position, parent=parentSymbiont)
      newOrg.symbionts.append(newSymbiont)

    self.orgs[position] = newOrg

  def update(self):
    '''A function that runs a single update'''
    self.currentUpdate+=1
    num_symbionts = 0
    for org in self.orgs:
      if not org.empty:
        result = org.update()
        num_symbionts += self.checkSymbionts(org)
        if org.fitness >= len(org.genome):
          self.reproduceOrg(org)


  def findBest(self, orgs_to_eval=False):
    '''Finds the best of the population or a provided subset'''
    if not orgs_to_eval:
      orgs_to_eval = self.orgs
    highest_fitness = 0
    fittest_org = False
    for org in orgs_to_eval:
      if org.fitness > highest_fitness:
        highest_fitness = org.fitness
        fittest_org = org
  
    if not fittest_org:
      print "Error! No Org selected!"
    return fittest_org

  def findBestSymbiont(self):
    orgs_to_eval = self.orgs
    highest_fitness = 0
    fittest_org = False
    for org in orgs_to_eval:
      if org.symbionts[0].fitness > highest_fitness:
        highest_fitness = org.symbionts[0].fitness
        fittest_symbiont = org
  
    if not fittest_symbiont:
      print "Error! No Org selected!"
    return fittest_symbiont

  def findAverageMutualism(self):
    sum_mutualism = 0.0
    sum_fitness = 0.0
    org_count = 0
    for org in self.orgs:
      sum_mutualism += org.vert_trans
      sum_fitness += org.fitness
      org_count += 1

    mutualism_avg = sum_mutualism/len(self.orgs)
    fitness_avg = sum_fitness/len(self.orgs)
    return mutualism_avg, fitness_avg

  def SymbiontStats(self):
    bonus_sum = 0.0
    symbiont_count = 0

    for org in self.orgs:
      for symbiont in org.symbionts:
        symbiont_count += 1
        if symbiont.symbiont_bonus != 0:
          bonus_sum += float(symbiont.symbiont_bonus)/float(org.fitness)

    if symbiont_count == 0:
      bonus_avg = 0
    else:
      bonus_avg = bonus_sum/symbiont_count

    return (bonus_avg, symbiont_count)

def parseArgs(parser):
  '''A helper function for parsing the possible command line arguments.'''
  parser.add_argument('pop_x', help="the x value for the population size", type=int)
  parser.add_argument('pop_y', help='the y value for the population size', type=int)
  parser.add_argument('seed', help='the seed for the random number generator', type=int)
  parser.add_argument('--donation', help="the value the host donation rate should be set to, -1 for evolvable, default is 0.6", type=float, default=0.6)
  parser.add_argument("--vert_trans", help="the value the vertical transmission rate should be set to, -1 for host evolvable, default is 0.75", type=float, default=0.75)
  parser.add_argument("--sym_mut", help="how mutualistic the symbiont should start, ie how many sites at the beginning of the genome should be forced to 1 to help the host. Negative values make it more parasitic instead. Default is 0", type=int, default=0)
    #might change these to mutation rates in the future:
  parser.add_argument("--host_evolve", help="Should the host's genome be able to evolve? 0 for no, 1 for yes, default is 1", type=int, default =1)
  parser.add_argument('--sym_evolve', help="Should the symbiont's genome be able to evolve? 0 for no, 1 for yes, default is 1", type=int, default=1) 
  parser.add_argument("--num_updates", help="the number of updates the experiment should run for, default is 20000", type=int, default=20000)
  parser.add_argument("filename", help="the prefix of the filename the data should be saved to")

  return parser.parse_args()

parser = argparse.ArgumentParser(description="Simulation the evolution of hosts and symbionts")
args = parseArgs(parser)

seed = args.seed
random.seed(seed)
numpy.random.seed(seed)
num_updates = args.num_updates
pop_x = args.pop_x
pop_y = args.pop_y
pop_size = pop_x*pop_y
donation = args.donation
if donation == -1:
  donation = 0.6
  donation_sd = 0.1
else:
  donation_sd = 0.0
  vertical_transmission = args.vert_trans
if vertical_transmission == -1:
  vertical_transmission = 0.75
  vert_sd = 0.1
else:
  vert_sd = 0.0
filename = args.filename
#This is how many initial positions will be set to 1, can be negative to indicate set to 0 instead!
symbiont_mutualism = args.sym_mut

starting_symbiont_proportion = 0.5
data_file = open(filename+str(seed)+".dat", 'w')
data_file.write("Update Vert_Avg Num_Symbionts Num_Hosts Symbiont_Bonus_Prop Host_Fitness Num_Vert_Trans Num_Hor_Trans\n")

population_orgs = Population(pop_size)
for i in range(num_updates):
  population_orgs.update()
  if i%100 == 0:
    #      print "Update: ", i
    avg_mut, avg_fit = population_orgs.findAverageMutualism()
    bonus_avg, symbiont_count = population_orgs.SymbiontStats()
    data_file.write('{} {} {} {} {} {} {} {}\n'.format(i,avg_mut, symbiont_count,len(population_orgs.orgs), bonus_avg, avg_fit, population_orgs.vert_trans_count, population_orgs.hor_trans_count))
    population_orgs.vert_trans_count = 0
    population_orgs.hor_trans_count = 0

data_file.close()
