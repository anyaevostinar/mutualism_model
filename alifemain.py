#!/usr/bin/env python

import random
import numpy
import sys

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
    if not self.empty:
      if len(genome):
        self.genome = genome
      elif parent:
        newGenome =[]
        for i in range(len(parent.genome)):
          newGenome.append(parent.genome[i])
        self.genome = newGenome
        self.mutate()
        parent.mutate()
        parent.fitness = 0
        parent.age = 0
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
      if len(self.symbionts) > 0 and random.random() < self.donation:
        #symbiont gets this update
        lucky_sym = random.choice(self.symbionts)
        self.fitness += lucky_sym.update()
        
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
    if donation_sd:
      self.donation = numpy.random.normal(self.donation, donation_sd*donation_sd)
    
      
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
    self.genome = genome
    if not self.empty:
      if len(genome):
        self.genome = genome
      elif parent:
        newGenome =[]
        for i in range(len(parent.genome)):
          newGenome.append(parent.genome[i])
        self.genome = newGenome
        self.mutate()
        parent.mutate()
        parent.fitness = 0
        parent.age = 0
      else:
        print "fail"

  def __repr__(self):
    info = "empty: " + str(self.empty) + ", ID: " + str(self.ID) + ", genome: " + str(self.genome) + ", fitness: " + str(self.fitness)
    return info

  def update(self):
    '''Updates the organism's fitness based on its age'''
    if not self.empty:
      self.age += 1
      cur_gene = self.genome[self.age%len(self.genome)]
      self.fitness += cur_gene
      ones = 0
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
    for i in range(popsize):
      self.orgs.append(self.makeOrg())

  def makeOrg(self):
    '''A function to make a new organism randomly'''
    randomBitArray = numpy.random.randint(2, size=(100,))
    newOrg = Organism(len(self.orgs), genome=list(randomBitArray))
    if random.random() < starting_symbiont_proportion:
      #Add a symbiont!
      randomBitArray2 = numpy.random.randint(2, size=(50,))
      newSymbiont = Symbiont(newOrg.ID, genome=list(randomBitArray2))
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
        if not newHost.empty and random.random < (1-vertical_transmission):
          newSymbiont = Symbiont(newHost.ID, parent = symbiont)
          ##currently only letting there be one symbiont per host, but will change this
          if len(newHost.symbionts) >0:
            newHost.symbionts[0] = newSymbiont
          else:
            newHost.symbionts.append(newSymbiont)
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
    if len(org.symbionts) > 0 and random.random() < vertical_transmission:
      parentSymbiont = random.choice(org.symbionts)
      newSymbiont = Symbiont(position, parent=parentSymbiont)
      newOrg.symbionts.append(newSymbiont)
    self.orgs.pop(position)
    self.orgs.insert(position, newOrg)

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
    #print "Symbionts: ", num_symbionts


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
    for org in self.orgs:
      sum_mutualism += org.donation

    return sum_mutualism/len(self.orgs)

  def SymbiontStats(self):
    zeroes_count = 0.0
    symbiont_count = 0
    for org in self.orgs:
      for symbiont in org.symbionts:
        symbiont_count += 1
        for gene in symbiont.genome:
          if gene == 0:
            zeroes_count += 1

    zeroes_avg = zeroes_count/symbiont_count

    return (zeroes_avg, symbiont_count)


if len(sys.argv) == 1 or sys.argv[1] == "--help":
  print "usage: pop_x pop_y seed vert_trans"
else:

  random.seed(int(sys.argv[3]))

  num_updates = 1000
  pop_x = int(sys.argv[1])
  pop_y = int(sys.argv[2])
  pop_size = pop_x*pop_y
  donation = 0.6
  donation_sd = 0.1
  vertical_transmission = float(sys.argv[4])
  starting_symbiont_proportion = 0.5
  data_file = open("mutualism.dat", 'w')
  data_file.write("Update Donation_Avg Num_Symbionts Num_Hosts Avg_Zeroes\n")

  population_orgs = Population(pop_size)
  for i in range(num_updates):
    population_orgs.update()
    if i%100 == 0:
      avg_mut = population_orgs.findAverageMutualism()
      zeroes_avg, symbiont_count = population_orgs.SymbiontStats()
      data_file.write(str(i)+" "+str(avg_mut)+" "+str(symbiont_count)+" "+str(len(population_orgs.orgs)) +" "+str(zeroes_avg)+"\n")


  data_file.close()
