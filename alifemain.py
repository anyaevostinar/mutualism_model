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
      if cur_gene == 1 or cur_gene == 0:
        self.fitness += cur_gene
        return True
      

  def mutate(self):
    if random.random() < .02:
      newGenome = self.genome
      flipBit = random.randint(0, len(newGenome)-1)
      newGenome[flipBit] = random.randint(0,1)
      self.genome = newGenome
    
      
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
    return newOrg

  def update(self):
    '''A function that runs a single update'''
    self.currentUpdate+=1
#    current_loc = self.currentUpdate%len(self.orgs[0].genome)
    for org in self.orgs:
      if not org.empty:
        result = org.update()
        if org.fitness >= len(org.genome):
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
            position = random.choice(org.findNeighbors())
          newOrg = Organism(position, parent = org)
          self.orgs.pop(position)
          self.orgs.insert(position, newOrg)


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


random.seed(int(sys.argv[3]))

num_updates = 1000
pop_x = int(sys.argv[1])
pop_y = int(sys.argv[2])
pop_size = pop_x*pop_y

population_orgs = Population(pop_size)
for i in range(num_updates):
  population_orgs.update()

print population_orgs.findBest()
