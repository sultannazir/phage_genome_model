class Genome {

  constructor(init_nc, b0, seg_cost) {
    this.uid = genomeIds.next()
    //initialise a genome with init_nc many non-coding segments and one attachment region
    this.chromosome = []
    this.chromosome.push(new Gene("att"))
    for (let i = 0; i < init_nc; i++) this.chromosome.push(new Gene("."))
    this.dupSites = Array(this.chromosome.length).fill(0)
    this.phages = []
    this.genome_length = this.chromosome.length
    this.fitness = b0 - this.genome_length*seg_cost

    this.susceptible = 1
  }

  copy() {
    let child = new Genome()
    child.chromosome = []
    for (let i = 0; i < this.chromosome.length; i++) child.chromosome.push(this.chromosome[i].copy())
    child.phages = []
    for (let i = 0; i < this.phages.length; i++){
      child.phages.push([])
      for (let j = 0; j < this.phages[i].length; j++) child.phages[i].push(this.phages[i][j].copy())
    }
    child.genome_length = this.genome_length
    child.dupSites = this.dupSites
    child.fitness = this.fitness
    child.susceptible = this.susceptible

    child.mutate()

    return child
  }

  calculate_fitness() {
    this.susceptible = 0
    for (let i = 0; i < this.chromosome.length; i++){
      if (this.chromosome[i].type == "att") this.susceptible = 1
    }
    if (this.phages.length == 0 && this.susceptible == 0) this.fitness = 0
    else{
      this.genome_length = this.chromosome.length
      this.dupSites = Array(this.chromosome.length).fill(0)
      for (let i = 0; i < this.phages.length; i++){
        this.genome_length = this.genome_length + this.phages[i].length
        this.dupSites.push(0)
        this.dupSites.concat(Array(this.phages[i].length+1).fill(i+1))
      }
      this.fitness = b0 - this.genome_length*seg_cost
    }
  }

  mutate() {
    var mutation = false
    var frac = (this.chromosome.length + this.phages.length)/this.num_dupSites
    for (let i = 0; i < this.chromosome.length; i++){
      var randomnr = sim.rng.genrand_real1()
      if (randomnr < gene_deletion_rate) {
        this.chromosome.splice(i, 1)
        mutation = true
      }
      else if (randomnr < gene_deletion_rate + gene_inactivation_rate) {
        this.chromosome[i].type = "."
        mutation = true
      }
      else if (randomnr < gene_deletion_rate + gene_inactivation_rate + gene_duplication_rate) {
        var insertion = sim.rng.genrand_int(0,this.dupSites.length-1)
        if (this.dupSites[insertion] == 0) this.chromosome.push(this.chromosome[i].copy())
        else if (this.chromosome[i].type != "att") this.phages[this.dupSites[insertion]-1].push(this.chromosome[i].copy())
        mutation = true
      }
    }
    var frac = (this.chromosome.length + this.phages.length)/this.num_dupSites
    for (let i = 0; i < this.phages.length; i++){
      for (let j = 0; j < this.phages[i].length; j++){
        var randomnr = sim.rng.genrand_real1()
        if (randomnr < gene_deletion_rate) {
          this.phages[i].splice(j, 1)
          mutation = true
        }
        else if (randomnr < gene_deletion_rate + gene_inactivation_rate) {
          this.phages[i][j].type = "."
          mutation = true
        }
        else if (this.phages[i][j].type != "L" && randomnr < gene_deletion_rate + gene_inactivation_rate + gene_duplication_rate*frac) {
          var insertion = sim.rng.genrand_int(0,this.dupSites.length-1)
          if (this.dupSites[insertion] == 0) this.chromosome.push(this.phages[i][j].copy())
          else this.phages[this.dupSites[insertion]-1].push(this.phages[i][j].copy())
          mutation = true
        }
      }
      var randomnr = sim.rng.genrand_real1()
      if (randomnr < phage_deletion_rate){
        this.phages.splice(i, 1)
        mutation = true
      }
      else if (randomnr < phage_duplication_rate){
        let new_phage = []
        for (j = 0; j < this.phages[i].length; j++) new_phage.push(this.phages[i][j].copy())
        this.phages.push(new_phage)
        mutation = true
      }
    }
    if (mutation) this.calculate_fitness()
  }

}

class Virion {

  constructor() {
    this.uid = virionIds.next()
    this.genome = []
    this.genome.push(new Gene("L"))
    this.genome.push(new Gene("R"))
  }
}

class Gene {

  constructor(type) {
    this.uid = geneIds.next()
    this.type = type
    if (type == "att") this.lock = sim.rng.genrand_real1()
    else if (type == "L") this.key = sim.rng.genrand_real1()
    else if (type == "R") {
      this.highX_state = 1
      this.lowX_state = 1
    }
  }

  copy() {
    let new_gene = new Gene(this.type)
    new_gene.mutate()
    return new_gene
  }

  mutate() {
    if (this.type == "att") {
      var rnr = sim.genrand_real1()
      if (rnr < coevolution_rate) new_gene.lock = ((this.lock + coevolution_step*(sim.rng.genrand_real1()*2-1)) + 1) % 1
      else new_gene.lock = this.lock
    }

    else if (this.type == "L") {
      var rnr = sim.genrand_real1()
      if (rnr < coevolution_rate) new_gene.key = ((this.key + coevolution_step*(sim.rng.genrand_real1()*2-1)) + 1) % 1
      else new_gene.key = this.key
    }

    else if (this.type == "R") {
      var rnr = sim.genrand_real1()
      if (rnr < 0.5*gene_regulation_mutation_rate) {
        new_gene.lowX_state = this.lowX_state
        new_gene.highX_state = 1 - this.highX_state
      }
      else if (rnr < gene_regulation_mutation_rate) {
        new_gene.lowX_state = 1 - this.lowX_state
        new_gene.highX_state = this.highX_state
      }
      else {
        new_gene.lowX_state = this.lowX_state
        new_gene.highX_state = this.highX_state
      }
    }
  }

}

function shuffle(a) {
  var j, x, i;
  for (i = a.length - 1; i > 0; i--) {
    j = Math.floor(sim.rng.genrand_real2() * (i + 1));
    x = a[i];
    a[i] = a[j];
    a[j] = x;
  }
  return a;
}

function* idGenerator() {
  let id = 1;
  while (true) {
    yield id
    id ++
  }
}

const genomeIds = idGenerator()
const virionIds = idGenerator()
const geneIds = idGenerator()
