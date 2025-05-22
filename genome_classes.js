class Genome {

  constructor() {
    this.uid = genomeIds.next()
    this.chromosome = []
    shuffle(hk_genes)
    for (let i = 0; i < init_hk; i++)this.chromosome.push(new Gene(hk_genes[i%hk_genes.length]))
    for (let i = 0; i < init_nc; i++) this.chromosome.push(new Gene('.'))
    shuffle(this.chromosome)
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
    child.calculate_fitness()
    child.mutate()

    return child
  }

  calculate_fitness() {
    this.fitness = 1
    let check = JSON.parse(JSON.stringify(hk_genes))

    for (let i = 0; i < this.chromosome.length; i++){
      var idx = check.indexOf(this.chromosome[i].type)
      if (idx >= 0) check.splice(idx, 1)
    }
    this.genome_length = this.chromosome.length
    this.dupSites = Array(this.chromosome.length).fill(0)
    for (let i = 0; i < this.phages.length; i++){
      this.genome_length = this.genome_length + this.phages[i].length
      this.dupSites.push(0)
      this.dupSites.concat(Array(this.phages[i].length+1).fill(i+1))
      for (let j = 0; j < this.phages[i].length; j++) {
        var idx = check.indexOf(this.phages[i][j].type)
        if (idx >= 0) check.splice(idx, 1)
      }
    }

    if (check.length == 0 && this.susceptible == 1) this.fitness = b0 - this.genome_length*seg_cost
    else this.fitness = 0

  }

  mutate() {
    for (let i = 0; i < this.chromosome.length; i++){
      var randomnr = sim.rng.genrand_real1()
      if (randomnr < gene_deletion_rate) {
        this.chromosome.splice(i, 1)
      }
      else if (randomnr < gene_deletion_rate + gene_inactivation_rate) {
        this.chromosome[i].type = "."
      }
      else if (randomnr < gene_deletion_rate + gene_inactivation_rate + gene_duplication_rate) {
        var insertion = sim.rng.genrand_int(0,this.dupSites.length-1)
        if (this.dupSites[insertion] == 0) this.chromosome.push(this.chromosome[i].copy())
        else this.phages[this.dupSites[insertion]-1].push(this.chromosome[i].copy())
      }
    }
    for (let i = 0; i < this.phages.length; i++){
      for (let j = 0; j < this.phages[i].length; j++){
        var randomnr = sim.rng.genrand_real1()
        if (this.phages[i][j].type != "L" && randomnr < gene_deletion_rate) {
          this.phages[i].splice(j, 1)
        }
        else if (this.phages[i][j].type != "L" && randomnr < gene_deletion_rate + gene_inactivation_rate) {
          this.phages[i][j].type = "."
        }
        else if (this.phages[i][j].type != "L" && randomnr < gene_deletion_rate + gene_inactivation_rate + gene_duplication_rate) {
          var insertion = sim.rng.genrand_int(0,this.dupSites.length-1)
          if (this.dupSites[insertion] == 0) this.chromosome.push(this.phages[i][j].copy())
          else this.phages[this.dupSites[insertion]-1].push(this.phages[i][j].copy())
        }
      }
      var randomnr = sim.rng.genrand_real1()
      if (randomnr < phage_deletion_rate){
        this.phages.splice(i, 1)
      }
      else if (randomnr < phage_deletion_rate + phage_duplication_rate){
        let new_phage = []
        for (let j = 0; j < this.phages[i].length; j++) new_phage.push(this.phages[i][j].copy())
        this.phages.push(new_phage)
      }
    }
    this.calculate_fitness()
  }

}

class Virion {

  constructor(type, parent) {
    this.genome = []
    if (type) {
      this.genome.push(new Gene("L"))
      this.genome.push(new Gene("R"))
    }
    else for (let i = 0; i < parent.length; i++) this.genome.push(parent[i].copy())
  }
}

class Gene {

  constructor(type) {
    this.type = type
    if (['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'R', '.'].includes(type)) {
      this.lock = []
      for (let i = 0; i < bit_string_length; i++) this.lock.push(Math.floor(sim.rng.genrand_real1()*2))
    }
    else if (type == "L") {
      this.key = []
      for (let i = 0; i < bit_string_length; i++) this.key.push(Math.floor(sim.rng.genrand_real1()*2))
      this.lysis = sim.rng.genrand_real1() * max_induction_rate
    }
    if (type == "R") {
      this.highX_state = 1
      this.lowX_state = 1
    }
  }

  copy() {
    let new_gene = new Gene(this.type)

    if (['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'R', '.'].includes(this.type)) {
      new_gene.lock = []
      for(let i = 0; i < bit_string_length; i++) {
        if (sim.rng.genrand_real1() < coevolution_rate) new_gene.lock.push(1 - this.lock[i])
        else new_gene.lock.push(this.lock[i])
      }
    }

    else if (this.type == "L") {
      new_gene.key = []
      for(let i = 0; i < bit_string_length; i++) {
        if (sim.rng.genrand_real1() < coevolution_rate) new_gene.key.push(1 - this.key[i])
        else new_gene.key.push(this.key[i])
      }

      if (sim.rng.genrand_real1() < lysis_mutation_rate) new_gene.lysis = sim.rng.genrand_real1() * max_induction_rate
      else new_gene.lysis = this.lysis
    }

    if (this.type == "R") {
      var rnr = sim.rng.genrand_real1()
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

    return new_gene
  }

}

function shuffle(a) {
  if (a.length > 1) {
    var j, x, i;
    for (i = a.length - 1; i > 0; i--) {
      j = Math.floor(sim.rng.genrand_real2() * (i + 1));
      x = a[i];
      a[i] = a[j];
      a[j] = x;
    }
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
