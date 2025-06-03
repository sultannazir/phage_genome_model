class Genome {

  constructor() {
    this.uid = genomeIds.next()
    this.chromosome = []
    shuffle(hk_genes)
    for (let i = 0; i < init_hk; i++)this.chromosome.push(new Gene(hk_genes[i%hk_genes.length]))
    for (let i = 0; i < init_nc; i++) this.chromosome.push(new Gene('.'))
    var rnr = sim.rng.genrand_real1()
    if (rnr < init_R_prob) this.chromosome.push(new Gene('R'))
    shuffle(this.chromosome)
    this.phages = []
    this.fitness = b0
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
    child.mutate()

    return child
  }

  calculate_fitness() {
    let check = JSON.parse(JSON.stringify(hk_genes))
    for (let i = 0; i < this.chromosome.length; i++){
      var idx = check.indexOf(this.chromosome[i].type)
      if (idx >= 0) check.splice(idx, 1)
    }
    for (let i = 0; i < this.phages.length; i++){
      for (let j = 0; j < this.phages[i].length; j++) {
        var idx = check.indexOf(this.phages[i][j].type)
        if (idx >= 0) check.splice(idx, 1)
      }
    }
    if (check.length > 0) this.fitness = 0
  }

  mutate() {
    for (let i = 0; i < this.chromosome.length; i++){
      var randomnr = sim.rng.genrand_real1()
      if (randomnr < att_deletion_rate) {
        this.chromosome.splice(i, 1)
      }
      else if (randomnr < att_deletion_rate + att_inactivation_rate) {
        this.chromosome[i].type = "."
      }
      else if (randomnr < att_deletion_rate + att_inactivation_rate + att_duplication_rate) {
        this.chromosome.push(this.chromosome[i].copy())
      }
    }
    for (let i = 0; i < this.phages.length; i++){
      for (let j = 0; j < this.phages[i].length; j++){
        var randomnr = sim.rng.genrand_real1()
        if (this.phages[i][j].type != "L" && randomnr < att_deletion_rate) {
          this.phages[i].splice(j, 1)
        }
        else if (this.phages[i][j].type != "L" && randomnr < att_deletion_rate + att_inactivation_rate) {
          this.phages[i][j].type = "."
        }
        else if (this.phages[i][j].type != "L" && randomnr < att_deletion_rate + att_inactivation_rate + att_duplication_rate) {
          this.phages[i].push(this.phages[i][j].copy())
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

  constructor(type, parent, chromosome) {
    this.genome = []
    if (type) {
      this.genome.push(new Gene("L"))
      if (sim.rng.genrand_real1() < init_R_prob) this.genome.push(new Gene("R"))
    }
    else {
      for (let i = 0; i < parent.length; i++) this.genome.push(parent[i].copy())
      if (sim.rng.genrand_real1() < excision_error_rate) {
        var rand_att = sim.rng.genrand_int(0,chromosome.length-1)
        this.genome.push(chromosome[rand_att].copy())
      }
    }
  }

}

class Gene {

  constructor(type) {
    this.type = type
    if (['A', 'B', 'C', 'R', '.'].includes(type)) {
      this.lock = []
      for (let i = 0; i < bit_string_length; i++) this.lock.push(Math.floor(sim.rng.genrand_real1()*2))
    }
    else if (type == "L") {
      this.key = []
      for (let i = 0; i < bit_string_length; i++) this.key.push(Math.floor(sim.rng.genrand_real1()*2))
      this.lysis = sim.rng.genrand_real1() * max_induction_rate
    }
    if (type == "R") this.Xhf = sim.rng.genrand_real1()
  }

  copy() {
    let new_gene = new Gene(this.type)

    if (['A', 'B', 'C', 'R', '.'].includes(this.type)) {
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
      var rnr = sim.rng.genrand_real1()
      if (rnr < lysis_mutation_rate) new_gene.lysis = sim.rng.genrand_real1() * max_induction_rate
      else new_gene.lysis = this.lysis
    }

    if (this.type == "R") {
      var rnr = sim.rng.genrand_real1()
      if (rnr < Xhf_mutation_rate) new_gene.Xhf = sim.rng.genrand_real1()
      else new_gene.Xhf = this.Xhf
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
