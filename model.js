Simulation = require('./cacatoo.js')
let yargs = require('yargs')
let cmd_params = yargs.argv

let iter = cmd_params.seed
let VD = cmd_params.VD

let speed = 1
// Parameters
let num_start_cells = 1000
let max_cells = 50000
let init_nc = 5
let init_hk = 5
let init_R_prob = 0.1
let init_Xhf = 0.5
let b0 = 0.0125
let repression_cost = 0.0025

let X_uptake = 0.1
let Y_uptake = 0.1
let X_max_influx = 0.01
let Y_max_influx = 0.01
let X_decay = 0.01
let Y_decay = 0.01

let d0 = 0.0025

let att_deletion_rate = 0.0001
let att_duplication_rate = 0.0001
let att_inactivation_rate = 0.001
let phage_deletion_rate = 0.0001
let phage_duplication_rate = 0.0001
let Xhf_mutation_rate = 0.0001
let lysis_mutation_rate = 0.0001
let excision_error_rate = 0.001
let coevolution_rate = 0.001

let bit_string_length = 20

let env_fluctuation_freq = 0.001

let specificity = 10
let virion_decay_rate = 0.001
let infection_rate = 0.1
let virion_diffusion_rate = 10**(-1*VD)

let bacteria_influx_rate = 0.001
let virion_influx_rate = 1

let max_induction_rate = 0.1
let burst_size = 30
let vir_cost = 1

let hk_genes = ['A', 'B', 'C']

var size = 80

// Configuration
let config = {
    title: "Evolution of prophage genome",
    description: "Virtual experiment of genome evolution in a bacteria-prophage system",
    maxtime: 1000000,
    seed: iter,
    skip: 0,
    ncol: size,
    nrow: size,		            // dimensions of the grid to build
    wrap: [false, false],       // Wrap boundary [COLS, ROWS]
    graph_interval: 1,
    fpsmeter: true
}

let flockconfig = {
  // Flock parameters
  num_boids: 0,            // Starting number of boids (flocking individuals)
  shape: 'dot',            // Shape of the boids drawn (options: bird, arrow, line, rect, dot, ant)
  click: 'none',          // Clicking the boids pushes them away from the mouse
  max_speed: 1.0,          // Maximum velocity of boids
  max_force: 1.0,          // Maximum steering force applied to boids (separation/cohesion/alignment rules)
  init_velocity:0,
  friction: 0.999,
  brownian: 0.05,
  // Mouse parameters
  mouse_radius: 20,        // Radius of boids captured by the mouse overlay
  draw_mouse_radius: 'false',    // Show a circle where the mouse is
  draw_mouse_colour: 'white',
  // Collision behaviour
  collision_force: 0.1,
  cohesion: {strength:0.0, radius:4.5},
  size: 1,                // Size of the boids (scales drawing and colision detection)

  // Optimalisation (speed) parameters
  //qt_colour: "white",       // Show quadtree (optimalisation by automatically tessalating the space)
  qt_capacity: 3        // How many boids can be in one subspace of the quadtree before it is divided further
}

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



sim = new Simulation(config)

sim.makeGridmodel("environment");

sim.initialGrid(sim.environment, "X", 0, 1)
sim.initialGrid(sim.environment, "Y", 0, 1)

for (let x = 0; x < sim.ncol; x++) {
  for (let y = 0; y < sim.nrow; y++) {
    sim.environment.grid[x][y].X = x/size
    sim.environment.grid[x][y].Y = 0
    sim.environment.grid[x][y].virions = []
    sim.environment.grid[x][y].density = 0
  }
}

sim.makeFlockmodel("flock", flockconfig)

sim.flock.boids = []
sim.flock.populateSpot(num_start_cells,sim.nr/2,sim.nc/2)
for (let boid of sim.flock.boids) {
  boid.genome = new Genome()
  boid.X = sim.environment.grid[Math.floor(boid.position.x)][Math.floor(boid.position.y)].X
  boid.Y = sim.environment.grid[Math.floor(boid.position.x)][Math.floor(boid.position.y)].Y
  boid.fitness = boid.genome.fitness
  boid.repression = false
  for (let i = 0; i < boid.genome.chromosome.length; i++) {
    if (boid.genome.chromosome[i].type == "R-") {
      if (boid.X < boid.genome.chromosome[i].Xhf) boid.repression = true
    }
  }
}

sim.steps = 0

sim.flock.update = function() {
  for (let i = 0; i < speed; i++) {
    sim.steps ++
    shuffle(this.boids)
    let newboids = []
    for (let boid of this.boids) {
      let death = false

      let x = Math.floor(boid.position.x);
      let y = Math.floor(boid.position.y);

      // substrate uptake
      for (let gp of sim.flock.getNearbyGridpoints(boid, sim.environment, boid.size + 2,)) {
        let v = gp.Y * Y_uptake
        gp.Y -= v
        boid.Y += v
      }
      boid.Y *= 1 - Y_decay
      boid.X = sim.environment.grid[Math.floor(boid.position.x)][Math.floor(boid.position.y)].X


      // infection
      if (sim.environment.grid[x][y].virions.length > 0) {
        shuffle(sim.environment.grid[x][y].virions)
        for (let vir = 0; vir < sim.environment.grid[x][y].virions.length; vir++) {
          if (sim.rng.genrand_real1() < infection_rate) {
            let integration = false
            shuffle(boid.genome.chromosome)
            for (let gnum = 0; gnum < boid.genome.chromosome.length && !integration; gnum++) {
              let key = sim.environment.grid[x][y].virions[vir].genome[0].key
              let lock = boid.genome.chromosome[gnum].lock
              let diff = 0
              for (let aa = 0; aa < bit_string_length; aa++) diff += Math.abs(lock[aa]-key[aa])
              let probability = (1 - diff/bit_string_length)**specificity
              if (sim.rng.genrand_real1() < probability) {
                integration = true
                boid.genome.chromosome.splice(gnum, 1)
                boid.genome.phages.push(sim.environment.grid[x][y].virions[vir].genome)
                boid.genome.calculate_fitness()
              }
            }
            sim.environment.grid[x][y].virions.splice(vir, 1)
          }
        }
      }

      // calculate fitness and repression state

      boid.repression = false
      boid.fitness = boid.genome.fitness
      for (let gnum = 0; gnum < boid.genome.chromosome.length && boid.fitness > 0; gnum++) {
        if (boid.genome.chromosome[gnum].type == "R"){
          if (boid.X <= boid.genome.chromosome[gnum].Xhf) {
            boid.repression = true
            boid.fitness -= repression_cost
          }
        }
      }
      for (let phg = 0; phg < boid.genome.phages.length; phg++) {
        for (let gnum = 0; gnum < boid.genome.phages[phg].length; gnum++) {
          if (boid.genome.phages[phg][gnum].type == "R"){
            if (boid.X <= boid.genome.phages[phg][gnum].Xhf) {
              boid.repression = true
              boid.fitness -= repression_cost
            }
          }
        }
      }

      // lysis

      shuffle(boid.genome.phages)
      if (boid.genome.fitness != 0 && !boid.repression) {
        for (let phg = 0; phg < boid.genome.phages.length; phg++) {
          if (!death && boid.genome.phages[phg].length > 0 && boid.genome.phages[phg][0].type == "L") {
            if (sim.rng.genrand_real1() < boid.genome.phages[phg][0].lysis) {
              for (let offnum = 0; offnum < Math.floor((burst_size - vir_cost*boid.genome.phages[phg].length)*boid.genome.fitness/b0); offnum++){
                let new_virion = new Virion(false, boid.genome.phages[phg], boid.genome.chromosome)
                sim.environment.grid[x][y].virions.push(new_virion)
              }
              death = true
            }
          }
        }
      }

      // birth
      if (this.boids.length < max_cells && sim.rng.random() < boid.genome.fitness*boid.Y) {
        let newboid = this.copyBoid(boid);
        let angle = sim.rng.random() * Math.PI * 2;
        newboid.position.x += 0.5 * boid.size * Math.cos(angle);
        newboid.position.y += 0.5 * boid.size * Math.sin(angle);
        newboid.genome = boid.genome.copy()
        newboid.X = boid.X
        newboid.Y /= 2
        boid.Y /= 2
        newboids.push(newboid);

      }

      //death
      if (sim.rng.genrand_real1() < d0 || death) {
        let idx = this.boids.indexOf(boid)
        this.boids.splice(idx, 1)
      }

    }
    for (let off = 0; off < newboids.length; off++) this.boids.push(newboids[off])
  }
}

sim.environment.nextState = function(x, y) {
  //this.grid[x][y].X = x/size
  this.grid[x][y].Y += Y_max_influx * y/size
  this.grid[x][y].Y *= 1 - Y_decay

  for (let vir = 0; vir < this.grid[x][y].virions.length; vir++){
    let rnr = sim.rng.genrand_real1()
    if (rnr < virion_decay_rate) {
      this.grid[x][y].virions.splice(vir, 1)
      this.grid[x][y].density--
    }
    else if (rnr < virion_decay_rate + virion_diffusion_rate) {
      let nbr = this.randomMoore8(this, x, y)
      nbr.virions.push(this.grid[x][y].virions[vir])
      this.grid[x][y].virions.splice(vir, 1)
    }
  }
}

sim.environment.update = function() {

  if (sim.rng.genrand_real1() < env_fluctuation_freq){
    for (let x = 0; x < sim.ncol; x++) {
      for (let y = 0; y < sim.nrow; y++) {
        this.grid[x][y].X = 1 - this.grid[x][y].X
      }
    }
  }

  if (sim.rng.genrand_real1() < bacteria_influx_rate) {
    let angle = this.random() * 2 * Math.PI
    let new_boid = {
            position: { x: sim.rng.genrand_real1()*sim.nrow , y: sim.rng.genrand_real1()*sim.ncol },
            velocity: { x: sim.flock.init_velocity*Math.cos(angle) * sim.flock.max_speed, y: sim.flock.init_velocity*Math.sin(angle) * sim.flock.max_speed },
            acceleration: { x: 0, y: 0 },
            size: sim.flock.config.size
    }
    new_boid.genome = new Genome()
    new_boid.X = sim.environment.grid[Math.floor(new_boid.position.x)][Math.floor(new_boid.position.y)].X
    new_boid.Y = 0
    new_boid.fitness = new_boid.genome.fitness
    new_boid.repression = false
    for (let i = 0; i < new_boid.genome.chromosome.length; i++) {
      if (new_boid.genome.chromosome[i].type == "R") {
        if (new_boid.X <= new_boid.genome.chromosome[i].Xhf) {
          new_boid.repression = true
          new_boid.fitness -= repression_cost
        }
      }
    }
    sim.flock.boids.push(new_boid)
  }

  if (sim.rng.genrand_real1() < virion_influx_rate) {
    let x = sim.rng.genrand_int(0, sim.nrow - 1)
    let y = sim.rng.genrand_int(0, sim.ncol - 1)
    this.grid[x][y].virions.push(new Virion(true))
  }

  this.asynchronous()
  let num_prophages = 0
  let avg_lysis = 0
  let num_Ls = 0
  for (let boid of sim.flock.boids) {
    num_prophages = num_prophages + boid.genome.phages.length
    for (let phg = 0; phg < boid.genome.phages.length; phg++) for (let gnum = 0; gnum < boid.genome.phages[phg].length; gnum++) {
      if (boid.genome.phages[phg][gnum].type == "L") {
        num_Ls++
        avg_lysis += boid.genome.phages[phg][gnum].lysis
      }
    }
  }

  let num_virions = 0
  for (let x = 0; x < sim.ncol; x++) {
    for (let y = 0; y < sim.nrow; y++) {
      this.grid[x][y].density = this.grid[x][y].virions.length
      num_virions += this.grid[x][y].virions.length
    }
  }

  if (sim.time%1000 == 1) {
    if (sim.flock.boids.length > 0) {
      shuffle(sim.flock.boids)
      for (let ind = 0; ind < Math.min(sim.flock.boids.length, 100); ind++) sim.write_append(JSON.stringify(sim.flock.boids[ind].position) + '\t' + JSON.stringify(sim.flock.boids[ind].genome.chromosome) + '\t' + JSON.stringify(sim.flock.boids[ind].genome.phages) + '\n', 'test_data/VD'+VD+'iter'+iter+'time'+sim.time+'.dat')
      sim.write_append(sim.time + '\t' + sim.flock.boids.length + '\t' + avg_lysis/num_prophages + '\t' + num_prophages + '\t' + num_virions + '\n', 'test_data/timeseries_VD'+VD+'iter'+iter+'.dat')
    }
  }
}


sim.start()
