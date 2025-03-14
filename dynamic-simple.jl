using DifferentialEquations, Plots, Statistics, Distributions

duration = 1000.0
# Parameters from the Cyton2 paper
λ_firstDivision = LogNormal(log(39.89), 0.28)
λ_subsequentDivision = 9.21
λ_divisionDestiny = LogNormal(log(71.86), 0.11)
λ_lifetime = LogNormal(log(116.8), 0.85)


struct Cell
  time_to_die::Float64
  time_to_divide::Float64
end
Cell() = Cell(rand(λ_lifetime), rand(λ_firstDivision))
Cell(c::Cell, t::Float64) = Cell(c.time_to_die, t+rand(λ_firstDivision))

population = [Cell() for i in 1:10]
tstops = vcat(map(c -> c.time_to_die, population), map(c -> c.time_to_divide, population))

function division(u, t, integrator)
  minimum(map(c -> c.time_to_divide-t, population))
end

function divide!(integrator)
  t = integrator.t
  indx = findmin(map(c -> c.time_to_divide-t, population))[2]
  cell = population[indx]
  deleteat!(population, indx)
  for _ in 1:2
    c = Cell(cell, t)
    push!(population, c)
    add_tstop!(integrator, c.time_to_die)
    add_tstop!(integrator, c.time_to_divide)
  end
end

division_cb = ContinuousCallback(division, divide!)

function death(u, t, integrator)
  minimum(map(c -> c.time_to_die-t, population))
end

function die!(integrator)
  t = integrator.t
  deads = findall(c -> c.time_to_die-t == 0, population)
  deleteat!(population, deads)

  @info "$(length(deads)) dead cells at $(integrator.t)"

  if length(population) == 0
    @info "Stopping at $(integrator.t) - no cells"
    terminate!(integrator)
  end
end

death_cb = ContinuousCallback(death, die!)

struct PopulationData
  count::Int
end

pop_data = SavedValues(Float64, PopulationData)
function pop_count_saver(u, t, integrator)
  n = length(population)
  return PopulationData(n)
end
pop_count_cb = SavingCallback(pop_count_saver, pop_data, saveat=0:0.1:duration)

cbs = CallbackSet(division_cb, death_cb, pop_count_cb)
u0 = [0.0]
tspan = (0.0, duration)
prob = ODEProblem((du, u, p, t)->du[1]=1, u0, tspan)
sol = solve(prob, Tsit5(), callback = cbs, tstops = tstops)

# h = plot(sol.t, map((x) -> length(x)/2, sol[:]), lw = 3, ylabel = "Number of Cells", xlabel = "Time")
cell_count = map((x) -> x.count, pop_data.saveval)

h = plot(pop_data.t, cell_count, lw = 3, ylabel = "Number of Cells", xlabel = "Time")
display(h)
population