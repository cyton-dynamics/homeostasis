using DifferentialEquations, Plots, Statistics

duration = 20
α = 0.3
β = 0.2

function f(du, u, p, t)
  n = length(u)
  for i = 1:2:n
    # division protein
    du[i] = α * u[i]
    # death protein
    du[i+1] = - β * u[i+1]
  end
end

division(u, t, integrator) = 1 - maximum(u[begin:2:end])

function divide!(integrator)
  u = integrator.u
  maxidx = findmax(u[begin:2:end])[2] * 2 - 1
  resize!(integrator, length(u) + 2)

  # Division protein
  Θ = rand()
  u[maxidx] = Θ
  u[end-1] = 1 - Θ

  # Inherit death protein
  u[end] = u[maxidx+1] + rand()*0.01
end

division_cb = ContinuousCallback(division, divide!)

T = 0.1
death(u, t, integrator) = minimum(u[begin+1:2:end]) - T

function die!(integrator)
  u = integrator.u
  minidx = findmin(u[begin+1:2:end] .- 0.1)[2] * 2 
  deleteat!(integrator, [minidx-1, minidx])

  @info "Cell died at $(integrator.t)!"

  if length(u) == 0
    @info "Stopping at $(integrator.t) - no cells"
    terminate!(integrator)
  end
end

death_cb = ContinuousCallback(death, die!)

struct PopulationData
  count::Int
  death_protein::Float64
  division_protein::Float64
end

pop_data = SavedValues(Float64, PopulationData)
function pop_count_saver(u, t, integrator)
  n = length(u)
  death_protein = mean(u[begin+1:2:end])
  division_protein = mean(u[begin:2:end])
  return PopulationData(n, death_protein, division_protein)
end
pop_count_cb = SavingCallback(pop_count_saver, pop_data, saveat=0:0.1:duration)

cbs = CallbackSet(division_cb, death_cb, pop_count_cb)
u0 = rand(10)
tspan = (0.0, duration)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), callback = cbs)

# h = plot(sol.t, map((x) -> length(x)/2, sol[:]), lw = 3, ylabel = "Number of Cells", xlabel = "Time")
cell_count = map((x) -> x.count, pop_data.saveval)
death_protein = map((x) -> x.death_protein, pop_data.saveval)
division_protein = map((x) -> x.division_protein, pop_data.saveval)

h = plot(pop_data.t, cell_count, lw = 3, ylabel = "Number of Cells", xlabel = "Time")
display(h)
h = plot(pop_data.t, [death_protein, division_protein], lw = 3, ylabel = "Average protein level", xlabel = "Time", label=["death" "division"])
display(h)
 