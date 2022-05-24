using homeostasis

using Cyton
import Cyton: shouldDie, shouldDivide, inherit, step, stimulate, FateTimer, Stimulus

using DataFrames, Gadfly
# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

# Parameters from fitting MR-70 with cyton solver. (σ for subsequent division is a guess)
λ_subsequentDivision = LogNormalParms(log(11.1), 0.08)
λ_death = LogNormalParms(log(175.8), 0.34)
# Made up Parameters
λ_mycDecay = LogNormalParms(log(log(2)/24), 0.34)
mycThreshold = LogNormalParms(log(1), 0.2)
mycInitial = LogNormalParms(log(2), 0.2)

"This function creates cells at the beginning of the simulation"
function cellFactory(birth::Time=0.0 ;parms::Parameters=NoParms(), cellType::T=GenericCell()) where T <: CellType
  cell = Cell(birth, cellType)
  
  myc = MycTimer(λ_mycDecay, mycInitial, mycThreshold)
  addTimer(cell, myc)

  divisionTimer = DivisionTimer(λ_subsequentDivision)
  addTimer(cell, divisionTimer)

  deathTimer = DeathTimer(λ_death)
  addTimer(cell, deathTimer)

  addObserver(DivisionDestiny(), cell, destinyReached)

  return cell
end

function run(model::CellPopulation, runDuration::Time)
  print("Time to run:")
  @time begin
    counts = DataFrame(time=Float64[], count=Int[])
    Δt = modelTimeStep(model)

    for tm in 1:Δt:runDuration
      step(model)
      push!(counts, (tm, length(model.cells)))
    end
  end

  h = plot(counts, x=:time, y=:count, Geom.line)
  display(h)

  println("Done at model time=$(modelTime(model))")
end


@info rpad(lpad(" start ", 30, "-"), 55, "-")
@info "Creating population"
population = createPopulation(100, (birth) -> cellFactory(birth))
@info "Running model"
run(population, 100.0)
@info "done"
@info rpad(lpad(" end ", 31, "-"), 56, "-")
