module homeostasis

using Cyton
import Cyton: shouldDie, shouldDivide, inherit, step, stimulate, FateTimer, Stimulus

#----------------------- Parameters -----------------------
abstract type Parameters end
struct NoParms <: Parameters end
#----------------------------------------------------------

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

#------------------- Myc destiny timer --------------------
"The event that indicates that the cell has reached division destiny"
struct DivisionDestiny <: CellEvent end

"The state of the Myc timer"
mutable struct MycTimer <: FateTimer
  # Myc decay rate
  λ::Float64
  # Current level of myc
  myc::Float64
  # If Myc drops below this threshold, no more dividing!
  threshold::Float64
end
MycTimer(λ::DistributionParmSet, myc::DistributionParmSet, threshold::DistributionParmSet) = MycTimer(draw(λ), draw(myc), draw(threshold))

"At each time step Myc decays but is also driven by constant exogenous stimulus"
function step(myc::MycTimer, time::Time, Δt::Duration)::Union{CellEvent, Nothing}
  myc.myc *= exp(-myc.λ*Δt)
  if myc.myc < myc.threshold
    return DivisionDestiny()
  else
    return nothing
  end
end

"Daughter cells inherit the mother's Myc timer"
inherit(myc::MycTimer, ::Time) = myc

"Once destiny is reached we need to tell the division timer to stop dividing"
function destinyReached(::DivisionDestiny, cell::Cell, ::Time)
  for timer in cell.timers
    if typeof(timer) == DivisionTimer
      timer.reachedDestiny = true
    end
  end
end
#----------------------------------------------------------

#------------------ Division machinery --------------------
"Time to divide drawn from distribution"
mutable struct DivisionTimer <: FateTimer
  nextDivision::Float64
  reachedDestiny::Bool
end
"Constructor for new cells"
DivisionTimer(division::DistributionParmSet) = DivisionTimer(draw(division), false)

function step(timer::DivisionTimer, time::Float64, ::Float64) 
  if time>timer.nextDivision
    return Division()
  else
    return nothing
  end
end

"Daughter cells get a new division timer"
inherit(::DivisionTimer, time::Time) = DivisionTimer(λ_subsequentDivision, time)
DivisionTimer(r::DistributionParmSet, start::Time) = DivisionTimer(draw(r) + start, false)

"Indicate the cell will divide. Must be earlier than destiny and after the next division time"
shouldDivide(division::DivisionTimer, time::Time) = !division.reachedDestiny && time > division.nextDivision
#----------------------------------------------------------

#--------------------- Death machinery --------------------
"The death timer"
struct DeathTimer <: FateTimer
  timeToDeath::Float64
  deathTimeDistribution::DistributionParmSet
end
"DeathTimer constructor for initial cells"
function DeathTimer(r::DistributionParmSet)
  DeathTimer(draw(r), r)
end
"DeathTimer constructor for division"
function DeathTimer(death::DeathTimer, time::Time)
  DeathTimer(draw(death.deathTimeDistribution)+time, death.deathTimeDistribution)
end

"On division, daughter cells inherit the death timer"
inherit(timer::DeathTimer, time::Time) = DeathTimer(timer, time)
function step(timer::DeathTimer, time::Time, ::Duration)
  if time > timer.timeToDeath
    return Death()
  else
    return nothing
  end
end
#----------------------------------------------------------

#------------------- Stimulus machinery -------------------
struct ExogeneousStimulus <: Stimulus
  strength::Float64
end
ExogeneousStimulus(d::DistributionParmSet) = ExogeneousStimulus(draw(d))

function stimulate(::FateTimer, ::Stimulus, ::Time, ::Duration) end

function stimulate(myc::MycTimer, stim::ExogeneousStimulus, ::Time, Δt::Duration)
  myc.myc += stim.strength * Δt 
end

function stimulate(cell::Cell{GenericCell}, stim::ExogeneousStimulus, time::Time, Δt::Duration)
  for timer in cell.timers
    stimulate(timer, stim, time, Δt)
  end
end
#----------------------------------------------------------

end # module homeostatis
