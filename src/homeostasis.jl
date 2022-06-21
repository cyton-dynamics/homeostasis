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



λ_dproDecay = LogNormalParms(log(log(2)/24), 0.34)
dproThreshold = LogNormalParms(log(1), 0.2)
dproInitial = LogNormalParms(log(2), 0.2)



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
MycTimer(λ::DistributionParmSet, myc::DistributionParmSet, threshold::DistributionParmSet) = MycTimer(sample(λ), sample(myc), sample(threshold))

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
#--------------------------Death Protien-------------------------------
"The state of the Dpro timer"
mutable struct DproTimer <: FateTimer
  # Dpro decay rate
  λ::Float64
  # Current level of Dpro
  dpro::Float64
  # If Dpro drops below this threshold, then the cell dies
  threshold::Float64
end
DproTimer(λ::DistributionParmSet, dpro::DistributionParmSet, threshold::DistributionParmSet) = DproTimer(sample(λ), sample(myc), sample(threshold))

"At each time step Dpro decays but is also driven by constant exogenous stimulus"
function step(dpro::DproTimer, time::Float64, Δt::Duration)::Union{CellEvent, Nothing}
  dpro.dpro *= exp(-dpro.λ*Δt)
  if dpro.dpro < dpro.threshold
    return Death()
  else
    return nothing
  end
end

"Daughter cells inherit the mother's Death protien timer"
inherit(dpro::DproTimer, ::Time) = dpro

#----------------------------------------------------------------------
#------------------ Division machinery --------------------
"Time to divide drawn from distribution"
mutable struct DivisionTimer <: FateTimer
  nextDivision::Float64
  reachedDestiny::Bool
end
"Constructor for new cells"
DivisionTimer(division::DistributionParmSet) = DivisionTimer(sample(division), false)

function step(timer::DivisionTimer, time::Float64, ::Float64) 
  if time>timer.nextDivision
    return Division()
  else
    return nothing
  end
end

"Daughter cells get a new division timer"
inherit(::DivisionTimer, time::Time) = DivisionTimer(λ_subsequentDivision, time)
DivisionTimer(r::DistributionParmSet, start::Time) = DivisionTimer(sample(r) + start, false)

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
  DeathTimer(sample(r), r)
end
"DeathTimer constructor for division"
function DeathTimer(death::DeathTimer, time::Time)
  DeathTimer(sample(death.deathTimeDistribution)+time, death.deathTimeDistribution)
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
ExogeneousStimulus(d::DistributionParmSet) = ExogeneousStimulus(sample(d))
ExogeneousStimulus(d::DistributionParmSet,fractional_occupancy::Float64) = ExogeneousStimulus(sample(d)*fractional_occupancy)
function stimulate(::FateTimer, ::Stimulus, ::Time, ::Duration) end

function stimulate(myc::MycTimer, stim::ExogeneousStimulus, ::Time, Δt::Duration)
  myc.myc +=stim.strength*(1- exp(-myc.λ*Δt))/myc.λ
end

function stimulate(dpro::DproTimer, stim::ExogeneousStimulus, ::Time, Δt::Duration)
  dpro.dpro +=stim.strength*(1- exp(-dpro.λ*Δt))/dpro.λ
end

function stimulate(cell::Cell{GenericCell}, stim::ExogeneousStimulus, time::Time, Δt::Duration)
  for timer in cell.timers
    stimulate(timer, stim, time, Δt)
  end
end
#----------------------------------------------------------
#------------------------------Agent based modelling of IL7------

"""
The current free level of the IL7 in the system 

"""

affinity_to_cell=4
production_rate=10
absorption_rate_IL7=2
#---------Environment type and population creation-----
mutable struct IL7<: EnvironmentalAgent
  concentration::Float64
end

function IL7()::IL7
  return IL7(100)
end

function IL7(conc::Float64)::IL7
  return IL7(conc)
end

"""
environmentFactory()
nEnvAgents:: Number of environment agents you want 
Function to create a bunch of environment agents at the start of the simulation
"""

function environmentFactory()::Vector{EnvironmentalAgent}
  environmentAgents=[]
  IL7=IL7(20)
  push!(environmentAgents,IL7)
  return environmentAgents
end

function frac_ocu(concentration::Float64,affinity_to_cell::Int)
  frac_ocu=log(concentration)^affinity_to_cell/(1+log(concentration)^affinity_to_cell)
  return frac_ocu
end

function step(IL7::IL7,time::Time,Δt::Duration,model::CytonModel )
  number_of_cells=length(model.cells)
  IL7.concentration+=production_rate*Δt 
  IL7.concentration-=number_of_cells*absorption_rate_IL7*Δt*frac_ocu(IL7.concentration,affinity_to_cell)

end
#------------------------------------------------------


#--------Defining the interact() funciton--------------
n_bound_to_cell::Float64=0.0

function interact(IL7::IL7, cell::Cell, time::Time, Δt::Duration)

  α_myc=LogNormalParms(log(log(2)/24), 0.34)
  myc_stim=ExogeneousStimulus(sample(α_myc)*frac_ocu(IL7.concentration,affinity_to_cell))
  α_dpro=LogNormalParms(log(log(2)/24), 0.34)
  dpro_stim=ExogeneousStimulus(sample(α_dpro)*frac_ocu(IL7.concentration,affinity_to_cell))
  
  myctimer=filter(x->typeof(x)==MycTimer,cell.timers)
  dprotimer=filter(x->typeof(x)==DproTimer,cell.timers)

  stimulate(myctimer,myc_stim,time,Δt)
  stimulate(dprotimer,dpro_stim,time,Δt) 

end


end # module homeostasis
