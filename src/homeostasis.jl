
module homeostasis

using Cyton
import Cyton: inherit, step, interact
using Parameters: @with_kw

export environmentFactory, cellFactory, DivisionTimer, SurvivalTimer, phase
export Phase, G1, S, G2M, timeInG1, divisionTimer, survivalTimer

@enum Phase G1 S G2M

#----------------------------------------------------------
"Parameters of the Cell Cycle"
minG1Duration = LogNormalParms(log(3), 0.3)
sDuration = LogNormalParms(log(5), 0.5)
g2mDuration = LogNormalParms(log(2), 0.2)

"Myc parameters"
tHalf = 24
λ = log(2)/tHalf
λ_mycDecay = LogNormalParms(log(λ), 1)
mycThreshold = LogNormalParms(log(1), 0.2)
mycInitial = FixedDistributionParms(0.0)# LogNormalParms(log(10), 0.2)
α_myc = LogNormalParms(log(1), 0.34)

"Survival protein parameters"
tHalf = 24
λ = log(2)/tHalf
λ_survivalProteinDecay = LogNormalParms(log(λ), 1)
#survivalProteinThreshold = LogNormalParms(log(1), 0.2)
survivalProteinInitial = LogNormalParms(log(1), 0.2)
α_survivalProtein = LogNormalParms(log(1), 0.34)

"IL7 parameters"
affinity_to_cell = 5
production_rate = 0.1
absorption_rate_IL7 = 0.1

function cellFactory(birth::Time, parms::TrialParameters)
  cell = Cell(birth)

  divisionTimer = DivisionTimer(minG1Duration, 
    sDuration, 
    g2mDuration, 
    λ_mycDecay, 
    α_myc, 
    mycInitial, 
    mycThreshold)
  addTimer(cell, divisionTimer)

  survivalProteinThreshold = LogNormalParms(parms.threshold, 0.2)
  st = SurvivalTimer(λ_survivalProteinDecay, α_survivalProtein, survivalProteinInitial, survivalProteinThreshold)
  addTimer(cell, st)

  return cell
end
#----------------------------------------------------------

#-------------------- Survival Protien -----------------------
"The state of the Dpro timer"
@with_kw mutable struct SurvivalTimer <: FateTimer
  # Dpro decay rate
  λ::Float64
  # Sensitivity to IL7
  α::Float64
  # Current level of Dpro
  survivalProtein::Float64
  # If Dpro drops below this threshold, then the cell dies
  threshold::Float64
end
function SurvivalTimer(λ::DistributionParmSet, α::DistributionParmSet, dpro::DistributionParmSet, threshold::DistributionParmSet) 
  SurvivalTimer(; λ=sample(λ), α=sample(α), survivalProtein=sample(dpro)+sample(threshold), threshold=sample(threshold))
end

function step(st::SurvivalTimer, time::Float64, Δt::Duration)::Union{CellEvent, Nothing}
  λ = st.λ
  threshold = st.threshold
  st.survivalProtein *= exp(-λ*Δt)
  if st.survivalProtein < threshold
    return Death()
  else
    return nothing
  end
end

"Daughter cells inherit the mother's Survival protein timer"
inherit(st::SurvivalTimer, ::Time) = SurvivalTimer(st)

function update(st::SurvivalTimer, ::Time, Δt::Duration, strength::Float64)
  st.survivalProtein += strength * st.α * Δt
end
#----------------------------------------------------------------------

#------------------ Division machinery --------------------
@with_kw mutable struct DivisionTimer <: FateTimer
  # This records the time that the cell commits to dividing
  committedToDivideAt::Union{Nothing, Time}

  # Parameters controlling the
  minG1Duration::Duration
  sDuration::Duration
  g2mDuration::Duration

  # Myc decay rate
  λ::Float64
  # Sensitivity to IL7
  α::Float64
  # Current level of myc
  myc::Float64
  # If Myc drops below this threshold, no more dividing!
  threshold::Float64
end
function DivisionTimer(minG1Duration::DistributionParmSet, 
  sDuration::DistributionParmSet, 
  g2mDuration::DistributionParmSet,
  λ::DistributionParmSet, 
  α::DistributionParmSet, 
  myc::DistributionParmSet, 
  threshold::DistributionParmSet) 

  DivisionTimer(;
  committedToDivideAt=nothing, 
  minG1Duration=sample(minG1Duration),
  sDuration=sample(sDuration),
  g2mDuration=sample(g2mDuration),
  λ=sample(λ), 
  α=sample(α), 
  myc=sample(myc), 
  threshold=sample(threshold))

end

timerFor(T::Type, cell::Cell) = first(filter(x -> x isa T, cell.timers))
divisionTimer(cell::Cell) = timerFor(DivisionTimer, cell)
survivalTimer(cell::Cell) = timerFor(SurvivalTimer, cell)

function update(myc::DivisionTimer, ::Time, Δt::Duration, strength::Float64)
  myc.myc += strength * myc.α  * Δt
end

function phase(cell::Cell, time::Time)
  divtimer = divisionTimer(cell)
  return phase(divtimer, time)
end

"Returns the phase of the cycle and the proportion of time spent in that phase"
function phase(cycle::DivisionTimer, time::Time)::Tuple{Phase, Union{Nothing, Duration}}
  t = cycle.committedToDivideAt

  if t === nothing
    return (G1, nothing)
  end
  
  δt = time - t

  if δt < cycle.minG1Duration
    return (G1, δt)
  end

  if δt < cycle.minG1Duration + cycle.sDuration
    return (S, δt)
  end

  return (G2M, δt)
end

function timeInG1(cell::Cell, time::Time)::Duration
  dt = divisionTimer(cell)
  cd = dt.committedToDivideAt
  if cd === nothing
    return time - cell.birth
  end

  if time - cd < dt.minG1Duration
    return time - cell.birth
  end

  return cd + dt.minG1Duration - cell.birth
end

cycleTime(myc::DivisionTimer) = myc.minG1Duration + myc.sDuration + myc.g2mDuration

"Step function that will help the cell to commit to division"
function step(myc::DivisionTimer, time::Time, Δt::Duration)::Union{CellEvent, Nothing}
  myc.myc *= exp(-myc.λ*Δt)

  (p, δt) = phase(myc, time)

  if p == G1 && δt === nothing && myc.myc > myc.threshold
    myc.committedToDivideAt = time
    return nothing
  end

  if p == G2M && δt >= cycleTime(myc)
    return Division()
  else
    return nothing
  end
end

"Daughter cells get a new division timer, inherits parents properties"
function inherit(d::DivisionTimer, ::Time) 
  DivisionTimer(;
    committedToDivideAt=nothing, 
    minG1Duration=d.minG1Duration,
    sDuration=d.sDuration,
    g2mDuration=d.g2mDuration,
    λ=d.λ, 
    α=d.α, 
    myc=d.myc,
    threshold=d.threshold)
end
#----------------------------------------------------------

#------------------------------Agent based modelling of IL7------

#---------Environment type and population creation-----
mutable struct IL7 <: EnvironmentalAgent
  concentration::Float64
end

makeIL7(conc::Float64=5.0) = IL7(conc)

"""
environmentFactory()

Function to create a bunch of environment agents at the start of the simulation
"""

environmentFactory()::Vector{EnvironmentalAgent} = [makeIL7(0.0)]

function frac_ocu(concentration::Float64, affinity_to_cell::Int)
  frac_ocu = concentration^affinity_to_cell/(1+concentration^affinity_to_cell)
  return frac_ocu
end

function step(iL7::IL7, time::Time, Δt::Duration, model::CytonModel)
  number_of_cells = length(model.cells)
  concentration = iL7.concentration
  concentration += production_rate * Δt
  concentration -= number_of_cells * absorption_rate_IL7 * Δt * frac_ocu(iL7.concentration, affinity_to_cell)
  if concentration < 0
    s = "CONCENTRATION=$(concentration) GOING BELOW 0!!!\n"
    s *= "Forcing concentration to 0."
    s *= "This suggests numerical instability. You should fix this!"
    @warn s
    concentration = 0
  end
  iL7.concentration = concentration
  return nothing
end
#--------------------------------------------------------]


#---------- Defining the interact() funciton ------------
function interact(iL7::IL7, cell::Cell, time::Time, Δt::Duration)
  o = frac_ocu(iL7.concentration, affinity_to_cell)
  update(divisionTimer(cell), time, Δt, o)
  update(survivalTimer(cell), time, Δt, o) 
end

end 

