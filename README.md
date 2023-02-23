# An agent based model of cell population homeostasis

*Note:* This is a rapidly evolving project and it is likely this README is out of date.

## Introduction
The model is intended to illustrate and explore the behaviour of a cell system undergoing cytokine directed homeostasis. In such a model, the cell number is maintained at a set point (or other attractor) by balancing the rates of survival and division. There are many examples of homeostatic systems in biology, but we will attempt to most closely align our model to T cell homeostasis as regulated by the cytokine IL7. However, as IL-7 has a complex receptor, complex signalling pathways and intersection with multiple intracellular fate processes we prefer to first simplify and experiment with key principles of the modular (Cyton) philosophy for model development. It is therefore important to understand that the IL7, as used in this model, is really and idealised and simplified cytokine. Similarly, the cytokine receptor, IL7R, is also simplified. We are particularly interested in agent based versions that can recreate and explore the effect of drug therapies (e.g. venetoclax) on the proportion of cells in cell cycle phases. 

## Model overview
We assume there is an occupancy function, $f_o([IL7])$, which determines the proportion amount of bound receptor as a function the concentration of IL7, $[IL7]$. This will be something like a Hill function which starts at 0 and saturates at 1 when $[IL7]=\infty$. Without loss of generality, we can assume the 0.5 level of $f_o$ occurs at $[IL7]=1$. This leaves a single free parameter which is essentially the slope at this point.

We model the $[IL7]$ economy by assuming it is produced at a constant, basal rate. Degradation is assumed to occur at a rate proportional to the amount of IL7 bound to receptor, e.g. through receptor internalisation. Thus:

$$
\frac{d[IL7]}{dt} = p_{[IL7]} - \sum_{i}\,\alpha_i\,f_o([IL7])
$$

where $p_{[IL7]}$ is the constant production rate, $\alpha_i$ is the absorption rate of the ith cell and the sum is over cells.

Each cell has two internal modules, or timers in Cyton terminology.

1. [A Myc driven destiny timer](https://pubmed.ncbi.nlm.nih.gov/27820810/). When Myc is above threshold cell division is initiated. The cell them moves though G1/S/G2M phases and then divides. If the Myc level in the daughter cells is above threshold cell cycle will be initiated again.
2. A survival timer. This is based on (as yet unpublished) work that shows that intrinsic apoptosis is triggered when an ensemble of apoptotic proteins falls below a threshold. In this simulation, the collection of proteins is modeled as a single quantity.

Protein levels providing these functions are modelled as

$$
\frac{dq_{i,p}}{dt} =  \alpha_i \, f_o([IL7]) - \lambda_i q_{i,p}
$$

where $q_{i,p}$ is level for protein $p$ in cell $i$. There is only one absorption for each cell because the IL7R will drive multiple internal pathways. The parameters $\alpha_i$ and $\lambda_i$ are typically per cell and drawn from distributions.

Cells initiate cell cycle when $q_{Myc} > T_i$, where $T_i$ is also a per cell threshold drawn from a distribution. (There is a redundancy between the $T_i$, $\alpha_i$ and $\lambda_i$ because they can be scaled by the same amount and reproduce the same behaviour of the model). Once cell cycle is initiated there is a per cell delay (drawn from a distribution) until the cell divides.

On division two daughter cells are created:
* Their initial protein values are those of the mother
* *All* Parameter values, for both the survival protein and Myc are copied from the mother

Cells die when $q_{survival} < T_i$, where $T_i$ is also a per cell threshold drawn from a distribution.

## Code issues
The Cyton module stepper updates the model by calling a `step` function to update each cell and IL7 autonomously and then an `interact` function to model the interaction between IL7 and the cells. (This is a poor design and it will change.)


The `step` function for the `Survival` timer looks like:
```
function step(st::SurvivalTimer, time::Float64, Δt::Duration)::Union{CellEvent, Nothing}
```
It can return a `Death()` event when the cell dies or `nothing` otherwise.


The step function for the `DivisionTimer` looks like
```
function step(myc::DivisionTimer, time::Time, Δt::Duration)::Union{CellEvent, Nothing}
```
It can return a `Division()` event when the cell dies or `nothing` otherwise. When `Division()` is returned the cell divides, the cell factory is called to create a new cell and the `inherit` methods are called for the fate timers.

IL7 `step` function looks like:
```
function step(iL7::IL7, time::Time, Δt::Duration, model::CytonModel)
```

The `interact` function looks like:
```
function interact(iL7::IL7, cell::Cell, time::Time, Δt::Duration)
```
