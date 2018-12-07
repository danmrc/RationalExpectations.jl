# Sims Solver

```@docs

sims(G0::AbstractArray,G1::AbstractArray,Pi::AbstractArray,Psi::AbstractArray)

irf(Theta1::AbstractArray,Theta2::AbstractArray,t::Int,shock)
```

## Important

When you add an autocorrelated shock to the system, the value on the lhs **must be on t+1**. Otherwise, the IRF will be wrongly computed.
