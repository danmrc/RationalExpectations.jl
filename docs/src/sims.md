# Sims Solver

```@docs

sims(G0::AbstractArray,G1::AbstractArray,Pi::AbstractArray,Psi::AbstractArray)

irf(Theta1::AbstractArray,Theta2::AbstractArray,t::Int,shock)
```

## Note

We don't force Julia to use the complex Schur decomposition. If you want it, just use `sims(complex(G0),complex(G1),Pi,Psi)`. However, you will not be able to plot (yet) simply passing the model to `irf`
