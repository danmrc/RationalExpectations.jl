# Klein Solver

```@docs
klein(A::AbstractArray,B::AbstractArray,C::AbstractArray,t::Int,k0::AbstractArray,shock_exp::AbstractArray,jumps::AbstractArray)
```

## Important

When you add an autocorrelated shock to the system, the value on the lhs **must be on t+1**. Otherwise, the IRF will be wrongly computed.
