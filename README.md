# NGSpiceTools.jl

Developing a toolset for simulating circuits, using NGSpice as the backend. In a perfect world, I would then combine this with [circuit-designer](git@github.com:pkeshava/circuit-designer.git) as the frontend for the schematic capture... let's see. 

## Docs

```zsh
cd docs
julia --project=. make.jl
```

## Testing

```zsh
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)
```

## Todo

- [ ] AC analysis of simple amplifier stages
- [ ] Incorporate more spice models in examples
- [ ] Work on automatic testing...
