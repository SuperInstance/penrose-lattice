# penrose-lattice

**Penrose tilings as spectral graphs. Fibonacci substitution. Inflation/deflation symmetry. Farey sequences.**

Pure Rust, zero dependencies. A library for computing spectral properties of Penrose tilings — non-periodic graphs whose edge ratios converge to the golden ratio φ.

## What it does

Penrose tilings are aperiodic tilings built from two rhombus types (thick and thin). Their ratio of thick-to-thin tiles converges to φ ≈ 1.618. This library models Penrose tilings as weighted graphs and computes their spectral properties — eigenvalues, conservation ratios, Laplacians — to explore how aperiodic structure affects spectral behavior.

## Core concepts

### Substitution Matrix

The Fibonacci substitution matrix `[[1,1],[1,0]]` drives the tiling. Its eigenvalues are φ and -1/φ. Powers of this matrix generate Fibonacci numbers, and the ratio of consecutive Fibonacci numbers converges to φ.

```rust
let (e1, e2) = SubstitutionMatrix::eigenvalues();
// e1 = φ ≈ 1.618, e2 = -1/φ ≈ -0.618

let ratios = SubstitutionMatrix::golden_ratio_convergence(10);
// Ratios converge to φ
```

### Penrose Graph

Inflation starts from a seed thick rhombus (4 vertices) and repeatedly subdivides edges using Fibonacci-like splitting. Each inflation adds bridge nodes at edge midpoints, splitting edges with weights proportional to φ.

```rust
let pg = PenroseGraph::from_inflation(4);
println!("{} nodes, {} edges", pg.nodes.len(), pg.edges.len());
println!("Tile ratio: {:.6} (approaches φ)", pg.tile_ratio());
println!("CR: {:.6}", pg.conservation());
```

### Inflation/Deflation

- **Inflation**: Each edge gets a bridge node, split into two edges weighted by φ/(φ+1) and 1/(φ+1).
- **Deflation**: Reverse — remove bridge nodes and reconnect.

```rust
let crs = PenroseGraph::cr_trajectory(5);
// CR values across inflation levels

let inflated = InflationDeflation::inflate_graph(&adjacency);
let repeated = InflationDeflation::repeated_inflation(&adjacency, 3);
```

### Farey Sequences

Farey sequences of order n enumerate reduced fractions with denominators ≤ n. Adjacent fractions satisfy |bc - ad| = 1. The Farey graph connects adjacent fractions and has interesting spectral properties.

```rust
let farey = FareySequence::of_order(5);
// 11 fractions: 0/1, 1/5, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 1/1

let adj = farey.farey_graph();
let cr = farey.conservation();
```

### Stern-Brocot Tree

The Stern-Brocot tree enumerates all positive rationals via mediants. Each level doubles the number of fractions.

```rust
let sb = FareySequence::stern_brocot_tree(4);
```

## Architecture

```
SubstitutionMatrix  — Fibonacci matrix powers, golden ratio convergence
PenroseGraph        — inflation, deflation, adjacency, Laplacian, CR
InflationDeflation  — graph-level inflation with CR tracking
FareySequence       — Farey sequences, adjacency graphs, Stern-Brocot tree
jacobi_eigenvalues  — Jacobi iteration for symmetric eigenvalue problems
```

## Running

```bash
cargo run
```

Prints substitution eigenvalues, golden ratio convergence, Penrose graph stats, CR trajectory, Farey sequence analysis, and inflation experiments.

## Testing

```bash
cargo test
```

21 tests covering: Fibonacci substitution, golden ratio convergence, Penrose inflation/deflation, tile ratio convergence to φ, Laplacian row sums, CR computation, Farey sequence counts and adjacency, Stern-Brocot tree, and integration tests.

## API

- `SubstitutionMatrix::eigenvalues()` → `(φ, -1/φ)`
- `SubstitutionMatrix::golden_ratio_convergence(n)` → `Vec<f64>`
- `SubstitutionMatrix::apply_n(n)` → Fibonacci matrix power
- `PenroseGraph::from_inflation(levels)` → graph with nodes, edges, tile counts
- `PenroseGraph::conservation()` → CR of the graph
- `PenroseGraph::cr_trajectory(max_levels)` → CR at each inflation level
- `PenroseGraph::tile_ratio()` → thick/thin ratio (→ φ)
- `InflationDeflation::inflate_graph(adj)` → inflated adjacency
- `InflationDeflation::repeated_inflation(adj, levels)` → CR trajectory
- `FareySequence::of_order(n)` → Farey sequence with graph and CR
- `FareySequence::stern_brocot_tree(depth)` → Stern-Brocot enumeration

Part of the [SuperInstance OpenConstruct](https://github.com/SuperInstance/OpenConstruct) ecosystem.

## License

MIT
