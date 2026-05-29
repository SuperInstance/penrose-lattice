# penrose-lattice

**Penrose tilings and aperiodic lattices — Fibonacci substitution, golden-ratio geometry, and spectral structure.**

Pure Rust. Generates Penrose tilings via substitution rules, constructs the corresponding graph, and analyzes spectral properties that emerge from the aperiodic structure. The golden ratio φ appears everywhere because the substitution matrix IS the Fibonacci matrix.

## What This Gives You

- **Substitution matrix** — Fibonacci recurrence as a 2×2 matrix with eigenvalues φ and −1/φ
- **Penrose graph construction** — nodes and edges from thick/thin rhomb tiles
- **Golden ratio convergence** — watch thick/thin ratios converge to φ as substitution depth increases
- **Spectral analysis** — eigenvalue structure of the Penrose graph
- **Zero dependencies** — all computation from scratch

## The Core Idea

Penrose tilings are aperiodic: they never exactly repeat, but they have long-range order. The thick and thin rhombs tile the plane according to substitution rules governed by the Fibonacci sequence. The substitution matrix `[[1,1],[1,0]]` has eigenvalues φ and −1/φ — the golden ratio is baked into the geometry at every scale.

## Quick Start

```rust
use penrose_lattice::SubstitutionMatrix;

// The Fibonacci substitution matrix
let eigenvalues = SubstitutionMatrix::eigenvalues();
// eigenvalues = (φ, -1/φ) = (1.618..., -0.618...)

// Apply substitution n times
let m = SubstitutionMatrix::apply_n(10);
// m = [[F(11), F(10)], [F(10), F(9)]]

// Watch the thick/thin ratio converge to φ
let ratios = SubstitutionMatrix::golden_ratio_convergence(20);
for (i, r) in ratios.iter().enumerate() {
    println!("Depth {}: thick/thin = {:.6}", i + 1, r);
}
```

## How It Fits

Part of the SuperInstance mathematical ecosystem:

- **[spectral-graph-core](https://github.com/SuperInstance/spectral-graph-core)** — Spectral analysis of the resulting graphs
- **[tropical-algebra](https://github.com/SuperInstance/tropical-algebra)** — Max-plus algebra connections
- **penrose-lattice** — Aperiodic geometry (this repo)

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
penrose-lattice = { git = "https://github.com/SuperInstance/penrose-lattice" }
```

## License

MIT

Part of the [SuperInstance](https://github.com/SuperInstance) ecosystem.
