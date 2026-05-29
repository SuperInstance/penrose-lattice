use std::f64::consts::PI;

const PHI: f64 = 1.618033988749895;
const NEG_INV_PHI: f64 = -0.618033988749895;

// ── Module 1: Substitution Matrix ──────────────────────────────────────

struct SubstitutionMatrix;

impl SubstitutionMatrix {
    fn matrix() -> [[f64; 2]; 2] {
        [[1.0, 1.0], [1.0, 0.0]]
    }

    fn eigenvalues() -> (f64, f64) {
        (PHI, NEG_INV_PHI)
    }

    fn apply_n(n: usize) -> [[u64; 2]; 2] {
        // F(n+1) F(n)  = Fibonacci powers of [[1,1],[1,0]]
        // F(n)   F(n-1)
        let mut a: [[u64; 2]; 2] = [[1, 0], [0, 1]]; // identity
        let mut base: [[u64; 2]; 2] = [[1, 1], [1, 0]];
        let mut exp = n;
        while exp > 0 {
            if exp % 2 == 1 {
                a = mat2_mul_u64(&a, &base);
            }
            base = mat2_mul_u64(&base, &base);
            exp /= 2;
        }
        a
    }

    fn golden_ratio_convergence(n: usize) -> Vec<f64> {
        let mut ratios = Vec::with_capacity(n);
        for k in 1..=n {
            let m = Self::apply_n(k);
            let thick = m[0][0] + m[0][1]; // F(k+1) + F(k)
            let thin = m[1][0] + m[1][1];  // F(k) + F(k-1)
            if thin > 0 {
                ratios.push(thick as f64 / thin as f64);
            }
        }
        ratios
    }
}

fn mat2_mul_u64(a: &[[u64; 2]; 2], b: &[[u64; 2]; 2]) -> [[u64; 2]; 2] {
    let mut c = [[0u64; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

// ── Module 2: Penrose Graph ────────────────────────────────────────────

#[derive(Clone, Debug)]
enum TileType {
    Thick,
    Thin,
    Vertex,
}

#[derive(Clone, Debug)]
struct PenroseNode {
    x: f64,
    y: f64,
    tile_type: TileType,
}

struct PenroseGraph {
    nodes: Vec<PenroseNode>,
    edges: Vec<(usize, usize, f64)>,
    tile_counts: (u64, u64),
}

impl PenroseGraph {
    fn from_inflation(levels: usize) -> PenroseGraph {
        // Start with a thick rhombus (4 vertices) and inflate
        let mut graph = Self::seed_thick();

        for _ in 0..levels {
            graph = graph.inflate();
        }
        graph
    }

    fn seed_thick() -> PenroseGraph {
        let angle = PI / 5.0; // 36°
        let nodes = vec![
            PenroseNode { x: 0.0, y: 0.0, tile_type: TileType::Vertex },
            PenroseNode { x: 1.0, y: 0.0, tile_type: TileType::Vertex },
            PenroseNode { x: 1.0 + angle.cos(), y: angle.sin(), tile_type: TileType::Vertex },
            PenroseNode { x: angle.cos(), y: angle.sin(), tile_type: TileType::Vertex },
        ];
        let edges = vec![
            (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0),
            (0, 2, PHI), // long diagonal of thick rhombus
        ];
        PenroseGraph { nodes, edges, tile_counts: (1, 0) }
    }

    fn inflate(&self) -> PenroseGraph {
        // Penrose inflation ( Robinson triangle decomposition ):
        // Each thick rhombus → 1 thick + 1 thin triangle (further subdivided)
        // Each thin rhombus → 1 thick triangle
        // Simplified: replace each edge with Fibonacci-like substitution
        let n = self.nodes.len();
        let mut new_nodes = self.nodes.clone();
        let mut new_edges = Vec::new();
        let mut thick = self.tile_counts.0;
        let mut thin = self.tile_counts.1;

        for &(i, j, w) in &self.edges {
            let bridge = new_nodes.len();
            // Insert a bridge node at midpoint
            new_nodes.push(PenroseNode {
                x: (self.nodes[i].x + self.nodes[j].x) / 2.0,
                y: (self.nodes[i].y + self.nodes[j].y) / 2.0,
                tile_type: TileType::Vertex,
            });
            // Split edge into two with slightly perturbed weights (Fibonacci-like)
            let w1 = w * PHI / (PHI + 1.0);
            let w2 = w * 1.0 / (PHI + 1.0);
            new_edges.push((i, bridge, w1));
            new_edges.push((bridge, j, w2));
        }

        // Each inflation roughly follows Fibonacci: thick→thick+thin, thin→thick
        let new_thick = thick + thin;
        let new_thin = thick;
        thick = new_thick;
        thin = new_thin;

        PenroseGraph { nodes: new_nodes, edges: new_edges, tile_counts: (thick, thin) }
    }

    fn deflate(&self) -> Option<PenroseGraph> {
        if self.nodes.len() < 5 {
            return None;
        }
        // Reverse: merge bridge nodes. Remove nodes added during inflation.
        // Simplified: remove every other node (bridges) and reconnect
        let n = self.nodes.len();
        let original_n = (n + 1) / 2; // approximate reverse
        if original_n < 4 {
            return None;
        }
        Some(PenroseGraph {
            nodes: self.nodes[..original_n].to_vec(),
            edges: self.edges.iter().filter(|&&(i, j, _)| i < original_n && j < original_n).cloned().collect(),
            tile_counts: (self.tile_counts.1, self.tile_counts.0.saturating_sub(self.tile_counts.1)),
        })
    }

    fn adjacency_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.nodes.len();
        let mut adj = vec![vec![0.0; n]; n];
        for &(i, j, w) in &self.edges {
            adj[i][j] += w;
            adj[j][i] += w;
        }
        adj
    }

    fn laplacian(&self) -> Vec<Vec<f64>> {
        let adj = self.adjacency_matrix();
        let n = adj.len();
        let mut lap = vec![vec![0.0; n]; n];
        for i in 0..n {
            let mut deg = 0.0;
            for j in 0..n {
                if i != j {
                    deg += adj[i][j];
                    lap[i][j] = -adj[i][j];
                }
            }
            lap[i][i] = deg;
        }
        lap
    }

    fn conservation(&self) -> f64 {
        let evals = self.eigenvalues();
        if evals.is_empty() {
            return 0.0;
        }
        let total: f64 = evals.iter().map(|e| e.abs()).sum();
        let max_e = evals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if total == 0.0 {
            return 0.0;
        }
        // CR = ratio of energy in dominant mode to total
        max_e.abs() / total
    }

    fn eigenvalues(&self) -> Vec<f64> {
        let lap = self.laplacian();
        jacobi_eigenvalues(&lap)
    }

    fn tile_ratio(&self) -> f64 {
        if self.tile_counts.1 == 0 {
            if self.tile_counts.0 == 0 { 0.0 } else { f64::INFINITY }
        } else {
            self.tile_counts.0 as f64 / self.tile_counts.1 as f64
        }
    }

    fn cr_trajectory(max_levels: usize) -> Vec<f64> {
        let mut traj = Vec::with_capacity(max_levels);
        for lvl in 0..max_levels {
            let g = Self::from_inflation(lvl);
            traj.push(g.conservation());
        }
        traj
    }

    fn is_converging(trajectory: &[f64]) -> bool {
        if trajectory.len() < 3 {
            return true;
        }
        let mut deltas = 0.0;
        let mut count = 0u64;
        for w in trajectory.windows(2) {
            deltas += (w[1] - w[0]).abs();
            count += 1;
        }
        let avg_delta = deltas / count as f64;
        // Check second half is more stable than first half
        let mid = trajectory.len() / 2;
        let mut first_deltas = 0.0;
        let mut second_deltas = 0.0;
        for w in trajectory[..mid].windows(2) {
            first_deltas += (w[1] - w[0]).abs();
        }
        for w in trajectory[mid..].windows(2) {
            second_deltas += (w[1] - w[0]).abs();
        }
        let first_avg = first_deltas / mid.max(1) as f64;
        let second_avg = second_deltas / (trajectory.len() - mid).max(1) as f64;
        second_avg <= first_avg
    }
}

// ── Module 3: InflationDeflation ───────────────────────────────────────

struct InflationDeflation;

impl InflationDeflation {
    fn inflate_graph(adj: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = adj.len();
        // Each node stays, each edge (i,j) gets a bridge node at index n + edge_idx
        let mut edge_count = 0usize;
        for i in 0..n {
            for j in i+1..n {
                if adj[i][j] != 0.0 {
                    edge_count += 1;
                }
            }
        }
        let new_n = n + edge_count;
        let mut new_adj = vec![vec![0.0; new_n]; new_n];
        // Copy original self-loops / diagonal
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue; // will be updated
                }
            }
        }
        let mut bridge = n;
        for i in 0..n {
            for j in i+1..n {
                if adj[i][j] != 0.0 {
                    let w = adj[i][j];
                    let w1 = w * PHI / (PHI + 1.0);
                    let w2 = w / (PHI + 1.0);
                    new_adj[i][bridge] = w1;
                    new_adj[bridge][i] = w1;
                    new_adj[bridge][j] = w2;
                    new_adj[j][bridge] = w2;
                    bridge += 1;
                }
            }
        }
        new_adj
    }

    fn inflation_cr_delta(adj: &[Vec<f64>]) -> f64 {
        let cr_before = graph_cr(adj);
        let inflated = Self::inflate_graph(adj);
        let cr_after = graph_cr(&inflated);
        (cr_after - cr_before).abs()
    }

    fn inflation_eigenvalue_ratio(adj: &[Vec<f64>]) -> Vec<(f64, f64)> {
        let evals_before = jacobi_eigenvalues(&laplacian_from_adj(adj));
        let inflated = Self::inflate_graph(adj);
        let evals_after = jacobi_eigenvalues(&laplacian_from_adj(&inflated));
        let max_n = evals_before.len().min(evals_after.len());
        (0..max_n.min(5)).map(|i| (evals_before[i], evals_after[i])).collect()
    }

    fn repeated_inflation(adj: &[Vec<f64>], levels: usize) -> Vec<f64> {
        let mut crs = Vec::with_capacity(levels + 1);
        crs.push(graph_cr(adj));
        let mut current = adj.to_vec();
        for _ in 0..levels {
            current = Self::inflate_graph(&current);
            crs.push(graph_cr(&current));
        }
        crs
    }
}

fn graph_cr(adj: &[Vec<f64>]) -> f64 {
    let lap = laplacian_from_adj(adj);
    let evals = jacobi_eigenvalues(&lap);
    if evals.is_empty() { return 0.0; }
    let total: f64 = evals.iter().map(|e| e.abs()).sum();
    let max_e = evals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if total == 0.0 { return 0.0; }
    max_e.abs() / total
}

fn laplacian_from_adj(adj: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = adj.len();
    let mut lap = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut deg = 0.0;
        for j in 0..n {
            if i != j {
                deg += adj[i][j];
                lap[i][j] = -adj[i][j];
            }
        }
        lap[i][i] = deg;
    }
    lap
}

// ── Module 4: Farey Sequence ───────────────────────────────────────────

struct FareySequence {
    n: usize,
    fractions: Vec<(u64, u64)>,
}

impl FareySequence {
    fn of_order(n: usize) -> FareySequence {
        let mut fracs = vec![(0u64, 1u64)];
        for denom in 1..=n as u64 {
            for numer in 1..=denom {
                if gcd(numer, denom) == 1 {
                    fracs.push((numer, denom));
                }
            }
        }
        fracs.push((1, 1));
        fracs.sort_by(|a, b| (a.0 as f64 / a.1 as f64).partial_cmp(&(b.0 as f64 / b.1 as f64)).unwrap());
        fracs.dedup();
        FareySequence { n, fractions: fracs }
    }

    fn farey_graph(&self) -> Vec<Vec<f64>> {
        let m = self.fractions.len();
        let mut adj = vec![vec![0.0; m]; m];
        for i in 0..m {
            for j in i+1..m {
                let (a, b) = self.fractions[i];
                let (c, d) = self.fractions[j];
                // Adjacent if |bc - ad| == 1
                if (b * c as u64).abs_diff(a * d as u64) == 1 {
                    adj[i][j] = 1.0;
                    adj[j][i] = 1.0;
                }
            }
        }
        adj
    }

    fn conservation(&self) -> f64 {
        let adj = self.farey_graph();
        graph_cr(&adj)
    }

    fn stern_brocot_tree(depth: usize) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        let mut stack = vec![(0u64, 1u64, 1u64, 0u64, 0usize)]; // (a, b, c, d, level)
        while let Some((a, b, c, d, lvl)) = stack.pop() {
            let (m_num, m_den) = Self::mediant(a, b, c, d);
            result.push((m_num, m_den));
            if lvl + 1 < depth {
                stack.push((m_num, m_den, c, d, lvl + 1));
                stack.push((a, b, m_num, m_den, lvl + 1));
            }
        }
        result
    }

    fn mediant(a: u64, b: u64, c: u64, d: u64) -> (u64, u64) {
        (a + c, b + d)
    }
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

// ── Numerics: Jacobi eigenvalue iteration ──────────────────────────────

fn jacobi_eigenvalues(mat: &[Vec<f64>]) -> Vec<f64> {
    let n = mat.len();
    if n == 0 { return vec![]; }
    let mut a = mat.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n { v[i][i] = 1.0; }

    for _ in 0..50 * n {
        // Find largest off-diagonal
        let (mut pi, mut pj) = (0, 1);
        let mut max_val = 0.0;
        for i in 0..n {
            for j in i+1..n {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    pi = i;
                    pj = j;
                }
            }
        }
        if max_val < 1e-12 { break; }

        let theta = if (a[pi][pi] - a[pj][pj]).abs() < 1e-15 {
            PI / 4.0
        } else {
            0.5 * (2.0 * a[pi][pj] / (a[pi][pi] - a[pj][pj])).atan()
        };
        let c = theta.cos();
        let s = theta.sin();

        // Rotate rows/cols pi, pj
        let mut new_a = a.clone();
        for k in 0..n {
            if k != pi && k != pj {
                new_a[pi][k] = c * a[pi][k] + s * a[pj][k];
                new_a[k][pi] = new_a[pi][k];
                new_a[pj][k] = -s * a[pi][k] + c * a[pj][k];
                new_a[k][pj] = new_a[pj][k];
            }
        }
        new_a[pi][pi] = c*c*a[pi][pi] + 2.0*s*c*a[pi][pj] + s*s*a[pj][pj];
        new_a[pj][pj] = s*s*a[pi][pi] - 2.0*s*c*a[pi][pj] + c*c*a[pj][pj];
        new_a[pi][pj] = 0.0;
        new_a[pj][pi] = 0.0;
        a = new_a;

        // Update eigenvectors
        for k in 0..n {
            let vki = v[k][pi];
            let vkj = v[k][pj];
            v[k][pi] = c * vki + s * vkj;
            v[k][pj] = -s * vki + c * vkj;
        }
    }

    let mut eigenvalues: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    eigenvalues.sort_by(|a, b| b.partial_cmp(a).unwrap()); // descending
    eigenvalues
}

// ── Main ───────────────────────────────────────────────────────────────

fn main() {
    println!("=== Penrose Lattice ===\n");

    // Substitution
    let (e1, e2) = SubstitutionMatrix::eigenvalues();
    println!("Substitution eigenvalues: φ={}, -1/φ={}", e1, e2);
    println!("Golden ratio convergence: {:?}", SubstitutionMatrix::golden_ratio_convergence(10));

    // Penrose graph
    let pg = PenroseGraph::from_inflation(4);
    println!("\nPenrose graph: {} nodes, {} edges", pg.nodes.len(), pg.edges.len());
    println!("Tile counts: thick={}, thin={}, ratio={:.6}", pg.tile_counts.0, pg.tile_counts.1, pg.tile_ratio());
    println!("CR: {:.6}", pg.conservation());

    // CR trajectory
    let traj = PenroseGraph::cr_trajectory(5);
    println!("\nCR trajectory: {:?}", traj);
    println!("Converging: {}", PenroseGraph::is_converging(&traj));

    // Farey
    let farey = FareySequence::of_order(5);
    println!("\nFarey(5): {} fractions", farey.fractions.len());
    println!("Farey CR: {:.6}", farey.conservation());

    // Stern-Brocot
    let sb = FareySequence::stern_brocot_tree(4);
    println!("Stern-Brocot depth 4: {} fractions", sb.len());

    // Inflation on a simple graph
    let simple = vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]];
    println!("\nInflation CR delta: {:.6}", InflationDeflation::inflation_cr_delta(&simple));
    let rep = InflationDeflation::repeated_inflation(&simple, 3);
    println!("Repeated inflation CRs: {:?}", rep);
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eigenvalues_are_phi() {
        let (e1, e2) = SubstitutionMatrix::eigenvalues();
        assert!((e1 - PHI).abs() < 1e-10, "e1 should be φ");
        assert!((e2 - NEG_INV_PHI).abs() < 1e-10, "e2 should be -1/φ");
    }

    #[test]
    fn test_substitution_matrix_structure() {
        let m = SubstitutionMatrix::matrix();
        assert_eq!(m[0][0], 1.0);
        assert_eq!(m[0][1], 1.0);
        assert_eq!(m[1][0], 1.0);
        assert_eq!(m[1][1], 0.0);
    }

    #[test]
    fn test_golden_ratio_convergence() {
        let ratios = SubstitutionMatrix::golden_ratio_convergence(10);
        assert!(!ratios.is_empty());
        // Last value should be close to φ
        let last = *ratios.last().unwrap();
        assert!((last - PHI).abs() < 0.05, "ratio {} should be close to φ", last);
    }

    #[test]
    fn test_apply_n_fibonacci() {
        let m = SubstitutionMatrix::apply_n(1);
        assert_eq!(m[0][0], 1);
        assert_eq!(m[0][1], 1);
        assert_eq!(m[1][0], 1);
        assert_eq!(m[1][1], 0);
        // F(5) = 5
        let m5 = SubstitutionMatrix::apply_n(5);
        assert_eq!(m5[0][0], 8); // F(6)
        assert_eq!(m5[0][1], 5); // F(5)
    }

    #[test]
    fn test_penrose_from_inflation() {
        let pg = PenroseGraph::from_inflation(2);
        assert!(pg.nodes.len() > 4, "should have more than seed nodes");
        assert!(!pg.edges.is_empty());
    }

    #[test]
    fn test_tile_ratio_approaches_phi() {
        let pg = PenroseGraph::from_inflation(8);
        let ratio = pg.tile_ratio();
        assert!((ratio - PHI).abs() < 0.1, "tile ratio {} should approach φ", ratio);
    }

    #[test]
    fn test_conservation_computed() {
        let pg = PenroseGraph::from_inflation(3);
        let cr = pg.conservation();
        assert!(cr > 0.0, "CR should be positive");
        assert!(cr <= 1.0, "CR should be ≤ 1");
    }

    #[test]
    fn test_cr_trajectory() {
        let traj = PenroseGraph::cr_trajectory(4);
        assert_eq!(traj.len(), 4);
        for cr in &traj {
            assert!(*cr >= 0.0 && *cr <= 1.0);
        }
    }

    #[test]
    fn test_laplacian_rows_sum_zero() {
        let pg = PenroseGraph::from_inflation(2);
        let lap = pg.laplacian();
        for row in &lap {
            let sum: f64 = row.iter().sum();
            assert!(sum.abs() < 1e-10, "Laplacian row should sum to 0, got {}", sum);
        }
    }

    #[test]
    fn test_inflate_graph_doubles_edges() {
        let simple = vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]];
        let inflated = InflationDeflation::inflate_graph(&simple);
        assert_eq!(inflated.len(), 6); // 3 original + 3 bridges
        let mut edges = 0;
        for i in 0..inflated.len() {
            for j in i+1..inflated.len() {
                if inflated[i][j] != 0.0 { edges += 1; }
            }
        }
        assert_eq!(edges, 6); // 3 original edges → 6 half-edges
    }

    #[test]
    fn test_inflation_cr_delta() {
        let simple = vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]];
        let delta = InflationDeflation::inflation_cr_delta(&simple);
        assert!(delta >= 0.0, "CR delta should be non-negative");
    }

    #[test]
    fn test_repeated_inflation() {
        let simple = vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]];
        let crs = InflationDeflation::repeated_inflation(&simple, 3);
        assert_eq!(crs.len(), 4);
    }

    #[test]
    fn test_farey_order_5_count() {
        let f = FareySequence::of_order(5);
        // |F(5)| = 1 + sum_{k=1}^{5} phi(k) = 1 + 1 + 1 + 2 + 2 + 4 = 11
        assert_eq!(f.fractions.len(), 11);
    }

    #[test]
    fn test_farey_adjacent_bc_ad_equals_1() {
        let f = FareySequence::of_order(5);
        for w in f.fractions.windows(2) {
            let (a, b) = w[0];
            let (c, d) = w[1];
            assert_eq!((b * c).abs_diff(a * d), 1, "Farey adjacency: bc-ad=1");
        }
    }

    #[test]
    fn test_mediant() {
        let (n, d) = FareySequence::mediant(1, 3, 2, 5);
        assert_eq!(n, 3);
        assert_eq!(d, 8);
    }

    #[test]
    fn test_farey_graph_laplacian() {
        let f = FareySequence::of_order(3);
        let adj = f.farey_graph();
        let lap = laplacian_from_adj(&adj);
        for row in &lap {
            let sum: f64 = row.iter().sum();
            assert!(sum.abs() < 1e-10);
        }
    }

    #[test]
    fn test_farey_conservation() {
        let f = FareySequence::of_order(5);
        let cr = f.conservation();
        assert!(cr > 0.0 && cr <= 1.0);
    }

    #[test]
    fn test_stern_brocot_tree() {
        let sb = FareySequence::stern_brocot_tree(3);
        assert!(!sb.is_empty());
        // All should be valid fractions
        for &(n, d) in &sb {
            assert!(d > 0);
            assert!(gcd(n, d) > 0);
        }
    }

    #[test]
    fn test_integration_penrose_cr_convergence() {
        let traj = PenroseGraph::cr_trajectory(6);
        // Trajectory should exist and be valid
        assert_eq!(traj.len(), 6);
        for cr in &traj {
            assert!(*cr > 0.0);
        }
    }

    #[test]
    fn test_integration_farey_fibonacci_indexing() {
        // Farey sequence of increasing order should grow ~ 3n²/π²
        let f3 = FareySequence::of_order(3);
        let f5 = FareySequence::of_order(5);
        assert!(f5.fractions.len() > f3.fractions.len());
    }

    #[test]
    fn test_deflation() {
        let pg = PenroseGraph::from_inflation(3);
        let deflated = pg.deflate();
        assert!(deflated.is_some());
        let d = deflated.unwrap();
        assert!(d.nodes.len() < pg.nodes.len());
    }
}
