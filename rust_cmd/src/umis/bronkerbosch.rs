use std::collections::HashSet;

use petgraph::graphmap::{GraphMap, NodeTrait};
use petgraph::Undirected;

/// Implementation according to "Algorithm 457: Finding All Cliques of an Undirected Graph"
/// by Bron and Kerbosch; http://doi.acm.org/10.1145/362342.362367
///
/// connected is a symmetrical boolean matrix, N the number of nodes in the graph,
/// values of the diagonal should be true.
#[allow(dead_code)]
pub struct BronKerbosch<N: NodeTrait, E> {
    pub graph: GraphMap<N, E, Undirected>,
    pub max_cliques: Vec<HashSet<N>>,
}

#[allow(dead_code)]
impl<N: NodeTrait, E> BronKerbosch<N, E> {
    pub fn new(graphmap: GraphMap<N, E, Undirected>) -> BronKerbosch<N, E> {
        BronKerbosch {
            graph: graphmap,
            max_cliques: Vec::new(),
        }
    }

    pub fn compute(&mut self) {
        let p = self.graph.nodes().collect::<HashSet<N>>();
        let r = HashSet::new();
        let x = HashSet::new();
        self.bronkerbosch(p, r, x);
    }

    #[allow(dead_code)]
    pub fn cliques(&self) -> &Vec<HashSet<N>> {
        &self.max_cliques
    }


    fn bronkerbosch(&mut self, p: HashSet<N>, r: HashSet<N>, x: HashSet<N>) {
        let mut p_fp = p.clone();
        let mut x_fp = x.clone();

        if p.is_empty() {
            if x.is_empty() {
                self.max_cliques.push(r.clone());
            }
            return;
        }

        for v in p.iter() {
            let v_neighbours = self.graph.neighbors(v.clone()).collect::<HashSet<N>>();

            let p_intersect_v_neighbors = p_fp.intersection(&v_neighbours).cloned().collect();
            let mut r_union_v = r.clone();
            r_union_v.insert(v.clone());
            let x_intersect_v_neighbors = x_fp.intersection(&v_neighbours).cloned().collect();

            self.bronkerbosch(p_intersect_v_neighbors, r_union_v, x_intersect_v_neighbors);

            p_fp.remove(v);
            x_fp.insert(*v);
        }
    }
}


#[cfg(test)]
mod tests {
    use petgraph::graphmap::UnGraphMap;
    use super::*;
    #[test]
    fn simple_wikipedia_test() {

        let mut g = UnGraphMap::new();

        // example from https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        g.add_edge("6", "4", 1);
        g.add_edge("4", "5", 1);
        g.add_edge("4", "3", 1);
        g.add_edge("3", "2", 1);
        g.add_edge("5", "2", 1);
        g.add_edge("5", "1", 1);
        g.add_edge("2", "1", 1);

        let mut bk = BronKerbosch::new(g);
        bk.compute();

        // Wikipedia example maximal cliques: {1,2,5}, {2,3}, {3,4}, {4,5}, {4,6}
        let cliques = bk.cliques();
        assert_eq!(cliques.len(), 5);

        let has_125 = cliques.iter().any(|c| c.contains(&"1") && c.contains(&"2") && c.contains(&"5") && c.len() == 3);
        assert!(has_125, "Expected clique {{1, 2, 5}}");

        let has_23 = cliques.iter().any(|c| c.contains(&"2") && c.contains(&"3") && c.len() == 2);
        assert!(has_23, "Expected clique {{2, 3}}");

        let has_34 = cliques.iter().any(|c| c.contains(&"3") && c.contains(&"4") && c.len() == 2);
        assert!(has_34, "Expected clique {{3, 4}}");

        let has_45 = cliques.iter().any(|c| c.contains(&"4") && c.contains(&"5") && c.len() == 2);
        assert!(has_45, "Expected clique {{4, 5}}");

        let has_46 = cliques.iter().any(|c| c.contains(&"4") && c.contains(&"6") && c.len() == 2);
        assert!(has_46, "Expected clique {{4, 6}}");
    }

    #[test]
    fn test_empty_graph() {
        // Empty graph: P=empty, X=empty → algorithm pushes empty R as a clique
        let g: UnGraphMap<u32, u32> = UnGraphMap::new();
        let mut bk = BronKerbosch::new(g);
        bk.compute();
        assert_eq!(bk.cliques().len(), 1);
        assert!(bk.cliques()[0].is_empty());
    }

    #[test]
    fn test_single_node() {
        let mut g = UnGraphMap::new();
        g.add_node(1u32);
        let mut bk: BronKerbosch<u32, u32> = BronKerbosch::new(g);
        bk.compute();
        assert_eq!(bk.cliques().len(), 1);
        assert!(bk.cliques()[0].contains(&1u32));
    }

    #[test]
    fn test_single_edge() {
        let mut g = UnGraphMap::new();
        g.add_edge(1u32, 2u32, 1);
        let mut bk = BronKerbosch::new(g);
        bk.compute();
        assert_eq!(bk.cliques().len(), 1);
        assert!(bk.cliques()[0].contains(&1u32));
        assert!(bk.cliques()[0].contains(&2u32));
    }

    #[test]
    fn test_complete_graph_k4() {
        let mut g = UnGraphMap::new();
        // Complete graph on 4 vertices
        for i in 0u32..4 {
            for j in (i + 1)..4 {
                g.add_edge(i, j, 1);
            }
        }
        let mut bk = BronKerbosch::new(g);
        bk.compute();
        // K4 has exactly one maximal clique of size 4
        assert_eq!(bk.cliques().len(), 1);
        assert_eq!(bk.cliques()[0].len(), 4);
    }

    #[test]
    fn test_triangle() {
        let mut g = UnGraphMap::new();
        g.add_edge(1u32, 2u32, 1);
        g.add_edge(2u32, 3u32, 1);
        g.add_edge(1u32, 3u32, 1);
        let mut bk = BronKerbosch::new(g);
        bk.compute();
        assert_eq!(bk.cliques().len(), 1);
        assert_eq!(bk.cliques()[0].len(), 3);
    }

    #[test]
    fn test_disconnected_edges() {
        let mut g = UnGraphMap::new();
        g.add_edge(1u32, 2u32, 1);
        g.add_edge(3u32, 4u32, 1);
        let mut bk = BronKerbosch::new(g);
        bk.compute();
        assert_eq!(bk.cliques().len(), 2);
        // Each clique should be size 2
        for c in bk.cliques() {
            assert_eq!(c.len(), 2);
        }
    }

    #[test]
    fn test_star_graph() {
        // Star: center connected to 4 outer nodes, no edges between outer nodes
        let mut g = UnGraphMap::new();
        for i in 1u32..=4 {
            g.add_edge(0, i, 1);
        }
        let mut bk = BronKerbosch::new(g);
        bk.compute();
        // Each edge {0, i} forms a maximal clique
        assert_eq!(bk.cliques().len(), 4);
        for c in bk.cliques() {
            assert_eq!(c.len(), 2);
            assert!(c.contains(&0));
        }
    }

}