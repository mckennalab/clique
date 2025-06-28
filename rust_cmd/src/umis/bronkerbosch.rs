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

        for clique in bk.cliques() {
            println!("Clique {:?}", clique);
        }
    }

}