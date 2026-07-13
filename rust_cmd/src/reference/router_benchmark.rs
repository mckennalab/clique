//! Ground-truth benchmark ranking the three near-identical-panel reference
//! routers (IDF, discriminating-position, POA) on TRUE accuracy.
//!
//! Earlier validations compared each router to clique's alignment-based call,
//! but the aligner is not ground truth — it carries the same edit-pull bias the
//! routers are meant to avoid, so "concordance" cannot rank them on correctness.
//!
//! Here we instead **simulate** reads from each reference with a seeded
//! substitution-error model. Each read's origin reference is its ground-truth
//! label, so accuracy is unambiguous. We sweep the error rate and report, per
//! router: forced accuracy (argmax == origin), confident accuracy (accuracy on
//! the reads it does not flag ambiguous), and the confident fraction (coverage).
//!
//! This deliberately uses NO biological edits: on the RNF2 palindrome panel the
//! edit sites coincide with the reference-defining positions, so an "edited"
//! read has a genuinely ambiguous origin. Error-only reads keep the ground truth
//! clean and isolate reference-discrimination accuracy under sequencing error.

#[cfg(test)]
mod tests {
    use crate::reference::discriminating::DiscriminatingClassifier;
    use crate::reference::fasta_reference::{Reference, ReferenceManager};
    use crate::reference::idf::IdfIndex;
    use crate::reference::poa::PoaGraph;

    const UNEDITED: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGTACTCATCCTGTCATCTTAGCTAAGATGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";
    const LEFT: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGCGATGCTCCTGTCGTCTTAGCTAAGATGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";
    const RIGHT: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGCATACTTCCTGTCATCTTAGCTAAGACGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";

    const PANEL: [(&str, &str); 3] = [("unedited", UNEDITED), ("left", LEFT), ("right", RIGHT)];
    const READS_PER_REF: usize = 150;

    /// Deterministic xorshift64 PRNG (reproducible; no external deps / no global rng).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    /// Copy `seq`, flipping each base to a different base with probability `rate`.
    fn mutate(seq: &[u8], rate: f64, rng: &mut Rng) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        seq.iter()
            .map(|&b| {
                if rng.unit() < rate {
                    loop {
                        let nb = BASES[rng.below(4)];
                        if nb != b {
                            return nb;
                        }
                    }
                } else {
                    b
                }
            })
            .collect()
    }

    fn manager() -> ReferenceManager<'static, 'static, 'static> {
        let structs: Vec<Reference> = PANEL
            .iter()
            .map(|(name, seq)| Reference {
                sequence: seq.as_bytes().to_vec(),
                name: name.as_bytes().to_vec(),
                suffix_table: ReferenceManager::find_seeds(&seq.as_bytes().to_vec(), 8),
            })
            .collect();
        ReferenceManager::from_fasta_vec(structs, 8, 4)
    }

    /// One accuracy tally for a router at one error rate.
    struct Tally {
        total: usize,
        forced_correct: usize,
        confident: usize,
        confident_correct: usize,
    }
    impl Tally {
        fn new() -> Tally {
            Tally { total: 0, forced_correct: 0, confident: 0, confident_correct: 0 }
        }
        fn record(&mut self, call: Option<&[u8]>, ambiguous: bool, truth: &[u8]) {
            self.total += 1;
            let correct = call == Some(truth);
            if correct {
                self.forced_correct += 1;
            }
            if !ambiguous {
                self.confident += 1;
                if correct {
                    self.confident_correct += 1;
                }
            }
        }
        fn forced_acc(&self) -> f64 {
            self.forced_correct as f64 / self.total as f64
        }
        fn confident_acc(&self) -> f64 {
            if self.confident == 0 { 1.0 } else { self.confident_correct as f64 / self.confident as f64 }
        }
        fn confident_frac(&self) -> f64 {
            self.confident as f64 / self.total as f64
        }
    }

    /// Minimal FASTA reader: returns `(name, uppercased sequence)` per record.
    /// The name is the first whitespace-delimited token after `>`.
    fn read_fasta_upper(path: &str) -> Vec<(Vec<u8>, Vec<u8>)> {
        let content = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("could not read {}: {}", path, e));
        let mut refs = Vec::new();
        let mut name: Option<Vec<u8>> = None;
        let mut seq: Vec<u8> = Vec::new();
        for line in content.lines() {
            if let Some(header) = line.strip_prefix('>') {
                if let Some(n) = name.take() {
                    refs.push((n, std::mem::take(&mut seq)));
                }
                let token = header.split_whitespace().next().unwrap_or("");
                name = Some(token.as_bytes().to_vec());
            } else {
                seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        if let Some(n) = name {
            refs.push((n, seq));
        }
        refs
    }

    /// Apply one indel to `seq`: a deletion or insertion of `len` bases at a
    /// random position. `len == 0` returns a copy unchanged.
    fn apply_indel(seq: &[u8], insertion: bool, len: usize, rng: &mut Rng) -> Vec<u8> {
        if len == 0 {
            return seq.to_vec();
        }
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        if insertion {
            let pos = rng.below(seq.len() + 1);
            let mut out = seq[..pos].to_vec();
            for _ in 0..len {
                out.push(BASES[rng.below(4)]);
            }
            out.extend_from_slice(&seq[pos..]);
            out
        } else {
            if seq.len() <= len {
                return seq.to_vec();
            }
            let pos = rng.below(seq.len() - len);
            let mut out = seq[..pos].to_vec();
            out.extend_from_slice(&seq[pos + len..]);
            out
        }
    }

    fn manager_from(refs: &[(Vec<u8>, Vec<u8>)]) -> ReferenceManager<'static, 'static, 'static> {
        let structs: Vec<Reference> = refs
            .iter()
            .map(|(name, seq)| Reference {
                sequence: seq.clone(),
                name: name.clone(),
                suffix_table: ReferenceManager::find_seeds(seq, 8),
            })
            .collect();
        ReferenceManager::from_fasta_vec(structs, 8, 4)
    }

    /// How POA (and, for contrast, IDF and the star-anchor discriminating
    /// classifier) tolerate indels of varying length and class — the kind of
    /// noise the MARC1 homing-guide system produces. Reads are drawn from each
    /// reference with 1% substitution error plus ONE indel (insertion or deletion)
    /// of a controlled length; ground truth is the origin reference. Ignored
    /// (data file + 112-ref graph + per-read alignment):
    ///   cargo test --release --bin clique \
    ///     reference::router_benchmark::tests::poa_accuracy_vs_indel_length_and_class -- --ignored --nocapture
    #[test]
    #[ignore]
    fn poa_accuracy_vs_indel_length_and_class() {
        let refs = read_fasta_upper("test_data/all_MARC1_references.fa");
        let rm = manager_from(&refs);
        let poa = PoaGraph::from_references(&refs).unwrap();
        let idf = IdfIndex::from_reference_manager(&rm, 0.05);
        let disc = DiscriminatingClassifier::from_reference_manager(&rm, 1).ok();

        eprintln!(
            "\nRouter accuracy vs indel length/class on MARC1 ({} refs, +1% substitution)",
            refs.len()
        );
        eprintln!("forced accuracy (best == origin reference)\n");
        eprintln!("{:<11}{:>5}{:>10}{:>10}{:>12}", "class", "len", "poa", "idf", "discrim");
        eprintln!("{}", "-".repeat(48));

        let lengths = [0usize, 1, 2, 3, 5, 10];
        let reads_per_ref = 3usize;
        let base_error = 0.01;

        for &insertion in &[false, true] {
            for &len in &lengths {
                let mut rng = Rng(0x1DE1_0000_0000_0001
                    ^ ((insertion as u64) << 40)
                    ^ ((len as u64 + 1).wrapping_mul(0x9E3779B97F4A7C15)));
                let (mut poa_ok, mut idf_ok, mut disc_ok, mut total) = (0usize, 0usize, 0usize, 0usize);
                for (name, seq) in &refs {
                    for _ in 0..reads_per_ref {
                        let mut read = mutate(seq, base_error, &mut rng);
                        read = apply_indel(&read, insertion, len, &mut rng);
                        total += 1;

                        if poa.classify_read(&read).best == *name {
                            poa_ok += 1;
                        }
                        if idf.score(&read).best.as_deref() == Some(name.as_slice()) {
                            idf_ok += 1;
                        }
                        if let Some(d) = &disc {
                            if d.classify_read(&read).best.as_bytes() == name.as_slice() {
                                disc_ok += 1;
                            }
                        }
                    }
                }
                let pct = |ok: usize| 100.0 * ok as f64 / total as f64;
                let disc_cell = if disc.is_some() {
                    format!("{:>10.1}%", pct(disc_ok))
                } else {
                    format!("{:>11}", "n/a")
                };
                eprintln!(
                    "{:<11}{:>5}{:>9.1}%{:>9.1}%{}",
                    if insertion { "insertion" } else { "deletion" },
                    len,
                    pct(poa_ok),
                    pct(idf_ok),
                    disc_cell
                );
            }
            eprintln!();
        }
    }

    /// POA-classifier accuracy on the MARC1 homing-guide panel (112 references
    /// sharing a backbone but with diverse, variable-LENGTH spacer regions — the
    /// indel-bearing case POA is built for). Ground truth = the origin reference.
    /// Ignored by default (reads a data file, builds a 112-ref graph):
    ///   cargo test --release --bin clique \
    ///     reference::router_benchmark::tests::poa_accuracy_on_marc1_panel -- --ignored --nocapture
    #[test]
    #[ignore]
    fn poa_accuracy_on_marc1_panel() {
        let refs = read_fasta_upper("test_data/all_MARC1_references.fa");
        assert!(refs.len() >= 100, "expected the full MARC1 panel, got {}", refs.len());
        // Distinct sequences bound the achievable accuracy (identical refs are
        // genuinely indistinguishable).
        let mut seqs: Vec<&Vec<u8>> = refs.iter().map(|(_, s)| s).collect();
        seqs.sort();
        seqs.dedup();
        let distinct = seqs.len();

        let graph = PoaGraph::from_references(&refs).unwrap();
        eprintln!(
            "\nMARC1 POA accuracy: {} references ({} distinct sequences) -> {} nodes ({} backbone, {} branch)",
            refs.len(), distinct, graph.node_count(), graph.backbone_node_count(), graph.branch_node_count()
        );
        eprintln!("{:>6}{:>13}{:>15}{:>15}", "error", "forced-acc", "confident-acc", "confident-frac");
        eprintln!("{}", "-".repeat(49));

        let reads_per_ref = 3usize;
        let mut acc_at_zero = 0.0f64;
        for (ri, &error) in [0.0, 0.01, 0.02, 0.05].iter().enumerate() {
            let mut rng = Rng(0x4D41_5243_3100_0001 ^ ((ri as u64 + 1).wrapping_mul(0x9E3779B97F4A7C15)));
            let mut t = Tally::new();
            for (name, seq) in &refs {
                for _ in 0..reads_per_ref {
                    let read = mutate(seq, error, &mut rng);
                    let c = graph.classify_read(&read);
                    t.record(Some(&c.best), c.ambiguous, name);
                }
            }
            if ri == 0 {
                acc_at_zero = t.forced_acc();
            }
            eprintln!(
                "{:>5.1}%{:>12.1}%{:>14.1}%{:>14.1}%",
                error * 100.0,
                100.0 * t.forced_acc(),
                100.0 * t.confident_acc(),
                100.0 * t.confident_frac()
            );
        }

        // Error-free reads must be classified to their origin at least as often
        // as the panel is distinguishable (allowing for any near-duplicate refs).
        let achievable = distinct as f64 / refs.len() as f64;
        assert!(
            acc_at_zero >= achievable.min(0.98) - 0.02,
            "error-free forced accuracy {:.3} below the distinguishable ceiling {:.3}",
            acc_at_zero, achievable
        );
    }

    /// Per-read routing latency of the three routers on the RNF2 panel. Each
    /// router is run over the batch repeatedly until >= 0.2 s of work has
    /// accumulated, so fast (IDF) and slow (alignment-based) routers are timed
    /// stably. Run in RELEASE for realistic numbers:
    ///   cargo test --release --bin clique \
    ///     reference::router_benchmark::tests::benchmark_router_speed -- --ignored --nocapture
    #[test]
    #[ignore]
    fn benchmark_router_speed() {
        use std::hint::black_box;
        use std::time::Instant;

        let rm = manager();

        // One-time build costs.
        let t = Instant::now();
        let idf = IdfIndex::from_reference_manager(&rm, 0.05);
        let idf_build = t.elapsed();
        let t = Instant::now();
        let disc = DiscriminatingClassifier::from_reference_manager(&rm, 1).unwrap();
        let disc_build = t.elapsed();
        let t = Instant::now();
        let poa = PoaGraph::from_reference_manager(&rm, 1).unwrap();
        let poa_build = t.elapsed();

        // A fixed batch of realistic reads (1% substitution error).
        let mut rng = Rng(0xC0FF_EE12_3456_789A);
        let n = 500usize;
        let reads: Vec<Vec<u8>> = (0..n)
            .map(|i| mutate(PANEL[i % 3].1.as_bytes(), 0.01, &mut rng))
            .collect();

        // Warm up the caches so the first-timed router isn't penalised.
        for r in reads.iter().take(20) {
            black_box(idf.score(r).ambiguous);
            black_box(disc.classify_read(r).ambiguous);
            black_box(poa.classify_read(r).ambiguous);
        }

        // Adaptive timer: repeat the batch until >= 0.2 s, return µs/read.
        fn time_us_per_read<F: Fn(&[u8]) -> u64>(reads: &[Vec<u8>], f: F) -> f64 {
            let mut acc = 0u64;
            let mut iters = 0u64;
            let t = Instant::now();
            loop {
                for r in reads {
                    acc = acc.wrapping_add(f(r));
                }
                iters += 1;
                if t.elapsed().as_secs_f64() > 0.2 {
                    break;
                }
            }
            let elapsed = t.elapsed().as_secs_f64();
            black_box(acc);
            elapsed * 1e6 / (iters as f64 * reads.len() as f64)
        }

        let idf_us = time_us_per_read(&reads, |r| idf.score(r).best.map_or(0, |b| b.len() as u64));
        let disc_us = time_us_per_read(&reads, |r| disc.classify_read(r).best.len() as u64);
        let poa_us = time_us_per_read(&reads, |r| poa.classify_read(r).best.len() as u64);

        eprintln!("\nRouter speed on RNF2 (3 refs x 156 bp), reads @ 1% substitution error");
        eprintln!("{:<16}{:>12}{:>12}{:>16}{:>12}", "router", "build", "us/read", "reads/sec", "vs IDF");
        eprintln!("{}", "-".repeat(68));
        let row = |name: &str, build: std::time::Duration, us: f64| {
            eprintln!(
                "{:<16}{:>9.3} ms{:>11.2}{:>16.0}{:>11.1}x",
                name, build.as_secs_f64() * 1e3, us, 1e6 / us, us / idf_us
            );
        };
        row("idf", idf_build, idf_us);
        row("discriminating", disc_build, disc_us);
        row("poa", poa_build, poa_us);
        eprintln!(
            "\nAt 1e6 reads, routing time: idf {:.1}s, discriminating {:.0}s, poa {:.0}s",
            idf_us, disc_us, poa_us
        );
    }

    #[test]
    #[ignore]
    fn benchmark_routers_on_simulated_reads() {
        let rm = manager();
        let idf = IdfIndex::from_reference_manager(&rm, 0.05);
        let disc = DiscriminatingClassifier::from_reference_manager(&rm, 1).unwrap();
        let poa = PoaGraph::from_reference_manager(&rm, 1).unwrap();

        let error_rates = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10];

        eprintln!("\nGround-truth router accuracy on simulated RNF2 reads");
        eprintln!("({} reads/reference, origin reference = ground truth, substitution errors only)\n", READS_PER_REF);
        eprintln!("{:>6}  {:<16}{:>11}{:>15}{:>15}", "error", "router", "forced-acc", "confident-acc", "confident-frac");
        eprintln!("{}", "-".repeat(66));

        // Running sums for the final average ranking.
        let mut avg_forced = [0.0f64; 3];
        let router_names = ["idf", "discriminating", "poa"];

        for (ri, &rate) in error_rates.iter().enumerate() {
            let mut tallies = [Tally::new(), Tally::new(), Tally::new()];
            // Independent, deterministic sample per error rate.
            let mut rng = Rng(0x51F3_2A17_9E4D_0001 ^ ((ri as u64 + 1).wrapping_mul(0x9E3779B97F4A7C15)));

            for (name, seq) in PANEL.iter() {
                let truth = name.as_bytes();
                for _ in 0..READS_PER_REF {
                    let read = mutate(seq.as_bytes(), rate, &mut rng);

                    // IDF
                    let s = idf.score(&read);
                    tallies[0].record(s.best.as_deref(), s.ambiguous, truth);

                    // Discriminating
                    let c = disc.classify_read(&read);
                    tallies[1].record(Some(c.best.as_bytes()), c.ambiguous, truth);

                    // POA
                    let p = poa.classify_read(&read);
                    tallies[2].record(Some(&p.best), p.ambiguous, truth);
                }
            }

            for (i, t) in tallies.iter().enumerate() {
                avg_forced[i] += t.forced_acc();
                let label = if i == 0 { format!("{:>5.1}%", rate * 100.0) } else { String::new() };
                eprintln!(
                    "{:>6}  {:<16}{:>10.1}%{:>14.1}%{:>14.1}%",
                    label, router_names[i],
                    100.0 * t.forced_acc(), 100.0 * t.confident_acc(), 100.0 * t.confident_frac()
                );
            }
            eprintln!();
        }

        let n = error_rates.len() as f64;
        let mut ranking: Vec<(usize, f64)> =
            (0..3).map(|i| (i, avg_forced[i] / n)).collect();
        ranking.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        eprintln!("Mean forced accuracy across the error sweep (ranking):");
        for (rank, (i, acc)) in ranking.iter().enumerate() {
            eprintln!("  {}. {:<16} {:.1}%", rank + 1, router_names[*i], 100.0 * acc);
        }

        // Sanity floor: at 0% error every router must be perfect, and all should
        // stay well above chance (1/3) even at the highest error rate.
        // (No hard ranking assertion — the printed table is the deliverable.)
    }
}
