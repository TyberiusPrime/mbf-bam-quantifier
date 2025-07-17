use bstr::BString;
use petgraph::graph::{Graph, UnGraph};
use std::collections::{HashMap, HashSet, VecDeque};

use rust_htslib::bam;

use crate::{bam_ext::BamRecordExtensions, engine};

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, Default)]
pub enum DeduplicationBucket {
    #[default]
    #[serde(alias = "per_position")]
    PerPosition,
    #[serde(alias = "per_reference")]
    PerReference,
}

#[derive(serde::Deserialize, Debug, Clone, Copy)]
#[serde(deny_unknown_fields)]
pub struct HammingDistance {
    #[serde(alias = "max_hamming")]
    max_distance: u64,
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, Copy)]
#[serde(tag = "kind")]
pub enum DeduplicationMode {
    #[serde(alias = "none")]
    NoDedup,

    #[serde(alias = "unique")]
    Unique,

    #[serde(alias = "percentile")]
    Percentile,

    #[serde(alias = "cluster")]
    Cluster(HammingDistance),

    #[serde(alias = "directional")]
    Directional(HammingDistance),

    #[serde(alias = "directional_starsolo")]
    #[serde(alias = "Directional_STARsolo")]
    DirectionalStarsolo(HammingDistance),
}

#[derive(PartialEq, Eq, Debug)]
pub struct MappingQuality {
    no_of_alignments: u8,
    mapq: u8,
    //consider deciding on cigar length as well?
}

impl Ord for MappingQuality {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.no_of_alignments
            .cmp(&other.no_of_alignments)
            .then(self.mapq.cmp(&other.mapq))
    }
}

impl PartialOrd for MappingQuality {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug)]
pub struct UmiGroupInfo {
    best_read_index: usize,
    best_mapping_quality: MappingQuality,
    count: u64, //todo: save bytes?
}

pub enum DedupPerBucket {
    None,
    Umi(HashMap<BString, UmiGroupInfo>),
    SingleCell(HashMap<BString, HashMap<BString, UmiGroupInfo>>),
}

pub enum AcceptReadResult {
    Duplicated,
    New,
    DuplicateButPrefered(usize),
}

impl DedupPerBucket {
    pub fn new(umi: bool, single_cell: bool) -> Self {
        match (umi, single_cell) {
            (true, true) => DedupPerBucket::SingleCell(HashMap::new()),
            (true, false) => DedupPerBucket::Umi(HashMap::new()),
            (false, true) => unimplemented!(), // TODO, I guess
            (false, false) => DedupPerBucket::None,
        }
    }

    pub fn accept_read_inner(
        map: &mut HashMap<BString, UmiGroupInfo>,
        this_index: usize,
        read: &bam::record::Record,
        umi: &BString,
    ) -> AcceptReadResult {
        let hit = map.get_mut(umi);
        let this_mq = MappingQuality {
            no_of_alignments: read.no_of_alignments().try_into().unwrap_or(255),
            mapq: read.mapq(),
        };
        match hit {
            Some(UmiGroupInfo {
                best_read_index: old_index,
                best_mapping_quality: mq,
                count,
            }) => {
                *count += 1;
                if this_mq > *mq {
                    *mq = this_mq;
                    let result = AcceptReadResult::DuplicateButPrefered(*old_index);
                    *old_index = this_index;
                    result
                } else {
                    AcceptReadResult::Duplicated
                }
            }
            None => {
                map.insert(
                    umi.clone(),
                    UmiGroupInfo {
                        best_read_index: this_index,
                        best_mapping_quality: this_mq,
                        count: 1,
                    },
                );
                AcceptReadResult::New
            }
        }
    }

    pub fn accept_read(
        &mut self,
        read: &bam::record::Record,
        this_index: usize,
        umi: Option<&BString>,
        barcode: Option<&BString>,
    ) -> AcceptReadResult {
        match self {
            DedupPerBucket::None => AcceptReadResult::New,
            DedupPerBucket::Umi(map) => Self::accept_read_inner(
                map,
                this_index,
                read,
                umi.expect("UMI should be extracted before deduplication"),
            ),
            DedupPerBucket::SingleCell(map) => {
                let barcode = barcode
                    .expect("Barcode should be extracted before deduplication")
                    .as_slice();

                let by_barcode = map.entry(barcode.into());
                let by_barcode = match by_barcode {
                    std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                    std::collections::hash_map::Entry::Vacant(e) => e.insert(HashMap::new()),
                };
                Self::accept_read_inner(
                    by_barcode,
                    this_index,
                    read,
                    umi.expect("UMI should be extracted before deduplication"),
                )
            }
        }
    }

    pub fn finish_bucket(
        &mut self,
        reads: &mut Vec<(engine::AnnotatedRead, usize)>,
        mode: DeduplicationMode,
    ) {
        match self {
            &mut DedupPerBucket::None => {}
            &mut DedupPerBucket::Umi(ref mut map) => {
                Self::finish_umi_subbucket(mode, map, reads);
            }
            &mut DedupPerBucket::SingleCell(ref mut by_barcode) => {
                for (_barcode, map) in by_barcode.iter_mut() {
                    Self::finish_umi_subbucket(mode, map, reads);
                }
            }
        }
    }

    fn finish_umi_subbucket(
        mode: DeduplicationMode,
        map: &mut HashMap<BString, UmiGroupInfo>,
        reads: &mut Vec<(engine::AnnotatedRead, usize)>,
    ) {
        match mode {
            DeduplicationMode::NoDedup | DeduplicationMode::Unique => {}
            DeduplicationMode::Percentile => {
                let filter = dedup_percentile(map); //does not change couns
                filter_reads_with_those_umis(reads, filter)
            }
            DeduplicationMode::Cluster(max_hamming) => {
                let (filter, add_counts) = dedup_cluster(map, max_hamming.max_distance);
                filter_reads_with_those_umis(reads, filter);
                update_umi_counts(map, add_counts);
            }
            DeduplicationMode::Directional(max_hamming) => {
                let (filter, add_counts) = dedup_directional(map, max_hamming.max_distance, -1);
                filter_reads_with_those_umis(reads, filter);
                update_umi_counts(map, add_counts);
            }
            DeduplicationMode::DirectionalStarsolo(max_hamming) => {
                let (filter, add_counts) = dedup_directional(map, max_hamming.max_distance, 0);
                filter_reads_with_those_umis(reads, filter);
                update_umi_counts(map, add_counts);
            }
        }
    }
}

fn filter_reads_with_those_umis(
    reads: &mut Vec<(engine::AnnotatedRead, usize)>,
    filter: HashSet<&BString>,
) {
    for (read, _org_index) in reads.iter_mut() {
        match read {
            engine::AnnotatedRead::Counted(info) => {
                if let Some(umi) = info.umi.as_ref() {
                    if filter.contains(umi) {
                        *read = engine::AnnotatedRead::new_approximate_duplicate(
                            info.corrected_position,
                        );
                    }
                }
            }
            _ => {}
        }
    }
}

fn dedup_percentile(map: &HashMap<BString, UmiGroupInfo>) -> HashSet<&BString> {
    {
        //UMI. UMIs with counts < 1% of the median counts for UMIs at the same position are ignored.
        const _: () = assert!(
            std::mem::size_of::<usize>() == 8,
            "usize isn't u64, this needs adjusting"
        );
        let mut counts: Vec<u64> = map
            .values()
            .map(
                |v| v.count as u64, //if usize > u64
                                    //the const assert above should have exploded
            )
            .collect();
        if counts.len() == 1 {
            return HashSet::new();
        }
        let medians = medians::medianu64(&mut counts[..]).unwrap_or(medians::Medians::Odd(&0));
        let median = match medians {
            medians::Medians::Odd(m) => *m,
            medians::Medians::Even((m1, m2)) => (m1 + m2) / 2,
        };

        let threshold = median / 100; //rounds towards 0, so should be fine.
        if counts.iter().any(|x| *x > 1) {
            dbg!(&counts);
            dbg!(median);
            dbg!(&threshold);
        }
        let filter: HashSet<_> = map
            .iter()
            .filter(|(_k, v)| (v.count as u64) <= threshold)
            .map(|(k, _v)| k)
            .collect();
        filter
    }
}

fn dedup_cluster(
    map: &HashMap<BString, UmiGroupInfo>,
    max_hamming: u64,
) -> (HashSet<&BString>, HashMap<BString, u64>) {
    if map.len() < 2 {
        return (HashSet::new(), HashMap::new()); //no need to cluster if there is only one UMI
                                                 //TODO: push upward.
    }
    let mut graph = UnGraph::<&BString, u64>::new_undirected(); //todo: do not store edge
                                                                //weights...
    let mut nodes = HashMap::new();
    for (umi, _info) in map.iter() {
        let node_index = graph.add_node(umi);
        nodes.insert(umi, node_index);
    }
    for (umi1, umi2) in map.keys().flat_map(|umi1| {
        map.keys()
            .filter(move |umi2| umi1 < *umi2)
            .map(move |umi2| (umi1, umi2))
    }) {
        let dist = bio::alignment::distance::hamming(umi1, umi2);
        if dist <= max_hamming {
            let node1 = nodes.get(umi1).expect("UMI should be in the graph");
            let node2 = nodes.get(umi2).expect("UMI should be in the graph");
            graph.add_edge(*node1, *node2, dist);
        }
    }
    connected_components_to_filter(connected_components_values_undirected(&graph), map)
}

fn dedup_directional(
    map: &HashMap<BString, UmiGroupInfo>,
    max_hamming: u64,
    lower_count_offset: isize,
) -> (HashSet<&BString>, HashMap<BString, u64>) {
    //create a undirected graph with edges where the hamming distance is leq than max_hamming
    //then find connected components
    //
    if map.len() < 2 {
        return (HashSet::new(), HashMap::new()); //no need to cluster if there is only one UMI
    }
    let mut graph = Graph::<(&BString, u64), ()>::new();
    let mut nodes = HashMap::new();
    for (umi, info) in map.iter() {
        let node_index = graph.add_node((umi, info.count as u64));
        nodes.insert(umi, node_index);
    }
    for (umi1, umi2) in map.keys().flat_map(|umi1| {
        map.keys()
            .filter(move |umi2| umi1 < *umi2)
            .map(move |umi2| (umi1, umi2))
    }) {
        let counts_1: isize = map
            .get(umi1)
            .map_or(0, |info| info.count)
            .try_into()
            .expect("counts exceeded isize");
        let counts_2: isize = map
            .get(umi2)
            .map_or(0, |info| info.count)
            .try_into()
            .expect("counts exceeded isize");
        let a_exceeds_b = counts_1 >= (2 * counts_2) + lower_count_offset;
        let b_exceeds_a = counts_2 >= (2 * counts_1) + lower_count_offset;

        if a_exceeds_b || b_exceeds_a {
            let dist = bio::alignment::distance::hamming(umi1, umi2);

            if a_exceeds_b && dist <= max_hamming {
                let node1 = nodes.get(umi1).expect("UMI should be in the graph");
                let node2 = nodes.get(umi2).expect("UMI should be in the graph");
                graph.add_edge(*node1, *node2, ());
            }
            //these can happen at the same time. for example with both counts being 1.
            if b_exceeds_a && dist <= max_hamming {
                let node1 = nodes.get(umi1).expect("UMI should be in the graph");
                let node2 = nodes.get(umi2).expect("UMI should be in the graph");
                graph.add_edge(*node2, *node1, ());
            }
        }
    }

    let mut umis_sorted_by_count = map
        .iter()
        .map(|(umi, info)| (umi, info.count, nodes.get(umi).unwrap()))
        .collect::<Vec<_>>();
    //we need to sort by count, and on ties by the umi.
    //(Can't use by_key for lifetime reasons)
    umis_sorted_by_count.sort_unstable_by(|(umi1, count1, _), (umi2, count2, _)| {
        count2.cmp(count1).then(umi1.cmp(umi2))
    });

    let mut filter = HashSet::new();
    let mut update_counts = HashMap::new();
    for (umi, _count, node_index) in umis_sorted_by_count {
        if !filter.contains(umi) {
            //identify and prune all reachable from this node.
            let mut add = 0;
            petgraph::visit::depth_first_search(&graph, Some(*node_index), |event| match event {
                petgraph::visit::DfsEvent::Discover(n, _) => {
                    let (connected_umi, umi_count) = graph[n];
                    if connected_umi != umi {
                        filter.insert(connected_umi);
                        add += umi_count;
                    }
                }
                _ => {}
            });
            if add > 0 {
                update_counts.insert(umi.clone(), add);
            }
        }
    }
    (filter, update_counts)
}

fn connected_components_to_filter<'a>(
    connected_components: Vec<Vec<&'a BString>>,
    map: &'a HashMap<BString, UmiGroupInfo>,
) -> (HashSet<&'a BString>, HashMap<BString, u64>) {
    let mut filter = HashSet::new();
    let mut count_updates = HashMap::new();

    for group in connected_components.iter() {
        if group.len() > 1 {
            //find the one with the highest count...
            let max_umi = group
                .iter()
                .max_by_key(|&&umi| {
                    map.get(umi).map_or(0, |info| info.count as u64) //if umi not in map, count is 0
                })
                .expect("There should be at least one UMI in the group");
            //and add all others to our filter set
            let mut add = 0;
            for umi in group.iter() {
                if umi != max_umi {
                    filter.insert(*umi);
                    add += map.get(*umi).map(|info| info.count as u64).unwrap();
                }
            }
            count_updates.insert((**max_umi).clone(), add);
        }
    }
    (filter, count_updates)
}

fn update_umi_counts(
    map: &mut HashMap<BString, UmiGroupInfo>,
    updated_counts: HashMap<BString, u64>,
) {
    for (umi, count) in updated_counts {
        let info = map.get_mut(&umi).unwrap();
        info.count += count;
    }
}

/// Returns a vector of components, each as a Vec<N> of node values.
pub fn connected_components_values_undirected<N: Clone, E>(graph: &UnGraph<N, E>) -> Vec<Vec<N>> {
    let mut visited = HashSet::new();
    let mut components = Vec::new();

    for start in graph.node_indices() {
        if visited.contains(&start) {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited.insert(start);

        while let Some(node) = queue.pop_front() {
            component.push(graph[node].clone());
            for neighbor in graph.neighbors(node) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }

        components.push(component);
    }

    components
}

pub fn connected_components_values_directed<N: Clone, E>(graph: &Graph<N, E>) -> Vec<Vec<N>> {
    let mut visited = HashSet::new();
    let mut components = Vec::new();

    for start in graph.node_indices() {
        if visited.contains(&start) {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited.insert(start);

        while let Some(node) = queue.pop_front() {
            component.push(graph[node].clone());
            for neighbor in graph.neighbors(node) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }

        components.push(component);
    }

    components
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use bstr::BString;

    use super::{dedup_directional, UmiGroupInfo};

    /// test that directional is +- equivalent to RSEC, based on the example given in the
    /// BD Rhapsody™ Sequence Analysis Pipeline User's Guide
    #[test]
    fn test_directional_is_rsec() {
        let mut map: HashMap<BString, UmiGroupInfo> = HashMap::new();
        let mut insert = |umi: &[u8], count: usize| {
            map.insert(
                umi.into(),
                UmiGroupInfo {
                    best_read_index: 0,
                    best_mapping_quality: super::MappingQuality {
                        no_of_alignments: 1,
                        mapq: 60,
                    },
                    count: count.try_into().expect("count should fit into u64"),
                },
            );
        };
        insert(b"GTCAAATT", 3);
        insert(b"GTCAAAAT", 24);
        insert(b"GTCAAAAA", 74);
        insert(b"TTCAAAAA", 153);
        insert(b"TTCAGAAA", 1);
        insert(b"CTCAAAAA", 2);
        insert(b"TTCAAAAT", 4);
        insert(b"TTCAAACT", 1);
        insert(b"TTCGGACA", 88);
        let (filter, update_counts) = dedup_directional(&map, 1, -1);
        //super::filter_reads_with_those_umis(reads, filter);
        dbg!(&update_counts);
        assert_eq!(update_counts.len(), 1);

        assert!(!filter.contains::<BString>(&(b"TTCAAAAA".into())));
        assert!(!filter.contains::<BString>(&(b"TTCGGACA".into())));
        assert_eq!(filter.len(), 7);
        super::update_umi_counts(&mut map, update_counts);
        assert_eq!(
            map.get::<BString>(&(b"TTCAAAAA".into())).unwrap().count,
            262 // - 153
        );

        // that one isn't being updated
        //assert_eq!(update_counts.get::<BString>(&(b"TTCGGACA".into())), 88 - 88);

        //todo check numbers
    }
}
