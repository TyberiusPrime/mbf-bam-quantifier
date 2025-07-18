use anyhow::{bail, Context, Result};
use bstr::BString;
use log::info;
use petgraph::graph::{Graph, UnGraph};
use std::{
    collections::{HashMap, HashSet, VecDeque},
    io::Write,
    path::PathBuf,
};

use rust_htslib::bam;

use crate::{
    bam_ext::BamRecordExtensions,
    engine::{self, AnnotatedRead, DuplicateReadInfo},
};

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, Default)]
pub enum DeduplicationBucket {
    #[default]
    #[serde(alias = "per_position")]
    PerPosition,
    #[serde(alias = "per_reference")]
    PerReference,
    #[serde(alias = "per_region")]
    PerRegion,
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
        info!("called finish_bucket");
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
        info!("called subbucket");
        match mode {
            DeduplicationMode::NoDedup | DeduplicationMode::Unique =>
                //we already have filtered unique in accept_read()
                {}
            DeduplicationMode::Percentile => {
                dedup_percentile(map, reads); //does not change couns
            }
            DeduplicationMode::Cluster(max_hamming) => {
                dedup_cluster(map, reads, max_hamming.max_distance);
            }
            DeduplicationMode::Directional(max_hamming) => {
                dedup_directional(map, reads, max_hamming.max_distance, -1);
            }
            DeduplicationMode::DirectionalStarsolo(max_hamming) => {
                dedup_directional(map, reads, max_hamming.max_distance, 0);
            }
        }
    }
}

fn dedup_percentile(
    map: &HashMap<BString, UmiGroupInfo>,
    reads: &mut Vec<(engine::AnnotatedRead, usize)>,
) {
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
        if counts.len() > 1 {
            let medians = medians::medianu64(&mut counts[..]).unwrap_or(medians::Medians::Odd(&0));
            let median = match medians {
                medians::Medians::Odd(m) => *m,
                medians::Medians::Even((m1, m2)) => (m1 + m2) / 2,
            };

            let threshold = median / 100; //rounds towards 0, so should be fine.
            for (_umi, info) in map.iter() {
                if info.count < threshold {
                    set_read_aproximate_duplicate(reads, info.best_read_index);
                }
            }
        }
    }
}

fn set_read_aproximate_duplicate(reads: &mut Vec<(engine::AnnotatedRead, usize)>, index: usize) {
    let old = &reads[index];
    match &old.0 {
        AnnotatedRead::Counted(old_info) => {
            reads[index] = (
                engine::AnnotatedRead::ApproximateUmiDuplicate(Box::new(DuplicateReadInfo {
                    corrected_position: old_info.corrected_position,
                })),
                reads[index].1,
            );
        }
        AnnotatedRead::ApproximateUmiDuplicate(_) => {
            // we can reach this if there are multiple ways in a directional graph
            // leading to this one.
        }
        _ => {
            unreachable!();
        }
    }
}

fn dedup_cluster(
    map: &mut HashMap<BString, UmiGroupInfo>,
    reads: &mut Vec<(engine::AnnotatedRead, usize)>,
    max_hamming: u64,
) {
    if map.len() < 2 {
        return;
    }
    let mut graph = UnGraph::<&BString, ()>::new_undirected(); //todo: do not store edge
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
            graph.add_edge(*node1, *node2, ());
        }
    }
    let connected_components = connected_components_values_undirected(&graph);
    for group in connected_components.iter() {
        if group.len() > 1 {
            //find the one with the highest count...
            let max_umi = group
                .iter()
                .max_by_key(|&umi| {
                    map.get(umi).map_or(0, |info| info.count as u64) //if umi not in map, count is 0
                })
                .expect("There should be at least one UMI in the group");
            //and clear all the others
            let mut add = 0;
            for umi in group.iter() {
                if umi != max_umi {
                    let this_info = &mut map.get_mut(umi).unwrap();
                    add += this_info.count;
                    this_info.count = 0;
                    set_read_aproximate_duplicate(reads, this_info.best_read_index);
                }
            }
            //don't froget to add to the max_umi
            map.get_mut(max_umi).unwrap().count += add;
        }
    }
}

fn dedup_directional(
    map: &mut HashMap<BString, UmiGroupInfo>,
    reads: &mut Vec<(engine::AnnotatedRead, usize)>,
    max_hamming: u64,
    lower_count_offset: isize,
) {
    //create a undirected graph with edges where the hamming distance is leq than max_hamming
    //then find connected components
    //
    if map.len() < 2 {
        return;
    }
    let mut graph = Graph::<(BString, u64), ()>::new();
    let mut nodes = HashMap::new();
    for (umi, info) in map.iter() {
        let node_index = graph.add_node((umi.clone(), info.count as u64));
        nodes.insert(umi.clone(), node_index);
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
        .map(|(umi, info)| (umi.clone(), info.count, nodes.get(umi).unwrap()))
        .collect::<Vec<_>>();
    //we need to sort by count, and on ties by the umi.
    //(Can't use by_key for lifetime reasons)
    umis_sorted_by_count.sort_unstable_by(|(umi1, count1, _), (umi2, count2, _)| {
        count2.cmp(count1).then(umi1.cmp(umi2))
    });

    let mut processed = HashSet::new();
    for (umi, _count, node_index) in umis_sorted_by_count {
        if !processed.contains(&umi) {
            //identify and prune all reachable from this node.
            let mut add = 0;
            let mut seen = false;
            petgraph::visit::depth_first_search(&graph, Some(*node_index), |event| match event {
                petgraph::visit::DfsEvent::Discover(n, _) => {
                    let (connected_umi, umi_count) = &graph[n];
                    if *connected_umi != umi {
                        processed.insert(connected_umi);
                        add += umi_count;
                        let info = map.get_mut(connected_umi).unwrap();
                        info.count = 0; //set to 0, so it will be removed.
                        set_read_aproximate_duplicate(reads, info.best_read_index);
                        seen = true;
                    }
                }
                _ => {}
            });
            if seen {
                map.get_mut(&umi).unwrap().count += add;
            }
        }
    }
}

fn update_umi_counts(
    map: &mut HashMap<BString, UmiGroupInfo>,
    updated_counts: HashMap<BString, u64>,
) {
    for (umi, count) in updated_counts {
        let info = map.get_mut(&umi).unwrap();
        info.count = count;
    }
}

/// Returns a vector of components, each as a Vec<N> of node values.
pub fn connected_components_values_undirected<E>(
    graph: &UnGraph<&BString, E>,
) -> Vec<Vec<BString>> {
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

/// the BD inspired 'distribution based error correction' for UMIs counts.
/// essentially: Per biomolecule (not per cell!):
/// If enough UMIs are present at high enough levels ((average depth of those not 0) >= 4 )
/// fit two binomial distributions. Calculate a threshold between them, and then remove all
/// molecules that are below the threshold.
pub fn apply_external_threshold(
    command: &Vec<String>,
    title: &str,
    block: &mut Vec<(engine::AnnotatedRead, usize)>,
    dedup_storage: &mut DedupPerBucket,
) -> Result<()> {
    let mut umi_counts: Vec<u64> = Vec::new();
    let mut per_umi_counts: HashMap<&BString, u64> = HashMap::new();
    match dedup_storage {
        DedupPerBucket::None => {
            unreachable!("prevented by config check")
        }
        DedupPerBucket::Umi(counts) => {
            umi_counts.extend(counts.values().map(|info| info.count));
            per_umi_counts.extend(counts.iter().map(|(umi, info)| (umi, info.count)));
        }
        DedupPerBucket::SingleCell(per_barcode_counts) => {
            for (_barcode, counts) in per_barcode_counts.iter_mut() {
                umi_counts.extend(counts.values().map(|info| info.count));
                for (k, info) in counts.iter() {
                    *(per_umi_counts.entry(k).or_insert(0)) += info.count;
                }
            }
        }
    }
    if umi_counts.is_empty() {
        //no UMIs, nothing to do.
        return Ok(());
    }
    let cmd = if command[0].starts_with('/') {
        PathBuf::from(&command[0])
    } else {
        std::env::current_dir().unwrap().join(&command[0])
    };
    let mut process = std::process::Command::new(cmd)
        .args(&command[1..])
        .arg(title)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .with_context(|| {
            format!(
                "Failed to start external threshold process: {}{:?}",
                std::env::current_dir().unwrap_or_default().display(),
                &command,
            )
        })?;
    //write comma sepearated umi_counts to stdin
    let str_counts = umi_counts
        .iter()
        .filter(|x| **x > 0) // we might have 'empty' umis that were collapsed onto others in there.
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let mut stdin = process.stdin.take().expect("failed to get stdin");
    //dbg!(&str_counts);
    stdin
        .write_all(str_counts.as_bytes())
        .context("Failed to write to stdin of external threshold process")?;
    drop(stdin);
    //collect stdout and stderr
    let output = process
        .wait_with_output()
        .context("Failed to wait for external threshold process")?;
    if !output.status.success() {
        bail!(
            "External threshold process failed with status: {}. Stderr: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let stdout = String::from_utf8(output.stdout)
        .context("Failed to read stdout of external threshold process, not utf-8")?;
    //parse the output as a single number
    let threshold: u64 = stdout.trim().parse().with_context(|| {
        format!(
            "Failed to parse threshold from external threshold process output. Stdout was {stdout}, stderr was:{}",
            String::from_utf8_lossy(&output.stderr)
        )
    })?;
    info!("The externally defined UMI min-coverage threshold for {title} is {threshold}");
    //dbg!("External threshold:", threshold);

    let umis_to_filter: HashSet<&BString> = per_umi_counts
        .iter()
        .filter(|(_umi, count)| **count < threshold)
        .map(|(umi, _count)| *umi)
        .collect();

    for read in block.iter_mut() {
        if let engine::AnnotatedRead::Counted(info) = &mut read.0 {
            if let Some(umi) = info.umi.as_ref() {
                if umis_to_filter.contains(umi) {
                    *read = (engine::AnnotatedRead::RemovedByExternalUmiTreshold, read.1)
                }
            }
        }
    }

    Ok(())
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
        let mut reads: Vec<(super::engine::AnnotatedRead, usize)> = Vec::new();
        let mut insert = |umi: &[u8], count: usize| {
            map.insert(
                umi.into(),
                UmiGroupInfo {
                    best_read_index: reads.len(),
                    best_mapping_quality: super::MappingQuality {
                        no_of_alignments: 1,
                        mapq: 60,
                    },
                    count: count.try_into().expect("count should fit into u64"),
                },
            );
            reads.push((
                super::engine::AnnotatedRead::Counted(Box::new(super::engine::AnnotatedReadInfo {
                    corrected_position: 0, // clipping corrected position. Samspec is limited to 0..2^31-1,
                    hits: crate::engine::Hits {
                        correct: Vec::new(),
                        reverse: Vec::new(),
                    },
                    umi: Some(umi.into()), // Optional: What's it's UMI. 24 bytes
                    barcode: None,         // Optional: What's it's cell-barcode 24 bytes
                    mapping_priority: (0, 0),
                    reverse: false,
                })),
                reads.len(),
            ));
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
        dedup_directional(&mut map, &mut reads, 1, -1);
        assert_eq!(
            map.get::<BString>(&(b"TTCAAAAA".into())).unwrap().count,
            262
        );
        assert_eq!(map.get::<BString>(&(b"TTCGGACA".into())).unwrap().count, 88);

        let pruned_read_count = reads
            .iter()
            .filter(|(read, _index)| {
                matches!(
                    read,
                    super::engine::AnnotatedRead::ApproximateUmiDuplicate(_)
                )
            })
            .count();
        assert_eq!(pruned_read_count, 9 - 2); //since I only add one read per umi here...

        assert_eq!(
            map.get::<BString>(&(b"TTCAAACT".into())).unwrap().count,
            0 // - 153
        );

        // that one isn't being updated
        //assert_eq!(update_counts.get::<BString>(&(b"TTCGGACA".into())), 88 - 88);

        //todo check numbers
    }
}
