use anyhow::{Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions as htslib_record_extensions;

pub trait BamRecordExtensions {
    fn blocks(&self) -> Vec<(u32, u32)>;
    ///find intron positions (start, stop)
    #[allow(dead_code)] //todo: check
    fn introns(&self) -> Vec<(u32, u32)>;

    ///return None on unaligned reads (pos = -1...)
    fn corrected_pos(&self, max_skip_len: u32) -> Option<i32>;

    fn no_of_alignments(&self) -> u32;
    fn replace_aux(&mut self, tag: &[u8], value: bam::record::Aux) -> Result<()>;
}

impl BamRecordExtensions for bam::Record {
    fn blocks(&self) -> Vec<(u32, u32)> {
        self.aligned_blocks()
            .map(|x| (x[0] as u32, x[1] as u32))
            .collect()
    }

    fn introns(&self) -> Vec<(u32, u32)> {
        htslib_record_extensions::introns(self)
            .map(|x| (x[0] as u32, x[1] as u32))
            .collect()
    }

    /// correct the position to what it would have if there was on clipping
    /// ie. if they were treated as mismatches
    /// i32 is ok, sam raneg is 0..2^31-1, and we can have negative corrected positions
    /// if it's aligned to the start
    fn corrected_pos(&self, max_skip_len: u32) -> Option<i32> {
        let p:i32 = self.pos().try_into().expect("bam pos exceeded i32?");
        if p < 0 {
            None
        } else {
            //it's always the leading ones... since the seq gets flipped

            let skip: i32 = self.cigar().leading_softclips().try_into().expect("softclip exceeded i64");
            if skip > max_skip_len.try_into().unwrap(){
                panic!("Your reads have skipped regions > max_skip_len ({skip}>{max_skip_len}). Increase the setting via input.max_skip_length. Or filter the reads?")
            }
            Some(p.saturating_sub(skip))
        }
    }
    /// try to retrieve the number of mapping coordinates
    /// for this read. Uses the NH tag. defaults to 1
    fn no_of_alignments(&self) -> u32 {
        //let's try the NH tag, as by the tag spec
        if let Ok(nh) = self.aux(b"NH") {
            match nh {
                bam::record::Aux::U8(x) => x as u32,
                bam::record::Aux::U16(x) => x as u32,
                bam::record::Aux::U32(x) => x,
                /* bam::record::Aux::I8(x) => x as u32,
                bam::record::Aux::I16(x) => x as u32,
                bam::record::Aux::I32(x) => x, */
                _ => {
                    panic!("Mapping coordinate tag NH wasn't an unsigned int.");
                }
            }
        } else {
            1 // we can't tell.
        }
    }

    fn replace_aux(&mut self, tag: &[u8], value: bam::record::Aux) -> Result<()> {
        let has_aux = self.aux_iter().any(|res| match res {
            Ok((present_tag, _value)) => present_tag == tag,
            Err(_) => false,
        });
        if has_aux {
            self.remove_aux(tag).context("failed to remove tag")?;
        }
        self.push_aux(tag, value)
            .context("failed to push aux tag (after removing)")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::bam_ext::BamRecordExtensions;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;

    #[test]
    fn spliced_reads() {
        let mut bam = bam::Reader::from_path("sample_data/spliced_reads.bam").unwrap();
        let mut it = bam.records();
        let blocks = it.next().expect("iter").unwrap().blocks();
        //6S45M - 0
        assert!(blocks[0] == (16050676, 16050721));

        let blocks = it.next().unwrap().unwrap().blocks();
        //7M2D44M - 1
        assert!(blocks[0] == (16096878, 16096885));
        //7M2D44M - 1
        assert!(blocks[1] == (16096887, 16096931));

        let blocks = it.next().unwrap().unwrap().blocks();
        //29M2D22M - 2
        assert!(blocks[0] == (16097145, 16097174));
        //29M2D22M - 2
        assert!(blocks[1] == (16097176, 16097198));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 3
        assert!(blocks[0] == (16117350, 16117401));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 4
        assert!(blocks[0] == (16118483, 16118534));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 5
        assert!(blocks[0] == (16118499, 16118550));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 6
        assert!(blocks[0] == (16118499, 16118550));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 7
        assert!(blocks[0] == (16118499, 16118550));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 8
        assert!(blocks[0] == (16123411, 16123462));

        let blocks = it.next().unwrap().unwrap().blocks();
        //6S45M - 9
        assert!(blocks[0] == (16123417, 16123462));

        let blocks = it.next().unwrap().unwrap().blocks();
        //41M10S - 10
        assert!(blocks[0] == (16165860, 16165901));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 11
        assert!(blocks[0] == (16180871, 16180922));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 12
        assert!(blocks[0] == (16189705, 16189756));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 13
        assert!(blocks[0] == (16231271, 16231322));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 14
        assert!(blocks[0] == (16237657, 16237708));

        let blocks = it.next().unwrap().unwrap().blocks();
        //9S42M - 15
        assert!(blocks[0] == (16255012, 16255054));

        let blocks = it.next().unwrap().unwrap().blocks();
        //51M - 16
        assert!(blocks[0] == (16255391, 16255442));

        let blocks = it.next().unwrap().unwrap().blocks();
        //50M1S - 17
        assert!(blocks[0] == (16255392, 16255442));

        let blocks = it.next().unwrap().unwrap().blocks();
        //45M6S - 18
        assert!(blocks[0] == (16256084, 16256129));

        let blocks = it.next().unwrap().unwrap().blocks();
        //3S48M - 19
        assert!(blocks[0] == (16256224, 16256272));

        let blocks = it.next().unwrap().unwrap().blocks();
        //42M9S - 20
        assert!(blocks[0] == (16325199, 16325241));

        let blocks = it.next().unwrap().unwrap().blocks();
        //13S38M - 21
        assert!(blocks[0] == (16352865, 16352903));

        let blocks = it.next().unwrap().unwrap().blocks();
        //44M7S - 22
        assert!(blocks[0] == (16352968, 16353012));

        let blocks = it.next().unwrap().unwrap().blocks();
        //5S46M - 23
        assert!(blocks[0] == (16414998, 16415044));

        let blocks = it.next().unwrap().unwrap().blocks();
        //23M4I24M - 24
        assert!(blocks[0] == (17031591, 17031614));
        //23M4I24M - 24
        assert!(blocks[1] == (17031614, 17031638));

        let blocks = it.next().unwrap().unwrap().blocks();
        //18M1I32M - 25
        assert!(blocks[0] == (17057382, 17057400));
        //18M1I32M - 25
        assert!(blocks[1] == (17057400, 17057432));

        let blocks = it.next().unwrap().unwrap().blocks();
        //17M2183N34M - 26
        assert!(blocks[0] == (17092766, 17092783));
        //17M2183N34M - 26
        assert!(blocks[1] == (17094966, 17095000));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1M2183N50M - 27
        assert!(blocks[0] == (17092782, 17092783));
        //1M2183N50M - 27
        assert!(blocks[1] == (17094966, 17095016));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1M2183N50M - 28
        assert!(blocks[0] == (17092782, 17092783));
        //1M2183N50M - 28
        assert!(blocks[1] == (17094966, 17095016));

        let blocks = it.next().unwrap().unwrap().blocks();
        //9S33M9S - 29
        assert!(blocks[0] == (17137287, 17137320));

        let blocks = it.next().unwrap().unwrap().blocks();
        //2S48M1S - 30
        assert!(blocks[0] == (17306238, 17306286));

        let blocks = it.next().unwrap().unwrap().blocks();
        //4S45M2S - 31
        assert!(blocks[0] == (17561868, 17561913));

        let blocks = it.next().unwrap().unwrap().blocks();
        //41M11832N10M - 32
        assert!(blocks[0] == (17566078, 17566119));
        //41M11832N10M - 32
        assert!(blocks[1] == (17577951, 17577961));

        let blocks = it.next().unwrap().unwrap().blocks();
        //11M11832N25M710N15M - 33
        assert!(blocks[0] == (17566108, 17566119));
        //11M11832N25M710N15M - 33
        assert!(blocks[1] == (17577951, 17577976));
        //11M11832N25M710N15M - 33
        assert!(blocks[2] == (17578686, 17578701));

        let blocks = it.next().unwrap().unwrap().blocks();
        //8M11832N25M710N18M - 34
        assert!(blocks[0] == (17566111, 17566119));
        //8M11832N25M710N18M - 34
        assert!(blocks[1] == (17577951, 17577976));
        //8M11832N25M710N18M - 34
        assert!(blocks[2] == (17578686, 17578704));

        let blocks = it.next().unwrap().unwrap().blocks();
        //8M11832N25M710N18M - 35
        assert!(blocks[0] == (17566111, 17566119));
        //8M11832N25M710N18M - 35
        assert!(blocks[1] == (17577951, 17577976));
        //8M11832N25M710N18M - 35
        assert!(blocks[2] == (17578686, 17578704));

        let blocks = it.next().unwrap().unwrap().blocks();
        //8M11832N25M710N18M - 36
        assert!(blocks[0] == (17566111, 17566119));
        //8M11832N25M710N18M - 36
        assert!(blocks[1] == (17577951, 17577976));
        //8M11832N25M710N18M - 36
        assert!(blocks[2] == (17578686, 17578704));

        let blocks = it.next().unwrap().unwrap().blocks();
        //8M11832N25M710N18M - 37
        assert!(blocks[0] == (17566111, 17566119));
        //8M11832N25M710N18M - 37
        assert!(blocks[1] == (17577951, 17577976));
        //8M11832N25M710N18M - 37
        assert!(blocks[2] == (17578686, 17578704));

        let blocks = it.next().unwrap().unwrap().blocks();
        //7M11832N25M710N19M - 38
        assert!(blocks[0] == (17566112, 17566119));
        //7M11832N25M710N19M - 38
        assert!(blocks[1] == (17577951, 17577976));
        //7M11832N25M710N19M - 38
        assert!(blocks[2] == (17578686, 17578705));

        let blocks = it.next().unwrap().unwrap().blocks();
        //6M11832N25M710N20M - 39
        assert!(blocks[0] == (17566113, 17566119));
        //6M11832N25M710N20M - 39
        assert!(blocks[1] == (17577951, 17577976));
        //6M11832N25M710N20M - 39
        assert!(blocks[2] == (17578686, 17578706));

        let blocks = it.next().unwrap().unwrap().blocks();
        //6M11832N25M710N20M - 40
        assert!(blocks[0] == (17566113, 17566119));
        //6M11832N25M710N20M - 40
        assert!(blocks[1] == (17577951, 17577976));
        //6M11832N25M710N20M - 40
        assert!(blocks[2] == (17578686, 17578706));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1S44M1467N6M - 41
        assert!(blocks[0] == (17579733, 17579777));
        //1S44M1467N6M - 41
        assert!(blocks[1] == (17581244, 17581250));

        let blocks = it.next().unwrap().unwrap().blocks();
        //2M1514N48M95N1M - 42
        assert!(blocks[0] == (17581369, 17581371));
        //2M1514N48M95N1M - 42
        assert!(blocks[1] == (17582885, 17582933));
        //2M1514N48M95N1M - 42
        assert!(blocks[2] == (17583028, 17583029));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1M1514N48M95N2M - 43
        assert!(blocks[0] == (17581370, 17581371));
        //1M1514N48M95N2M - 43
        assert!(blocks[1] == (17582885, 17582933));
        //1M1514N48M95N2M - 43
        assert!(blocks[2] == (17583028, 17583030));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1M1514N48M95N2M - 44
        assert!(blocks[0] == (17581370, 17581371));
        //1M1514N48M95N2M - 44
        assert!(blocks[1] == (17582885, 17582933));
        //1M1514N48M95N2M - 44
        assert!(blocks[2] == (17583028, 17583030));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1S22M95N28M - 45
        assert!(blocks[0] == (17582911, 17582933));
        //1S22M95N28M - 45
        assert!(blocks[1] == (17583028, 17583056));

        let blocks = it.next().unwrap().unwrap().blocks();
        //37M538N13M1S - 46
        assert!(blocks[0] == (17588621, 17588658));
        //37M538N13M1S - 46
        assert!(blocks[1] == (17589196, 17589209));

        let blocks = it.next().unwrap().unwrap().blocks();
        //37M538N13M1S - 47
        assert!(blocks[0] == (17588621, 17588658));
        //37M538N13M1S - 47
        assert!(blocks[1] == (17589196, 17589209));

        let blocks = it.next().unwrap().unwrap().blocks();
        //37M538N13M1S - 48
        assert!(blocks[0] == (17588621, 17588658));
        //37M538N13M1S - 48
        assert!(blocks[1] == (17589196, 17589209));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1S25M1D25M - 49
        assert!(blocks[0] == (17591770, 17591795));
        //1S25M1D25M - 49
        assert!(blocks[1] == (17591796, 17591821));

        let blocks = it.next().unwrap().unwrap().blocks();
        //24M1D24M3S - 50
        assert!(blocks[0] == (17593855, 17593879));
        //24M1D24M3S - 50
        assert!(blocks[1] == (17593880, 17593904));

        let blocks = it.next().unwrap().unwrap().blocks();
        //16M1D28M7S - 51
        assert!(blocks[0] == (17593863, 17593879));
        //16M1D28M7S - 51
        assert!(blocks[1] == (17593880, 17593908));

        let blocks = it.next().unwrap().unwrap().blocks();
        //11S7M1I32M - 52
        assert!(blocks[0] == (17596476, 17596483));
        //11S7M1I32M - 52
        assert!(blocks[1] == (17596483, 17596515));

        let blocks = it.next().unwrap().unwrap().blocks();
        //5S9M1892N37M - 53
        assert!(blocks[0] == (17624012, 17624021));
        //5S9M1892N37M - 53
        assert!(blocks[1] == (17625913, 17625950));

        let blocks = it.next().unwrap().unwrap().blocks();
        //2S9M1892N40M - 54
        assert!(blocks[0] == (17624012, 17624021));
        //2S9M1892N40M - 54
        assert!(blocks[1] == (17625913, 17625953));

        let blocks = it.next().unwrap().unwrap().blocks();
        //1S7M3D19M2285N24M - 55
        assert!(blocks[0] == (31796700, 31796707));
        //1S7M3D19M2285N24M - 55
        assert!(blocks[1] == (31796710, 31796729));
        //1S7M3D19M2285N24M - 55
        assert!(blocks[2] == (31799014, 31799038));

        let blocks = it.next().unwrap().unwrap().blocks();
        //14M799N28M13881N7M2S - 56
        assert!(blocks[0] == (36722692, 36722706));
        //14M799N28M13881N7M2S - 56
        assert!(blocks[1] == (36723505, 36723533));
        //14M799N28M13881N7M2S - 56
        assert!(blocks[2] == (36737414, 36737421));

        let blocks = it.next().unwrap().unwrap().blocks();
        //4S21M1696N23M2331N3M - 57
        assert!(blocks[0] == (44587963, 44587984));
        //4S21M1696N23M2331N3M - 57
        assert!(blocks[1] == (44589680, 44589703));
        //4S21M1696N23M2331N3M - 57
        assert!(blocks[2] == (44592034, 44592037));
    }

    #[test]
    fn test_introns() {
        let mut bam = bam::Reader::from_path("sample_data/spliced_reads.bam").unwrap();
        let mut it = bam.records();

        //6S45M - 0
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //7M2D44M - 1
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //29M2D22M - 2
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 3
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 4
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 5
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 6
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 7
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 8
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //6S45M - 9
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //41M10S - 10
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 11
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 12
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 13
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 14
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //9S42M - 15
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //51M - 16
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //50M1S - 17
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //45M6S - 18
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //3S48M - 19
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //42M9S - 20
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //13S38M - 21
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //44M7S - 22
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //5S46M - 23
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //23M4I24M - 24
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //18M1I32M - 25
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //17M2183N34M - 26
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17092783, 17094966));
        //1M2183N50M - 27
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17092783, 17094966));
        //1M2183N50M - 28
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17092783, 17094966));
        //9S33M9S - 29
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //2S48M1S - 30
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //4S45M2S - 31
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //41M11832N10M - 32
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17566119, 17577951));
        //11M11832N25M710N15M - 33
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //8M11832N25M710N18M - 34
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //8M11832N25M710N18M - 35
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //8M11832N25M710N18M - 36
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //8M11832N25M710N18M - 37
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //7M11832N25M710N19M - 38
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //6M11832N25M710N20M - 39
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //6M11832N25M710N20M - 40
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17566119, 17577951));
        assert_eq!(introns[1], (17577976, 17578686));
        //1S44M1467N6M - 41
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17579777, 17581244));
        //2M1514N48M95N1M - 42
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17581371, 17582885));
        assert_eq!(introns[1], (17582933, 17583028));
        //1M1514N48M95N2M - 43
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17581371, 17582885));
        assert_eq!(introns[1], (17582933, 17583028));
        //1M1514N48M95N2M - 44
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (17581371, 17582885));
        assert_eq!(introns[1], (17582933, 17583028));
        //1S22M95N28M - 45
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17582933, 17583028));
        //37M538N13M1S - 46
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17588658, 17589196));
        //37M538N13M1S - 47
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17588658, 17589196));
        //37M538N13M1S - 48
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17588658, 17589196));
        //1S25M1D25M - 49
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //24M1D24M3S - 50
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //16M1D28M7S - 51
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //11S7M1I32M - 52
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 0);
        //5S9M1892N37M - 53
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17624021, 17625913));
        //2S9M1892N40M - 54
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (17624021, 17625913));
        //1S7M3D19M2285N24M - 55
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (31796729, 31799014));
        //14M799N28M13881N7M2S - 56
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (36722706, 36723505));
        assert_eq!(introns[1], (36723533, 36737414));
        //4S21M1696N23M2331N3M - 57
        let introns = it.next().unwrap().unwrap().introns();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (44587984, 44589680));
        assert_eq!(introns[1], (44589703, 44592034));
    }
}
