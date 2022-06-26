use std::convert::TryInto;
use std::io::Cursor;

use twobit::convert::fasta::{to_fasta, FastaReader};
use twobit::convert::to_2bit;
use twobit::TwoBitFile;

fn u64_to_usize(v: u64) -> usize {
    v.try_into().expect("testing")
}

#[test]
fn fasta_to_twobit() {
    let seq1 = "NNNNNNAGTCAGT\nCGTCAGCTAGTacgtcagtctgacgtTACG\nTGCATGCGATNNNNACGTANNN";
    let seq2 = "gtca\ngtcgtagtgctaTGACG\nTAGCTGACGT\nagtcgtc\ngctaggcACGTCGTAGCTACGT";
    let seq3 = "AGTCGTAG\nNCTNGNNN\nNNACGTTAGANGTCGTCAGTNGTACNCGNTATCNGTG\nNCNNNN";
    let data = format!(">seq1\n{}\n>seq2\n{}\n>seq3\n{}", seq1, seq2, seq3);
    let fasta: Vec<u8> = data.as_bytes().into();
    let reader = FastaReader::mem_open(fasta).expect("testing");
    let mut out = Cursor::new(vec![]);
    to_2bit(&mut out, &reader).expect("testing");
    let twobit_file_data: Vec<u8> = out.into_inner();

    let twobit = TwoBitFile::from_buf(twobit_file_data).expect("testing");
    for (i, info) in twobit.sequence_info().iter().enumerate() {
        match i {
            0 => {
                assert_eq!(info.chr, "seq1");
                assert_eq!(info.length, 65);
                assert_eq!(info.hard_masks_total_length, 13);
                assert_eq!(info.soft_masks_total_length, 15);
            }
            1 => {
                assert_eq!(info.chr, "seq2");
                assert_eq!(info.length, 60);
                assert_eq!(info.hard_masks_total_length, 0);
                assert_eq!(info.soft_masks_total_length, 30);
            }
            2 => {
                assert_eq!(info.chr, "seq3");
                assert_eq!(info.length, 59);
                assert_eq!(info.hard_masks_total_length, 17);
                assert_eq!(info.soft_masks_total_length, 0);
            }
            _ => assert!(false),
        }
    }
    assert_eq!(twobit.sequence_info().len(), 3);
}

#[test]
fn twobit_fasta_roundtrip() {
    // generate Fasta file
    let mut reader = TwoBitFile::open("assets/foo.2bit")
        .expect("testing")
        .enable_softmask(true);
    let mut fasta_buf = Cursor::new(Vec::with_capacity(300));
    to_fasta(&mut reader, &mut fasta_buf).expect("testing");

    // generate 2bit file
    let mut fasta_reader = FastaReader::mem_open(fasta_buf.into_inner()).expect("testing");
    let mut twobit_buf = Cursor::new(Vec::with_capacity(u64_to_usize(
        reader.info().expect("testing").file_size,
    )));
    to_2bit(&mut twobit_buf, &mut fasta_reader).expect("testing");

    // read the original file and compare the generated in-memory file byte by byte
    let original = std::fs::read("assets/foo.2bit").expect("testing");
    assert_eq!(original.len(), twobit_buf.get_ref().len());
    for (old, new) in original.iter().zip(twobit_buf.get_ref()) {
        assert_eq!(old, new);
    }
}
