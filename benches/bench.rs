use criterion::{black_box, criterion_group, criterion_main, Criterion};

use twobit::TwoBitFile;

pub fn criterion_benchmark(c: &mut Criterion) {
    let open_tb = || TwoBitFile::open_and_read("assets/foo.2bit").unwrap();
    c.bench_function("sequence", |b| {
        let mut tb = open_tb();
        b.iter(|| black_box(tb.sequence("chr1", 24..74).unwrap().len()))
    });
    c.bench_function("full_sequence", |b| {
        let mut tb = open_tb();
        b.iter(|| black_box(tb.sequence("chr1", ..).unwrap().len()))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
