#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use vcf::*;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter, Write};
use rayon::prelude::*;
use clap::Parser;
use bgzip::BGZFWriter;

pub struct FileReader {
    _filename: String,
    reader : VCFReader<BufReader<MultiGzDecoder<File>>>,
    vcf_record: VCFRecord
}

impl FileReader {
    pub fn new(filename: &str) -> FileReader {
        let reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(
            filename,
        ).unwrap()))).unwrap();
        let vcf_record = VCFRecord::new(reader.header().to_owned());
        FileReader { 
            reader, 
            _filename:filename.to_string(), 
            vcf_record
        }
    }

    pub fn load_next(&mut self) -> usize {
        match self.reader.next_record(&mut self.vcf_record) {
            Ok(lines_read) => if lines_read {1} else {0}
            _ => 0
        }
    }

    pub fn check_meta(&self, origin: &Self) -> bool {
        (origin.vcf_record.chromosome == self.vcf_record.chromosome) &&
        (origin.vcf_record.position == self.vcf_record.position) &&
        (origin.vcf_record.id == self.vcf_record.id) &&
        (origin.vcf_record.reference == self.vcf_record.reference) &&
        (origin.vcf_record.alternative == self.vcf_record.alternative)
    }
}

#[derive(Parser, Debug)]
#[clap(author="Magnus Manske <magnusmanske@googlemail.com>", version="0.0.2", about, long_about = None)]
/// Merge large numbers of VCF files with identical (CROM/POS/ID/REF/ALT) variants.
/// Pipe in a list of vcf.gz files; umcompressed VCF will be written to STDOUT.
struct Args {
    /// Check for identical CHROM, POS, ID, REF, ALT in every file; ~5% slower
    #[clap(short, long)]
    check: bool,

    /// Serial not parallel file reading; faster if you have few input files, or few CPU cores
    #[clap(short, long)]
    serial: bool,

    /// Output in bgzip format
    #[clap(short, long)]
    bgzip: bool,
}

fn get_out_header(readers: &Vec<FileReader>) -> VCFHeader {
    let all_samples : Vec<U8Vec> = readers.iter().map(|reader|reader.reader.header().samples().to_owned()).flatten().collect();
    let mut items : Vec<VCFHeaderLine> = readers[0].reader.header().items().to_vec() ;
    items.pop();
    VCFHeader::new(items, all_samples)
}

fn read_one_line_from_every_file(readers: &mut Vec<FileReader>, row: usize, serial: bool, check_metadata: bool) -> bool {
    let total_read : usize = if serial {
        readers.iter_mut().map(|reader|reader.load_next()).sum()
    } else {
        readers.par_iter_mut().map(|reader|reader.load_next()).sum()
    } ;
    if total_read != readers.len() {
        if total_read > 0 { // Some VCF files were read but some were not, that's a problem
            panic!("Could only read from {} of {} files on data row {}",total_read,readers.len(),row);
        }
        // All VCF files have reached the end simultaneously, we're done!
        return false ;
    }
    if check_metadata {
        let has_problem = if serial {
            readers.iter().skip(1).any(|x|{!x.check_meta(&readers[0])})
        } else {
            readers.par_iter().skip(1).any(|x|{!x.check_meta(&readers[0])})
        } ;
        if has_problem {
            panic!("Row {} has a problem with CROM/POS/ID/REF/ALT.",row);
        }
    }
    true
}

fn join_vcf_records(readers: &mut Vec<FileReader>) -> VCFRecord {
    let mut joined_vcf_record = readers[0].vcf_record.clone();
    readers.iter_mut().skip(1).for_each(|record|{
        joined_vcf_record.genotype.append(&mut record.vcf_record.genotype);
    });
    joined_vcf_record
}

fn main() {
    let args = Args::parse();

    // Prepare output
    let stdout = stdout() ;
    let stdout_lock = stdout.lock();
    let buffered_writer = BufWriter::new(stdout_lock) ;
    let buffered_writer : Box<dyn Write> = if args.bgzip {
        Box::new(BGZFWriter::new(buffered_writer, flate2::Compression::default()))
    } else {
        Box::new(buffered_writer)
    } ;

    // Read files and write output
    let mut readers: Vec<FileReader> = stdin().lock().lines().map(|line|{FileReader::new(&line.unwrap())}).collect();
    let header = get_out_header(&readers) ;
    let mut vcf_writer = VCFWriter::new(buffered_writer,&header).unwrap();
    let mut row: usize = 1 ;
    while read_one_line_from_every_file(&mut readers, row, args.serial , args.check ) {
        vcf_writer.write_record(&join_vcf_records(&mut readers)).unwrap();
        row += 1 ;
    }
}
