#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use vcf::*;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{stdin, stdout, BufRead, BufReader, BufWriter};
use rayon::prelude::*;

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
            reader:reader, 
            _filename:filename.to_string(), 
            vcf_record:vcf_record
        }
    }

    pub fn load_next(&mut self) -> Result<bool,VCFError> {
        Ok(self.reader.next_record(&mut self.vcf_record)?)
    }

    pub fn check_meta(&self, origin: &Self) -> bool {
        (origin.vcf_record.chromosome == self.vcf_record.chromosome) &&
        (origin.vcf_record.position == self.vcf_record.position) &&
        (origin.vcf_record.reference == self.vcf_record.reference) &&
        (origin.vcf_record.alternative == self.vcf_record.alternative)
    }
}

fn main() {
    let mut readers : Vec<FileReader> = stdin().lock().lines().map(|line|{FileReader::new(&line.unwrap())}).collect();
    let all_samples : Vec<U8Vec> = readers.iter().map(|reader|reader.reader.header().samples().to_owned()).flatten().collect();
    let mut items : Vec<VCFHeaderLine> = readers[0].reader.header().items().to_vec() ;
    items.pop();
    let out_header = VCFHeader::new(items, all_samples);
    let s = stdout() ;
    let writer = BufWriter::new(s.lock());
    let mut out = VCFWriter::new(writer,&out_header).unwrap();

    let mut row = 0 ;
    loop {
        row += 1 ;
        let total_read : usize = readers
        .par_iter_mut()
        .map(|reader|{
            match reader.load_next() {
                Ok(true) => 1 ,
                _ => 0 ,
            }
        })
        .sum();
        if total_read != readers.len() {
            if total_read > 0 {
                println!("Could only read from {} of {} files on data row {}",total_read,readers.len(),row);
            }
            break
        }
        
        /*
        // Paranoia
        if readers.iter().any(|x|{!x.check_meta(&readers[0])}) {
            println!("Row {} has a problem!",row);
            break ;
        }
        */
        let mut joined_vcf_record = readers[0].vcf_record.clone();
        readers.iter_mut().skip(1).for_each(|record|{
            joined_vcf_record.genotype.append(&mut record.vcf_record.genotype);
        });
        out.write_record(&joined_vcf_record).unwrap();
    }
}
