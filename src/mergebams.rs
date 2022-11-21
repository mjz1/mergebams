
/**

*********create test files*********

cd ~/develop/mergebams
cargo build --release
ml SAMtools
samtools view -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D383_3/outs/per_sample_outs/ITS_D383_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam1.bam
samtools view -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam2.bam
samtools view -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam3.bam
samtools view -h -s 0.0001 /home/sfurlan/scratch/NB1/AML_601_34/outs/per_sample_outs/AML_601_34/count/sample_alignments.bam > ~/develop/mergebams/test/bam4.bam


head /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_barcodes.csv

************************************


*********installation of mergebams*************
* 
* 
* *********************************************

*********usage of mergebams*********

cd ~/develop/mergebams/test
cargo build --release
../target/release/mergebams -i bam1.bam,bam2.bam -l test1_,test2_ -b barcodes1.tsv.gz,barcodes2.tsv.gz -o .
zcat < out_barcodes.tsv.gz
samtools view out_bam.bam | head -n 200
samtools sort out_bam.bam > out_bam.sorted.bam
samtools view out_bam.sorted.bam | head -n 200

************************************


**/


extern crate simple_log;
extern crate clap;
extern crate bam;
extern crate csv;
extern crate flate2;
extern crate differ;

use differ::{Differ, Tag};
use flate2::GzBuilder;
use flate2::Compression;
use bam::RecordWriter;
use bam::record::tags::TagValue;
use clap::{App, load_yaml};
use std::str;
use flate2::{read};
use std::{
    // error::Error,
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, BufRead, Write},
    path::Path,
};


struct Params {
    inputs: String,
    labels: String,
    bcs: String,
    out: String,
    threads: usize,
}


fn main() {
    let params = load_params();
    let params = checkheaders(params);
    let params = addtags(params);
    if params.bcs == "0"{
        return;
    } else {
        let _result = lines_from_file(params);
        return;
    }
}


fn load_params() -> Params {
    let yaml = load_yaml!("params_mergebams.yml");
    let params = App::from_yaml(yaml).get_matches();
    let inputs = params.value_of("inputs").unwrap();
    let labels = params.value_of("labels").unwrap();
    let threads: usize = params.value_of("threads").unwrap_or("1").parse().unwrap();
    let bcs = params.value_of("bcs").unwrap_or("0");
    let out = params.value_of("out").unwrap();
    Params{
        inputs: inputs.to_string(),
        labels: labels.to_string(),
        bcs: bcs.to_string(),
        out: out.to_string(),
        threads: threads,
    }
}


fn lines_from_file(params: Params) -> io::Result<()>{
    let out_csv = params.out.to_string()+"/out_barcodes.tsv.gz";
    let s1 = params.bcs.to_string();
    let tsvvec = s1.split(",").collect::<Vec<&str>>();
    let labvec1 = params.labels.to_string();
    let labvec = labvec1.split(",").collect::<Vec<&str>>();
    let f = File::create(out_csv)?;
    let mut gz = GzBuilder::new()
                    .filename("/out_barcodes.tsv.gz")
                    .write(f, Compression::default());

    for (pos, filename) in tsvvec.iter().enumerate() {
        let path = Path::new(filename);
        let file = match File::open(&path) {
            Err(_why) => panic!("couldn't open {}", path.display()),
            Ok(file) => file,
        };
        if path.extension() == Some(OsStr::new("gz")){
            let buf = BufReader::new(read::GzDecoder::new(file));
            for line in buf.lines(){
                writeln!(&mut gz, "{}", labvec[pos].to_string() + &line.unwrap())?;
            }
        }else{
            let buf = BufReader::new(file);
            for line in buf.lines(){
                writeln!(&mut gz, "{}", labvec[pos].to_string() + &line.unwrap())?;
            }
        }
    }
    return Ok(());
}

// fn print_lines(c: char, lines: &[&str]) {
//     for line in lines {
//         println!("{} {}", c, line);
//     }
// }

fn checkheaders(params: Params) -> Params{
    let inputs = params.inputs.to_string();
    let bam_vec = inputs.split(",").collect::<Vec<&str>>();
    let mut header_vec = Vec::new();
    for inbam in bam_vec.iter() {
        let hreader = bam::BamReader::from_path(inbam, 0).unwrap();
        let header = hreader.header().clone();
        header_vec.push(header);
    }
    // let result = checkheaders(header_vec);
    // for header in header_vec.iter() {
    //     println!("{:?}", header.reference_names())
    // }
    let a = header_vec[0].reference_names();
    let b = header_vec[1].reference_names();
    let differ = Differ::new(&a, &b);
    for span in differ.spans() {
        match span.tag {
            Tag::Equal => println!("all headers are equal!"),
            Tag::Insert => println!("{:?}", span),
            Tag::Delete => println!("{:?}", span),
            Tag::Replace => println!("{:?}", span),
        }
    }
    return params;
}

// fn checkheaders(header_vec: Vec<bam::Header>) -> io::Result<()>{

//     return Ok(());
// }

fn addtags(params: Params) -> Params{
    let out_bam = params.out.to_string()+"/out_bam.bam";
    let inputs = params.inputs.to_string();
    let bam_vec = inputs.split(",").collect::<Vec<&str>>();
    let labels = params.labels.to_string();
    let lab_vec = labels.split(",").collect::<Vec<&str>>();
    let (read_threads, write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let hreader = bam::BamReader::from_path(bam_vec[0], read_threads).unwrap();
    let mut writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(write_threads)
        .from_path(out_bam, hreader.header().clone()).unwrap();
    for (pos, inbam) in bam_vec.iter().enumerate() {
        eprintln!("Opening: {:?}", inbam);
        let reader = bam::BamReader::from_path(inbam.to_string(), read_threads).unwrap();
        for record in reader {
            let mut newrecord = record.as_ref().unwrap().clone();
            match record.unwrap().tags().get(b"CB") {
                Some(TagValue::String(array_view, _)) => {
                    let labstr = lab_vec[pos];
                    let preftag = labstr.as_bytes().to_vec();
                    let oldtag = array_view.to_vec();
                    let new_cb = &[preftag, oldtag].concat();
                    newrecord.tags_mut().remove(b"CB");
                    newrecord.tags_mut().push_string(b"CB", &new_cb);
                    writer.write(&newrecord).unwrap();
                },
                Some(TagValue::Char(value)) => println!("Char = {}", value),
                _ => panic!("'CB' does not appear to have string type"),
            }
        }
    }
    return params;
}
    

