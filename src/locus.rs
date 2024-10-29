use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
};

#[derive(Debug)]
pub struct Locus {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
}

pub fn load_loci(path: PathBuf) -> Result<Vec<Locus>, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let reader = BufReader::new(file);
    let mut loci = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let rec: Vec<&str> = line.split_whitespace().collect();
        if rec.len() < 4 {
            return Err(format!("Bad input line {line}"));
        }
        let (chrom, start, end, name) = (rec[0].to_string(), rec[1], rec[2], rec[3]);
        let start = start
            .parse::<i64>()
            .map_err(|_| format!("Bad input line {line}"))?;
        let end = end
            .parse::<i64>()
            .map_err(|_| format!("Bad input line {line}"))?;
        let name = name.to_string();
        loci.push(Locus {
            chrom,
            start,
            end,
            name,
        });
    }

    Ok(loci)
}
