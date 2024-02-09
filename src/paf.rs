use std::io::BufRead;

pub struct PafRecord {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize, 
    pub query_end: usize,
    pub target_name: String,
    pub target_length: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub matches: usize,
    pub mapq: u8,
}

impl PafRecord {

    pub fn parse(line: &str) -> Result<Self, ParseErr> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(ParseErr::NotEnoughFields);
        }
        
        let query_name = fields[0].to_string();
        let query_length = fields[1].parse::<usize>()?;
        let query_start = fields[2].parse::<usize>()?;
        let query_end = fields[3].parse::<usize>()?;
        let target_name = fields[5].to_string();
        let target_length = fields[6].parse::<usize>()?;
        let target_start = fields[7].parse::<usize>()?;
        let target_end = fields[8].parse::<usize>()?;
        let matches = fields[9].parse::<usize>()?;
        let mapq = fields[11].parse::<u8>()?;

        Ok(Self {
            query_name,
            query_length,
            query_start,
            query_end,
            target_name,
            target_length,
            target_start,
            target_end, 
            matches,
            mapq
        })
    }

}

#[derive(Debug)]
pub enum ParseErr {
    NotEnoughFields,
    IoError(std::io::Error),
    InvalidField(std::num::ParseIntError),
}

impl From<std::num::ParseIntError> for ParseErr {
    fn from(err: std::num::ParseIntError) -> Self {
        ParseErr::InvalidField(err)
    }
}

pub fn parse_paf<R: BufRead>(reader: R) -> Result<Vec<PafRecord>, ParseErr> {
    let mut records = Vec::new();

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| ParseErr::IoError(e))?; // Convert io::Error to ParseErr
        let record = PafRecord::parse(&line)?;
        records.push(record);
    }

    Ok(records)
}


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse() {
        let line = "q1\t100\t10\t20\t+\tt1\t200\t30\t40\t10\t20\t255";
        let record = PafRecord::parse(line).unwrap();
        assert_eq!(record.query_name, "q1");
        assert_eq!(record.query_length, 100);
        assert_eq!(record.query_start, 10);
        assert_eq!(record.query_end, 20);
        assert_eq!(record.target_name, "t1");
        assert_eq!(record.target_length, 200);
        assert_eq!(record.target_start, 30);
        assert_eq!(record.target_end, 40);
        assert_eq!(record.matches, 10);
        assert_eq!(record.mapq, 255);

        let empty = "";
        assert!(PafRecord::parse(empty).is_err());
    }
}
