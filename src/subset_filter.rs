use std::collections::HashSet;
use std::io;

/// Filter for matching sequence names against a subset list
#[derive(Default)]
pub struct SubsetFilter {
    exact: HashSet<String>,
    normalized: HashSet<String>,
    sample_ids: HashSet<String>,
    sample_haps: HashSet<(String, String)>,
}

impl SubsetFilter {
    /// Returns the number of entries in the filter
    pub fn entry_count(&self) -> usize {
        self.exact.len()
    }

    /// Check if a sequence name matches the filter
    pub fn matches(&self, seq_name: &str) -> bool {
        if self.exact.contains(seq_name) {
            return true;
        }

        let no_coords = seq_name.split(':').next().unwrap_or(seq_name);
        if seq_name != no_coords && self.exact.contains(no_coords) {
            return true;
        }
        if self.normalized.contains(no_coords) {
            return true;
        }

        if self.matches_sample_keys(no_coords) {
            return true;
        }

        // As a fallback, try again on the raw sequence in case coordinates were not present
        self.matches_sample_keys(seq_name)
    }

    fn matches_sample_keys(&self, seq_name: &str) -> bool {
        if let Some((sample, hap)) = extract_sample_and_hap(seq_name) {
            if let Some(hap) = hap {
                if self.sample_haps.contains(&(sample.clone(), hap.clone())) {
                    return true;
                }
            }

            if self.sample_ids.contains(&sample) {
                return true;
            }
        }

        false
    }
}

/// Load a subset filter from a file
pub fn load_subset_filter(path: &str) -> io::Result<SubsetFilter> {
    let contents = std::fs::read_to_string(path).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Failed to read subset sequence list '{path}': {e}"),
        )
    })?;

    let filter = parse_subset_filter(&contents);
    if filter.entry_count() == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Subset sequence list '{path}' did not contain any sequence names"),
        ));
    }

    Ok(filter)
}

fn parse_subset_filter(contents: &str) -> SubsetFilter {
    let mut filter = SubsetFilter::default();

    for line in contents.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        filter.exact.insert(trimmed.to_string());

        let no_coords = trimmed.split(':').next().unwrap_or(trimmed);
        filter.normalized.insert(no_coords.to_string());

        if let Some((sample, hap)) = extract_sample_and_hap(no_coords) {
            if let Some(hap) = hap {
                filter.sample_haps.insert((sample, hap));
            } else {
                filter.sample_ids.insert(sample);
            }
        }
    }

    filter
}

fn extract_sample_and_hap(name: &str) -> Option<(String, Option<String>)> {
    if let Some(idx) = name.find("_hap") {
        let sample = name[..idx].to_string();
        let hap_digits: String = name[idx + 4..]
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        let hap = if hap_digits.is_empty() {
            None
        } else {
            Some(hap_digits)
        };
        return Some((sample, hap));
    }

    if let Some((sample, rest)) = name.split_once('#') {
        let hap_fragment = rest.split('#').next().unwrap_or(rest);
        let hap_digits: String = hap_fragment
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        let hap = if hap_digits.is_empty() {
            None
        } else {
            Some(hap_digits)
        };
        return Some((sample.to_string(), hap));
    }

    if !name.contains(':') && !name.trim().is_empty() {
        return Some((name.to_string(), None));
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subset_filter_matches_variants() {
        let contents = "# comment\nchr1\nchr2\n\nchr1\t\n  chr3  \nHG00097_hap1_hprc_r2_v1.0.1\nHG00098#2#chr5\n";
        let filter = parse_subset_filter(contents);

        // Basic names and coordinate variants
        assert!(filter.matches("chr1"));
        assert!(filter.matches("chr1:10-20"));
        assert!(filter.matches("chr3"));

        // Sample + hap conversions
        assert!(filter.matches("HG00097#1#chr7"));
        assert!(filter.matches("HG00097#1"));

        // Exact hashed names
        assert!(filter.matches("HG00098#2#chr5"));
        assert!(!filter.matches("HG00098#1#chr5"));
    }
}
