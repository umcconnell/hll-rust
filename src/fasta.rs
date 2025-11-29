use std::collections::VecDeque;
use std::io::{self, BufRead};

/// A simple FASTA reader that reads records one by one.
///
/// It provides methods to iterate over kmers and canonical kmers of the current record.
pub struct FastaReader<R: BufRead> {
    reader: R,
    line: String,
    finished: bool,
    pub id: Option<Vec<u8>>,
}

impl<R: BufRead> FastaReader<R> {
    /// Creates a new `FastaReader` from a type implementing `BufRead`.
    pub fn new(reader: R) -> Self {
        FastaReader {
            reader,
            line: String::new(),
            finished: false,
            id: None,
        }
    }

    /// Advances the reader to the next record.
    ///
    /// Returns `Ok(true)` if a record was found, `Ok(false)` if EOF was reached.
    /// The record ID is stored in `self.id`.
    pub fn next_record(&mut self) -> io::Result<bool> {
        if self.finished {
            return Ok(false);
        }

        if self.line.is_empty() {
            self.line.clear();
            if self.reader.read_line(&mut self.line)? == 0 {
                self.finished = true;
                return Ok(false);
            }
        }

        if !self.line.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Expected '>' at the start of a fasta record.",
            ));
        }

        self.id = Some(
            self.line
                .trim_start_matches('>')
                .trim_end()
                .as_bytes()
                .to_vec(),
        );
        self.line.clear();

        Ok(true)
    }

    /// Returns an iterator over the kmers of the current record.
    pub fn kmers<'a>(&'a mut self, k: usize) -> KmerStream<'a, R> {
        KmerStream::new(self, k)
    }

    /// Returns an iterator over the canonical kmers of the current record.
    ///
    /// A canonical kmer is the lexicographically smaller of the kmer and its reverse complement.
    pub fn canonical_kmers<'a>(&'a mut self, k: usize) -> CanonicalKmerStream<KmerStream<'a, R>> {
        CanonicalKmerStream::new(self.kmers(k))
    }

    /// Reads the full sequence of the current record.
    ///
    /// This consumes the rest of the current record.
    pub fn read_sequence(&mut self) -> io::Result<Vec<u8>> {
        let mut sequence = Vec::new();
        loop {
            self.line.clear();
            let bytes_read = self.reader.read_line(&mut self.line)?;
            if bytes_read == 0 {
                self.finished = true;
                break;
            }
            if self.line.starts_with('>') {
                break;
            }
            sequence.extend_from_slice(self.line.trim().as_bytes());
        }
        Ok(sequence)
    }
}

/// An iterator over the kmers of a FASTA record.
pub struct KmerStream<'a, R: BufRead> {
    reader: &'a mut FastaReader<R>,
    k: usize,
    buffer: VecDeque<u8>,
    stream_finished: bool,
}

impl<'a, R: BufRead> KmerStream<'a, R> {
    fn new(reader: &'a mut FastaReader<R>, k: usize) -> Self {
        KmerStream {
            reader,
            k,
            buffer: VecDeque::with_capacity(k * 2),
            stream_finished: false,
        }
    }

    fn fill_buffer(&mut self) -> io::Result<()> {
        while self.buffer.len() < self.k && !self.stream_finished {
            self.reader.line.clear();
            let bytes_read = self.reader.reader.read_line(&mut self.reader.line)?;

            if bytes_read == 0 || self.reader.line.starts_with('>') {
                self.stream_finished = true;
                if bytes_read == 0 {
                    self.reader.finished = true;
                }
                break;
            }

            self.buffer.extend(self.reader.line.trim().as_bytes());
        }
        Ok(())
    }
}

impl<'a, R: BufRead> Drop for KmerStream<'a, R> {
    fn drop(&mut self) {
        if self.stream_finished {
            return;
        }

        // Consume the rest of the lines of the current sequence until the next record or EOF
        loop {
            self.reader.line.clear();
            if let Ok(bytes_read) = self.reader.reader.read_line(&mut self.reader.line) {
                if bytes_read == 0 {
                    self.reader.finished = true;
                    break;
                }
                if self.reader.line.starts_with('>') {
                    break;
                }
            } else {
                // On an IO error, we can't do much but stop.
                self.reader.finished = true;
                break;
            }
        }
    }
}

impl<'a, R: BufRead> Iterator for KmerStream<'a, R> {
    type Item = io::Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.stream_finished && self.buffer.len() < self.k {
            return None;
        }

        if let Err(e) = self.fill_buffer() {
            return Some(Err(e));
        }

        if self.buffer.len() < self.k {
            return None;
        }

        let kmer: Vec<u8> = self.buffer.iter().take(self.k).cloned().collect();
        self.buffer.pop_front();

        Some(Ok(kmer))
    }
}

/// An iterator over the canonical kmers of a FASTA record.
///
/// Wraps another iterator yielding kmers and converts them to canonical form.
pub struct CanonicalKmerStream<I> {
    iter: I,
}

impl<I> CanonicalKmerStream<I> {
    pub fn new(iter: I) -> Self {
        CanonicalKmerStream { iter }
    }
}

impl<I> Iterator for CanonicalKmerStream<I>
where
    I: Iterator<Item = io::Result<Vec<u8>>>,
{
    type Item = io::Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.next() {
            Some(Ok(kmer)) => Some(Ok(get_canonical(&kmer))),
            other => other,
        }
    }
}

pub fn get_canonical(kmer: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(kmer);
    if kmer <= &rc { kmer.to_vec() } else { rc }
}

pub fn get_canonical_into<'a>(kmer: &'a [u8], buffer: &'a mut [u8]) -> &'a [u8] {
    reverse_complement_into(kmer, buffer);
    if kmer <= buffer { kmer } else { buffer }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            x => x,
        })
        .collect()
}

pub fn reverse_complement_into(seq: &[u8], out: &mut [u8]) {
    assert_eq!(seq.len(), out.len());
    for (i, &base) in seq.iter().enumerate() {
        out[out.len() - 1 - i] = match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            x => x,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_canonical_kmers() {
        let data = b">seq1\nATCG\n";
        let cursor = Cursor::new(data);
        let mut reader = FastaReader::new(cursor);
        reader.next_record().unwrap();

        let kmers: Vec<Vec<u8>> = reader.canonical_kmers(3).map(|r| r.unwrap()).collect();

        // Expected:
        // ATC -> rev_comp: GAT. min: ATC
        // TCG -> rev_comp: CGA. min: CGA

        assert_eq!(kmers, vec![b"ATC".to_vec(), b"CGA".to_vec()]);
    }

    #[test]
    fn test_canonical_kmers_lowercase_and_n() {
        let data = b">seq1\natcn\n";
        let cursor = Cursor::new(data);
        let mut reader = FastaReader::new(cursor);
        reader.next_record().unwrap();

        let kmers: Vec<Vec<u8>> = reader.canonical_kmers(3).map(|r| r.unwrap()).collect();

        // atc -> rev_comp: gat. min: atc
        // tcn -> rev_comp: nga. min: nga (since n < t)

        assert_eq!(kmers, vec![b"atc".to_vec(), b"nga".to_vec()]);
    }

    #[test]
    fn test_canonical_kmers_palindromes() {
        let data = b">seq1\nGCGC\n";
        let cursor = Cursor::new(data);
        let mut reader = FastaReader::new(cursor);
        reader.next_record().unwrap();

        // k=4: GCGC -> rev_comp: GCGC. min: GCGC
        let kmers: Vec<Vec<u8>> = reader.canonical_kmers(4).map(|r| r.unwrap()).collect();
        assert_eq!(kmers, vec![b"GCGC".to_vec()]);

        // k=2: GC, CG, GC
        // GC -> GC
        // CG -> CG
        // GC -> GC
        let data = b">seq1\nGCGC\n";
        let cursor = Cursor::new(data);
        let mut reader = FastaReader::new(cursor);
        reader.next_record().unwrap();
        let kmers: Vec<Vec<u8>> = reader.canonical_kmers(2).map(|r| r.unwrap()).collect();
        assert_eq!(kmers, vec![b"GC".to_vec(), b"CG".to_vec(), b"GC".to_vec()]);
    }

    #[test]
    fn test_canonical_kmers_multiple_records() {
        let data = b">seq1\nAAA\n>seq2\nTTT\n";
        let cursor = Cursor::new(data);
        let mut reader = FastaReader::new(cursor);

        reader.next_record().unwrap();
        let kmers1: Vec<Vec<u8>> = reader.canonical_kmers(3).map(|r| r.unwrap()).collect();
        // AAA -> TTT. min: AAA
        assert_eq!(kmers1, vec![b"AAA".to_vec()]);

        reader.next_record().unwrap();
        let kmers2: Vec<Vec<u8>> = reader.canonical_kmers(3).map(|r| r.unwrap()).collect();
        // TTT -> AAA. min: AAA
        assert_eq!(kmers2, vec![b"AAA".to_vec()]);
    }

    #[test]
    fn test_short_sequence() {
        let data = b">seq1\nAT\n";
        let cursor = Cursor::new(data);
        let mut reader = FastaReader::new(cursor);
        reader.next_record().unwrap();

        let kmers: Vec<Vec<u8>> = reader.canonical_kmers(3).map(|r| r.unwrap()).collect();
        assert!(kmers.is_empty());
    }
}
