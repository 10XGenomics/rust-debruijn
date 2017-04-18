use Kmer;

pub struct PathCompression<T: Kmer> {
    switch_strands: bool,
    k: T,
}

impl<T: Kmer> PathCompression<T> {
    pub fn go(&self) {
        println!("{}", self.switch_strands);
    }
}