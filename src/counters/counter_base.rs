pub trait Counter {
    fn new(size: usize) -> Self;
    fn add(&mut self, item: &[u8]);
    fn estimate(&self) -> f64;
}
