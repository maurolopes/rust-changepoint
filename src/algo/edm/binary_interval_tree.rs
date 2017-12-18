use num::{One, Num, Signed};
use std::ops::Sub;

pub trait TreeNum: PartialOrd + Num + One + Clone + Signed + Sub + From<f64> + Copy {}

#[derive(PartialEq)]
pub struct BinaryIntervalTree<T: TreeNum> {
    pub counter: usize,
    lower: T,
    upper: T,
    children: Option<Box<(BinaryIntervalTree<T>, BinaryIntervalTree<T>)>>
}

impl <T> BinaryIntervalTree<T>
    where T: TreeNum,
{
    pub fn new(depth: usize) -> BinaryIntervalTree<T> {
        BinaryIntervalTree::new3(depth, T::from(0.0), T::from(1.0))
    }

    fn new3(depth: usize, lower: T, upper: T) -> BinaryIntervalTree<T>
    {
        if depth == 0 {
            BinaryIntervalTree{ counter: 0,
                                lower: lower,
                                upper: upper,
                                children: None }
        }
        else {
            let mid = (lower + upper) / T::from(2.0);
            let left = BinaryIntervalTree::new3(depth - 1, lower, mid);
            let right = BinaryIntervalTree::new3(depth - 1, mid, upper);
            BinaryIntervalTree { counter: 0,
                                 lower: lower,
                                 upper: upper,
                                 children: Some(Box::new((left,right))) }
        }
    }

    pub fn insert(&mut self, value: T)
        where T: TreeNum
    {
        if value >= self.lower && value <= self.upper {
            self.counter += 1;
            if let Some(ref mut c) = self.children {
                c.0.insert(value);
                c.1.insert(value);
            };
        }
    }

    pub fn approximate_median(&self) -> T {
        let mut K = ((self.counter as f64) / 2.0f64).ceil() as usize;
        let &mut node = self;

        while let Some(ref c) = node.children {
            let a = c.0.counter;
            if a > K { node = c.0; }
            else if a < K { K -= a; node = c.1; }
            else {
                let b = c.1.counter;
                let x = (c.0.lower + c.0.upper) / 2.0;
                let y = (c.1.lower + c.1.upper) / 2.0;
                return (a*x + b*y) / (a+b);
            }

        }
        (node.lower + node.upper) / 2.0
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::ops::RangeFull;

    #[test]
    fn approximate_median() {
        let mut source_numbers = vec![0.09, 0.42, 0.99, 0.3];
        let mut tree: BinaryIntervalTree<f64> = BinaryIntervalTree::new();
        for source_number in source_numbers.drain(RangeFull) {
            tree.insert(source_number.into());
        }
        assert_eq!(tree.approximate_median(), 0.375);
    }
}
