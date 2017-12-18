use algo::best_candidate::BestCandidate;
use algo::changepoint::ChangePointDetector;
use algo::edm::binary_interval_tree::{BinaryIntervalTree,TreeNum};

use errors::*;

fn calculate_stat<T>(tau: usize, kappa: usize, m1: T, m2: T, m3: T) -> T
    where T: TreeNum + From<f64>
{
    let ftau1 = T::from(tau as f64);
    let ftau2 = T::from(tau as f64);
    let fkappa1 = T::from(kappa as f64);
    let fkappa2 = T::from(kappa as f64);
    ftau1 * (fkappa1 - ftau2) / fkappa2 * (m1.clone() + m1 - m2 - m3)
}

fn approx_median<T>(tree: &BinaryIntervalTree<T>) -> T
    where T: TreeNum + From<f64>
{
    T::from(0.0)
}

fn forward_update<T>(z: &[T],
                     delta: usize,
                     t_a: &mut BinaryIntervalTree<T>,
                     t_b: &mut BinaryIntervalTree<T>,
                     t_ab: &mut BinaryIntervalTree<T>,
                     tau: &mut usize,
                     best_stat: &mut T,
                     best_loc: &mut usize)
    where
    T: TreeNum + From<f64>,
{
    let n = z.len();
    *tau += 1;
    // update counts in t_a, t_b and t_ab resulting from new tau value
    for kappa in (*tau+delta)..n {
        t_b.insert((z[kappa]-z[kappa-1]).abs());
        let m1 = approx_median(&t_ab);
        let m2 = approx_median(&t_a);
        let m3 = approx_median(&t_b);
        let stat = calculate_stat(t_a.counter, t_b.counter, m1, m2, m3);
        if stat > *best_stat {
            *best_stat = stat;
            *best_loc = *tau;
        }
    }
}

fn backward_update<T>(z: &[T],
                      delta: usize,
                      t_a: &mut BinaryIntervalTree<T>,
                      t_b: &mut BinaryIntervalTree<T>,
                      t_ab: &mut BinaryIntervalTree<T>,
                      tau: &mut usize,
                      best_stat: &mut T,
                      best_loc: &mut usize)
    where
    T: TreeNum + From<f64>,
{
    let n = z.len();
    *tau += 1;
    // update counts in t_a, t_b and t_ab resulting from new tau value
    let mut kappa = n;
    while kappa >= *tau+delta {
        t_b.insert((z[kappa] - z[kappa-1]).abs());
        let m1 = t_ab.approximate_median();
        let m2 = t_a.approximate_median();
        let m3 = t_b.approximate_median();
        let stat = calculate_stat(t_a.counter, t_b.counter, m1, m2, m3);
        if stat > *best_stat {
            *best_stat = stat;
            *best_loc = *tau;
        }
        kappa -= 1;
    }
}

fn edm<T>(z: &[T], delta: usize, depth: usize) -> BestCandidate<T>
where
    T: TreeNum + From<f64>,
{
    let n = z.len();
    let mut t_a = BinaryIntervalTree::new(depth);
    let mut t_b = BinaryIntervalTree::new(depth);
    let mut t_ab = BinaryIntervalTree::new(depth);

    for i in 1..delta {
        for j in (i+1)..delta {
            t_a.insert((z[i]-z[j]).abs());
            t_b.insert((z[i+delta] - z[j+delta]).abs());
        }
    }
    for i in 1..delta {
        for j in 1..delta {
            t_ab.insert((z[i] - z[j+delta]).abs());
        }
    }
    let m1 = t_ab.approximate_median();
    let m2 = t_a.approximate_median();
    let m3 = t_b.approximate_median();

    let mut best_stat = calculate_stat(t_a.counter, t_b.counter, m1, m2, m3);
    let mut best_loc = delta;

    let mut tau = delta;
    let mut forward_move = 0;

    while tau <= n-delta {
        if forward_move == 1 {
            forward_update(&z,
                           delta,
                           &mut t_a,
                           &mut t_b,
                           &mut t_ab,
                           &mut tau,
                           &mut best_stat,
                           &mut best_loc);
        } else {
            backward_update(&z,
                            delta,
                            &mut t_a,
                            &mut t_b,
                            &mut t_ab,
                            &mut tau,
                            &mut best_stat,
                            &mut best_loc);
        }
        forward_move = 1-forward_move;
    }
    return BestCandidate { statistic: best_stat,
                           location: best_loc }
}



#[derive(Clone, Debug)]
pub struct EDM {
    delta: usize,
    depth: usize
}

impl EDM {
    pub fn new(delta: usize, depth: usize) -> Self {
        EDM { delta: delta, depth: depth }
    }
}

impl<T> ChangePointDetector<T> for EDM
where T: TreeNum,
{
    fn find_candidate(&self, observations: &[T]) -> Result<BestCandidate<T>> {
        if observations.len() < self.delta * 2 {
            Err(
                ErrorKind::NotEnoughValues(observations.len(), self.delta).into(),
            )
        } else {
            Ok(edm(observations, self.delta, self.depth))
        }
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use algo::non_nan::NonNaN;
//     use rand::SeedableRng;
//     use rand::distributions::{IndependentSample, Normal};
//     use mersenne_twister::MersenneTwister;
//     use num::abs;

//     #[test]
//     fn heaps_find_the_median() {
//         let initial_number: NonNaN<f32> = NonNaN::new(1.0).unwrap();
//         let mut heaps: Heaps<NonNaN<f32>> = Heaps::new();
//         heaps.add_to_heaps(initial_number.clone());
//         assert_eq!(heaps.get_median(), initial_number.clone());
//         heaps.add_to_heaps(NonNaN::new(2.0).unwrap());
//         heaps.add_to_heaps(NonNaN::new(3.0).unwrap());
//         heaps.add_to_heaps(NonNaN::new(4.0).unwrap());
//         heaps.add_to_heaps(NonNaN::new(5.0).unwrap());
//         heaps.add_to_heaps(NonNaN::new(6.0).unwrap());
//         heaps.add_to_heaps(NonNaN::new(7.0).unwrap());
//         heaps.add_to_heaps(NonNaN::new(8.0).unwrap());
//         assert_eq!(heaps.get_median(), NonNaN::new(4.0).unwrap());
//     }

//     #[test]
//     fn edm_on_central_tendency() {
//         let before_change_count = 100;
//         let after_change_count = 400;
//         let delta = 10;
//         let tolerance = 50;
//         let mut rng: MersenneTwister = SeedableRng::from_seed(0x1234);
//         let mut input: Vec<NonNaN<f64>> = Vec::new();
//         let before_change_dist = Normal::new(10.0, 5.0);
//         for _ in 0..before_change_count {
//             input.push(
//                 NonNaN::new(before_change_dist.ind_sample(&mut rng)).unwrap(),
//             );
//         }
//         let after_change_dist = Normal::new(30.0, 5.0);
//         for _ in 0..after_change_count {
//             input.push(NonNaN::new(after_change_dist.ind_sample(&mut rng)).unwrap());
//         }
//         let best_candidate = edm(&input, delta);
//         let abs_loc_diff = abs(best_candidate.location as i64 - before_change_count as i64);
//         assert!(abs_loc_diff < tolerance);
//     }
// }
