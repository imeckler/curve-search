use rug::Integer;
use rug::integer::IsPrime;
use rug::Assign;
use std::thread;
use primal::Sieve;
use quadratic;
use integer_sqrt::IntegerSquareRoot;
use rayon::prelude::*;
use rayon::range::*;

/*
squarefree check:
check all
12.5 * 10^6 primes  < 10^16 ?
*/

/*
want to solve
X^2 - 3DY^2 = -8

write
D' = 3D

X^2 - D' Y^2 = -8

Fix 0 < D' < 3 * 10^16, D' / 3 squarefree, D' = 9 mod 24, -2 is square mod D'.

fn is_square_free(s : &Sieve, n : u64) -> bool {
    s.factor(n as usize).unwrap().iter().all(|(_, n)| *n < 2)
}

*/

fn is_suitable(sieve : &Sieve, d : u64) -> bool {
    let d3_factorization = sieve.factor((d/3) as usize).unwrap();
    let d3_square_free = d3_factorization.iter().all(|(_, n)| *n < 2);

    if d3_square_free {
        let a = (d - 2) as isize;
        if quadratic::jacobi(a, 3).unwrap() == 1 {
            let a_is_square_mod_d : bool = d3_factorization.iter().all(|(p, _)| {
                *p == 3 || quadratic::jacobi(a, *p as isize).unwrap() == 1
            });

            return a_is_square_mod_d
        }
    }
    return false
}

fn pre_pell(dint : u64) -> Vec<(Integer, Integer)> {
    let floor_sqrt_d = Integer::from(dint.integer_sqrt());

    let mut b_i_sub_1 = Integer::from(0);
    let mut g_i_sub_1 = Integer::from(1);

    let mut p_i = Integer::from(0);
    let mut q_i = Integer::from(1);
    let mut a_i = Integer::from(&floor_sqrt_d);

    let mut b_i = Integer::from(1);
    let mut g_i = Integer::from(&a_i);

    let mut p_i_squared = Integer::from(&p_i * &p_i);

    let mut i = 0;

    let d = Integer::from(dint);

    let mut gbs : Vec<(Integer, Integer)> = Vec::new();
    gbs.push((Integer::from(&g_i), Integer::from(&b_i)));

    loop {
        i += 1;
        p_i.assign(Integer::from(&a_i * &q_i) - &p_i);
        p_i_squared.assign(p_i.square_ref());

        q_i.assign( Integer::from(&d - &p_i_squared).div_exact(&q_i) );
        /*
        let (qquot, qrem) = <(Integer, Integer)>::from(Integer::from(&d - &p_i_squared).div_rem_ref(&q_i));
        assert!(qrem == 0);
        q_i.assign(qquot); */

        let (quot, _) = 
            <(Integer, Integer)>::from(
            Integer::from(&p_i + &floor_sqrt_d).div_rem_floor_ref(&q_i));

        a_i.assign(quot);

        let b_new = Integer::from(&a_i * &b_i + &b_i_sub_1);
        b_i_sub_1.assign(&b_i);
        b_i.assign(b_new);

        let g_new = Integer::from(&a_i * &g_i + &g_i_sub_1);
        g_i_sub_1.assign(&g_i);
        g_i.assign(g_new);

        gbs.push((Integer::from(&g_i), Integer::from(&b_i)));

        if q_i == 1 && (i % 2) == 0 {
            break;
        }
    }

    return gbs
}

fn pell1(d : u64, m : i64) -> (Option<(Integer, Integer)>, Vec<(Integer, Integer)>) {
    let gbs = pre_pell(d);
    let d = Integer::from(d);

    /* All solutions to X^2 - DY^2 = 1 are in the same class due to
    * proposition 2 of the paper.
    * Here we find the minimal element of the class.
    * */

    /* A minimal solution (x, y) the is one for which x is non-negative and
     * y is as small as possible */
    let zero = Integer::from(0);
    /* Could squash this into the pre_pell loop */
    let result1 : Option<(&Integer, &Integer)> = gbs.iter().fold(None, |acc, (g, b)| {
        if Integer::from(g.square_ref()) - &d & Integer::from(b.square_ref()) == 1 && g.ge(&zero) {
            match acc {
                None => Some ((g, b)),
                Some ((_g_acc, b_acc)) =>
                    if b.lt(b_acc) {
                        Some ((g, b))
                    } else {
                        acc
                    }
            }
        } else {
            acc
        }
    });
    let result1 = result1.map(|(g, b)| (Integer::from(g), Integer::from(b)));

    let results_m : Vec<(Integer, Integer)> = gbs.iter().filter_map(|(g, b)| {
        let x = Integer::from(g.square_ref()) - &d * Integer::from(b.square_ref());
        let (quot, rem) = Integer::from(m).div_rem_floor(x);

        if rem == 0 && quot > 0 && quot.is_perfect_square() {
            let f = quot.sqrt();
            Some ((Integer::from(&f * g), Integer::from(&f * b)))
        } else {
            None
        }
    }).collect();

    return (result1, results_m)
}

/*
fn pell2(d : u64, m : i64) {
} */

fn pow2(n : u64) -> Integer {
    let mut res = Integer::from(1);
    for _ in 0..n {
        res.assign(Integer::from(&res + &res));
    }
    return res
}

fn smooth_part(n : &Integer) -> Integer {
    let small_primes : [u32; 7] = [ 2, 3, 5, 7, 11, 13, 17 ];
    let mut x = Integer::from(n);
    let mut res = Integer::from(1);
    for p in small_primes.iter() {
        while x.mod_u(*p) == 0 {
            x.div_exact_u_mut(*p);
            res *= *p;
        }
    }
    res
}

fn is_potentially_prime(n : &Integer) -> bool {
    /*
    let small_primes : [u32; 10] = [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 ];
    small_primes.iter().all(|p| n.mod_u(*p) != 0) */

    match n.is_probably_prime(0) {
        IsPrime::No => false,
        _ => true
    }
}

fn process_pell_solution(
    max_x : &Integer,
    min_x : &Integer,
    d_prime : u64,
    (u, v) : (&Integer, &Integer),
    (x0, y0) : (Integer, Integer)) {

    let v_d_prime = Integer::from(d_prime) * v;
    let mut x = Integer::from(x0);
    let mut y = Integer::from(y0);
    let mod6 = x.mod_u(6);
    let sign = if mod6 == 5 { 1 } else { -1 };

    if mod6 == 1 || mod6 == 5 {
        while x.as_abs().lt(max_x) {
            // l = (x -+ 1)/6
            let mut l = Integer::from(&x + sign);
            l.div_exact_u_mut(6);

            let s = smooth_part(&l).to_u64();
            let other = Integer::from(2* &l) + 1;
            let s_n = smooth_part(&Integer::from(other)).to_u64();

            // TODO: Check if l is in the right range
            if  l.as_abs().ge(min_x) &&
                match (s, s_n) {
                    (None, _) | (_, None) => true,
                    (Some (s), Some (s_n)) =>  
                        4 * s > 100000 / s && 2 * s_n > 10000 / s
                } {
                let l2 = Integer::from(l.square_ref());

                // q = 4l^2 + 1
                let mut q = Integer::from(4 * &l2);
                q += 1;

                // n = 4l^2 -+ 2l + 1 = 2l (2l + 1) + 1
                let mut n = Integer::from(4 * &l2);
                n += 2 * &l;
                n += 1;

                if is_potentially_prime(&q) && is_potentially_prime(&n) {
                    println!("q, n = {} , {} , {}", q, n, d_prime)
                }
                /*
                match (q.is_probably_prime(1), n.is_probably_prime(1)) {
                    (IsPrime::No, _) => (),
                    (_, IsPrime::No) => (),
                    _ => {
                        println!("q, n = {} , {} , {}", q, n, s)
                    }
                } */
            }

            let prev_x = Integer::from(&x);
            x *= u;
            x += &v_d_prime * &y;
            y *= u;
            y += prev_x * v;
        }
    }
}



fn main() {
    /*
    static NTHREADS:i64 = 8;
    let n = NTHREADS * 10000000000;
    let chunk_size = n / NTHREADS;
    */

    let max_determinant_log10 = 6 + 3;
    let upper_bound : u64 = 3 * u64::pow(10, max_determinant_log10);
    let square_root_upper_bound : u64 = 1 + upper_bound.integer_sqrt();

    let sieve = primal::Sieve::new(square_root_upper_bound as usize);

    let mut results = 0;
    let mut i = 0;
    let max = upper_bound / 24;

    let max_bits = 1500;
    let min_bits = 400;

    let max_x = pow2((max_bits / 2) + 1);
    let min_x = pow2(min_bits / 2);

    // (9..upper_bound).step_by(24).for_each(|d_prime| {

    let chunk_size = 10000;
    (0..(upper_bound / 24 / chunk_size)).into_par_iter().for_each(|chunk| {
        for idx in 0..chunk_size {
            let d_pre = chunk * chunk_size + idx;
            let d_prime : u64 = 9 + 24  * d_pre;

            if is_suitable(&sieve, d_prime) {
                if d_prime > 64 {
                    let (result1, results8) = pell1(d_prime, -8);
                    // TODO: Filter out duplicate classes from results8
                    if results8.len() > 0 && result1.is_some() {
                        let (u, v) = result1.unwrap();
                        let (x0, y0) = &results8[0];
                        process_pell_solution(
                            &max_x, &min_x, d_prime, (&u, &v),
                            (Integer::from(x0), Integer::from(y0)));

                        let x = Integer::from(x0 * &u) - Integer::from(y0 * &v) * Integer::from(d_prime);
                        let y = Integer::from(y0 * &u) - Integer::from(x0 * &v);
                        process_pell_solution(
                            &max_x, &min_x, d_prime,
                            (&u, &Integer::from(-v)),
                            (x, y));
                    }
                }
            }
        }
    });

    println!("RESULTS {}, {}", results, (results as f64) / (upper_bound / 24) as f64);
}

// where do we go from here
// 1. speed up this program
// 2a. Figure out the BLS equations for lower embedding degrees than 12
// 2b. solve the BLS quartic to go from MNT to BLS?
