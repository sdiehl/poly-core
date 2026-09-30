use num_rational::BigRational;
use polycore::{Echelon, Lead};

fn main() {
    let q = |n: i64| BigRational::from_integer(n.into());
    let species = ["C3H8", "O2", "CO2", "H2O"];
    let mut ech = Echelon::new(Lead::Low);
    ech.insert(&[(0, q(3)), (2, q(-1))]);
    ech.insert(&[(0, q(8)), (3, q(-2))]);
    ech.insert(&[(1, q(2)), (2, q(-2)), (3, q(-1))]);
    let v = &ech.nullspace(4)[0];
    let scale = v[0].1.clone();
    let coef: Vec<String> = v.iter().map(|(_, c)| (c / &scale).to_string()).collect();
    println!(
        "{} {} + {} {} -> {} {} + {} {}",
        coef[0], species[0], coef[1], species[1], coef[2], species[2], coef[3], species[3]
    );
}
