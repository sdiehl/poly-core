use polycore::{Order, Ring};

fn main() -> Result<(), polycore::ParseError> {
    let ring = Ring::new(["x", "y"], Order::Lex);
    let f = ring.parse("x^2*y + x*y^2 + y^2")?;
    let gens = [ring.parse("x*y - 1")?, ring.parse("y^2 - 1")?];
    let (quots, rem) = f.divide(&gens);
    println!("{} =", ring.show(&f));
    println!("  ({}) * (x*y - 1)", ring.show(&quots[0]));
    println!("+ ({}) * (y^2 - 1)", ring.show(&quots[1]));
    println!("+ {}", ring.show(&rem));
    Ok(())
}
