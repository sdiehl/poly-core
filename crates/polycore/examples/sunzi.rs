use polycore::crt;

fn main() {
    let (x, m) = crt::crt(&[(2, 3), (3, 5), (2, 7)]);
    println!("x = {x} mod {m}");
}
