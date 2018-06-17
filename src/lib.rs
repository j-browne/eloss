#![feature(proc_macro, specialization)]

#[macro_use]
extern crate lazy_static;
extern crate pyo3;
use pyo3::prelude::*;
use pyo3::py::modinit;
use interpolation::interpolate;
use std::collections::HashMap;

mod interpolation;

lazy_static! {
    static ref STOPPING_POWERS: HashMap<String, (Vec<f64>, Vec<f64>)> = {
        let mut map = HashMap::new();
        for (proj, targ, input) in [
            ("34S", "Butane", include_str!("data/34S_butane.txt")),
            ("34Cl", "Butane", include_str!("data/34Cl_butane.txt")),
            ("34Ar", "Butane", include_str!("data/34Ar_butane.txt")),
            ("37Cl", "Butane", include_str!("data/37Cl_butane.txt")),
            ("37Ar", "Butane", include_str!("data/37Ar_butane.txt")),
            ("37K", "Butane", include_str!("data/37K_butane.txt")),
            ("34S", "Mylar", include_str!("data/34S_mylar.txt")),
            ("34Cl", "Mylar", include_str!("data/34Cl_mylar.txt")),
            ("34Ar", "Mylar", include_str!("data/34Ar_mylar.txt")),
            ("37Cl", "Mylar", include_str!("data/37Cl_mylar.txt")),
            ("37Ar", "Mylar", include_str!("data/37Ar_mylar.txt")),
            ("37K", "Mylar", include_str!("data/37K_mylar.txt")),
            ("34S", "He", include_str!("data/34S_he.txt")),
            ("34Cl", "He", include_str!("data/34Cl_he.txt")),
            ("34Ar", "He", include_str!("data/34Ar_he.txt")),
            ("37Cl", "He", include_str!("data/37Cl_he.txt")),
            ("37Ar", "He", include_str!("data/37Ar_he.txt")),
            ("37K", "He", include_str!("data/37K_he.txt")),
        ].into_iter()
        {
            let mut xs = Vec::new();
            let mut ys = Vec::new();
            for line in input.lines().skip(1) {
                let mut line = line.split_whitespace().skip(2);
                xs.push(line.next().unwrap().parse().unwrap());
                ys.push(line.next().unwrap().parse().unwrap());
            }
            map.insert(format!("{}\u{31}{}", proj, targ), (xs, ys));
        }
        map
    };
    static ref MASSES: HashMap<String, f64> = {
        let mut map = HashMap::new();
        map.insert("34S".to_string(), 33.967867012);
        map.insert("34Cl".to_string(), 33.973762491);
        map.insert("34Ar".to_string(), 33.980270093);
        map.insert("37Cl".to_string(), 36.965902584);
        map.insert("37Ar".to_string(), 36.966776314);
        map.insert("37K".to_string(), 36.973375889);
        map
    };
}

/// Calculate the energy loss of a projectile in a target.
///
/// * proj is the name of the projectile (`"34S"`, `"34Cl"`, `"34Ar"`, `"37Cl"`, `"37Ar"`, `"37K"`)
/// * e is the total kinetic energy of the projectile in MeV
/// * targ is the name of the target (`"Butane"`, `"Mylar"`, or `"He"`)
/// * thick is the thickness of the target in mg/cm^2
pub fn eloss(proj: &str, e: f64, targ: &str, thick: f64) -> f64 {
    let stop = &STOPPING_POWERS[&format!("{}\u{31}{}", proj, targ)];
    let mass = MASSES[proj];
    let step_size = 1e-5;

    let mut energy_u = e / mass;
    let mut rem_thick = thick;
    let d_thick = rem_thick * step_size;
    while rem_thick > 0.0 && energy_u > 0.0 {
        let s = interpolate(energy_u, &stop.0, &stop.1).to_value().unwrap();
        let eloss = s * d_thick / mass;
        energy_u -= eloss;
        rem_thick -= d_thick;
    }

    e - energy_u * mass
}

#[modinit(eloss)]
fn init_mod(py: Python, m: &PyModule) -> PyResult<()> {

    #[pyfn(m, "eloss")]
    // ``#[pyfn()]` converts the arguments from Python objects to Rust values
    // and the Rust return value back into a Python object.
    fn eloss_py(proj: &str, e: f64, targ: &str, thick: f64) -> PyResult<f64> {
       let out = eloss(proj, e, targ, thick);
       Ok(out)
    }

    Ok(())
}
