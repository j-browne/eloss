extern crate eloss;

use eloss::eloss;
use std::iter::repeat;
use eloss::MOLAR_MASSES;

const IC_TEMP: f64 = 300.0; // K
const JET_DIST: f64 = 0.3; // cm
const AVOGADRO_CONSTANT: f64 = 6.022140857e23; // 1/mol
const GAS_CONSTANT: f64 = 8.3144598; // J/mol/K
const MYLAR_DENSITY: f64 = 1.39; // g/cm^3
const REACTION_LOCATION: f64 = 0.5;

#[derive(Debug, Clone)]
enum Error {}


#[derive(Debug, Clone)]
pub struct Projectile {
    nuc: &'static str,
    energy: f64,
}

impl Projectile {
    pub fn new(nuc: &'static str, energy: f64) -> Self {
        Self { nuc, energy }
    }

    pub fn nuc(&self) -> &'static str {
        &self.nuc
    }

    pub fn energy(&self) -> f64 {
        self.energy
    }

    pub fn set_energy(&mut self, energy: f64) {
        self.energy = energy;
    }
}

#[derive(Debug, Clone)]
pub struct Target {
    material: &'static str,
    /// mg/cm^2
    thickness: f64,
    density: f64,
}

impl Target {
    pub fn new(material: &'static str) -> Self {
        Self {
            material,
            thickness: 0.0,
            density: 0.0,
        }
    }

    pub fn material(&self) -> &'static str {
        &self.material
    }

    /// thickness: mg/cm^2
    pub fn thickness(&self) -> f64 {
        self.thickness
    }

    /// density: g/cm^3
    pub fn density(&self) -> f64 {
        self.density
    }

    /// distance: cm
    pub fn distance(&self) -> f64 {
        self.thickness() / self.density() / 1000.0
    }

    /// molar_mass: g/mol
    pub fn molar_mass(&self) -> f64 {
        MOLAR_MASSES[self.material]
    }

    /// density: g/cm^3
    pub fn set_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// rhoa: atoms/cm^2
    /// distance: cm
    /// thickness: mg/cm^2
    /// density: g/cm^3
    pub fn set_density_thickness_with_rhoa_distance(mut self, rhoa: f64, distance: f64) -> Self {
        self.thickness = rhoa / AVOGADRO_CONSTANT * (1000.0 * self.molar_mass());
        self.density = (self.thickness / 1000.0) / distance;
        self
    }

    /// press: torr
    /// temp: K
    /// density: g/cm^3
    pub fn set_density_with_press_temp(mut self, press: f64, temp: f64) -> Self {
        self.density = ((press * 133.322) * (self.molar_mass()) / GAS_CONSTANT / temp) / 1000000.0;
        self
    }

    /// distance: cm
    /// thickness: mg/cm^2
    pub fn set_thickness_with_distance(mut self, distance: f64) -> Self {
        self.thickness = 1000.0 * self.density * distance;
        self
    }
}


struct Setup {
    proj_1: Projectile,
    proj_2: Projectile,
    jet_targs_1: Vec<Target>,
    jet_targs_2: Vec<Target>,
    window_targs: Vec<Target>,
    ic_targs: Vec<Target>,
}

impl Setup {
    pub fn new(proj_1: Projectile, proj_2: Projectile, ic_press: f64, rhoa: f64) -> Self {
        let jet_targs_1 = vec![
            Target::new("He")
                .set_density_thickness_with_rhoa_distance(rhoa, REACTION_LOCATION * JET_DIST),
        ];
        let jet_targs_2 = vec![
            Target::new("He")
                .set_density_thickness_with_rhoa_distance(rhoa, (1.0 - REACTION_LOCATION) * JET_DIST),
        ];
        let window_targs = vec![
            Target::new("Mylar")
                .set_density(MYLAR_DENSITY)
                .set_thickness_with_distance(3e-4),
        ];
        let ic_targs = vec![
            Target::new("Butane")
                .set_density_with_press_temp(ic_press, IC_TEMP)
                .set_thickness_with_distance(2.0),
            Target::new("Butane")
                .set_density_with_press_temp(ic_press, IC_TEMP)
                .set_thickness_with_distance(3.66),
            Target::new("Butane")
                .set_density_with_press_temp(ic_press, IC_TEMP)
                .set_thickness_with_distance(3.66),
            Target::new("Butane")
                .set_density_with_press_temp(ic_press, IC_TEMP)
                .set_thickness_with_distance(7.32),
            Target::new("Butane")
                .set_density_with_press_temp(ic_press, IC_TEMP)
                .set_thickness_with_distance(18.3),
        ];
        Self {
            proj_1,
            proj_2,
            jet_targs_1,
            jet_targs_2,
            window_targs,
            ic_targs,
        }
    }

    fn set_jet_rhoa(&mut self, rhoa: f64) {
        let mut old_targs = Vec::new();
        std::mem::swap(&mut old_targs, &mut self.jet_targs_1);
        for t in old_targs {
            let t = t.set_density_thickness_with_rhoa_distance(rhoa, REACTION_LOCATION * JET_DIST);
            self.jet_targs_1.push(t);
        }
        let mut old_targs = Vec::new();
        std::mem::swap(&mut old_targs, &mut self.jet_targs_2);
        for t in old_targs {
            let t = t.set_density_thickness_with_rhoa_distance(rhoa, (1.0 - REACTION_LOCATION) * JET_DIST);
            self.jet_targs_2.push(t);
        }
    }

    fn set_ic_press(&mut self, ic_press: f64) {
        let mut old_targs = Vec::new();
        std::mem::swap(&mut old_targs, &mut self.ic_targs);
        for t in old_targs {
            let d = t.distance();
            let t = t.set_density_with_press_temp(ic_press, IC_TEMP)
                .set_thickness_with_distance(d);
            self.ic_targs.push(t);
        }
    }

    fn calculate(&self) -> Vec<f64> {
        let mut e_losses = vec![];
        let mut e_diff = 0.0;
        for (mut p, t) in repeat(self.proj_1.clone()).zip(self.jet_targs_1.iter())
                    .chain(repeat(self.proj_2.clone()).zip(self.jet_targs_2.iter()))
                    .chain(repeat(self.proj_2.clone()).zip(self.window_targs.iter()))
                    .chain(repeat(self.proj_2.clone()).zip(self.ic_targs.iter()))
        {
            let e_curr = p.energy() - e_diff;
            p.set_energy(e_curr);
            let e_loss = eloss(p.nuc(), p.energy(), t.material(), t.thickness());
            e_diff += e_loss;
            e_losses.push(e_loss);
        }
        e_losses
    }
}

fn main() -> Result<(), Error> {
    let mut setup = Setup::new(Projectile::new("34Ar", 55.4), Projectile::new("34Ar", 55.4), 15.0, 1e19);

    for rhoa in &[5e18, 1e19] {
        setup.set_jet_rhoa(*rhoa);
        for ic_press in &[14.0, 15.0, 16.0] {
            setup.set_ic_press(*ic_press);
            println!("{:?}", setup.calculate());
        }
    }

    Ok(())
}
