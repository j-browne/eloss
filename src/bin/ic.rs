extern crate eloss;

use eloss::eloss;
use std::iter::repeat;
use eloss::MOLAR_MASSES;
use std::{fs, io, num};
use std::path::Path;
use std::collections::HashMap;
use std::str::FromStr;

const IC_TEMP: f64 = 300.0; // K
const JET_DIST: f64 = 0.3; // cm
const AVOGADRO_CONSTANT: f64 = 6.022140857e23; // 1/mol
const GAS_CONSTANT: f64 = 8.3144598; // J/mol/K
const MYLAR_DENSITY: f64 = 1.39; // g/cm^3
const REACTION_LOCATION: f64 = 0.5;

#[derive(Debug, Clone)]
enum Error {
    IO,
    ParseFloatError(num::ParseFloatError),
    ParseRunTypeError,
}

impl From<io::Error> for Error {
    fn from(_e: io::Error) -> Self {
        Error::IO
    }
}

impl From<num::ParseFloatError> for Error {
    fn from(e: num::ParseFloatError) -> Self {
        Error::ParseFloatError(e)
    }
}

#[allow(dead_code)]
#[derive(Debug)]
struct ValUnc {
    val: f64,
    unc_stat: f64,
    unc_sys: f64,
}

#[allow(dead_code)]
struct RunInfo {
    run_type: RunType,
    start_time: f64,
    stop_time: f64,
    cap_ic: Option<ValUnc>,
    cap_in: Option<ValUnc>,
    rhoa: Option<ValUnc>,
}

enum RunType {
    Run(f64),
    NozTest,
}

impl FromStr for RunType {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1.625" | "1.71" => Ok(RunType::Run(s.parse().or(Err(Error::ParseRunTypeError))?)),
            "noz_test" => Ok(RunType::NozTest),
            _ => Err(Error::ParseRunTypeError),
        }
    }
}

fn get_run_info<P: AsRef<Path>>(filename: P) -> Result<HashMap<String, RunInfo>, Error> {
    let mut run_info = HashMap::new();
    let data = fs::read_to_string(filename)?;
    for line in data.lines() {
        let x: Vec<_> = line.split_whitespace().collect();

        if x.len() == 0 || x[0].starts_with("#") {
            continue;
        }

        let run_name = x[0].to_string();
        let run_type = x[1].parse()?;
        let start_time = x[2].parse()?;
        let stop_time = x[3].parse()?;
        let cap_ic = if let (Ok(val), Ok(unc_stat), Ok(unc_sys)) = (x[4].parse(), x[5].parse(), x[6].parse()) {
            Some(ValUnc{ val, unc_stat, unc_sys })
        } else {
            None
        };
        let cap_in = if let (Ok(val), Ok(unc_stat), Ok(unc_sys)) = (x[7].parse(), x[8].parse(), x[9].parse()) {
            Some(ValUnc{ val, unc_stat, unc_sys })
        } else {
            None
        };
        let rhoa = if let (Ok(val), Ok(unc_stat), Ok(unc_sys)) = (x[10].parse(), x[11].parse(), x[12].parse()) {
            Some(ValUnc{ val, unc_stat, unc_sys })
        } else {
            None
        };
        run_info.insert(run_name, RunInfo{run_type, start_time, stop_time, cap_ic, cap_in, rhoa});
    }
    Ok(run_info)
}

#[derive(Debug, Clone)]
struct Projectile {
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
struct Target {
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

    pub fn set_proj_1(&mut self, p: Projectile) {
        self.proj_1 = p
    }

    pub fn set_proj_2(&mut self, p: Projectile) {
        self.proj_2 = p
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
    let run_info = get_run_info("run_info.txt")?;
    for (name, info) in run_info {
        for proj in &[
        Projectile::new("34S", 54.170),
        Projectile::new("34Cl", 54.179),
        Projectile::new("34Ar", 54.190),
        ] {
            setup.set_proj_1(proj.clone());
            setup.set_proj_2(proj.clone());
            if let (Some(rhoa_val_unc), Some(ic_press_val_unc)) = (info.rhoa.as_ref(), info.cap_ic.as_ref()) {
                let mut xs = Vec::new();
                let mut ys = Vec::new();
                let mut des = Vec::new();
                let mut es = Vec::new();
                for (rhoa, ic_press) in &[
                    (rhoa_val_unc.val, ic_press_val_unc.val),
                    (rhoa_val_unc.val + rhoa_val_unc.unc_sys, ic_press_val_unc.val + ic_press_val_unc.unc_sys),
                    (rhoa_val_unc.val - rhoa_val_unc.unc_sys, ic_press_val_unc.val - ic_press_val_unc.unc_sys)
                ] {
                    setup.set_jet_rhoa(*rhoa);
                    setup.set_ic_press(*ic_press);
                    let elosses = setup.calculate();
                    xs.push(elosses[4]);
                    ys.push(elosses[5]);
                    des.push(elosses[6]);
                    es.push(elosses[7]);
                }
                let x = ValUnc { val: xs[0], unc_stat: 0.0, unc_sys: f64::max(f64::abs(xs[1] - xs[0]), f64::abs(xs[2] - xs[0]))};
                let y = ValUnc { val: ys[0], unc_stat: 0.0, unc_sys: f64::max(f64::abs(ys[1] - ys[0]), f64::abs(ys[2] - ys[0]))};
                let de = ValUnc { val: des[0], unc_stat: 0.0, unc_sys: f64::max(f64::abs(des[1] - des[0]), f64::abs(des[2] - des[0]))};
                let e = ValUnc { val: es[0], unc_stat: 0.0, unc_sys: f64::max(f64::abs(es[1] - es[0]), f64::abs(es[2] - es[0]))};

                for chan in 0..32 {
                    println!("{}\tX\t{}\t{}\t{}\t{}", name, chan, proj.nuc(), x.val * 1000.0, x.unc_sys * 1000.0);
                    println!("{}\tY\t{}\t{}\t{}\t{}", name, chan, proj.nuc(), y.val * 1000.0, y.unc_sys * 1000.0);
                }
                println!("{}\tdE\t0\t{}\t{}\t{}", name, proj.nuc(), de.val * 1000.0, de.unc_sys * 1000.0);
                println!("{}\tE\t0\t{}\t{}\t{}", name, proj.nuc(), e.val * 1000.0, e.unc_sys * 1000.0);
            }
        }
    }

    Ok(())
}
