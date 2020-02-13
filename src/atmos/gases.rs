pub struct GasDatum {
    pub id: &'static str,
    pub name: &'static str,
    pub specific_heat: f32,
    pub fusion_power: f32,
}

impl GasDatum {
    pub fn heat_capacity(&self,amount: f32) -> f32 {
        return amount * self.specific_heat;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_heat_cap() {
        assert_eq!(OXYGEN.heat_capacity(100.0),2_000.0);
    }
}

pub const OXYGEN : GasDatum = GasDatum {
    id: "o2",
    name: "Oxygen",
    specific_heat: 20.0,
    fusion_power: 0.0,
};

pub const NITROGEN : GasDatum = GasDatum {
    id: "n2",
    name: "Nitrogen",
    specific_heat: 20.0,
    fusion_power: 0.0,
};

pub const CARBON_DIOXIDE : GasDatum = GasDatum {
    id: "co2",
    name: "Carbon Dioxide",
    specific_heat: 30.0,
    fusion_power: 3.0,
};

pub const PLASMA : GasDatum = GasDatum {
    id: "plasma",
    name: "Plasma",
    specific_heat: 200.0,
    fusion_power: 0.0,
};

pub const WATER_VAPOR : GasDatum = GasDatum {
    id: "water_vapor",
    name: "Water Vapor",
    specific_heat: 40.0,
    fusion_power: 8.0,
};

pub const HYPERNOBLIUM : GasDatum = GasDatum {
    id: "nob",
    name: "Hyper-noblium",
    specific_heat: 2000.0,
    fusion_power: 0.0,
};

pub const NITROUS_OXIDE : GasDatum = GasDatum {
    id: "n2o",
    name: "Nitrous Oxide",
    specific_heat: 40.0,
    fusion_power: 0.0,
};

pub const NITRYL : GasDatum = GasDatum {
    id: "no2",
    name: "Nitryl",
    specific_heat: 20.0,
    fusion_power: 15.0,
};

pub const TRITIUM : GasDatum = GasDatum {
    id: "tritium",
    name: "Tritium",
    specific_heat: 10.0,
    fusion_power: 1.0,
};

pub const BZ : GasDatum = GasDatum {
    id: "bz",
    name: "BZ",
    specific_heat: 20.0,
    fusion_power: 8.0,
};

pub const STIMULUM : GasDatum = GasDatum {
    id: "stim",
    name: "Stimulum",
    specific_heat: 5.0,
    fusion_power: 7.0,
};

pub const PLUOXIUM : GasDatum = GasDatum {
    id: "pluox",
    name: "Pluoxium",
    specific_heat: 80.0,
    fusion_power: 10.0,
};

pub const MIASMA : GasDatum = GasDatum {
    id: "miasma",
    name: "Miasma",
    specific_heat: 0.00001,
    fusion_power: 50.0,
};

// this list is hardcoded sorted by rarity, greatest to smallest, so we don't need to store rarity
pub const GASES: &'static [GasDatum] = &[ 
    NITROGEN, // 1000
    OXYGEN, // 900
    PLASMA, // 800
    CARBON_DIOXIDE, // 700
    NITROUS_OXIDE,  // 600
    WATER_VAPOR, // 500
    BZ, // 400
    TRITIUM, // 300
    MIASMA, // 250
    PLUOXIUM, // 200
    NITRYL, // 100
    HYPERNOBLIUM, // 50
    STIMULUM // 1
    ];

use super::phf::phf_map;

pub const GASES_BY_ID: phf::Map<&'static str,usize> = phf_map! {
    "n2" => 0,
    "o2" => 1,
    "plasma" => 2,
    "co2" => 3,
    "n2o" => 4,
    "water_vapor" => 5,
    "bz" => 6,
    "tritium" => 7,
    "miasma" => 8,
    "pluox" => 9,
    "no2" => 10,
    "nob" => 11,
    "stim" => 12
};

