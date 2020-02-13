use atmos::mix::*;

use atmos::gases::*;

use atmos::constants::*;

#[derive(PartialEq)]
pub enum ReactionResult {
    NoReaction,
    Reacting,
    StopReactions
}

fn nobliumsupression(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.gases[GASES_BY_ID["nob"]] >= 5.0 {
        return (ReactionResult::StopReactions,[].to_vec())
    }
    return (ReactionResult::NoReaction,[].to_vec());
}

fn tritfire(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    let initial_oxy = air.gases[GASES_BY_ID["o2"]];
    let initial_trit = air.gases[GASES_BY_ID["tritium"]];
    if air.temperature < FIRE_MINIMUM_TEMPERATURE_TO_EXIST || initial_trit < MINIMUM_MOLE_COUNT  || initial_oxy < MINIMUM_MOLE_COUNT {
        return (ReactionResult::NoReaction,[].to_vec());
    } else {
        let tritium_burn_oxy_factor = 100.0;
        let tritium_burn_trit_factor = 10.0;
        let fire_hydrogen_energy_released = 560000.0;
        let old_energy = air.thermal_energy();
        let mut energy_released = 0.0;
        let mut burned_fuel;
        if initial_oxy < initial_trit {
            burned_fuel = air.gases[GASES_BY_ID["o2"]] / tritium_burn_oxy_factor;
            if burned_fuel > initial_trit {
                burned_fuel = initial_trit;
            }
            air.gases[GASES_BY_ID["tritium"]] -= burned_fuel;
        } else {
            burned_fuel = initial_trit;
            air.gases[GASES_BY_ID["tritium"]] *= 1.0 - 1.0/tritium_burn_trit_factor;
            air.gases[GASES_BY_ID["o2"]] -= air.gases[GASES_BY_ID["tritium"]];
            energy_released += fire_hydrogen_energy_released * burned_fuel * (tritium_burn_trit_factor - 1.0);
        }
        let mut strings = [].to_vec();
        if burned_fuel > 0.0 {
            energy_released += fire_hydrogen_energy_released * burned_fuel;
            air.gases[GASES_BY_ID["water_vapor"]] += burned_fuel;
            strings = ["tritfire".to_string(),burned_fuel.to_string(),energy_released.to_string()].to_vec();
        }
        if energy_released > 0.0 {
            air.temperature = (old_energy+energy_released)/air.heat_capacity();
        }
        return (ReactionResult::Reacting,strings);
    }
}

fn plasmafire(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.temperature < FIRE_MINIMUM_TEMPERATURE_TO_EXIST || air.gases[GASES_BY_ID["plasma"]] < MINIMUM_MOLE_COUNT  || air.gases[GASES_BY_ID["o2"]] < MINIMUM_MOLE_COUNT
    {
        return (ReactionResult::NoReaction,[].to_vec());
    } else {
        use std::cmp::min;
        use super::float_ord::FloatOrd;
        let plasma_upper_temperature = 1390.0+T0C;
        let oxygen_burn_rate_base = 1.4;
        let super_saturation_threshold = 96.0;
        let plasma_oxygen_fullburn = 10.0;
        let plasma_burn_rate_delta = 9.0;
        let fire_plasma_energy_released = 3000000.0;
        let super_saturation;
        let temperature_scale;
        let mut plasma_burn_rate;
        let oxygen_burn_rate;
        let old_energy = air.thermal_energy();
        let energy_released;
        if air.temperature>plasma_upper_temperature {
            temperature_scale = 1.0;
        } else {
            temperature_scale = (air.temperature-FIRE_MINIMUM_TEMPERATURE_TO_EXIST)/(plasma_upper_temperature-FIRE_MINIMUM_TEMPERATURE_TO_EXIST);
        }
        let mut strings = [].to_vec();
        if temperature_scale > 0.0
        {
            oxygen_burn_rate = oxygen_burn_rate_base-temperature_scale;
            super_saturation = air.gases[GASES_BY_ID["o2"]]/air.gases[GASES_BY_ID["plasma"]]>super_saturation_threshold;
            if air.gases[GASES_BY_ID["o2"]]>air.gases[GASES_BY_ID["plasma"]]*plasma_oxygen_fullburn
            {
                plasma_burn_rate = air.gases[GASES_BY_ID["plasma"]]*temperature_scale/plasma_burn_rate_delta;
            }
            else
            {
                plasma_burn_rate = (temperature_scale*(air.gases[GASES_BY_ID["o2"]]/plasma_oxygen_fullburn))/plasma_burn_rate_delta;
            }
            plasma_burn_rate = min(min(FloatOrd(plasma_burn_rate),FloatOrd(air.gases[GASES_BY_ID["plasma"]])),FloatOrd(air.gases[GASES_BY_ID["o2"]]/oxygen_burn_rate)).0;
            air.gases[GASES_BY_ID["plasma"]] -= plasma_burn_rate;
            air.gases[GASES_BY_ID["o2"]] -= plasma_burn_rate*oxygen_burn_rate;
            if super_saturation
            {
                air.gases[GASES_BY_ID["tritium"]]+=plasma_burn_rate;
            }
            else
            {
                air.gases[GASES_BY_ID["co2"]]+=plasma_burn_rate;
            }
            energy_released = fire_plasma_energy_released * plasma_burn_rate;
            let reacted = plasma_burn_rate*(1.0+oxygen_burn_rate)>0.0;
            if energy_released > 0.0
            {
                air.temperature = (old_energy+energy_released)/air.heat_capacity();
            }
            if reacted
            {
                strings = ["plasfire".to_string(),(plasma_burn_rate*(1.0+oxygen_burn_rate)).to_string()].to_vec();
            }
        }
        return (ReactionResult::Reacting,strings);
    }
}

fn fusion(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.temperature < 10000.0 || air.gases[GASES_BY_ID["tritium"]] < 1.0 || air.gases[GASES_BY_ID["co2"]] < 250.0 || air.gases[GASES_BY_ID["plasma"]] < 250.0 {
        return (ReactionResult::NoReaction,[].to_vec())
    } else {
        use std::f32::consts::PI;
        let toroid_volume_breakeven = 1000.0;
        let instability_gas_factor = 0.003;
        let plasma_binding_energy = 20_000_000.0;
        let fusion_tritium_moles_used = 1.0;
        let fusion_instability_endothermality = 2.0;
        let fusion_tritium_conversion_coefficient = 1e-10;
        let mut reaction_energy = 0.0 ;
        let initial_energy = air.thermal_energy();
        let initial_plasma = air.gases[GASES_BY_ID["plasma"]];
        let initial_carbon = air.gases[GASES_BY_ID["co2"]];
        let scale_factor = air.volume / PI;
        //The size of the phase space hypertorus
        let toroidal_size = (3.0*PI)+((air.volume-toroid_volume_breakeven)/toroid_volume_breakeven).atan();
        //3*PI above rather than 2.0*PI because atan can return -pi
        let mut gas_power = 0.0; // uh??
        for (i, amount) in air.gases.iter().enumerate()
        {
            gas_power+=GASES[i].fusion_power*amount;
        }
        let instability = (gas_power*instability_gas_factor).powf(2.0)%toroidal_size;
        let mut plasma = (initial_plasma - FUSION_MOLE_THRESHOLD) / scale_factor;
        let mut carbon = (initial_carbon - FUSION_MOLE_THRESHOLD) / scale_factor;
        plasma = ((plasma - (instability*(carbon.sin()))%toroidal_size)).abs();
        //count the rings. ss13's modulus is positive, this ain't, who knew
        carbon = ((carbon - plasma)%toroidal_size).abs();
        air.gases[GASES_BY_ID["plasma"]] = plasma*scale_factor + FUSION_MOLE_THRESHOLD;
        air.gases[GASES_BY_ID["co2"]] = carbon*scale_factor + FUSION_MOLE_THRESHOLD;
        let delta_plasma = initial_plasma - air.gases[GASES_BY_ID["plasma"]];
        reaction_energy += delta_plasma*plasma_binding_energy;
        if instability < fusion_instability_endothermality && reaction_energy < 0.0
        {
            reaction_energy = 0.0;
        }
        else if reaction_energy < 0.0
        {
            reaction_energy *= (instability-fusion_instability_endothermality).sqrt();
        }
        if air.thermal_energy() + reaction_energy < 0.0
        {
            air.gases[GASES_BY_ID["plasma"]] = initial_plasma;
            air.gases[GASES_BY_ID["co2"]] = initial_carbon;
            return (ReactionResult::NoReaction,[].to_vec());
        }
        air.gases[GASES_BY_ID["tritium"]] -= fusion_tritium_moles_used;
        if reaction_energy > 0.0
        {
            air.gases[GASES_BY_ID["o2"]] += fusion_tritium_moles_used*(reaction_energy*fusion_tritium_conversion_coefficient);
            air.gases[GASES_BY_ID["n2o"]] += fusion_tritium_moles_used*(reaction_energy*fusion_tritium_conversion_coefficient);
        }
        else
        {
            air.gases[GASES_BY_ID["bz"]] += fusion_tritium_moles_used*(reaction_energy*-fusion_tritium_conversion_coefficient);
            air.gases[GASES_BY_ID["no2"]] += fusion_tritium_moles_used*(reaction_energy*-fusion_tritium_conversion_coefficient);
        }
        if reaction_energy != 0.0 {
            use std::cmp::max;
            use super::float_ord::FloatOrd;
            let particle_chance_constant = -20_000_000.0;
            let particle_chance = (particle_chance_constant/(reaction_energy-particle_chance_constant)) + 1.0;
            let FloatOrd(rad_power) = max(FloatOrd((-2000.0/instability)+1000.0),FloatOrd(0.0));
            air.temperature = (initial_energy+reaction_energy)/air.heat_capacity();
            return (ReactionResult::Reacting,[
                "fusion".to_string(),
                particle_chance.to_string(),
                rad_power.to_string()].to_vec());
        }
        return (ReactionResult::NoReaction,[].to_vec());
    }
}

fn nitrylformation(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.temperature < FIRE_MINIMUM_TEMPERATURE_TO_EXIST * 400.0 || air.gases[GASES_BY_ID["n2o"]] < 0.5 || air.gases[GASES_BY_ID["o2"]] < 20.0 || air.gases[GASES_BY_ID["n2"]] < 20.0 {
        return (ReactionResult::NoReaction,[].to_vec())
    } else {
        use std::cmp::min;
        use super::float_ord::FloatOrd;
        let nitryl_formation_energy = 100000.0;
        let old_energy = air.thermal_energy();
        let FloatOrd(reaction_efficiency) = min(FloatOrd(air.temperature/(FIRE_MINIMUM_TEMPERATURE_TO_EXIST*100.0)),min(FloatOrd(air.gases[GASES_BY_ID["o2"]]),FloatOrd(air.gases[GASES_BY_ID["n2"]])));
        let energy_used = reaction_efficiency*nitryl_formation_energy;
        if (air.gases[GASES_BY_ID["o2"]] < reaction_efficiency ) || (air.gases[GASES_BY_ID["n2"]] < reaction_efficiency) //Shouldn't produce gas from nothing.
        {
            return (ReactionResult::NoReaction,[].to_vec());
        }
        air.gases[GASES_BY_ID["o2"]] -= reaction_efficiency;
        air.gases[GASES_BY_ID["n2"]] -= reaction_efficiency;
        air.gases[GASES_BY_ID["no2"]] += reaction_efficiency*2.0;
        air.temperature = (old_energy-energy_used)/air.heat_capacity();
        return (ReactionResult::Reacting,[].to_vec());
    }
}

fn bzformation(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.gases[GASES_BY_ID["n2o"]] < 10.0 || air.gases[GASES_BY_ID["plasma"]] < 10.0 {
        return (ReactionResult::NoReaction,[].to_vec())
    } else {
        let fire_carbon_energy_released = 100_000.0;
        let old_energy = air.thermal_energy();
        let old_pressure = air.pressure();
        use std::cmp::min;
        use std::cmp::max;
        use super::float_ord::FloatOrd;
        //nobody ever said this code was good. let nobody ever try to tell you this code is good
        let FloatOrd(reaction_efficiency) = min(
                FloatOrd(1.0/(old_pressure/(0.1*ATMOSPHERE))*(max(
                        FloatOrd(air.gases[GASES_BY_ID["plasma"]]/air.gases[GASES_BY_ID["n2o"]]),
                        FloatOrd(1.0))).0),
                min(
                    FloatOrd(air.gases[GASES_BY_ID["n2o"]]),
                    FloatOrd(air.gases[GASES_BY_ID["plasma"]]/2.0)));
        let energy_released = 2.0*reaction_efficiency*fire_carbon_energy_released;
        if (air.gases[GASES_BY_ID["n2o"]] < reaction_efficiency )|| (air.gases[GASES_BY_ID["plasma"]] < (2.0*reaction_efficiency) || energy_released <= 0.0 ) //Shouldn't produce gas from nothing.
        {
            return (ReactionResult::NoReaction,[].to_vec());
        }
        air.gases[GASES_BY_ID["bz"]] += reaction_efficiency;
        if reaction_efficiency == air.gases[GASES_BY_ID["n2o"]]
        {
            let FloatOrd(nitrous_balance_change) = min(FloatOrd(old_pressure),FloatOrd(1.0));
            air.gases[GASES_BY_ID["bz"]] -= nitrous_balance_change;
            air.gases[GASES_BY_ID["o2"]] += nitrous_balance_change;
        }
        air.gases[GASES_BY_ID["n2o"]] -= reaction_efficiency;
        air.gases[GASES_BY_ID["plasma"]]  -= 2.0*reaction_efficiency;
        //SSresearch.science_tech.add_point_type(TECHWEB_POINT_TYPE_DEFAULT, min((reaction_efficency**2)*BZ_RESEARCH_SCALE),BZ_RESEARCH_MAX_AMOUNT)
        //we already returned if energy_released is non-positive
        air.temperature = (old_energy+energy_released)/air.heat_capacity();
        return (ReactionResult::Reacting,
        ["bzFormation".to_string(),
         reaction_efficiency.to_string()].to_vec());
    }
}

fn stimformation(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.temperature < STIMULUM_HEAT_SCALE / 2.0 || air.gases[GASES_BY_ID["no2"]] < 30.0 || air.gases[GASES_BY_ID["bz"]] < 20.0 || air.gases[GASES_BY_ID["tritium"]] < 30.0 || air.gases[GASES_BY_ID["plasma"]] < 10.0 {
        return (ReactionResult::NoReaction,[].to_vec())
    } else {
        use std::cmp::min;
        use std::cmp::max;
        use super::float_ord::FloatOrd;
        let stimulum_first_rise = 0.65;
        let stimulum_first_drop = 0.065;
        let stimulum_second_rise = 0.0009;
        let stimulum_absolute_drop = 0.00000335;
        let old_energy = air.thermal_energy();
        let FloatOrd(heat_scale) = min(FloatOrd(air.temperature/STIMULUM_HEAT_SCALE),min(FloatOrd(air.gases[GASES_BY_ID["tritium"]]),min(FloatOrd(air.gases[GASES_BY_ID["plasma"]]),FloatOrd(air.gases[GASES_BY_ID["no2"]]))));
        let stim_energy_change = heat_scale + stimulum_first_rise*(heat_scale.powf(2.0)) - stimulum_first_drop*(heat_scale.powf(3.0)) + stimulum_second_rise*(heat_scale.powf(4.0)) - stimulum_absolute_drop*(heat_scale.powf(5.0));
        //i mean it's not THAT odd it's O(-T^5) which you might think doesn't matter much except fusion tends to get to levels where that dominates
        if air.gases[GASES_BY_ID["no2"]] < heat_scale || air.gases[GASES_BY_ID["tritium"]] < heat_scale || air.gases[GASES_BY_ID["plasma"]] < heat_scale
        {
            return (ReactionResult::NoReaction,[].to_vec());
        }
        air.gases[GASES_BY_ID["stim"]] += heat_scale/10.0;
        air.gases[GASES_BY_ID["tritium"]] -= heat_scale;
        air.gases[GASES_BY_ID["plasma"]] -= heat_scale;
        air.gases[GASES_BY_ID["no2"]] -= heat_scale;
        if stim_energy_change != 0.0 
        {
            air.temperature = (old_energy+stim_energy_change)/air.heat_capacity();
        }
        return (ReactionResult::Reacting,[
        ("stimformation".to_string()),(max(FloatOrd(0.0),FloatOrd(stim_energy_change)).0).to_string()].to_vec());
    }
}

fn nobliumformation(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.temperature < 5000000.0 || air.gases[GASES_BY_ID["tritium"]] < 5.0 || air.gases[GASES_BY_ID["n2"]] < 10.0 {
        return (ReactionResult::NoReaction,[].to_vec())
    } else {
        use std::cmp::min;
        use std::cmp::max;
        use super::float_ord::FloatOrd;
        let noblium_formation_energy = 2e9 * 1.0;
        let old_energy = air.thermal_energy();
        let FloatOrd(nob_formed) = min(FloatOrd(air.gases[GASES_BY_ID["n2"]]+air.gases[GASES_BY_ID["tritium"]]/100.0),min(FloatOrd(air.gases[GASES_BY_ID["tritium"]]/10.0),FloatOrd(air.gases[GASES_BY_ID["n2"]]/20.0)));
        let energy_taken = (nob_formed * noblium_formation_energy) / max(FloatOrd(air.gases[GASES_BY_ID["bz"]]),FloatOrd(1.0)).0;
        if (air.gases[GASES_BY_ID["tritium"]] < 10.0*nob_formed) || (air.gases[GASES_BY_ID["n2"]] < 20.0*nob_formed)
        {
            return (ReactionResult::NoReaction,[].to_vec());
        }
        air.gases[GASES_BY_ID["tritium"]] -= 10.0*nob_formed;
        air.gases[GASES_BY_ID["n2"]] -= 20.0*nob_formed;
        air.gases[GASES_BY_ID["nob"]] += nob_formed;
        if energy_taken != 0.0
        {
            air.temperature = (old_energy-energy_taken)/air.heat_capacity();
        }
        return (ReactionResult::Reacting,[
        "nobliumformation".to_string(),nob_formed.to_string()].to_vec());
    }
}

fn miaster(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>) {
    if air.temperature < FIRE_MINIMUM_TEMPERATURE_TO_EXIST+70.0 || air.gases[GASES_BY_ID["miasma"]] < MINIMUM_MOLE_COUNT || air.gases[GASES_BY_ID["water_vapor"]]/air.total_moles() > 0.1 {
        return (ReactionResult::NoReaction,[].to_vec())
    } else {
        use std::cmp::min;
        use super::float_ord::FloatOrd;
        let FloatOrd(cleaned_air) = min(FloatOrd(air.gases[GASES_BY_ID["miasma"]]),
        FloatOrd(20.0 + (air.temperature - FIRE_MINIMUM_TEMPERATURE_TO_EXIST - 70.0) / 20.0));
        air.gases[GASES_BY_ID["miasma"]] -= cleaned_air;
        air.gases[GASES_BY_ID["o2"]] += cleaned_air;
        air.temperature += cleaned_air * (0.002);
        return (ReactionResult::Reacting,[
        "miaster".to_string(),cleaned_air.to_string()].to_vec());
    }
}

pub const REACTIONS: &'static [fn(air: &mut AtmosMixture) -> (ReactionResult,Vec<String>)] = &[
    nobliumsupression, // infinity
    nobliumformation, // 6
    stimformation, // 5
    bzformation, // 4
    nitrylformation, // 3
    fusion, // 2
    tritfire, // -1
    plasmafire, // -2
    miaster // -10
];


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_reactions() {
        let mut mix = AtmosMixture {
            gases: [2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 0.0, 2500.0],
            temperature: 1000000.0,
            volume: 2500.0};
        let reaction_results = mix.react();
        let mut ret_string : String = mix.to_params();
        for result in reaction_results {
            if(result.len() > 0) {
                ret_string = format!("{}$",ret_string);
                for part_of_result in result {
                    ret_string = format!("{}{};",ret_string,part_of_result);
                    println!("{}",part_of_result);
                }
            }
            println!("{}",ret_string);
        }
    }
}
