use atmos::gases::*;

use atmos::constants::*;

pub struct AtmosMixture {
    pub gases: [f32;GASES.len()],
    pub temperature: f32,
    pub volume: f32,
}

impl Default for AtmosMixture {
    fn default() -> AtmosMixture {
        AtmosMixture {
            gases: [0.0;GASES.len()], 
            temperature: 293.15, 
            volume: 2500.0 
        }
    }
}

impl AtmosMixture {
    pub fn empty(&mut self) {
        self.gases.iter_mut().for_each(|x| *x=0.0);
        self.temperature = 293.15;
    }
    pub fn parse_gas_string(&mut self,gas_string: &str) {
        self.empty();
        let arg_split = gas_string.split(";");
        for entry in arg_split {
            let mut indiv_split = entry.split("=");
            let name = indiv_split.next();
            let value = indiv_split.next();
            if name.is_some() && value.is_some() {
                if name == Some("TEMP") {
                    self.temperature = value.unwrap().parse().unwrap();
                } else if GASES_BY_ID.contains_key(name.unwrap()) {
                    self.gases[GASES_BY_ID[name.unwrap()]] = value.unwrap().parse().unwrap();
                }
            }
        }
    }
    pub fn heat_capacity(&self) -> f32 {
        return self.gases.iter().enumerate().fold(0.0,|x,(i,amount)| x + GASES[i].heat_capacity(*amount));
    }
    pub fn total_moles(&self) -> f32 {
        return self.gases.iter().fold(0.0,|x,amount| x + amount);
    }
    pub fn pressure(&self) -> f32 {
        return (self.total_moles()*self.temperature*GAS_CONSTANT)/self.volume;
    }
    pub fn thermal_energy(&self) -> f32 {
        return self.temperature * self.heat_capacity()
    }
    pub fn merge(&mut self, other: AtmosMixture) {
        self.temperature = (self.thermal_energy() + other.thermal_energy()) / (self.heat_capacity() + other.heat_capacity());
        for (i,amount) in self.gases.iter_mut().enumerate() {
            *amount = *amount + other.gases[i];
        }
    }
    pub fn remove(&mut self, amount: f32) -> AtmosMixture {
        use std::cmp::min;
        use super::float_ord::FloatOrd;
        let mut removed : AtmosMixture = Default::default();
        let cached_total = self.total_moles();
        let FloatOrd(real_amount) = min(FloatOrd(cached_total),FloatOrd(amount));
        for (i,this_amount) in self.gases.iter_mut().enumerate() {
            removed.gases[i] = (*this_amount / cached_total)*real_amount;
            *this_amount -= removed.gases[i];
        }
        return removed;
    }
    pub fn remove_ratio(&mut self, ratio: f32) -> AtmosMixture {
        use std::cmp::min;
        use super::float_ord::FloatOrd;
        let mut removed : AtmosMixture = Default::default();
        let FloatOrd(real_ratio) = min(FloatOrd(1.0),FloatOrd(ratio));
        for (i,this_amount) in self.gases.iter_mut().enumerate() {
            removed.gases[i] = *this_amount * real_ratio;
            *this_amount -= removed.gases[i];
        }
        return removed;
    }
    pub fn share(&mut self,other : &mut AtmosMixture, adjacent_turfs: i16) -> (f32,f32) { // returns pressure, moles
        let old_self_heat_capacity = self.heat_capacity();
        let old_other_heat_capacity = other.heat_capacity();
        let mut heat_capacity_self_to_other = 0.0;
        let mut heat_capacity_other_to_self = 0.0;
        let mut moved_moles = 0.0;
        let mut abs_moved_moles = 0.0;
        for (i,this_amount) in self.gases.iter_mut().enumerate() {
            let delta = (*this_amount - other.gases[i])/(adjacent_turfs+1) as f32;
            let heat_cap_delta = GASES[i].heat_capacity(delta);
            if delta > 0.0
            {
                heat_capacity_self_to_other += heat_cap_delta;
            }
            else
            {
                heat_capacity_other_to_self -= heat_cap_delta;
            }
            let abs_moles_this_time=delta.abs();
            if abs_moles_this_time>0.1
            {
                *this_amount -= delta;
                other.gases[i] += delta;
                moved_moles += delta;
                abs_moved_moles += abs_moles_this_time;
            }
        }
        if abs_moved_moles > 0.1 {
            let new_self_heat_capacity = old_self_heat_capacity + heat_capacity_other_to_self - heat_capacity_self_to_other;
            let new_other_heat_capacity = old_other_heat_capacity + heat_capacity_self_to_other - heat_capacity_other_to_self;
            self.temperature = (old_self_heat_capacity * self.temperature - heat_capacity_self_to_other * self.temperature + heat_capacity_other_to_self * self.temperature) / new_self_heat_capacity;
            other.temperature = (old_other_heat_capacity * other.temperature - heat_capacity_other_to_self * other.temperature + heat_capacity_self_to_other * other.temperature) / new_other_heat_capacity;
            self.temperature_share(other,OPEN_HEAT_TRANSFER_COEFFICIENT);
            return ((self.temperature*(self.total_moles() + moved_moles)-other.temperature*(other.total_moles() - moved_moles) * GAS_CONSTANT) / self.volume,abs_moved_moles);
        }
        else
        {
            self.temperature_share(other,OPEN_HEAT_TRANSFER_COEFFICIENT);
            return (0.0, 0.0);
        }
    }
    pub fn temperature_share(&mut self, other: &mut AtmosMixture, conduction_coefficient: f32) -> f32 {
        use std::cmp::max;
        use super::float_ord::FloatOrd;
        let temperature_delta = self.temperature-other.temperature;
        if temperature_delta.abs() < 0.1
        {
            return other.temperature;
        }
        let self_heat_capacity = self.heat_capacity();
        let other_heat_capacity = other.heat_capacity();
        let heat = conduction_coefficient*temperature_delta*(self_heat_capacity*other_heat_capacity/(self_heat_capacity+other_heat_capacity));
        self.temperature = max(FloatOrd(self.temperature - heat/self_heat_capacity), FloatOrd(CMB_TEMP)).0;
        other.temperature = max(FloatOrd(other.temperature + heat/other_heat_capacity), FloatOrd(CMB_TEMP)).0;
        return other.temperature;
    }
    pub fn temperature_share_turf(&mut self, conduction_coefficient: f32, other_temperature: f32, other_heat_capacity: f32) -> f32 {
        use std::cmp::max;
        use super::float_ord::FloatOrd;
        let temperature_delta = self.temperature-other_temperature;
        let self_heat_capacity = self.heat_capacity();
        let heat = conduction_coefficient*temperature_delta*(self_heat_capacity*other_heat_capacity/(self_heat_capacity+other_heat_capacity));
        self.temperature = max(FloatOrd(self.temperature - heat/self_heat_capacity), FloatOrd(CMB_TEMP)).0;
        return max(FloatOrd(other_temperature + heat/other_heat_capacity), FloatOrd(CMB_TEMP)).0;
    }
    pub fn react(&mut self) -> Vec<Vec<String>>
    {
        use atmos::reaction::*;
        let mut ret = [].to_vec();
        for reaction in REACTIONS {
            let (reaction_return, reaction_vec) = reaction(self);
            if reaction_return == ReactionResult::StopReactions {
                return [].to_vec();
            }
            else {
                ret.push(reaction_vec);
            }
        }
        return ret;
    }
    pub fn to_params(&mut self) -> String
    {
        let mut params : String = format!("TEMP={};",self.temperature);
        for (i,this_amount) in self.gases.iter_mut().enumerate()
        {
            if this_amount > &mut 0.0 {
                params = format!("{}{}={};",params,GASES[i].id,this_amount);
            }
        }
        return params;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gas_parsing() {
        let mut mix : AtmosMixture = Default::default();
        mix.parse_gas_string("TEMP=293.15;o2=22;n2=82;");
        assert!(mix.gases[GASES_BY_ID["o2"]] == 22.0);
        mix.parse_gas_string("TEMP=304.3;plasma=40;tritium=40;miasma=30;no2=40;");
        assert_eq!(mix.to_params(),"TEMP=304.3;plasma=40;tritium=40;miasma=30;no2=40;");
        mix.parse_gas_string("TEMP=303.4;palsma=30;tritium=30;"); // sic
        assert_eq!(mix.to_params(),"TEMP=303.4;tritium=30;");
    }
    #[test]
    fn test_sharing() {
        let mut mix1 : AtmosMixture = Default::default();
        let mut mix2 : AtmosMixture = Default::default();
        mix2.parse_gas_string("TEMP=304.3;plasma=40;tritium=40;miasma=30;no2=40;");
        let initial_total_moles = mix1.total_moles() + mix2.total_moles();
        mix1.share(&mut mix2,1);
        assert_eq!(mix1.total_moles()+mix2.total_moles(),initial_total_moles);
    }
}