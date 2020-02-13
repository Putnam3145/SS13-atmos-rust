#[macro_use]
mod atmos;

#[macro_use]
mod byond;

byond_fn! { react(atmos_string) {
    use atmos::mix::AtmosMixture;
    let mut mix : AtmosMixture = Default::default();
    mix.parse_gas_string(atmos_string);
    let reaction_results = mix.react();
    let mut ret_string : String = mix.to_params();
    for result in reaction_results {
        if(result.len() > 0) {
            ret_string = format!("{}$",ret_string);
            for part_of_result in result {
                ret_string = format!("{}{};",ret_string,part_of_result);
            }
        }
    }
    Some(ret_string)
} }

byond_fn! { share(atmos_string,sharer_atmos_string,adjacent_turfs_string) {
    use atmos::mix::AtmosMixture;
    let mut mix1 : AtmosMixture = Default::default();
    let mut mix2 : AtmosMixture = Default::default();
    mix1.parse_gas_string(atmos_string);
    mix2.parse_gas_string(sharer_atmos_string);
    let (pressure_change,last_share) = mix1.share(&mut mix2,adjacent_turfs_string.parse().unwrap());
    Some(format!("{}${}${}${}",mix1.to_params(),mix2.to_params(),last_share,pressure_change))
} }