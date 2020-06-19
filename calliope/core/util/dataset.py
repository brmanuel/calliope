"""
Copyright (C) 2013-2019 Calliope contributors listed in AUTHORS.
Licensed under the Apache 2.0 License (see LICENSE file).

"""

from calliope import exceptions
import xarray as xr
import pandas as pd


def get_loc_techs(loc_techs, tech=None, loc=None):
    """
    Get a list of loc_techs associated with the given technology and/or location.
    If multiple of both loc and tech are given, the function will return any
    combination of members of loc and tech lists found in loc_techs.

    Parameters
    ----------
    loc_techs : list
        set of loc_techs to search for the relevant tech and/or loc
    tech : string or list of strings, default None
        technology/technologies to search for in the set of location:technology
    loc : string or list of strings, default None
        location(s) to search for in the set of location:technology

    Returns
    -------
    relevant_loc_techs : list of strings

    """
    # If both are strings, there is only one loc:tech possibility to look for
    if (
        isinstance(tech, str)
        and isinstance(loc, str)
        and "::".join([loc, tech]) in loc_techs
    ):
        relevant_loc_techs = ["::".join([loc, tech])]

    tech = [tech] if tech is not None and isinstance(tech, str) else tech
    loc = [loc] if loc is not None and isinstance(loc, str) else loc

    if tech and not loc:
        relevant_loc_techs = [i for i in loc_techs if i.split("::")[1] in tech]
    elif loc and not tech:
        relevant_loc_techs = [i for i in loc_techs if i.split("::")[0] in loc]
    elif loc and tech:
        loc_techs_set = set(tuple(i.split("::")) for i in loc_techs)
        possible_loc_techs = set((l, t) for l in loc for t in tech)
        relevant_loc_techs = [
            "::".join(i) for i in possible_loc_techs.intersection(loc_techs_set)
        ]
    else:
        relevant_loc_techs = [None]

    return relevant_loc_techs


def reorganise_xarray_dimensions(data):
    """
    Reorganise Dataset or DataArray dimensions to be alphabetical *except*
    `timesteps`, which must always come last in any DataArray's dimensions
    """

    if not (isinstance(data, xr.Dataset) or isinstance(data, xr.DataArray)):
        raise TypeError(
            "Must provide either xarray Dataset or DataArray to be reorganised"
        )

    steps = [i for i in ["datesteps", "timesteps"] if i in data.dims]

    if isinstance(data, xr.Dataset):
        new_dims = (sorted(list(set(data.dims.keys()) - set(steps)))) + steps
    elif isinstance(data, xr.DataArray):
        new_dims = (sorted(list(set(data.dims) - set(steps)))) + steps

    updated_data = data.transpose(*new_dims).reindex({k: data[k] for k in new_dims})

    return updated_data


def split_loc_techs(data_var, return_as="DataArray"):
    """
    Get a DataArray with locations technologies, and possibly carriers
    split into separate coordinates.

    Parameters
    ----------
    data_var : xarray DataArray
        Variable from Calliope model_data, to split loc_techs dimension
    return_as : string
        'DataArray' to return xarray DataArray, 'MultiIndex DataArray' to return
        xarray DataArray with loc_techs as a MultiIndex,
        or 'Series' to return pandas Series with dimensions as a MultiIndex

    Returns
    -------
    updated_data_var : xarray DataArray of pandas Series
    """

    # Separately find the loc_techs(_carriers) dimension and all other dimensions
    loc_tech_dim = [i for i in data_var.dims if "loc_tech" in i]
    if not loc_tech_dim:
        loc_tech_dim = [i for i in data_var.dims if "loc_carrier" in i]

    if not loc_tech_dim:
        if return_as == "Series":
            return data_var.to_series()
        elif return_as in ["DataArray", "MultiIndex DataArray"]:
            return data_var
        else:
            raise ValueError(
                "`return_as` must be `DataArray`, `Series`, or "
                "`MultiIndex DataArray`, but `{}` given".format(return_as)
            )

    elif len(loc_tech_dim) > 1:
        e = exceptions.ModelError
        raise e(
            "Cannot split loc_techs or loc_tech_carriers dimension "
            "for DataArray {}".format(data_var.name)
        )

    loc_tech_dim = loc_tech_dim[0]
    # xr.Datarray -> pd.Series allows for string operations
    data_var_idx = data_var[loc_tech_dim].to_index()
    index_list = data_var_idx.str.split("::").tolist()

    # carrier_prod, carrier_con, and carrier_export will return an index_list
    # of size 3, all others will be an index list of size 2
    possible_names = ["loc", "tech", "carrier"]
    names = [i + "s" for i in possible_names if i in loc_tech_dim]

    data_var_midx = pd.MultiIndex.from_tuples(index_list, names=names)

    # Replace the Datarray loc_tech_dim with this new MultiIndex
    updated_data_var = data_var.copy()
    updated_data_var.coords[loc_tech_dim] = data_var_midx

    if return_as == "MultiIndex DataArray":
        return updated_data_var

    elif return_as == "Series":
        return reorganise_xarray_dimensions(updated_data_var.unstack()).to_series()

    elif return_as == "DataArray":
        return reorganise_xarray_dimensions(updated_data_var.unstack())

    else:
        raise ValueError(
            "`return_as` must be `DataArray`, `Series`, or "
            "`MultiIndex DataArray`, but `{}` given".format(return_as)
        )


def best_effort_factors(data_ranges_per_unit):
    factors = {'const': 1} # cannot scale const

    unscaled = {k: v for k,v in data_ranges_per_unit.items() if v['num'] not in factors or v['den'] not in factors}
    
    #while len(unscaled):
     #   pass        
    
    return factors
    
    

def scale(data, transform=lambda x: x):

    # will hold max and min value information for each unit 
    data_ranges_per_unit = {}
    data_ranges_per_unit_scaled = {}
    def update_ranges(ranges, num, denom, val):
        tag = '{}/{}'.format(num, denom)
        valabs = xr.ufuncs.fabs(val) # absolute value
        valfin = valabs.where(xr.ufuncs.isfinite(valabs), drop=True) # finite values (drops nan and inf)
        valpos = valfin.where(valfin > 0.0, drop=True) # drop zeros
        if valpos.size > 0:
            minval = valpos.min().values
            maxval = valpos.max().values
            if tag in ranges:
                minval = min(minval, ranges[tag]["min"])
                maxval = max(maxval, ranges[tag]["max"])
            ranges[tag] = {"num": num, "den": denom, "min": minval, "max": maxval}
            

    # extract scaling factors from data for easier accessing
    factors = {data.scale.unit.values[i] : data.scale.values[i] for i in range(0, len(data.scale))}
    
    # TODO: need a better way to distinguish between costs and non-costs!
    for key, val in data.data_vars.items():
        if 'cost' in key.split('_'): # scale cost in all cost classes (also group constraints with costs)
            for i in range(0, len(data.costs)):
                costclass = data.costs[i].values.item(0)
                factor = get_cost_scaling_factor(key, costclass, factors)
                if not factor is None:
                    update_ranges(data_ranges_per_unit, factor[1], factor[2], data[key].loc[dict(costs=costclass)])
                    data[key].loc[dict(costs=costclass)] *= transform(factor[0])
                    update_ranges(data_ranges_per_unit_scaled, factor[1], factor[2], data[key].loc[dict(costs=costclass)])

        else: # scale constraint
            factor = get_scaling_factor(key, factors)
            if not factor is None:
                update_ranges(data_ranges_per_unit, factor[1], factor[2], data[key])
                data[key] = transform(factor[0]) * val
                update_ranges(data_ranges_per_unit_scaled, factor[1], factor[2], data[key])

    # scale all resources according to respective unit
    # resources need to be handled specially because they can have units in
    # {kWh, kWh/m2, kWh/kW}
    if 'loc_techs_finite_resource' in data:
        for res in data.loc_techs_finite_resource:
            resource_unit = data.resource_unit.sel(loc_techs_finite_resource=res).values.item(0)
            factor = 1
            if resource_unit in ['energy', 'energy_per_area']:
                # do some scaling
                num = 'power'
                if resource_unit == 'energy':
                    factor = factors.get('power', 1)
                    den = 'const'
                elif resource_unit == 'energy_per_area':
                    factor = factors.get('power', 1)/factors.get('area', 1)
                    den = 'area'
                update_ranges(data_ranges_per_unit, num, den, data["resource"].loc[dict(loc_techs_finite_resource=res)])
                data["resource"].loc[dict(loc_techs_finite_resource=res)] *= transform(factor)
                update_ranges(data_ranges_per_unit_scaled, num, den, data["resource"].loc[dict(loc_techs_finite_resource=res)])
            else:
                # don't need to scale kWh/kW = h
                assert(resource_unit == 'energy_per_power' and 'sanity check on my naming')

                
    # print data ranges for inspection purposes

    print('\nUnit ranges before scaling')
    print('{:20} {:12} {:12}'.format("unit", "min", "max"))
    for k,v in data_ranges_per_unit.items():
        print('{:20} {:<12.8f} {:<12.8f}'.format(k, v["min"], v["max"]))
        
    print('\nUnit ranges after scaling')
    print('{:20} {:12} {:12}'.format("unit", "min", "max"))    
    for k,v in data_ranges_per_unit_scaled.items():
        print('{:20} {:<12.8f} {:<12.8f}'.format(k, v["min"], v["max"]))

    factors = best_effort_factors(data_ranges_per_unit)
    print('\nSuggesting the following factors')
    for k,v in factors.items():
        print('{:12} {:<12.8f}'.format(k, v))
        

    print('\n')
    return data



def get_cost_scaling_factor(cost_name, cost_class, scaling_factors):
    per_power = [
        "cost_energy_cap",
        "cost_resource_cap",
        "cost_om_annual",
        "cost_export",
        "cost_om_con",
        "cost_om_prod",
        "cost_storage_cap"
    ]
    per_power_distance = [
        "cost_energy_cap_per_distance"
    ]
    per_area = [
        "cost_resource_area"
    ]
    cost = [
        "cost_purchase_unit",
        "cost",
        "cost_var",
        "cost_investment",
        "group_cost_max",
        "group_cost_min",
        "group_cost_equals",
        "group_cost_var_max",
        "group_cost_var_min",
        "group_cost_var_equals",
        "group_cost_investment_max",
        "group_cost_investment_min",
        "group_cost_investment_equals"
    ]
    per_cost = [
        "cost_om_annual_investment_fraction"
    ]

    fac = scaling_factors.get(cost_class, 1)
    if cost_name in per_power:
        return (
            fac/scaling_factors.get("power", 1),
            cost_class, "power"
        )
    elif cost_name in per_power_distance:
        return (
            fac/(scaling_factors.get("power", 1) * scaling_factors.get("distance", 1)),
            cost_class, "power_distance"
        )
    elif cost_name in per_area:
        return (
            fac/scaling_factors.get("area", 1),
            cost_class, "area"
        )
    elif cost_name in cost: 
        return (
            fac,
            cost_class, "const"
        )
    elif cost_name in per_cost:
        return (
            1/fac,
            "const", cost_class
        )
    else:
        print("returning none for ", cost_name, cost_class)
        return None
    

def get_scaling_factor(variable_name, scaling_factors):
    """
    - in order to scale time we also would need to modify the timeseries data

    # hour^-1
    if variable_name in ["charge_rate", "energy_cap_per_storage_cap_min", "energy_cap_per_storage_cap_max", "storage_loss"]:
        return units["per_hour"]
    # fraction/hour
    elif variable_name in ["energy_ramping"]:
        return units["fraction_per_hour"]
    # years
    elif variable_name in ["lifetime"]:
        return units["time"]

    - scaling integers, floats and fractions doesn't make sense 
    - is energy_cap_min/max/equals really kWh?
    """
    power = [
        "resource_cap_equals",
        "resource_cap_max",
        "resource_cap_min",
        "energy_cap_equals",
        "energy_cap_equals_systemwide",
        "energy_cap_max",
        "energy_cap_max_systemwide",
        "energy_cap_min",
        "energy_cap_per_unit",
        "export_cap",
        "units_max_systemwide",
        "units_equals_systemwide",
        "energy_cap",
        "resource_cap",
        "storage_cap_per_unit",
        "storage_cap_equals",
        "storage_cap_min",
        "storage_cap_max",
        "carrier_prod_min",
        "carrier_prod_max",
        "carrier_prod_equals",
        "carrier_con",
        "carrier_prod",
        "carrier_export",
        "storage_cap",
        "storage",
        "resource_con",
        "unmet_demand",
        "unused_supply",
    ]
    distance_inv = [
        "energy_eff_per_distance"
    ]
    distance = [
        "distance"
    ]
    area = [
        "resource_area_equals",
        "resource_area_max",
        "resource_area_min",
        "available_area",
        "resource_area"
    ]
    area_per_power = [
        "resource_area_per_energy_cap"
    ]
    non_scalable = [
        "units_min", # integer
        "units_equals", # integer
        "units_max", # integer
        "lifetime", # years
        "carrier_ratios", # fraction
        "parasitic_eff", # fraction
        "energy_eff", # fraction
        "energy_cap_min_use", # fraction
        "resource_min_use", # fraction
        "resource_eff", # fraction
        "resource_scale", # fraction
        "storage_discharge_depth", # fraction
        "storage_initial", # fraction
        "cost_depreciation_rate", # fraction
        "interest_rate", # fraction
        "energy_ramping", # fraction / hour
        "storage_loss", # fraction/hour
        "resource_cap_equals_energy_cap", # boolean
        "force_asynchronous_prod_con", # boolean
        "force_resource", # boolean
        "one_way", # boolean 
        "energy_con", # boolean
        "energy_prod", # boolean
        "energy_cap_scale", # float
        "energy_cap_per_storage_cap_min", # hour -1
        "energy_cap_per_storage_cap_max", # hour -1
        "energy_cap_per_storage_cap_equals", # hour -1
        "charge_rate", # hour -1
        "purchased", # boolean
        "units", #integer
        "operating_units", # integer
        "max_demand_timesteps", # timestamps
        "resource_unit", # N/A
        "objective_cost_class", # N/A
        "lookup_loc_techs_export",# N/A
        "lookup_loc_techs",# N/A
        "lookup_loc_techs_conversion_plus",# N/A
        "lookup_loc_techs_conversion",# N/A
        "lookup_remotes",# N/A
        "lookup_loc_carriers",# N/A
        "lookup_loc_techs_area",# N/A
        "lookup_primary_loc_tech_carriers_in",# N/A
        "lookup_primary_loc_tech_carriers_out",# N/A
        "resource_unit",# N/A
        "inheritance",# N/A
        "scale",# N/A
        "timestep_resolution",# N/A
        "timestep_weights",# N/A
        "export_carrier",# N/A
        "colors",# N/A
        "names",# N/A
        "loc_coordinates",# N/A
    ]

    #CONSTRAINTS
    # kW
    if variable_name in power:
        return (
            scaling_factors.get("power", 1),
            "power", "const"
        )
    elif variable_name in distance_inv:
        return (
            1/scaling_factors.get("distance", 1),
            "const", "distance"
        )
    elif variable_name in distance:
        return (
            scaling_factors.get("distance", 1),
            "distance", "const"
        )
    elif variable_name in area:
        return (
            scaling_factors.get("area", 1),
            "area", "const"
        )
    elif variable_name in area_per_power:
        return (
            scaling_factors.get("area", 1) / scaling_factors.get("power", 1),
            "area", "power"
        )
    elif variable_name in non_scalable:
        return None
    else:
        print("returning none for variable ", variable_name)
        return None
