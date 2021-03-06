"""
Copyright (C) 2013-2019 Calliope contributors listed in AUTHORS.
Licensed under the Apache 2.0 License (see LICENSE file).

"""

from calliope import exceptions
import xarray as xr
import pandas as pd

# includes for the scaling
import math
import pyomo.core as po
from pyomo.opt import SolverFactory 
import pyomo.environ


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



def bound_rule(model, num1, den1, v1, num2, den2, v2):
    def g(acc):
        return model.x[acc] if acc != 'const' else 0
    
    return g(num1) - g(den1) + v1 - g(num2) + g(den2) - v2 <= model.r

def lower_limits_rule(model, num1, den1, v):
    def g(acc):
        return model.x[acc] if acc != 'const' else 0
    
    return g(num1) - g(den1) >= v


'''
Solve an auxiliary LP to find optimal scaling factors for given unit ranges:
Given a set of units {u1, u2, ..., un}
and the range of absolute values of each unit 
{ u1 -> [l1, r1], ..., un -> [ln, rn]}
we want to find scaling factors f1, ..., fn that minimize
max(fi*ri)/min(fj*rj)

A difficulty arises because each unit is a fraction of two base units: ui = bj/bk
and we need to retain consistency between scaling factors. 
Thus we optimize the scaling factors Fj of base units bj and compute factors of "composite" units, i.e.

ui = bj/bk --> fi = Fj/Fk
'''
def lp_unit_factors(ranges, solver_name, solver_io, solver_tolerance):
    model = po.ConcreteModel()

    '''
    variables
    (a) We create one variable for the log of the scaling factor of each unit
    (b) And one variable r to be minimized 
    '''
    unitvars = [unit for unit in filter(
        lambda u: u != 'const',
        list(set(
            list(map(lambda r: r['num'], ranges)) +
            list(map(lambda r: r['den'], ranges))
        )))
    ]

    model.x = po.Var(unitvars, domain=po.Reals)
    model.r = po.Var(domain=po.Reals)

    '''
    set objective to minimize r
    '''
    model.cost = po.Objective(expr = model.r)
    
    '''
    Constraints:
    limit si - sj + (u_hi){i/j} - sk + sl - (u_lo){k/l} <= r 
    for all pairs of units u{i/j}, u{k/l}
    to find max_{{i/j}, {k/l} \in units}(si*sl/sj*sk (u_hi){i/j} / (u_lo){k/l})
    '''
    ranges_wo_const = [rng for rng in ranges if rng['num'] != 'const' or rng['den'] != 'const']
    boundvals = [
        (r1['num'], r1['den'], math.log(r1['max'], 10), r2['num'], r2['den'], math.log(r2['min'], 10))
        for r1 in ranges for r2 in ranges
    ]
    
    model.bounds = po.Constraint(boundvals, rule=bound_rule)

    solver = SolverFactory(solver_name, solver_io=solver_io)
    solver.solve(model)
    temp_facs = {k: 10**model.x[k]() for k in unitvars}
    maxs = [r['max'] * temp_facs.get(r['num'], 1) / temp_facs.get(r['den'], 1) for r in ranges]
    mins = [r['min'] * temp_facs.get(r['num'], 1) / temp_facs.get(r['den'], 1) for r in ranges]
    best_range = max(maxs)/min(mins)
    print('best range: {}'.format(best_range))
    limit = math.ceil(math.log10(best_range))/2

    '''
    ensure that absolute values in model are not scaled below a certain threshold.
    this generally limits the objective function -> need to find good tradeoff.

    Practical Guidelines for Solving Difficult Linear Programs suggests we should ensure that user input consists of values larger than the solver tolerances tol, thus limit values to scaling_tolerance_threshold*tol. Or center the values around 0 if this gives an even higher threshold.

    note: our coefficient rounding below may lead to values 2 times smaller than the set limit, thus limit by 2*scaling_tolerance_threshold*tol 
    '''
    lower_limit = max(10*solver_tolerance, 10*10**(-limit))
    print('setting limit {}'.format(lower_limit))

    '''
    limit_vals = [
        (r['num'], r['den'], math.log(lower_limit/r['min'], 2))
        for r in ranges
    ]
    model.lower_limits = po.Constraint(limit_vals, rule=lower_limits_rule)
    '''

    violating_mins = set()
    while True:
        changed = False
        factors = {k: 10**model.x[k]() for k in unitvars}
        for rng in ranges_wo_const: 
            num = factors[rng['num']] if rng['num'] != 'const' else 1
            den = factors[rng['den']] if rng['den'] != 'const' else 1
            if rng['min']*num/den < lower_limit/10:
                print('violating {}/{}: {} -> add constr'.format(rng['num'], rng['den'], rng['min']*num/den))
                violating_mins.add((rng['num'], rng['den'], math.log(lower_limit/rng['min'], 10)))
                changed = True
            
        if not changed:
            break
                   
        if not model.component('lower_limits') is None:
            model.del_component(model.lower_limits)
            model.del_component(model.lower_limits_index)
                   
        model.add_component('lower_limits', po.Constraint(list(violating_mins), rule=lower_limits_rule))

        solver.solve(model)    

    '''
    we want all factors to be an exponent of 2 in order not to tamper with precision of user values (c.f. tomlin - on scaling linear programming problems)
    we achieve this by rounding the optimal values we just computed to integers before exponentiating them
    '''
    facs = {k: 10**math.floor(model.x[k]()) for k in unitvars}
    return facs



    


def lp_unit_factors_balanced(ranges, solver_name, solver_io, center, r_bound):
    

    def g(acc):
        return model.x[acc] if acc != 'const' else 0

    def bound_rule_lower(model, num, den, v, center):        
        return g(num) - g(den) + v >= center-model.r

    def bound_rule_upper(model, num, den, v, center):
        return g(num) - g(den) + v <= center+model.r

    # anti rules are for experimentation purposes only
    def anti_lower(model, num, den, v, center):
        return g(num) - g(den) + v <= -r_bound    

    def anti_upper(model, num, den, v, center):
        return g(num) - g(den) + v >= r_bound

    model = po.ConcreteModel()

    '''
    variables
    (a) We create one variable for the log of the scaling factor of each unit
    (b) And one variable r to be minimized 
    '''
    unitvars = [unit for unit in filter(
        lambda u: u != 'const',
        list(set(
            list(map(lambda r: r['num'], ranges)) +
            list(map(lambda r: r['den'], ranges))
        )))
    ]

    model.x = po.Var(unitvars, domain=po.Reals)
    model.r = po.Var(domain=po.Reals)

    '''
    set objective to minimize r
    '''
    model.cost = po.Objective(expr = model.r)

    '''
    Constraints:
    limit si - sj + (u_hi){i/j} <= center + r
          si - sj + (u_lo){i/j} >= center - r
    for all units u{i/j}. Then u scaled with si/sj will be between 2^{center-r} and 2^{center+r} 
    '''    
    ranges_wo_const = [rng for rng in ranges if rng['num'] != 'const' or rng['den'] != 'const']

    boundvals_lower = [
        (r['num'], r['den'], math.log(r['min'], 2), center)
        for r in ranges_wo_const 
    ]

    boundvals_upper = [
        (r['num'], r['den'], math.log(r['max'], 2), center)
         for r in ranges_wo_const
    ]
    
    model.bounds_lower = po.Constraint(boundvals_lower, rule=bound_rule_lower)
    model.bounds_upper = po.Constraint(boundvals_upper, rule=bound_rule_upper)

    model.anti_bounds_lower = po.Constraint(boundvals_lower, rule=anti_lower)
    model.anti_bounds_upper = po.Constraint(boundvals_upper, rule=anti_upper)

    #model.pprint()
    print('Scaling into [2^({}-r), 2^({}+r)]\nr >= {}'.format(
        center,center,r_bound
    ))
    
    solver = SolverFactory(solver_name, solver_io=solver_io)
    solver.solve(model)

    print('Achieved r = {}'.format(model.cost()))

    '''
    we want all factors to be an exponent of 2 in order not to tamper with precision of user values (c.f. tomlin - on scaling linear programming problems)
    we achieve this by rounding the optimal values we just computed to integers before exponentiating them
    '''
    facs = {k: 2**math.floor(model.x[k]()) for k in unitvars}
    return facs



'''
Iterate over all the nested data
Retrieve the unit of each variable in the data using 'get_unit' function
And return a dictionary mapping each unit to its min and max abs value
'''
def compute_unit_ranges(data):
    # will hold max and min value information for each unit 
    data_ranges_per_unit = {}
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

    # TODO: need a better way to distinguish between costs and non-costs!
    for key, val in data.data_vars.items():
        if 'cost' in key.split('_'): 
            for i in range(0, len(data.costs)):
                costclass = data.costs[i].values.item(0)
                unit = get_unit(key, costclass)
                if not unit is None:
                    update_ranges(data_ranges_per_unit, unit[0], unit[1], val.loc[dict(costs=costclass)])

        else: 
            unit = get_unit(key)
            if not unit is None:
                update_ranges(data_ranges_per_unit, unit[0], unit[1], val)

    # special case of resources. here we need to use the mapping resource_unit to get the unit of a resource
    if 'loc_techs_finite_resource' in data:
        for res in data.loc_techs_finite_resource:
            resource_unit = data.resource_unit.sel(loc_techs_finite_resource=res).values.item(0)
            elems = data["resource"].loc[dict(loc_techs_finite_resource=res)]
            if resource_unit == 'energy':
                update_ranges(data_ranges_per_unit, 'power', 'const', elems)
            elif resource_unit == 'energy_per_area':
                update_ranges(data_ranges_per_unit, 'power', 'area', elems)
            elif resource_unit == 'energy_per_cap':
                # this has unit hours, which we don't report!
                # new: we report this as const because we cannot scale it and stay consistent but we still want to consider it in optimization
                update_ranges(data_ranges_per_unit, 'const', 'const', elems)
            else:
                assert(False and 'there shouldnt be a resource of this type')

    return data_ranges_per_unit
    


'''
compute an auxiliary LP to get optimal scaling factors given the ranges of each unit and an lp solver
'''
def get_scale(data_ranges_per_unit, solver, solver_io, tolerance, center, r_bound):
    factors = lp_unit_factors_balanced(list(data_ranges_per_unit.values()), solver, solver_io, center, r_bound)
    return factors

'''
apply a scaling to the data
the scaling factors are assumed to be stored as part of the data at data.scale
'''
def scale(data, transform=lambda x: x):
    # extract scaling factors from data for easier accessing
    factors = {data.scale.unit.values[i] : data.scale.values[i] for i in range(0, len(data.scale))}
    factors['const'] = 1
    
    # TODO: need a better way to distinguish between costs and non-costs!
    for key, val in data.data_vars.items():
        if 'cost' in key.split('_'): # scale cost in all cost classes (also group constraints with costs)
            for i in range(0, len(data.costs)):
                costclass = data.costs[i].values.item(0)
                factor = get_scaling_factor(factors, key, costclass)
                if not factor is None:
                    data[key].loc[dict(costs=costclass)] *= transform(factor)

        else: # scale constraint
            factor = get_scaling_factor(factors, key)
            if not factor is None:
                data[key] = transform(factor) * val

    # scale all resources according to respective unit
    # resources need to be handled specially because they can have units in
    # {kWh, kWh/m2, kWh/kW}
    if 'loc_techs_finite_resource' in data:
        for res in data.loc_techs_finite_resource:
            resource_unit = data.resource_unit.sel(loc_techs_finite_resource=res).values.item(0)
            factor = 1
            if resource_unit in ['energy', 'energy_per_area', 'energy_per_cap']:
                # do some scaling
                if resource_unit == 'energy':
                    factor = factors.get('power', 1)
                elif resource_unit == 'energy_per_area':
                    factor = factors.get('power', 1)/factors.get('area', 1)
                else:
                    #factor = factors.get('energy', 1)/factors.get('power', 1)
                    pass
                data["resource"].loc[dict(loc_techs_finite_resource=res)] *= transform(factor)
            else:
                assert(False and 'sanity check on my naming failed')

    print('done')
    return data

    


units_to_names = {
    'power': [
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
        "group_carrier_prod_min", 
        "group_carrier_prod_max", 
        "group_carrier_prod_equals",
        "group_energy_cap_min", 
        "group_energy_cap_max", 
        "group_energy_cap_equals", 
    ],

    'energy': [
        "storage_cap_equals",
        "storage_cap_min",
        "storage_cap_max",
        "storage_cap_per_unit",        
    ],
    
    'distance_inv': [
        "energy_eff_per_distance"
    ],

    'distance': [
        "distance"
    ],

    'area': [
        "resource_area_equals",
        "resource_area_max",
        "resource_area_min",
        "available_area",
        "resource_area",
        "group_resource_area_min",
        "group_resource_area_max",
        "group_resource_area_equals",
    ],

    'area_per_power': [
        "resource_area_per_energy_cap"
    ],

    'power_per_energy': [
        "energy_ramping", # fraction / hour
        "storage_loss", # fraction/hour
        "energy_cap_per_storage_cap_min", # hour -1
        "energy_cap_per_storage_cap_max", # hour -1
        "energy_cap_per_storage_cap_equals", # hour -1
        "charge_rate", # hour -1
    ],
    
    'non_scalable': [
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
        "energy_cap_scale", # float
        "units", #integer
        "operating_units", # integer
        "group_demand_share_min", # fraction
        "group_demand_share_max", # fraction
        "group_demand_share_equals", # fraction
        "group_demand_share_per_timestep_min", # fraction
        "group_demand_share_per_timestep_max", # fraction
        "group_demand_share_per_timestep_equals", # fraction
        "group_demand_share_per_timestep_decision", # fraction
        "group_carrier_prod_share_min", # fraction
        "group_carrier_prod_share_max", # fraction
        "group_carrier_prod_share_equals", # fraction
        "group_carrier_prod_share_per_timestep_min", # fraction
        "group_carrier_prod_share_per_timestep_max", # fraction
        "group_carrier_prod_share_per_timestep_equals", # fraction
        "group_net_import_share_min", # fraction
        "group_net_import_share_max", # fraction
        "group_net_import_share_equals", # fraction
        "group_energy_cap_share_min", # fraction
        "group_energy_cap_share_max", # fraction
        "group_energy_cap_share_equals", # fraction
    ],

    'non_numeric': [
        "resource_cap_equals_energy_cap", # boolean
        "force_asynchronous_prod_con", # boolean
        "force_resource", # boolean
        "one_way", # boolean 
        "energy_con", # boolean
        "energy_prod", # boolean
        "purchased", # boolean
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
    ],
    
    'cost_per_power': [
        "cost_energy_cap",
        "cost_resource_cap",
        "cost_om_annual",
    ],

    'cost_per_energy': [
        "cost_export",
        "cost_om_con",
        "cost_om_prod",
        "cost_storage_cap"        
    ],
    
    'cost_per_power_distance': [
        "cost_energy_cap_per_distance"
    ],

    'cost_per_area': [
        "cost_resource_area"
    ],

    'cost': [
        "cost_purchase_unit",
        "cost_purchase",
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
        "group_cost_max",
        "group_cost_min", 
        "group_cost_equals", 
        "group_cost_var_max",
        "group_cost_var_min",
        "group_cost_var_equals",
        "group_cost_investment_max",
        "group_cost_investment_min",
        "group_cost_investment_equals",
    ],

    'per_cost': [
        "cost_om_annual_investment_fraction"
    ]
}


def get_unit(variable_name, cost_class=''):    
    if variable_name in units_to_names['power']:
        return ("power", "const")
    elif variable_name in units_to_names['energy']:
        return ("power", "const") #return ("energy", "const")
    elif variable_name in units_to_names['distance_inv']:
        return ("const", "distance")
    elif variable_name in units_to_names['distance']:
        return ("distance", "const")
    elif variable_name in units_to_names['area']:
        return ("area", "const")
    elif variable_name in units_to_names['area_per_power']:
        return ("area", "power")
    elif variable_name in units_to_names['power_per_energy']:
        return ("const", "const") # return ("power", "energy")
    elif variable_name in units_to_names['cost_per_power']:
        return (cost_class, "power")
    elif variable_name in units_to_names['cost_per_energy']:
        return (cost_class, "power") # return (cost_class, "energy")
    elif variable_name in units_to_names['cost_per_power_distance']:
        return (cost_class, "power_distance")
    elif variable_name in units_to_names['cost_per_area']:
        return (cost_class, "area")
    elif variable_name in units_to_names['cost']: 
        return (cost_class, "const")
    elif variable_name in units_to_names['per_cost']:
        return ("const", cost_class)
    elif variable_name in units_to_names['non_scalable']:
        return ("const", "const")
    elif variable_name in units_to_names['non_numeric']:
        return None
    else:
        # this is intended for "resource" but for any other variables this shouldn't happen
        # if this happens for a variable other than resource, this variable needs to be added to the units_to_names dict
        print("returning none for variable ", variable_name)
        return None


    
def get_scaling_factor(scaling_factors, variable_name, cost_class=''):

    unit = get_unit(variable_name, cost_class)
    if unit is None:
        return None
    else:
        # FIXME: this assert will fail for cost_energy_cap_per_distance but so far we haven't encountered this
        # the simple fix is treating it as a special case in the scaling
        # but with the current implementation we cannot incorporate this variable in the optimization problem because the optimization assumes all units are a fraction of exactly 2 basic units u1/u2...
        # cost_energy_cap_per_distance is cost/(energy*distance)...  
        assert(unit[0] in scaling_factors and unit[1] in scaling_factors and 'wot the heck')
        return scaling_factors[unit[0]]/scaling_factors[unit[1]]
