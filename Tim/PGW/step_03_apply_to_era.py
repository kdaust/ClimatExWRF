#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
description     PGW for ERA5 main routine to update ERA5 files with climate
                deltas to transition from ERA climate to PGW climate.
authors		Before 2022:    original developments by Roman Brogli
                Since 2022:     upgrade to PGW for ERA5 by Christoph Heim 
                2022:           udpates by Jonas Mensch
"""
##############################################################################
import argparse, os
import xarray as xr
import numpy as np
from argparse import RawDescriptionHelpFormatter
from pathlib import Path
from datetime import datetime, timedelta
from functions import (
    specific_to_relative_humidity,
    relative_to_specific_humidity,
    load_delta,
    load_delta_interp,
    integ_geopot,
    interp_logp_4d,
    determine_p_ref,
    integrate_tos
    )
from constants import CON_G, CON_RD
from parallel import IterMP
from settings import (
    i_debug,
    era5_file_name_base,
    var_name_map,
    TIME_ERA, LEV_ERA, HLEV_ERA, LON_ERA, LAT_ERA, SOIL_HLEV_ERA,
    TIME_GCM, PLEV_GCM,
    i_reinterp,
    p_ref_inp,
    thresh_phi_ref_max_error,
    max_n_iter,
    adj_factor,
    file_name_bases,
    )
##############################################################################

def pgw_for_era5(inp_era_file_path, out_era_file_path,
                delta_input_dir, era_step_dt,
                ignore_top_pressure_error,
                debug_mode=None):
    if i_debug >= 0:
        print('Start working on input file {}'.format(inp_era_file_path))

    #########################################################################
    ### PREPARATION STEPS
    #########################################################################
    # containers for variable computation
    vars_era = {}
    vars_pgw = {}
    deltas = {}

    # open data set
    era_file = xr.open_dataset(inp_era_file_path, decode_cf=False)

    ## compute pressure on ERA5 full levels and half levels
    # pressure on half levels
    pa_hl_era = (era_file.ak + 
                era_file[var_name_map['ps']] * era_file.bk).transpose(
                TIME_ERA, HLEV_ERA, LAT_ERA, LON_ERA)
    # if akm and akb coefficients (for full levels) exist, use them
    if 'akm' in era_file:
        akm = era_file.akm
        bkm = era_file.bkm
    # if akm and abk  coefficients do not exist, computed them
    # with the average of the half-level coefficients above and below
    else:
        akm = (
            0.5 * era_file.ak.diff(
            dim=HLEV_ERA, 
            label='lower').rename({HLEV_ERA:LEV_ERA}) + 
            era_file.ak.isel({HLEV_ERA:range(len(era_file.level1)-1)}).values
        )
        bkm = (
            0.5 * era_file.bk.diff(
            dim=HLEV_ERA, 
            label='lower').rename({HLEV_ERA:LEV_ERA}) + 
            era_file.bk.isel({HLEV_ERA:range(len(era_file.level1)-1)}).values
        )
    # pressure on full levels
    pa_era = (akm + era_file[var_name_map['ps']] * bkm).transpose(
                TIME_ERA, LEV_ERA, LAT_ERA, LON_ERA)

    # compute relative humidity in ERA climate state
    era_file[var_name_map['hur']] = specific_to_relative_humidity(
                        era_file[var_name_map['hus']], 
                        pa_era, era_file[var_name_map['ta']]).transpose(
                        TIME_ERA, LEV_ERA, LAT_ERA, LON_ERA)

    #########################################################################
    ### UPDATE SURFACE AND SOIL TEMPERATURE
    #########################################################################
    # update surface skin temperature using SST delta over sea and and 
    # surface skin temperature delta over land and over sea ice
    if i_debug >= 2:
        print('update surface skin temperature (ts)')
    delta_siconc = load_delta(delta_input_dir, 'siconc',
                        era_file[TIME_ERA], era_step_dt) 
    era_file[var_name_map['sic']].values += delta_siconc.values/100
    era_file[var_name_map['sic']].values = np.clip(
                        era_file[var_name_map['sic']].values, 0, 1)
    print(np.nanmin(era_file[var_name_map['sic']].values))
    #deltas['siconc'] = delta_siconc
    # load surface temperature climate delta
    #(for grid points over land and sea ice)
    delta_ts = load_delta(delta_input_dir, 'ts',
                           era_file[TIME_ERA], era_step_dt)
    # load SST climate delta (for grid points over open water)
    delta_tos = load_delta(delta_input_dir, 'tos',
                           era_file[TIME_ERA], era_step_dt)
    # combine using land and sea-ice mask in ERA5
    delta_ts_combined = integrate_tos(
        delta_tos.values,
        delta_ts.values, 
        era_file[var_name_map['sftlf']].isel({TIME_ERA:0}).values, 
        era_file[var_name_map['sic']].isel({TIME_ERA:0}).values
    )
    era_file[var_name_map['ts']].values += delta_ts_combined
    delta_ts.values = delta_ts_combined
    # store delta for output in case of --debug_mode = interpolate_full
    deltas['ts'] = delta_ts
    print(np.nanmin(era_file[var_name_map['sic']].values))
    # update temperature of soil layers
    if i_debug >= 2:
        print('update soil layer temperature (st)')
    # set climatological lower soil temperature delta to annual mean
    # climate delta of surface skin temperature.
    delta_st_clim = load_delta(delta_input_dir, 'ts',
                            era_file[TIME_ERA], 
                            target_date_time=None).mean(dim=[TIME_GCM])
    # interpolate between surface temperature and deep soil temperature
    # using exponential decay of annual cycle signal
    delta_soilt = (
            delta_st_clim + np.exp(-era_file.soil1/2.8) * 
                    (delta_ts - delta_st_clim)
    )
    delta_soilt = delta_soilt.transpose(TIME_ERA, SOIL_HLEV_ERA, LAT_ERA, LON_ERA)
    era_file[var_name_map['st']].values += delta_soilt
    # store delta for output in case of --debug_mode = interpolate_full
    deltas['st'] = delta_soilt

    #########################################################################
    ### START UPDATING 3D FIELDS
    #########################################################################
    # If no re-interpolation is done, the final PGW climate state
    # variables can be computed already now, before updating the 
    # surface pressure. This means that the climate deltas or
    # interpolated on the ERA5 model levels of the ERA climate state.
    if not i_reinterp:

        ### interpolate climate deltas onto ERA5 grid
        for var_name in ['ta','hur','ua','va']:
            if i_debug >= 2:
                print('update {}'.format(var_name))

            ## interpolate climate deltas to ERA5 model levels
            ## use ERA climate state
            delta_var = load_delta_interp(delta_input_dir,
                    var_name, pa_era, era_file[TIME_ERA], era_step_dt,
                    ignore_top_pressure_error)
            deltas[var_name] = delta_var

            ## compute PGW climate state variables
            vars_pgw[var_name] = (
                    era_file[var_name_map[var_name]] + 
                    deltas[var_name]
            )


    #########################################################################
    ### UPDATE SURFACE PRESSURE USING ITERATIVE PROCEDURE
    #########################################################################
    if i_debug >= 2:
        print('###### Start with iterative surface pressure adjustment.')
    # change in surface pressure between ERA and PGW climate states
    delta_ps = xr.zeros_like(era_file[var_name_map['ps']])
    # increment to adjust delta_ps with each iteration
    adj_ps = xr.zeros_like(era_file[var_name_map['ps']])
    # maximum error in geopotential (used in iteration)
    phi_ref_max_error = np.inf

    it = 1
    while phi_ref_max_error > thresh_phi_ref_max_error:

        # update surface pressure
        delta_ps += adj_ps
        ps_pgw = era_file[var_name_map['ps']] + delta_ps

        # recompute pressure on full and half levels
        pa_pgw = (akm + ps_pgw * bkm).transpose(
                    TIME_ERA, LEV_ERA, LAT_ERA, LON_ERA)
        pa_hl_pgw = (era_file.ak + ps_pgw * era_file.bk).transpose(
                    TIME_ERA, HLEV_ERA, LAT_ERA, LON_ERA)


        if i_reinterp:
            # interpolate ERA climate state variables as well as
            # climate deltas onto updated model levels, and
            # compute PGW climate state variables
            if i_debug >= 2:
                print('reinterpolate ta and hur')
            for var_name in ['ta', 'hur']:
                vars_era[var_name] = interp_logp_4d(
                                era_file[var_name_map[var_name]], 
                                pa_era, pa_pgw, extrapolate='constant')
                deltas[var_name] = load_delta_interp(delta_input_dir,
                                                var_name, pa_pgw,
                                                era_file[TIME_ERA], era_step_dt,
                                                ignore_top_pressure_error)
                vars_pgw[var_name] = vars_era[var_name] + deltas[var_name]

        # Determine current reference pressure (p_ref)
        if p_ref_inp is None:
            # get GCM pressure levels as candidates for reference pressure
            p_ref_opts = load_delta(delta_input_dir, 'zg',
                                era_file[TIME_ERA], era_step_dt)[PLEV_GCM]
            # maximum reference pressure in ERA and PGW climate states
            # (take 95% of surface pressure to ensure that a few model
            # levels are located in between which makes the solution
            # smoother).
            p_min_era = pa_hl_era.isel(
                        {HLEV_ERA:len(pa_hl_era[HLEV_ERA])-1}) * 0.95
            p_min_pgw = pa_hl_pgw.isel(
                        {HLEV_ERA:len(pa_hl_era[HLEV_ERA])-1}) * 0.95
            # reference pressure from a former iteration already set?
            try:
                p_ref_last = p_ref
            except UnboundLocalError:
                p_ref_last = None
            # determine local reference pressure
            p_ref = xr.apply_ufunc(determine_p_ref, p_min_era, p_min_pgw, 
                    p_ref_opts, p_ref_last,
                    input_core_dims=[[],[],[PLEV_GCM],[]],
                    vectorize=True)
            if HLEV_ERA in p_ref.coords:
                del p_ref[HLEV_ERA]
            # make sure a reference pressure above the required model
            # level could be found everywhere
            if np.any(np.isnan(p_ref)):
                raise ValueError('No reference pressure level above the ' +
                        'required local minimum pressure level could not ' +
                        'be found everywhere. ' +
                        'This is likely the case because your geopotential ' +
                        'data set does not reach up high enough (e.g. only to ' +
                        '500 hPa instead of e.g. 300 hPa?)')
        else:
            p_ref = p_ref_inp

        #p_sfc_era.to_netcdf('psfc_era.nc')
        #p_ref.to_netcdf('pref.nc')
        #quit()

        # convert relative humidity to speicifc humidity in pgw
        # take PGW climate state temperature and relative humidity
        # and pressure of current iteration
        vars_pgw['hus'] = relative_to_specific_humidity(
            vars_pgw['hur'], 
            pa_pgw, 
            vars_pgw['ta']
        )

        # compute updated geopotential at reference pressure
        phi_ref_pgw = integ_geopot(
            pa_hl_pgw, 
            era_file[var_name_map['zgs']], 
            vars_pgw['ta'], 
            vars_pgw['hus'], 
            era_file[HLEV_ERA], 
            p_ref
        )

        # recompute original geopotential at currently used 
        # reference pressure level
        phi_ref_era = integ_geopot(
            pa_hl_era, 
            era_file[var_name_map['zgs']],
            era_file[var_name_map['ta']], 
            era_file[var_name_map['hus']], 
            era_file[HLEV_ERA],
            p_ref
        )

        delta_phi_ref = phi_ref_pgw - phi_ref_era

        ## load climate delta at currently used reference pressure level
        climate_delta_phi_ref = load_delta(delta_input_dir, 'zg',
                            era_file[TIME_ERA], era_step_dt) * CON_G
        climate_delta_phi_ref = climate_delta_phi_ref.sel({PLEV_GCM:p_ref})
        del climate_delta_phi_ref[PLEV_GCM]

        # error in future geopotential
        phi_ref_error = delta_phi_ref - climate_delta_phi_ref

        # adjust surface pressure by some amount in the right direction
        adj_ps = - adj_factor * ps_pgw / (
                CON_RD * 
                vars_pgw['ta'].sel({LEV_ERA:np.max(era_file[LEV_ERA])})
            ) * phi_ref_error
        if LEV_ERA in adj_ps.coords:
            del adj_ps[LEV_ERA]

        phi_ref_max_error = np.abs(phi_ref_error).max().values
        if i_debug >= 2:
            print('### iteration {:03d}, phi max error: {}'.
                            format(it, phi_ref_max_error))

        it += 1

        if it > max_n_iter:
            raise ValueError('ERROR! Pressure adjustment did not converge '+
                  'for file {}. '.format(inp_era_file_path) +
                  'Consider increasing the value for "max_n_iter" in ' +
                  'settings.py')

    #########################################################################
    ### FINISH UPDATING 3D FIELDS
    #########################################################################
    # store computed delta ps for output in case of 
    # --debug_mode = interpolate_full
    deltas['ps'] = ps_pgw - era_file.PS

    ## If re-interpolation is enabled, interpolate climate deltas for
    ## ua and va onto final PGW climate state ERA5 model levels.
    if i_reinterp:
        for var_name in ['ua', 'va']:
            if i_debug >= 2:
                print('add {}'.format(var_name))
            var_era = interp_logp_4d(era_file[var_name_map[var_name]], 
                            pa_era, pa_pgw, extrapolate='constant')
            delta_var = load_delta_interp(delta_input_dir,
                    var_name, pa_pgw,
                    era_file[TIME_ERA], era_step_dt,
                    ignore_top_pressure_error)
            vars_pgw[var_name] = var_era + delta_var
            # store delta for output in case of 
            # --debug_mode = interpolate_full
            deltas[var_name] = delta_var

    #########################################################################
    ### DEBUG MODE
    #########################################################################
    ## for debug_mode == interpolate_full, write final climate deltas
    ## to output directory
    if debug_mode == 'interpolate_full':
        var_names = ['ps','ta','hur','ua','va','st','ts']
        for var_name in var_names:
            print(var_name)
            # creat output file name
            out_file_path = os.path.join(Path(out_era_file_path).parents[0],
                                '{}_delta_{}'.format(var_name_map[var_name], 
                                            Path(out_era_file_path).name))
            # convert to dataset
            delta = deltas[var_name].to_dataset(name=var_name_map[var_name])
            # save climate delta
            delta.to_netcdf(out_file_path, mode='w')

    #########################################################################
    ### SAVE PGW ERA5 FILE
    #########################################################################
    ## for production mode, modify ERA5 file and save
    else:
        ## update 3D fields in ERA file
        era_file[var_name_map['ps']] = ps_pgw
        for var_name in ['ta','hus','ua','va']:
            era_file[var_name_map[var_name]] = vars_pgw[var_name]
        ## remove manually computed RH field in ERA5 file
        del era_file[var_name_map['hur']]


        ## save updated ERA5 file
        print(np.nanmin(era_file[var_name_map['sic']].values))
        era_file.to_netcdf(out_era_file_path, mode='w')
        era_file.close()
        if i_debug >= 1:
            print('Done. Saved to file {}.'.format(out_era_file_path))



##############################################################################

def debug_interpolate_time(
                inp_era_file_path, out_era_file_path,
                delta_input_dir, era_step_dt,
                ignore_top_pressure_error,
                debug_mode=None):
    """
    Debugging function to test time interpolation. Is called if input
    inputg argument --debug_mode is set to "interpolate_time".
    """
    # load input ERA5 file
    # in this debugging function, the only purpose of this is to obtain 
    # the time format of the ERA5 file
    era_file = xr.open_dataset(inp_era_file_path, decode_cf=False)

    var_names = ['tos','tas','hurs','ps','ta','hur','ua','va','zg']
    for var_name in var_names:
        print(var_name)
        # creat output file name
        out_file_path = os.path.join(Path(out_era_file_path).parents[0],
                                    '{}_{}_{}'.format("delta",var_name, 
                                        Path(out_era_file_path).name))
        # load climate delta interpolated in time only
        delta = load_delta(delta_input_dir, var_name, era_file[TIME_ERA], 
                       target_date_time=era_step_dt)
        # convert to dataset
        delta = delta.to_dataset(name=var_name)
        delta.to_netcdf(out_file_path, mode='w')
    era_file.close()






##############################################################################
if __name__ == "__main__":
    ## input arguments
    parser = argparse.ArgumentParser(description =
    """
    Perturb ERA5 with PGW climate deltas. Settings can be made in
    "settings.py".
    ##########################################################################

    Main function to update ERA5 files with the PGW signal.
    The terminology used is HIST referring to the historical (or reference)
    climatology, SCEN referring to the future (climate change scenario)
    climatology, and SCEN-HIST (a.k.a. climate delta) referring to the
    PGW signal which should be applied to the ERA5 files.
    The script computes and adds a climate change signal for:
        - ua
        - va
        - ta (including tas for interpolation near the surface)
        - hus (computed using a hur and hurs climate delta)
        - surface skin temperature (including SST) and soil temperature
    and consequently iteratively updates ps to maintain hydrostatic
    balance. During this, the climate delta for zg is additionally required.
    A list of all climate deltas required is shown in settings.py.

    ##########################################################################

    If the variable names in the ERA5 files to be processed deviate from
    the CMOR convention, the dict 'var_name_map' in the file 
    settings.py allows to map between the CMOR names and the names in the ERA5
    file. Also the coordinate names in the ERA5 or the GCM climate
    delta files can be changed in settings.py, if required.

    ##########################################################################

    The code can be run in parallel on multiple ERA5 files at the same time.
    See input arguments.

    ##########################################################################

    Note that epxloring the option --debug_mode (-D) can provide a lot of
    insight into what the code does and can help gain confidence using the
    code (see argument documentation below for more information).

    ##########################################################################

    Some more information about the iterative surface pressure
    adjustment:

    - The procedure requires a reference pressure level (e.g. 500 hPa) for
    which the geopotential is computed. Based on the deviation between the
    computed and the GCM reference pressure geopotential, the surface pressure
    is adjusted. Since the climate deltas may not always be available at 
    native vertical GCM resolution, but the climate delta for the geopotential
    on one specific pressure level itself is computed by the GCM using data
    from all GCM model levels, this introduces an error in the surface
    pressure adjustment used here. See publication for more details.
    The higher (in terms of altitdue) the reference pressure is chosen, 
    the larger this error may get. 
    Alternatively, the reference pressure can be determined locally 
    as the lowest possible pressure above the surface for which a climate 
    delta for the geopotential is available (see settings.py).

    - If the iteration does not converge, 'thresh_phi_ref_max_error' in
    the file settings.py may have to be raised a little bit. Setting
    i_debug = 2 may help to diagnose if this helps.

    - As a default option, the climate deltas are interpolated to
    the ERA5 model levels of the ERA climate state before the surface
    pressure is adjusted (i_reinterp = 0).
    There is an option implemented (i_reinterp = 1) in which the
    deltas are re-interpolated on the updated ERA5 model levels
    with each iteration of surface pressure adjustment. However, this
    implies that the ERA5 fields are extrapolated at the surface
    (if the surface pressure increases) the effect of which was not
    tested in detail. The extrapolation is done assuming that the
    boundary values are constant, which is not ideal for height-dependent
    variables like e.g. temperature. As a default, it is recommended to set
    i_reinterp = 0.

    ##########################################################################

    """, formatter_class=RawDescriptionHelpFormatter)

    # input era5 directory
    parser.add_argument('-i', '--input_dir', type=str, default=None,
            help='Directory with ERA5 input files to process. ' +
                 'These files are not overwritten but copies will ' +
                 'be save in --output_dir .')

    # output era5 directory
    parser.add_argument('-o', '--output_dir', type=str, default=None,
            help='Directory to store processed ERA5 files.')

    # first bc step to compute 
    parser.add_argument('-f', '--first_era_step', type=str,
            default='2006080200',
            help='Date of first ERA5 time step to process. Format should ' +
            'be YYYYMMDDHH.')

    # last bc step to compute 
    parser.add_argument('-l', '--last_era_step', type=str,
            default='2006080300',
            help='Date of last ERA5 time step to process. Format should ' +
            'be YYYYMMDDHH.')

    # delta hour increments
    parser.add_argument('-H', '--hour_inc_step', type=int, default=3,
            help='Hourly increment of the ERA5 time steps to process '+
            'between --first_era_step and --last_era_step. Default value ' +
            'is 3-hourly, i.e. (00, 03, 06, 09, 12, 15, 18, 21 UTC).')

    # climate delta directory (already remapped to ERA5 grid)
    parser.add_argument('-d', '--delta_input_dir', type=str, default=None,
            help='Directory with GCM climate deltas (SCEN-HIST) to be used. ' +
            'This directory should have a climate delta for ta,hur,' +
            'ua,va,zg,tas,hurs,ts,tos (e.g. ta_delta.nc), as well as the ' +
            'HIST climatology value for ps (e.g. ps_historical.nc). ' +
            'All files have to be horizontally remapped to the grid of ' +
            'the ERA5 files used (see step_02_preproc_deltas.py).')

    # number of parallel jobs
    parser.add_argument('-p', '--n_par', type=int, default=1,
            help='Number of parallel tasks. Parallelization is done ' +
            'on the level of individual ERA5 files being processed at ' +
            'the same time.')

    # flag to ignore the error from to pressure extrapolation at the model top
    parser.add_argument('-t', '--ignore_top_pressure_error',
            action='store_true',
            help='Flag to ignore an error due to pressure ' +
            'extrapolation at the model top if GCM climate deltas reach ' +
            'up less far than ERA5. This can only be done if ERA5 data ' +
            'is not used by the limited-area model '+
            'beyond the upper-most level of the GCM climate ' +
            'deltas!!')

    # input era5 directory
    parser.add_argument('-D', '--debug_mode', type=str, default=None,
            help='If this flag is set, the ERA5 files will not be ' +
                 'modified but instead the processed climate deltas '
                 'are written to the output directory. There are two ' +
                 'options: for "-D interpolate_time", the climate deltas ' +
                 'are only interpolated to the time of the ERA5 files ' +
                 'and then stored. for "-D interpolate_full", the ' +
                 'full routine is run but instead of the processed ERA5 ' +
                 'files, only the difference between the processed and ' +
                 'the unprocessed ERA5 files is store (i.e. the climate ' +
                 'deltas after full interpolation to the ERA5 grid).')


    args = parser.parse_args()
    ##########################################################################

    # make sure required input arguments are set.
    if args.input_dir is None:
        raise ValueError('Input directory (-i) is required.')
    if args.output_dir is None:
        raise ValueError('Output directory (-o) is required.')
    if args.delta_input_dir is None:
        raise ValueError('Delta input directory (-d) is required.')

    # check for debug mode
    if args.debug_mode is not None:
        if args.debug_mode not in ['interpolate_time', 'interpolate_full']:
            raise ValueError('Invalid input for argument --debug_mode! ' +
                            'Valid arguments are: ' +
                             '"interpolate_time" or "interpolate_full"')

    # first date and last date to datetime object
    first_era_step = datetime.strptime(args.first_era_step, '%Y%m%d%H')
    last_era_step = datetime.strptime(args.last_era_step, '%Y%m%d%H')

    # time steps to process
    era_step_dts = np.arange(first_era_step,
                        last_era_step+timedelta(hours=args.hour_inc_step),
                        timedelta(hours=args.hour_inc_step)).tolist()

    # if output directory doesn't exist create it
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    IMP = IterMP(njobs=args.n_par, run_async=True)
    fargs = dict(
        delta_input_dir = args.delta_input_dir,
        ignore_top_pressure_error = args.ignore_top_pressure_error,
        debug_mode = args.debug_mode,
    )
    step_args = []

    ##########################################################################
    # iterate over time step and prepare function arguments
    for era_step_dt in era_step_dts:
        print(era_step_dt)

        # set output and input ERA5 file
        inp_era_file_path = os.path.join(args.input_dir, 
                era5_file_name_base.format(era_step_dt))
        out_era_file_path = os.path.join(args.output_dir, 
                era5_file_name_base.format(era_step_dt))

        step_args.append(dict(
            inp_era_file_path = inp_era_file_path,
            out_era_file_path = out_era_file_path,
            era_step_dt = era_step_dt
            )
        )

    # choose either main function (pgw_for_era5) for production mode and
    # debug mode "interpolate_full", or function time_interpolation 
    # for debug mode "interpolate_time"
    if (args.debug_mode is None) or (args.debug_mode == 'interpolate_full'):
        run_function = pgw_for_era5
    elif args.debug_mode == 'interpolate_time':
        run_function = debug_interpolate_time
    else:
        raise NotImplementedError()

    # run in parallel if args.n_par > 1
    IMP.run(run_function, fargs, step_args)

