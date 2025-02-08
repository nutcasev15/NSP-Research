# OpenMC 0.15.0 Parametrised Hexagonal Fuel Element Model


############### Library Imports
from collections.abc import Callable
import openmc
from openmc.checkvalue import check_type, check_iterable_type


############### Define Function to Calculate Reactivity Coefficients
def calculate_reactivity_coef(model : openmc.Model, param_list,
                                model_updater, updater_static_args=None,
                                print_output=True, init_args=None,
                                started_initialized=False, **kwargs):
    """Perform an accelerated reactivity calculation by modifying the
    model in memory using the C API.


    Parameters
    ----------
    param_list : collections.Iterable
        Iterable of floats to pass to the `model_updater` method
        that serves as the independent variable during the reactivity
        calculation.
    model_updater : collections.Callable
        Callable function which updates an `openmc.Model` object passed
        to it using a single numerical input variable. Cell temperature,
        rotation and translation updates are allowed,
        as well as material density and composition updates.
    updater_static_args : dict, optional
        Keyword-based arguments to pass to the `model_updater` method
        which do not change during the reactivity calculation. Defaults
        to no arguments.
    print_output : bool
        Whether or not to output the results during the calculation
        process. Defaults to True.
    init_args : dict, optional
        Keyword arguments to pass to :meth:`openmc.Model.init_lib`. Defaults
        to no arguments.
    started_initialized : bool
        Whether or not to initialize and finalize the C API interface. Defaults
        to False, which assumes the C API interface was previously inactive.
    **kwargs
        All remaining keyword arguments are passed to :meth:`openmc.Model.run`

    Returns
    -------
    values : List of float
        Numerical outputs of the `model_updater` method. Defaults to a
        duplicate of the param_list input variable if the method returns
        None.
    coefficients : List of float
        List of the calculated reactivity coefficients

    """

    check_iterable_type('param_list', param_list, float)

    if updater_static_args is None:
        updater_static_args = {}
    else:
        check_type('updater_static_args', updater_static_args, dict)
    check_type('model_updater', model_updater, Callable)
    check_type('model_updater',
                model_updater(model, param_list[0], **updater_static_args),
                float, none_ok=True)

    if init_args is None:
        init_args = {}
    else:
        check_type('init_args', init_args, dict)

    check_type('print_output', print_output, bool)

    # Check the return of model_updater method. If there is output,
    # we need to save the output to return it to the user.
    if model_updater(model, param_list[0], **updater_static_args) is None:
        duplicate_input = True
        values = param_list
    else:
        duplicate_input = False
        values = []

    # Define container for the calculated data
    coefficients = []

    # Initialize C API and export initial model xml file
    if not started_initialized:
        model.init_lib(**init_args, output=print_output)

    # Update model and calculate raw keff for each input value
    for param in param_list:
        # Update model and record the input corresponding
        # to the keff simulation.
        if duplicate_input is True:
            input = param
            model_updater(model, input, **updater_static_args)
        else:
            input = model_updater(model, param, **updater_static_args)
        values.append(input)

        # Run using C API in memory and save raw simulation output
        sp_filepath = model.run(**kwargs, output=print_output)
        with openmc.StatePoint(sp_filepath) as sp:
            raw_output = sp.keff
            text = 'Input Variable : {:.3f}; Model produced ' + \
                'a keff of {:1.5f} +/- {:1.5f}'
            print(text.format(input,  raw_output.n, raw_output.s))
            coefficients.append(raw_output)

    # Forward difference raw keff values and calculate reactivity
    # coefficients.
    for i in range(0, len(values) - 1):
        coef = coefficients[i + 1] - coefficients[i]
        coef /= coefficients[i + 1] * coefficients[i]
        coef /= values[i + 1] - values[i]
        coefficients[i] = coef

    # Discard last data point due to forward differencing
    del values[-1]
    del coefficients[-1]

    # Deallocate resources
    if not started_initialized:
        model.finalize_lib()

    return values, coefficients
