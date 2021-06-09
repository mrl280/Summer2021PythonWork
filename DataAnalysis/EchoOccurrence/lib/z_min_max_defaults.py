

def z_min_max_defaults(parameter):
    """
    :param parameter: The fitACF parameter as a string (e.g. 'v')
    :return: (<int>, <int>): The default z-min and z-max values for the provided parameter
    """

    defaultzminmax = {'p_l': [0, 50], 'v': [-600, 600],
                      'w_l': [0, 250], 'elv': [0, 50]}

    try:
        return defaultzminmax[parameter][0], defaultzminmax[parameter][1]
    except KeyError:
        raise KeyError("z_min_max_defaults() does not recognize the parameter: " + str(parameter))
