

def centroid(vertexes):
    """
    Compute the centroid of a polygon
    :param vertexes: a tuple with coordinates of the polygon points
    :return: a tuple with the centroid coordinates
    """
    _x_list = [vertex[0] for vertex in vertexes]
    _y_list = [vertex[1] for vertex in vertexes]
    _len = len(vertexes)
    _x = sum(_x_list) / _len
    _y = sum(_y_list) / _len
    return _x, _y
