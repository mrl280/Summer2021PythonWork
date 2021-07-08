import pydarn

import datetime as datetime

from pydarn import radar_fov, SuperDARNRadars

try:
    from .add_mlt_to_df import centroid
except ImportError:
    # We are performing local testing
    from DataAnalysis.EchoOccurrence.lib.add_mlt_to_df import centroid


def get_middle_of_fov(station, beam_range, gate_range, date_time_est=None, coords='geo'):
    """

    Return the lon and lat coordinates for the middle of the fov (the centroid)

    :param station: str:
            The radar station to consider, a 3 character string (e.g. "rkn").
            For a complete listing of available stations, please see https://superdarn.ca/radar-info
    :param gate_range: (int, inti):
            Inclusive. The gate range to consider.  If omitted (or None), then all the gates will be considered.
            Note that gates start at 0, so gates (0, 3) is 4 gates.
    :param beam_range: (int, int):
            Inclusive. The beam range to consider.  If omitted (or None), then all beams will be considered.
            Note that beams start at 0, so beams (0, 3) is 4 beams.
    :param date_time_est: datetime.datetime (optional)
            A datetime to use for the magnetic field estimate.  Only required if coords='aacgm'.
    :param coords: str: (optional)
            'geo' for geographical coordinates (default) or 'aacgm' for altitude corrected geomagnetic coordinates
            The coordinates that you would like the middle of the field of view returned in
    :return: (float, float):
            The longitude and latitude of the middle of the field of view
    """

    if date_time_est is None:
        date_time_est = datetime.datetime.now()

    # Compute fov indexes for the extreme corners of the area of interest
    # Remember that a straight indexing gets the bottom left corner of the cell
    starting_gate_idx = gate_range[0]
    ending_gate_idx = gate_range[1] + 1
    starting_beam_idx = beam_range[0]
    ending_beam_idx = beam_range[1] + 1

    # Gets the radars field of view
    all_radars_info = SuperDARNRadars()
    this_radars_info = all_radars_info.radars[pydarn.read_hdw_file(station).stid]  # Grab radar info
    radar_id = this_radars_info.hardware_info.stid
    cell_corners_lats, cell_corners_lons = radar_fov(stid=radar_id, coords=coords, date=date_time_est)

    # Compute coordinate tuples for each of the extreme corners
    lower_left_corner = (cell_corners_lons[starting_gate_idx, starting_beam_idx],
                         cell_corners_lats[starting_gate_idx, starting_beam_idx])
    lower_right_corner = (cell_corners_lons[starting_gate_idx, ending_beam_idx],
                          cell_corners_lats[starting_gate_idx, ending_beam_idx])
    upper_left_corner = (cell_corners_lons[ending_gate_idx, starting_beam_idx],
                         cell_corners_lats[ending_gate_idx, starting_beam_idx])
    upper_right_corner = (cell_corners_lons[ending_gate_idx, ending_beam_idx],
                          cell_corners_lats[ending_gate_idx, ending_beam_idx])

    # And finally go ahead and compute the centroid
    cent_lon, cent_lat = centroid(vertexes=[lower_left_corner, lower_right_corner,
                                            upper_left_corner, upper_right_corner])

    return cent_lon, cent_lat
