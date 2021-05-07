RISR data is labeled with experiment start date, each experiment is one file, even if the experiment ran for a couple days.

RISR data is stored as hdf5 files.  Each experiment has a number of different files (file type checked /Metadata/Experiment Notes)
RISR‐N is typically scheduled to run for 5–6 days per month in a common mode WorldDay.

Sometimes there may be seemingly duplicate files.  If so, the files differ in integration times (/Data/Array Layout/1D Parameters/inttms).
The longer the integration time, the fewer data points held within a file.

Density (Ne) data is corrected.  I believe Te/Ti is used to correct.

Newer files provide data access through both a table layout and an array layout.
The array layout is more user friendly, but older data only supports the table layout

RISR file types:
    - Long Pulse (480)
        A long pulse (LP) with 72 km range resolution designed for F region studies.
        Contains 1D parameters that vary with time (e.g. beam id, frequency, and power info)
        Contains 2D parameters that vary with time and range (e.g. Ne, ion and electron temps, LOS velocities, and composition)

    - Alternating Code (AC16-30)
        An alternating code pulse with 4.5 km resolution used for E region studies.
        Same 1D and 2D parameters as Long Pulse, just different data resolution

    - Vector velocities
        Experiment data derived from Long Pulse measurements
        Contains 1D parameters that vary with time (e.g. altitude, integration time)
        Contains 2D parameters that vary with time and corrected geomagnetic latitude
            (e.g. vector velocity and electric field components in three dimensions)

    - Long Pulse Uncorrected Ne
        Experiment data derived from Long Pulse measurements
        Contains 1D parameters that vary with time (e.g. beam id, frequency, and power info)
        Contains 2D parameters that vary with time and corrected geomagnetic latitude (e.g. Ne and altitude)