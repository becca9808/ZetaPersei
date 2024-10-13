import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import Angle

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from astroplan import Observer, FixedTarget
from astropy.time import Time
from datetime import timezone  # Importing timezone

def query_csv(filename, dec_range = None, ra_range = None, 
    brightness_range = None, max_pol = None, max_e_pol = None, columns = None, 
    sort_by = None, ascending = True):
    # Read the CSV file
    df = pd.read_csv(filename)

    # Initialize filter with all True
    filters = [True] * len(df)

    # Convert RA and DEC columns to SkyCoord objects
    coords = SkyCoord(ra=df['RA (J2000)'], dec=df['Dec (J2000)'], unit=(u.hourangle, u.deg))

    if ra_range:
        # Convert RA range from hourangle to degrees if necessary and ensure it's a pure number
        ra_start, ra_end = ra_range
        ra_start = ra_start.to(u.deg).value if isinstance(ra_start, u.Quantity) else ra_start
        ra_end = ra_end.to(u.deg).value if isinstance(ra_end, u.Quantity) else ra_end
        filters &= (coords.ra.degree >= ra_start) & (coords.ra.degree <= ra_end)
    
    if dec_range:
        # Ensure dec_range is a pure number for comparison
        dec_start, dec_end = dec_range
        dec_start = dec_start.value if isinstance(dec_start, u.Quantity) else dec_start
        dec_end = dec_end.value if isinstance(dec_end, u.Quantity) else dec_end
        filters &= (coords.dec.degree >= dec_start) & (coords.dec.degree <= dec_end)
    
    if brightness_range:
        brightness_max, brightness_min = brightness_range
        brightness_filter = (df['Vmag'] >= brightness_max) & (df['Vmag'] <= brightness_min)
        # brightness_filter = (df['Vmag'] >= brightness_max)
        # brightness_filter = (df['Vmag'] <= brightness_min)
        filters &= brightness_filter
    
    if max_e_pol is not None:
        max_e_pol_filter = (df['e_Pol'] <= max_e_pol)
        filters &= max_e_pol_filter

    if max_pol is not None:
        max_pol_filter = (df['Pol'] <= max_pol)
        filters &= max_pol_filter

    # Apply the filters
    filtered_df = df[filters]

    # If columns are provided, filter the dataframe to only include those columns
    if columns is not None:
        if not all(col in df.columns for col in columns):
            raise ValueError("One or more columns specified do not exist in the dataframe")
        filtered_df = filtered_df[columns]

    # Sort the DataFrame if a sort column is specified
    if sort_by is not None:
        if sort_by not in filtered_df.columns:
            raise ValueError(f"Column '{sort_by}' does not exist in the dataframe")
        filtered_df = filtered_df.sort_values(by=sort_by, ascending=ascending)

    return filtered_df


def query_vizier(identifier):
    try:
        # Query Vizier for the identifier
        result = Vizier.query_object(identifier)
        
        # If the result is empty, print a message and return None
        if not result:
            print("No results found for the given identifier.")
            return None, None
        
        # Extract RA and Dec
        ra_deg = result[0]['RAJ2000'][0]
        dec_deg = result[0]['DEJ2000'][0]
        
        # Convert RA and Dec to SkyCoord object
        coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit='deg')
        
        ra_string = coord.ra.to_string(unit='hour', sep=':', precision=2)
        dec_string = coord.dec.to_string(unit='degree', sep=':', precision=2)
        
        # Print in the desired format
        print(f"RA: {ra_string}")
        print(f"Dec: {dec_string}")

        return ra_string, dec_string
    except Exception as e:
        print(f"Error encountered: {e}")
        return None, None

def query_simbad(identifier):
    """
    This function queries SIMBAD for a given identifier and returns the associated RA and Dec in J2000 coordinates.
    RA is returned in hour angle units, and Dec is returned in degrees.
    The function prints the RA and Dec in the original string format provided by SIMBAD.

    Parameters:
    identifier (str): The identifier name to search for.

    Returns:
    tuple: A tuple containing the RA in hour angle units and Dec in degrees.
    """
    try:
        # Query SIMBAD for the object
        result_table = Simbad.query_object(identifier)
        
        if result_table is None:
            raise ValueError(f"No results found for identifier {identifier}.")
        
        # Extracting the RA and Dec values in string format
        ra_string = result_table['RA'][0]  # RA in hours:minutes:seconds
        dec_string = result_table['DEC'][0]  # Dec in degrees:minutes:seconds
        
        # Print out the original RA and Dec strings
        print(f"RA: {ra_string}, Dec: {dec_string}")
        
        # Convert RA and Dec to Angle objects
        ra_angle = Angle(ra_string, unit=u.hourangle)
        dec_angle = Angle(dec_string, unit=u.deg)
        
        # Return RA as hour angle and Dec in degrees
        return ra_angle.hour, dec_angle.degree
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def plot_altitude_for_targets(observer_location, targets, start_time, end_time,
        filename = "altitude_plot", dpi = 300):
    # Create an Observer object for the observatory location
    observer = Observer(location=observer_location)
    
    # Convert string targets to FixedTarget objects
    target_objs = [FixedTarget.from_name(target) for target in targets]

    # Convert start_time and end_time to Time objects in UTC
    start_time_utc = Time(start_time).to_datetime(observer.timezone).astimezone(timezone.utc)
    end_time_utc = Time(end_time).to_datetime(observer.timezone).astimezone(timezone.utc)
    start_time_utc = Time(start_time_utc)
    end_time_utc = Time(end_time_utc)
    
    # Create a time array to sample the altitudes between start_time and end_time
    time_sampling = start_time_utc + (end_time_utc - start_time_utc) * np.linspace(0, 1, 300)
    
    # Plot the altitude curves for each target
    for target in target_objs:
        altaz = observer.altaz(time_sampling, target)
        altitude = altaz.alt
        plt.plot(time_sampling.datetime, altitude, label=target.name)

    # Set the x-axis to the defined range and format the labels
    plt.xlim(start_time_utc.datetime, end_time_utc.datetime)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=timezone.utc))
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=1))
    plt.gcf().autofmt_xdate()
    
    # Customize the plot
    plt.ylabel("Altitude [deg]")
    plt.xlabel("UTC Time")
    plt.legend(loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename + ".png", dpi = dpi)
    plt.show()

def rotation_mueller(angle_rad):
    cos2a = np.cos(2 * angle_rad)
    sin2a = np.sin(2 * angle_rad)
    return np.array([
        [1, 0,      0,       0],
        [0, cos2a,  sin2a,   0],
        [0, -sin2a, cos2a,   0],
        [0, 0,      0,       1]
    ])

def UV_flip_mueller():
    return np.array([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, -1, 0],
                    [0, 0, 0, -1]
    ])

def calculate_input_aolp(on_sky_aolp, parallactic_angle, altitude, 
    include_M3 = True, instrument = "VAMPIRES"):
    """
    Calculates the input Angle of Linear Polarization (AoLP) after propagating through specified rotations.
    
    Parameters:
        on_sky_aolp (float): The initial angle of linear polarization on the sky in degrees.
        parallactic_angle (float): The parallactic angle in degrees.
        altitude (float): The altitude angle in degrees.
        include_M3 (bool, optional): Whether to include the M3 UV flip. Default is True.
        instrument (str, optional): The instrument being used. Default is "VAMPIRES".
        
    Returns:
        float: The resulting angle of linear polarization in degrees after propagation.
    """
    # Convert degrees to radians
    on_sky_aolp_rad = np.deg2rad(on_sky_aolp)
    parallactic_angle_rad = np.deg2rad(parallactic_angle)
    altitude_rad = np.deg2rad(altitude)
    
    # Define initial Stokes vector for linear polarization
    # Assuming normalized intensity (I=1) and fully polarized light
    S = np.array([
        1,
        np.cos(2 * on_sky_aolp_rad),
        np.sin(2 * on_sky_aolp_rad),
        0
    ])
    
    # Apply rotations
    M_parallactic = rotation_mueller(parallactic_angle_rad)
    if include_M3:
        M_M3 = UV_flip_mueller()
    else:
        M_M3 = np.identity(4)
    M_altitude = rotation_mueller(-altitude_rad)
    
    # Propagate Stokes vector through the rotations
    S_prime = S @ M_parallactic @ M_M3 @ M_altitude
    
    # Calculate resulting AoLP in radians
    aolp_rad = 0.5 * np.arctan2(S_prime[2], S_prime[1])
    
    # Convert resulting AoLP to degrees
    aolp_deg = np.rad2deg(aolp_rad)
    
    # Ensure AoLP is within [0, 180) degrees
    aolp_deg = aolp_deg % 180
    
    return aolp_deg