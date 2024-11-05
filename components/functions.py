from pyopenms import MSExperiment, MzMLFile
import numpy as np
import plotly.express as px
import pandas as pd

def parse_mzml_to_dict(filename):
    # Load the mzML file
    exp = MSExperiment()
    MzMLFile().load(filename, exp)

    # Initialize the dictionary
    ms_dict = {"MS1": {}, "MS2": {}}

    # Iterate through spectra in the experiment
    for spectrum in exp:
        # Get the MS level
        ms_level = spectrum.getMSLevel()

        # Only process MS1 and MS2 spectra
        if ms_level in [1, 2]:
            # Get the retention time
            rt = spectrum.getRT()/60

            # Get m/z and intensity arrays
            mz_array, intensity_array = spectrum.get_peaks()

            # Store in the dictionary
            key = f"MS{ms_level}"
            ms_dict[key][rt] = {
                "m/z": np.array(mz_array),
                "intensity": np.array(intensity_array)
            }

    return ms_dict

def generate_tic(ms_dict):
    """
    Generates a Total Ion Chromatogram (TIC) that includes both MS1 and MS2 data.
    
    Parameters:
    ms_dict (dict): The dictionary containing MS1 and MS2 data with retention time and intensity values.
    
    Returns:
    tuple: Two arrays, one for retention times and one for the summed intensities at each retention time.
    """
    # Initialize dictionaries to store summed intensities
    retention_times = sorted(set(list(ms_dict["MS1"].keys()) + list(ms_dict["MS2"].keys())))
    total_intensities = []

    # Sum intensities for MS1 and MS2 at each unique retention time
    for rt in retention_times:
        ms1_intensity = np.sum(ms_dict["MS1"].get(rt, {"intensity": []})["intensity"])
        ms2_intensity = np.sum(ms_dict["MS2"].get(rt, {"intensity": []})["intensity"])
        total_intensity = ms1_intensity + ms2_intensity
        total_intensities.append(total_intensity)

    return np.array(retention_times), np.array(total_intensities)

def plot_tic(retention_times, total_intensities):
    """
    Plots the Total Ion Chromatogram (TIC) using Plotly, combining MS1 and MS2 data.
    
    Parameters:
    retention_times (array): Array of retention times.
    total_intensities (array): Array of total ion intensities (combined MS1 and MS2).
    """
    # Create a DataFrame for plotting
    tic_data = {"Retention Time": retention_times, "Total Intensity": total_intensities}
    
    # Plot TIC
    fig = px.line(tic_data, x="Retention Time", y="Total Intensity")
    fig.update_xaxes(title="Retention Time (min)")
    fig.update_yaxes(title="Total Intensity")
    fig.show()

def generate_bpc(ms_dict, ms_level="MS1"):
    """
    Generates the Base Peak Chromatogram (BPC) for the specified MS level.
    
    Parameters:
    ms_dict (dict): The dictionary generated from the mzML file with retention time and intensity data.
    ms_level (str): The MS level to generate BPC for (default is "MS1").
    
    Returns:
    tuple: Two arrays, one for retention times and one for the base peak intensities at each retention time.
    """
    # Extract the MS1 or MS2 data
    ms_data = ms_dict.get(ms_level, {})

    # Prepare arrays for retention time and base peak intensity
    retention_times = []
    base_peak_intensities = []

    # Find the base peak (max intensity) for each retention time
    for rt, data in ms_data.items():
        # Check if there are intensities in the current data
        if len(data["intensity"]) > 0:
            base_peak_intensity = np.max(data["intensity"])
        else:
            base_peak_intensity = 0  # No intensity data for this retention time

        retention_times.append(rt)
        base_peak_intensities.append(base_peak_intensity)

    return np.array(retention_times), np.array(base_peak_intensities)

def plot_bpc(retention_times, base_peak_intensities):
    """
    Plots the Base Peak Chromatogram (BPC) using Plotly.
    
    Parameters:
    retention_times (array): Array of retention times.
    base_peak_intensities (array): Array of base peak intensities.
    """
    # Create a DataFrame for plotting
    bpc_data = {"Retention Time": retention_times, "Base Peak Intensity": base_peak_intensities}
    
    # Plot BPC
    fig = px.line(bpc_data, x="Retention Time", y="Base Peak Intensity")
    fig.update_xaxes(title="Retention Time (min)")
    fig.update_yaxes(title="Base Peak Intensity")
    fig.show()

def generate_xic(result_dict, mz_range):
    """
    Generates an Extracted Ion Chromatogram (XIC) for a specified m/z range.
    
    Parameters:
    result_dict (dict): The dictionary with MS1 and MS2 data, containing retention time, m/z, and intensity values.
    mz_range (tuple): A tuple (min_mz, max_mz) defining the m/z range for the XIC.
    
    Returns:
    tuple: Two arrays, one for retention times and one for the summed intensities in the specified m/z range.
    """
    ms1_data = result_dict.get("MS1", {})

    # Prepare arrays for retention times and extracted intensities
    retention_times = []
    xic_intensities = []

    # Iterate over each retention time in MS1 data
    for rt, data in ms1_data.items():
        mz_values = data["m/z"]
        intensity_values = data["intensity"]

        # Find indices where m/z values are within the specified range
        mask = (mz_values >= mz_range[0]) & (mz_values <= mz_range[1])
        
        # Sum intensities within the m/z range for this retention time
        xic_intensity = np.sum(intensity_values[mask]) if np.any(mask) else 0

        # Append to result lists
        retention_times.append(rt)
        xic_intensities.append(xic_intensity)

    return np.array(retention_times), np.array(xic_intensities)

def plot_xic(retention_times, xic_intensities):
    """
    Plots the Extracted Ion Chromatogram (XIC) using Plotly.
    
    Parameters:
    retention_times (array): Array of retention times.
    xic_intensities (array): Array of intensities within the specified m/z range.
    """
    # Create a DataFrame for plotting
    xic_data = {"Retention Time": retention_times, "XIC Intensity": xic_intensities}
    
    # Plot XIC
    fig = px.line(xic_data, x="Retention Time", y="XIC Intensity")
    fig.update_xaxes(title="Retention Time (min)")
    fig.update_yaxes(title="Intensity")
    fig.show()


def get_ms_spectrum(result_dict, rt_value_or_range, mz_grid_resolution=0.01, ms_level="MS1"):
    """
    Finds and processes MS1 spectra for a given retention time range or a single retention time.
    
    Parameters:
    result_dict (dict): The dictionary with MS1 and MS2 data, containing retention time, m/z, and intensity values.
    rt_value_or_range (float or tuple): A single retention time (float) or a range (tuple of start and end retention times).
    
    Returns:
    dict: A dictionary with 'm/z' and 'intensity' arrays. If given a range, intensities are averaged.
    """
    ms_data = result_dict.get(ms_level, {})
    
    if isinstance(rt_value_or_range, tuple):
        # Extract the retention time range
        rt_start, rt_end = rt_value_or_range

        # Filter retention times within the range and accumulate intensities and m/z values
        mz_values = []
        intensity_values = []

        for rt, data in ms_data.items():
            if rt_start <= rt <= rt_end:
                mz_values.append(data["m/z"])
                intensity_values.append(data["intensity"])

        # If we have multiple spectra within the range, interpolate and average intensities
        if intensity_values:
            # Create a common m/z grid over the range of observed m/z values with specified resolution
            min_mz = min(mz.min() for mz in mz_values)
            max_mz = max(mz.max() for mz in mz_values)
            common_mz_grid = np.arange(min_mz, max_mz, mz_grid_resolution)

            # Interpolate each intensity array onto the common m/z grid
            interpolated_intensities = []
            for mz, intensity in zip(mz_values, intensity_values):
                interpolated_intensity = np.interp(common_mz_grid, mz, intensity, left=0, right=0)
                interpolated_intensities.append(interpolated_intensity)
            
            # Average the interpolated intensities
            avg_intensity = np.mean(interpolated_intensities, axis=0)

            return {"m/z": common_mz_grid, "intensity": avg_intensity}
        else:
            raise ValueError("No data found within the specified retention time range.")
    

    elif isinstance(rt_value_or_range, (float, int)):
        # Single retention time provided
        rt = rt_value_or_range
        ms_idx = (result_dict[ms_level][rt] if rt in result_dict[ms_level] else result_dict[ms_level][min(result_dict[ms_level].keys(), key=lambda x: abs(x-rt))])
        return ms_idx
    
    else:
        raise TypeError("Input should be a single float/int retention time or a tuple (start, end) for a range.")
