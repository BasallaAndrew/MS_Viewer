import pyopenms as oms
import plotly.express as px
import pandas as pd

def generate_TIC(file) -> px.line:
    exp = oms.MSExperiment()
    oms.MzMLFile().load(file, exp)
    retention_times = []
    intensities = []
    for spec in exp:
        if spec.getMSLevel() == 1:
            if sum(spec.get_peaks()[1]) != 0:
                retention_times.append(spec.getRT()/60)
                intensities.append(sum(spec.get_peaks()[1]))
    centroided_spectra = oms.MSExperiment()
    oms.PeakPickerHiRes().pickExperiment(
    exp, centroided_spectra, True
    )
    idx = 0 
    ms1_dict = {}
    ms2_dict = {}
    for spec in exp:
        if spec.getMSLevel() == 1:
            ms1_dict[spec.getRT()/60] = idx
        if spec.getMSLevel() == 2:
            ms2_dict[spec.getRT()/60] = idx
        idx += 1
    return retention_times, intensities, ms1_dict, ms2_dict


