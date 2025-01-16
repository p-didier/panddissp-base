import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from package.gui_utils import RIRg_GUI

def check_plot_tdoas(doaEstTarget, doaEstAll, asc: RIRg_GUI):
    """
    Generates a plot of the acoustic scenario to check the validity of TDOA
    estimation obtained, e.g., via the frequency-domain MUSIC algorithm.

    Paramters
    ---------
    doaEstTarget : list of floats
        Estimated DOA for target talker(s).
    doaEstAll : list of floats
        All estimated DOAs (target talker(s) and noise source(s)).
    asc : RIRg_GUI class instance
        Acoustic scenario parameters.
    """

    def _get_line(doaRadiants, startCoords, roomDim):
        """
        Get 2D line (x and y coordinates) from DOA.
        """
        x0, y0 = startCoords[0], startCoords[1]
        x1 = x0 + roomDim * np.cos(doaRadiants)
        y1 = y0 + roomDim * np.sin(doaRadiants)
        return [x0, x1], [y0, y1]

    # Check inputs
    if isinstance(doaEstAll, float):
        doaEstAll = [doaEstAll]  # convert to list
    # Plot room
    fig = asc.plot_asc()
    # Plot DOA estimates
    for ii in range(len(doaEstAll)):
        x, y = _get_line(
            doaRadiants=doaEstAll[ii] + np.pi / 2,  # align to space orientation
            startCoords=np.mean(asc.micCoords, axis=0),
            roomDim=asc.roomDim
        )
        if doaEstAll[ii] in doaEstTarget:
            rayColor = 'g'  # indicate selected (target) DOAs with green
        else:
            rayColor = '0.7'
        plt.plot(x, y, '--', color=rayColor)
        plt.text(
            np.mean(x),
            np.mean(y),
            f'$\\hat{{\\theta}}_{ii+1}$={np.round(doaEstAll[ii] * 180 / np.pi, 2)}$^\\circ$',
            c=rayColor
        )
    return fig


def select_latest_rir(path):
    """
    Returns the full path to the latest RIR that was generated (via the
    RIR-generating GUI) in the folder `path`.
    """

    p = Path(path).glob('**/*')
    files = [x for x in p if x.is_file()]
    if str(files[-1])[-2:] == 'gz':
        rirFile = files[-1]
    else:
        rirFile = files[-2]
    
    return rirFile


def auto_choice_doa(doaEsts, asc: RIRg_GUI):
    """
    Automatically selects the DOA(s) corresponding to the target speaker(s)
    among all estimated DOAs contained in `doaEsts`, based the information
    contained in the RIR-GUI output object `asc`.
    """

    coordinatesTalkers = asc.audioCoords
    coordinatesNoises = asc.noiseCoords
    coordinatesMicArray = np.mean(np.array(asc.micCoords), axis=0)

    oracleDOAtalkers = np.zeros(len(coordinatesTalkers))
    for ii in range(len(coordinatesTalkers)):
        oracleDOAtalkers[ii] = np.arctan2(
            coordinatesMicArray[0] - coordinatesTalkers[ii][0],
            coordinatesMicArray[1] - coordinatesTalkers[ii][1]
        )

    oracleDOAnoises = np.zeros(len(coordinatesNoises))
    for ii in range(len(coordinatesNoises)):
        oracleDOAnoises[ii] = np.arctan2(
            coordinatesMicArray[0] - coordinatesNoises[ii][0],
            coordinatesMicArray[1] - coordinatesNoises[ii][1]
        )

    chosenDOAs = np.zeros_like(oracleDOAtalkers)
    for ii in range(len(oracleDOAtalkers)):
        chosenDOAs[ii] = doaEsts[
            (np.abs(doaEsts - oracleDOAtalkers[ii])).argmin()
        ]

    return chosenDOAs, oracleDOAtalkers