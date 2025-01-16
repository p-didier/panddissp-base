import os
import copy
import gzip
import pickle
import datetime
import numpy as np
import PySimpleGUI as sg  # tip: run `psgdemos.exe` for GUI demos
from pathlib import Path
import pyroomacoustics as pra
import matplotlib.pyplot as plt

# In order: [audio, noise, mics]
GRAPHELEMENTSCOLORS = np.array(['red', 'black', 'blue']) 
RADIUS_DOTS = 3  # radius of dots used to place the elements in the canvas
DEFAULT_NMICSPERARRAY = 5
DEFAULT_DBWMICS = 5
DEFAULT_ROOMDIM = 5
DEFAULT_T60 = 0.
DEFAULT_RIRLENGTH = 22050
DEFAULT_FS = 44100
class RIRg_GUI:

    def __init__(
            self,
            theme='DarkAmber',
            c=340,
            nMicsPerArray=DEFAULT_NMICSPERARRAY,
            distBwMics=DEFAULT_DBWMICS,
            roomDim=DEFAULT_ROOMDIM,
            t60=DEFAULT_T60,
            rirLength=DEFAULT_RIRLENGTH,
            fs=DEFAULT_FS,
            exportFolder=os.getcwd(),
            outputRIRplot=False
        ) -> None:

        # Defaults
        self.theme = sg.theme(theme)   # Add a touch of color
        self.c = c    # speed of sound [m/s]
        self.nMicsPerArray = nMicsPerArray
        self.distBwMics = distBwMics
        self.roomDim = roomDim
        self.t60 = t60
        self.rirLength = rirLength
        if int(fs) < 8000:
            print(f'/!\ /!\ /!\ The desired sampling frequency ({int(fs)} Hz) is too low. Setting fs = 8 kHz.')
            fs = 8000
        self.fs = int(fs)
        self.exportFolder = exportFolder
        self.outputRIRplot = outputRIRplot

        # Run
        self.run_window()

    def __repr__(self) -> str:
        """On `print()`."""
        strout = "Acoustic scenario object containing the RIRs, with parameters:\n"
        # vvv Case-insensitive sorting of list of strings
        # https://stackoverflow.com/a/10269828
        fields = sorted(list(self.__dict__.keys()), key=str.casefold)
        for field in fields:
            if 'RIR' not in field:  # avoid printing the RIRs themselves
                strout += f'>> Field "{field}": {getattr(self, field)}\n'
        return strout
        
    def run_window(self):
        """Callback function for GUI window and user interaction."""

        exclusiveButtons = (
            '-AUDIO BUTTON-',
            '-NOISE BUTTON-',
            '-MIC BUTTON-'
        ) 
        # ^^^ graph scatter element colours
        np_exclusiveButtons = np.array(exclusiveButtons)

        prevButtonStates = np.array([False, False, False])

        window = self.make_window()

        # Lists
        self.micsIds = []
        self.audioIds = []
        self.noiseIds = []
        self.micCoords = []
        self.audioCoords = []
        self.noiseCoords = []
        self.lineIds = []
        self.lineTextIds = []

        # Initialise graph area
        previousRoomDim = float(window['-ROOMDIM INPUT-'].get())
        self.update_graph_area(window)

        # Event Loop 
        while True:
            event, values = window.read(timeout=100)
            currButtonStates = np.array(
                [window[b].get() for b in exclusiveButtons]
            )
            if event in ('Exit', sg.WIN_CLOSED, None):
                break

            # Ensure that only one of the three exclusive buttons
            # (audio, noise, mic) is ticked at the same time.
            if any(currButtonStates != prevButtonStates):
                if sum(currButtonStates) > 1:
                    newestTickedButton = np_exclusiveButtons[
                        currButtonStates != prevButtonStates
                    ]
                    newestTickedButton = newestTickedButton[0]
                    for b in exclusiveButtons:
                        if b != newestTickedButton:
                            window[b].update(value=False)
                currButtonStates = np.array([
                    window[b].get() for b in exclusiveButtons
                ])
                prevButtonStates = currButtonStates

            # Draw circles on graph
            elif event == "-GRAPH-":
                self.update_attributes_from_user_inputs(window)
                if any(currButtonStates):
                    graph = window['-GRAPH-']       # type: sg.Graph
                    color = GRAPHELEMENTSCOLORS[currButtonStates][0]
                    if np_exclusiveButtons[currButtonStates][0] ==\
                        '-MIC BUTTON-':
                        if values['-NMICS INPUT-'] == '':
                            print('!! Need a number of microphones !!')
                        else:
                            for ii in range(int(values['-NMICS INPUT-'])):
                                # Convert mic distance from cm to canvas
                                # dimensions.
                                delta = float(values['-MICDIST-']) * 1e-2 *\
                                    graph.CanvasSize[1] /\
                                    float(window['-ROOMDIM INPUT-'].get())
                                # Compute coordinates in canvas dimensions.
                                canvasCoords = (
                                    values['-GRAPH-'][0],
                                    values['-GRAPH-'][1] +\
                                        ii * delta
                                )
                                # Draw circle
                                id = graph.draw_circle(
                                    canvasCoords,
                                    fill_color=color,
                                    radius=RADIUS_DOTS
                                )
                                self.micsIds.append(str(id))  # store mic ID
                                # Store coordinates for RIRs computation
                                coords = list(np.array(canvasCoords) /\
                                    np.array(graph.CanvasSize) * self.roomDim)
                                self.micCoords.append(coords)
                    else:
                        id = graph.draw_circle(
                            values['-GRAPH-'],
                            fill_color=color,
                            radius=RADIUS_DOTS
                        )
                        if np_exclusiveButtons[currButtonStates][0] ==\
                            '-AUDIO BUTTON-':
                            self.audioIds.append(str(id))     # store source ID
                            # Store coordinates for RIRs computation
                            coords = list(np.array(values['-GRAPH-']) /\
                                np.array(graph.CanvasSize) * self.roomDim)
                            self.audioCoords.append(coords)
                        if np_exclusiveButtons[currButtonStates][0] ==\
                            '-NOISE BUTTON-':
                            self.noiseIds.append(str(id))     # store source ID
                            # Store coordinates for RIRs computation
                            coords = list(np.array(values['-GRAPH-']) /\
                                np.array(graph.CanvasSize) * self.roomDim)
                            self.noiseCoords.append(coords)
            # "Delete" buttons
            elif event == '-DEL AUDIO BUTTON-':
                graph = window['-GRAPH-']
                graph._TKCanvas2.delete(*self.audioIds)
                self.audioIds = []  # reset IDs list
                self.audioCoords = []
            elif event == '-DEL NOISE BUTTON-':
                graph = window['-GRAPH-']
                graph._TKCanvas2.delete(*self.noiseIds)
                self.noiseIds = []  # reset IDs list
                self.noiseCoords = []
            elif event == '-DEL MIC BUTTON-':
                graph = window['-GRAPH-']
                graph._TKCanvas2.delete(*self.micsIds)
                self.micsIds = []  # reset IDs list
                self.micCoords = []
            elif event == '-DEL ALL COMPONENTS-':
                graph = window['-GRAPH-']
                graph._TKCanvas2.delete('all')
                self.micsIds = []  # reset IDs list
                self.noiseIds = []  # reset IDs list
                self.audioIds = []  # reset IDs list
                self.micCoords = []
                self.audioCoords = []
                self.noiseCoords = []
                # Add grid
                self.draw_line_at_every_meter(window)
            elif event == '-RESET ALL-':
                output = sg.popup_yes_no(
                    'Reset everything?',
                    keep_on_top=True
                )
                if output == 'Yes':
                    # Reset everything
                    graph = window['-GRAPH-']
                    graph._TKCanvas2.delete('all')
                    self.micsIds = []  # reset IDs list
                    self.noiseIds = []  # reset IDs list
                    self.audioIds = []  # reset IDs list
                    self.micCoords = []
                    self.audioCoords = []
                    self.noiseCoords = []
                    # Add grid
                    self.draw_line_at_every_meter(window)
                    # Reset inputs
                    window['-NMICS INPUT-'].Update(DEFAULT_NMICSPERARRAY)
                    window['-MICDIST-'].Update(DEFAULT_DBWMICS)
                    window['-ROOMDIM INPUT-'].Update(DEFAULT_ROOMDIM)
                    window['-T60 INPUT-'].Update(DEFAULT_T60)
                    window['-RIRLEN INPUT-'].Update(DEFAULT_RIRLENGTH)
                    window['-FS CHOICE-'].Update(DEFAULT_FS)
            elif event == '-COMPUTE RIRS-':
                # Check for invalid inputs first
                if len(self.micCoords) == 0:
                    print('WARNING: no microphone array defined. Cannot compute RIRs.')
                elif len(self.audioCoords) == 0 and len(self.noiseCoords) == 0:
                    print('WARNING: no source defined. Cannot compute RIRs.')
                else:
                    # Read latest user inputs
                    self.update_attributes_from_user_inputs(window)
                    # ===================================
                    # If all good, compute and store RIRs
                    # ===================================
                    if len(self.micCoords) > 0:
                        if len(self.audioCoords) == 0 and\
                            len(self.noiseCoords) == 0:
                            print('No sound source present. Please place audio and/or noise source(s) and try again.')
                        else:
                            fname = f'{self.exportFolder}/rirs_{get_datetime()}'
                            self.RIRsAudio, self.RIRsNoise = compute_rirs(
                                micPos=np.array(self.micCoords),
                                audioPos=np.array(self.audioCoords),
                                noisePos=np.array(self.noiseCoords),
                                roomDim=self.roomDim,
                                t60=self.t60,
                                fsRIR=self.fs,
                                rirLength=self.rirLength,
                                outputRIRplot=self.outputRIRplot,
                                fname=fname,
                                c=self.c,
                            )
                            if self.RIRsAudio is not None and self.RIRsNoise is not None:
                                print('RIRs computed successfully. Saving as Pickle archive.')
                                fullFname = fname + '.pkl.gz'
                                pickle.dump(self, gzip.open(fullFname, 'wb'))
                                print(f'RIRs saved in file: "{Path(fullFname).name}", in folder\n"{Path(fullFname).parent}"')
                                print('You may close the GUI if not needed anymore!')
                            else:
                                print('RIRs computation failed. Please try again.')
                    else:
                        print('No microphones are present. Please place microphones and try again.')
                    # ===================================
                    # ===================================
            elif event == 'Load layout':
                pathToFile = sg.popup_get_file(
                    'Choose your GUI layout file (.pkl.gz)', keep_on_top=True
                )
                if pathToFile is not None and pathToFile != '':
                    session = pickle.load(gzip.open(pathToFile, 'r'))
                    self.update_from_loaded_layout(window, session)
                else:
                    print('No file selected.')
                
            elif event == 'Save layout':
                # Read latest user inputs
                self.update_attributes_from_user_inputs(window)
                # Save layout as Pickle archive
                fname = f'GUIlayout_{get_datetime()}.pkl.gz'
                pickle.dump(self, gzip.open(fname, 'wb'))
                print(f'GUI layout saved in file: "{Path(fname).name}", in folder\n"{Path(fname).parent}"')
                
            elif event == 'Help content':
                sg.popup(
                    'Please refer to the PDF for Session 1',
                    keep_on_top=True
                )
            elif event == 'About':
                sg.popup(
                    'Speech and Audio Python GUI v1.0',
                    'Authors: Paul Didier',
                    'December 16, 2022',
                    'Original GUI (Matlab) by R. Ali, G. Bernardi, and A. Bertrand',
                    keep_on_top=True
                )
            
            # Other
            if window['-ROOMDIM INPUT-'].get() != '':
                lastInputRoomDim = float(window['-ROOMDIM INPUT-'].get())
                if lastInputRoomDim < 2:
                    print('[INVALID VALUE] Please enter a room dimension larger than 2 metres.')
                elif previousRoomDim != lastInputRoomDim:
                    self.update_graph_area(window)
                    previousRoomDim = lastInputRoomDim
        
        window.close()
        # exit(0)   # <- comment out to continue running after closing the GUI

    def update_graph_area(self, window):
        """
        Update the graph area according to the room dimensions.
        """

        def _update_points(coords, color):
            newCoords = []
            ids = []
            for ii in range(len(coords)):
                newCanvasCoords = list(np.array(coords[ii]) * multRatio)
                id = graph.draw_circle(
                    newCanvasCoords,
                    fill_color=color,
                    radius=RADIUS_DOTS
                )
                ids.append(id)
                newCoords.append(
                    list(np.array(newCanvasCoords) / multRatio)
                )
            return newCoords, ids

        graph = window['-GRAPH-']   # type: sg.Graph
        currRoomDim = float(window['-ROOMDIM INPUT-'].get())
        multRatio = graph.CanvasSize[1] / currRoomDim
        
        # Delete previous lines and texts
        graph._TKCanvas2.delete(*self.lineIds)
        graph._TKCanvas2.delete(*self.lineTextIds)

        # Delete previous points
        graph._TKCanvas2.delete(*self.audioIds)
        graph._TKCanvas2.delete(*self.noiseIds)
        graph._TKCanvas2.delete(*self.micsIds)
        self.audioIds = []
        self.noiseIds = []
        self.micsIds = []
        # Update audio points
        newCoords, newIds = _update_points(
            self.audioCoords,
            GRAPHELEMENTSCOLORS[0]
        )
        self.audioCoords = newCoords
        self.audioIds = newIds
        # Update noise points
        newCoords, newIds = _update_points(
            self.noiseCoords,
            GRAPHELEMENTSCOLORS[1]
        )
        self.noiseCoords = newCoords
        self.noiseIds = newIds
        # Update mic points
        newCoords, newIds = _update_points(
            self.micCoords,
            GRAPHELEMENTSCOLORS[2]
        )
        self.micCoords = newCoords
        self.micsIds = newIds

        # Add grid
        self.draw_line_at_every_meter(window)

    def draw_line_at_every_meter(self, window):
        """
        Draws a grid on the canvas, with a line marking every meter in both
        dimensions.
        """

        graph = window['-GRAPH-']   # type: sg.Graph
        currRoomDim = float(window['-ROOMDIM INPUT-'].get())

        # Draw lines at every meter
        deltaLine = graph.CanvasSize[1] / currRoomDim

        for ii in range(int(np.ceil(currRoomDim))):
            # Draw horizontal line
            id = graph.draw_line(
                point_from=(0.0, ii * deltaLine),
                point_to=(graph.CanvasSize[1], ii * deltaLine),
                color='gray',
                width=1
            )
            self.lineIds.append(id)
            # Add text
            id = graph.draw_text(
                text=f'{ii} m',
                location=(12, ii * deltaLine + 7),
                color='gray'
            )
            self.lineTextIds.append(id)
            # Draw vertical line
            id = graph.draw_line(
                point_from=(ii * deltaLine, 0.0),
                point_to=(ii * deltaLine, graph.CanvasSize[1]),
                color='gray',
                width=1
            )
            self.lineIds.append(id)
            # Add text
            if ii > 0:
                id = graph.draw_text(
                    text=f'{ii} m',
                    location=(ii * deltaLine + 12, 7),
                    color='gray',
                )
                self.lineTextIds.append(id)

    def plot_asc(self):
        """
        Plots a visual representation of the acoustic scenario (in 2D).
        """
        fig, axes = plt.subplots(1,1)
        fig.set_size_inches(8.5, 3.5)
        for ii in range(len(self.audioCoords)):
            axes.scatter(
                x=self.audioCoords[ii][0],
                y=self.audioCoords[ii][1],
                color='r',
                marker='x'
            )
            axes.text(
                x=self.audioCoords[ii][0],
                y=self.audioCoords[ii][1],
                s=f'a{ii+1}',
                color='r'
            )
        for ii in range(len(self.noiseCoords)):
            axes.scatter(
                x=self.noiseCoords[ii][0],
                y=self.noiseCoords[ii][1],
                color='k',
                marker='s'
            )
            axes.text(
                x=self.noiseCoords[ii][0],
                y=self.noiseCoords[ii][1],
                s=f'n{ii+1}',
                color='k'
            )
        for ii in range(len(self.micCoords)):
            axes.scatter(
                x=self.micCoords[ii][0],
                y=self.micCoords[ii][1],
                color='b',
                marker='o'
            )
        axes.grid()
        axes.set_xlim([0, self.roomDim])
        axes.set_ylim([0, self.roomDim])
        axes.set_aspect('equal', 'box')
        axes.set_xlabel('$x$ [m]')
        axes.set_ylabel('$y$ [m]')
        plt.tight_layout()


    def update_attributes_from_user_inputs(self, window):
        """Updates class parameters from user input in GUI window."""
        self.numMics = int(window['-NMICS INPUT-'].get())
        self.distBwMics = float(window['-MICDIST-'].get()) * 1e-2 # [m]
        self.roomDim = float(window['-ROOMDIM INPUT-'].get())
        self.t60 = float(window['-T60 INPUT-'].get())
        self.fs = int(window['-FS CHOICE-'].get() * 1e3) # [Hz]
        self.rirLength = int(window['-RIRLEN INPUT-'].get())


    def update_from_loaded_layout(self, window: sg.Window, layout):
        """Updates current window from loaded GUI layout."""
        # Update graph zone
        graph: sg.Graph = window['-GRAPH-']
        graph._TKCanvas2.delete(*self.audioIds)
        graph._TKCanvas2.delete(*self.noiseIds)
        graph._TKCanvas2.delete(*self.micsIds)

        multRatio = graph.CanvasSize[0] / layout.roomDim

        def _update_canvas(loadedCoords, color):
            """Helper function."""
            newIds = []
            for ii in range(len(loadedCoords)):
                canvasCoords = list(
                    np.array(loadedCoords[ii]) *  multRatio
                )
                id = graph.draw_circle(
                    canvasCoords,
                    fill_color=color,
                    radius=RADIUS_DOTS
                )
                newIds.append(id)
            return newIds

        # Plot audio sources
        newIds = _update_canvas(layout.audioCoords, GRAPHELEMENTSCOLORS[0])
        self.audioIds = newIds
        # Plot noise sources
        newIds = _update_canvas(layout.noiseCoords, GRAPHELEMENTSCOLORS[1])
        self.noiseIds = newIds
        # Plot mics
        newIds = _update_canvas(layout.micsCoords, GRAPHELEMENTSCOLORS[2])
        self.micsIds = newIds
        # Update coordinates
        self.audioCoords = layout.audioCoords
        self.micCoords = layout.micsCoords
        self.noiseCoords = layout.noiseCoords

        return None


    def make_window(self, theme=None):
        """Generate GUI window."""

        if theme is None:
            sg.theme(self.theme)
        else:
            sg.theme(theme)

        menu_def = [['&File', ['Load layout', 'Save layout']],
                    ['&Help', ['Help content', 'About']] ]

        layout_graphCol = [
            [sg.Text("Acoustic environment")],
            [sg.Graph(
                canvas_size=(200,200),
                graph_bottom_left=(0,0),
                graph_top_right=(200,200),
                background_color="white",
                key='-GRAPH-',
                enable_events=True
            )],
            # [sg.Image(
            #     data=sg.DEFAULT_BASE64_LOADING_GIF,
            #     enable_events=True,
            #     key='-GIF-IMAGE-'
            # )],
            [sg.Button("RESET ALL", key='-RESET ALL-')],
            [sg.Sizegrip()]
        ]

        layout_components = [
            [sg.Text('Add component', s=12)],
            [
                sg.Checkbox("Audio", key='-AUDIO BUTTON-'),
                sg.Button("Del", key='-DEL AUDIO BUTTON-'),
                sg.Text('', key='-N AUDIOSOURCES-')
            ],
            [
                sg.Checkbox("Noise", key='-NOISE BUTTON-'),
                sg.Button("Del", key='-DEL NOISE BUTTON-'),
                sg.Text('', key='-N NOISESOURCES-')
            ],
            [
                sg.Checkbox("Mic", key='-MIC BUTTON-'),
                sg.Button("Del", key='-DEL MIC BUTTON-'),
                sg.Text('', key='-NMICS DISPLAY-')
            ],
            [sg.Text('# mics'), sg.Input(
                s=15, key='-NMICS INPUT-', default_text=self.nMicsPerArray
            )],
            [sg.Text('d [cm]'), sg.Input(
                s=15, key='-MICDIST-', default_text=self.distBwMics
            )],
            [sg.Button("Reset Components", key='-DEL ALL COMPONENTS-')],
        ]

        layout_params = [
            [sg.Text('Parameters', s=12)],
            [sg.Text('Room\ndim. [m]'), sg.Input(
                s=15, key='-ROOMDIM INPUT-', default_text=self.roomDim
            )],
            [sg.Text('T60 [s]'), sg.Input(
                s=15, key='-T60 INPUT-', default_text=self.t60
            )],
            [sg.Text('Length RIR\n[samples]'), sg.Input(
                s=15, key='-RIRLEN INPUT-', default_text=self.rirLength
            )],
            [sg.Text('fs [kHz]'), sg.Combo(
                [8, 16, 44.1],
                default_value=self.fs / 1e3,
                s=(15,22),
                enable_events=True,
                readonly=True,
                k='-FS CHOICE-'
            )],
            [sg.Button("Create/Store RIR", key='-COMPUTE RIRS-')],
        ]
        
        layout = [
            [
                sg.MenubarCustom(
                    menu_def, key='-MENU-', tearoff=True
                )
            ],
            [
                sg.Col(layout_graphCol, p=0),
                sg.Col(layout_components, p=0),
                sg.Col(layout_params, p=0),        
            ]
        ]

        window = sg.Window(
            'Speech and Audio GUI',
            layout,
            right_click_menu_tearoff=True,
            grab_anywhere=False,
            resizable=True,
            margins=(0,0),
            # use_custom_titlebar=True,
            finalize=True,
            keep_on_top=True
        )
        window.set_min_size(window.size)
        return window


def compute_rirs(
        micPos,
        audioPos,
        noisePos,
        roomDim,
        t60,
        fsRIR,
        rirLength=None,
        outputRIRplot=False,
        fname='',
        c=340
    ):
    """
    Computes RIRs based on MATLAB function from the KUL course "P&D ISSP 2022"
    `create_rirs.m`.
    Major difference: not using the `simroommex.m` MATLAB MEX file.
    Instead: pyroomacoustics Python package.

    Parameters
    ----------
    -micPos : [Nm x 2] np.ndarray (float)
        Microphones coordinates [m].
    -audioPos : [Ns x 2] np.ndarray (float)
        Sources coordinates [m].
    -noisePos : [Nn x 2] np.ndarray (float)
        Noise sources coordinates [m].
    -roomDim : float
        Room dimension [m].
    -t60 : float
        T60 reverberation time [s].
    -fsRIR : int or float
        RIR sampling rate [samples/s].
    -rirLength : int
        Impulse response length [samples].
    -outputRIRplot : bool
        If True, outputs a plot of the RIRs.
    -fname : str
        Path (incl. name) of file for output RIR plot.
    -c : float
        Speed of sound [m/s].

    Returns
    -------
    RIRsAudio : [`rirLength` x Nm x Ns] np.ndarray (float)
        RIRs between each audio source and each microphone.
    RIRsNoise : [`rirLength` x Nm x Nn] np.ndarray (float)
        RIRs between each noiise source and each microphone.
    """

    rd = [roomDim, roomDim, 4]     # Room dimensions [x, y, z] (m)
    nMics = micPos.shape[0]
    nAudio = audioPos.shape[0]
    nNoise = noisePos.shape[0]
    micPos3D = np.concatenate(
        (micPos, np.full((nMics, 1), fill_value=2)), axis=1
    )
    if nAudio > 0:
        audioPos3D = np.concatenate(
            (audioPos, np.full((nAudio, 1), fill_value=2)), axis=1
        )
    if nNoise > 0:
        noisePos3D = np.concatenate(
            (noisePos, np.full((nNoise, 1), fill_value=2)), axis=1
        )

    if t60 != 0:  # check T60
        # Compute using Sabine's formula
        vol = np.prod(rd)
        surf = (rd[0] * rd[1] + rd[0] * rd[2] + rd[2] * rd[1]) * 2
        refl = np.amax((0, 1 - 0.161 * vol / (t60 * surf)))
        if refl > 0.92:  # hard threshold
            refl = 0.92
            t60 = 0.161 * vol / ((1 - refl) * surf)
            print(f'WARNING: too large reverberation time, RIRs created with reverberation time T60={np.round(t60, 3)} s')
            
    if rirLength is None:
        rirLength = np.round(np.amax(
            t60 * 1.2 * fsRIR,
            2 * fsRIR * np.amax(rd) / c
        ))  # Number of samples
        rirLength += rirLength % 2

    print(f'Note that the RIRs are sampled at {fsRIR} Hz.')

    # Invert Sabine's formula to obtain the parameters for the ISM simulator
    if t60 == 0:
        max_order = 0
        e_absorption = 0.5  # <-- arbitrary (unused)
    else:
        try:
            e_absorption, max_order = pra.inverse_sabine(t60, rd, c)
        except ValueError as err:
            print(f'/!\ /!\ PyRoomAcoustics failed to compute the RIRs. Error message:\n"{err}"')
            print('\nLarge room with small T60 do not work well together (too unrealistic).')
            return None, None
    # Create source list
    if nAudio > 0 and nNoise > 0:
        sources =\
            [pra.SoundSource(audioPos3D[ii, :]) for ii in range(nAudio)] +\
            [pra.SoundSource(noisePos3D[ii, :]) for ii in range(nNoise)]
    elif nAudio == 0 and nNoise > 0:
        sources =\
            [pra.SoundSource(noisePos3D[ii, :]) for ii in range(nNoise)]
    elif nAudio > 0 and nNoise == 0:
        sources =\
            [pra.SoundSource(audioPos3D[ii, :]) for ii in range(nAudio)]
    
    # Create acoustic environment (room)
    room = pra.ShoeBox(
        p=np.array(rd),   # room dimensions
        fs=fsRIR,
        t0=0,
        materials=pra.Material(e_absorption),
        max_order=max_order,
        sources=sources,
    )
    room.add_microphone_array(micPos3D.T)  # add microphones

    # Compute RIRs
    print(f'PyRoomAcoustics is computing the {int(nMics * (nAudio + nNoise))} RIRs...')
    room.compute_rir()

    # Post-process RIRs
    if nAudio > 0:
        RIRsAudio = np.zeros((rirLength, nMics, nAudio))
    else:
        RIRsAudio = None
    if nNoise > 0:
        RIRsNoise = np.zeros((rirLength, nMics, nNoise))
    else:
        RIRsNoise = None
    for ii in range(nMics):
        currRIRs = room.rir[ii]
        for jj in range(nAudio):
            if len(currRIRs[jj]) < rirLength:
                RIRsAudio[:, ii, jj] = np.concatenate(
                    (currRIRs[jj], np.zeros(rirLength - len(currRIRs[jj])))
                )
            else:
                RIRsAudio[:, ii, jj] = currRIRs[jj][:rirLength]
                print(f'/!\ Truncated PyRoomAcoustics RIR! (mic #{ii+1} - audio source #{jj+1})')
        for jj in range(nNoise):
            if len(currRIRs[nAudio + jj]) < rirLength:
                RIRsNoise[:, ii, jj] = np.concatenate(
                    (currRIRs[nAudio + jj],
                    np.zeros(rirLength - len(currRIRs[nAudio + jj])))
                )
            else:
                RIRsNoise[:, ii, jj] = currRIRs[nAudio + jj][:rirLength]
                print(f'/!\ Truncated PyRoomAcoustics RIR! (mic #{ii+1} - noise source #{jj+1})')

    if outputRIRplot:
        fig = plot_rirs(RIRsAudio, RIRsNoise)
        fig.savefig(fname + '.png')

    return RIRsAudio, RIRsNoise


def plot_rirs(RIRsAudio, RIRsNoise):
    """
    Output a plot of the RIRs.
    """
    nMics = RIRsAudio.shape[1]
    if RIRsAudio is not None:
        nAudio = RIRsAudio.shape[-1]
    else:
        nAudio = 0
    if RIRsNoise is not None:
        nNoise = RIRsNoise.shape[-1]
    else:
        nNoise = 0
    nCols = 1
    if nAudio > 0 and nNoise > 0:
        nCols = 2
    nRows = np.amax((nAudio, nNoise))

    fig, axes = plt.subplots(nRows, nCols)
    fig.set_size_inches(8.5, 3.5/3 * nMics * nRows)
    for ii in range(nRows):
        for jj in range(nCols):
            if jj == 0:
                if nAudio > 0:
                    if ii < nAudio:
                        RIRsCurr = RIRsAudio[:, :, ii]
                        title = f'Audio source {ii + 1}'
                    else:
                        RIRsCurr = None
                        title = ''
                elif nNoise > 0 and nAudio == 0:
                    if ii < nNoise:
                        RIRsCurr = RIRsNoise[:, :, ii]
                        title = f'Noise source {ii + 1}'
                    else:
                        RIRsCurr = None
                        title = ''
            elif jj == 1:
                if ii < nNoise:
                    RIRsCurr = RIRsNoise[:, :, ii]
                    title = f'Noise source {ii + 1}'
                else:
                    RIRsCurr = None
                    title = ''
            # One subplot per source
            if not isinstance(axes, np.ndarray):
                currAx = axes
            elif len(axes.shape) == 2:
                currAx = axes[ii, jj]
            elif nRows > 1:
                currAx = axes[ii]
            elif nCols > 1:
                currAx = axes[jj]
            # Plot
            if RIRsCurr is not None:
                for kk in range(nMics):
                        delta = np.amax(np.abs(RIRsCurr))
                        currAx.plot(
                            RIRsCurr[:, kk].T - kk * delta,
                            label=f'Mic. #{kk + 1}'
                        )
                currAx.grid()
                currAx.set_title(title)
                currAx.set_xlabel('Samples')
                if ii == 0 and jj == 0:
                    currAx.legend(loc='upper right')
    plt.tight_layout()

    return fig


def get_datetime():
    """
    Returns a formatted date-time string.
    """

    now = datetime.datetime.now()

    if len(str(now.month)) == 1:
        month = '0' + str(now.month)
    else:
        month = now.month
    if len(str(now.day)) == 1:
        day = '0' + str(now.day)
    else:
        day = now.day
    if len(str(now.minute)) == 1:
        minute = '0' + str(now.minute)
    else:
        minute = now.minute
    if len(str(now.hour)) == 1:
        hour = '0' + str(now.hour)
    else:
        hour = now.hour
    if len(str(now.second)) == 1:
        second = '0' + str(now.second)
    else:
        second = now.second
    date = f'{now.year}{month}{day}'
    time = f'{hour}{minute}{second}'

    return date + '_' + time


def load_rirs(path) -> RIRg_GUI:

    RIRGUIobject = pickle.load(gzip.open(path, 'r'))

    return RIRGUIobject