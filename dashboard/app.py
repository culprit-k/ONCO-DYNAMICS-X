"""Streamlit research console for the 2D and 3D dimensionless model."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import streamlit as st

from onco_dynamics_x.hardware.serial_protocol import (
    ControlMessage,
    apply_controls,
    read_control_message,
)
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import simulate
from onco_dynamics_x.state import initial_state

st.set_page_config(page_title="ONCO-DYNAMICS X", layout="wide")
st.title("ONCO-DYNAMICS X | Research simulation console")
st.caption("Educational, dimensionless model only — not a diagnostic, prognostic, or treatment tool.")

with st.sidebar:
    st.header("Physical console")
    serial_port = st.text_input("ESP32 serial port", value="/dev/ttyUSB0")
    if st.button("Read controls from ESP32"):
        try:
            message = read_control_message(serial_port)
            st.session_state["tumour_growth_knob"] = message.tumour_growth
            st.session_state["immune_activity_knob"] = message.immune_activity
            st.session_state["hardware_run"] = message.run
            st.session_state["hardware_reset"] = message.reset
            st.success("Validated console controls received.")
        except (ImportError, ValueError, TimeoutError, ConnectionError) as error:
            st.error(str(error))

    if st.session_state.pop("hardware_reset", False):
        st.session_state["tumour_growth_knob"] = 0.4
        st.session_state["immune_activity_knob"] = 0.5

    st.header("Simulation controls")
    if "tumour_growth_knob" not in st.session_state:
        st.session_state["tumour_growth_knob"] = 0.4
    if "immune_activity_knob" not in st.session_state:
        st.session_state["immune_activity_knob"] = 0.5

    tumour_growth = st.slider(
        "Tumour growth scenario control",
        0.0,
        1.0,
        key="tumour_growth_knob",
        help="Dimensionless scenario control; not a patient measurement.",
    )
    immune_activity = st.slider(
        "Immune activity scenario control",
        0.0,
        1.0,
        key="immune_activity_knob",
        help="Dimensionless scenario control; not a patient measurement.",
    )
    spatial_mode = st.selectbox("Spatial simulation", ["2D slice", "3D tissue volume"])
    if spatial_mode == "3D tissue volume":
        grid_size = st.slider("Voxel grid edge length", 11, 41, 21, step=10)
        steps = st.slider("Time steps", 20, 250, 100, step=20)
    else:
        grid_size = st.slider("Grid size", 21, 81, 41, step=10)
        steps = st.slider("Time steps", 20, 500, 200, step=20)
    dt = st.number_input(
        "Time increment (dimensionless)", value=0.0003, min_value=0.00001, format="%.5f"
    )
    infusion = st.slider("Illustrative drug-field input", 0.0, 0.5, 0.0, 0.01)
    manual_run = st.button("Run simulation", type="primary")
    hardware_run = st.session_state.pop("hardware_run", False)
    run = manual_run or hardware_run

if run:
    dimensions = 3 if spatial_mode == "3D tissue volume" else 2
    parameters = apply_controls(
        ModelParameters(),
        ControlMessage(tumour_growth, immune_activity),
    )
    dx = 2.0 / (grid_size - 1)

    def schedule(_time: float) -> float:
        return float(infusion)

    try:
        _, snapshots = simulate(
            initial_state(grid_size, spatial_dimensions=dimensions),
            parameters,
            dx=dx,
            dt=float(dt),
            steps=steps,
            infusion_schedule=schedule,
            save_every=steps,
        )
        final = snapshots[-1]
        fields = final.as_dict()

        if dimensions == 3:
            st.subheader("Computed cancer-cell field — 3D voxel volume")
            coordinates = np.linspace(-1.0, 1.0, grid_size)
            xx, yy, zz = np.meshgrid(
                coordinates, coordinates, coordinates, indexing="ij"
            )
            volume = go.Figure(
                data=[
                    go.Volume(
                        x=xx.ravel(),
                        y=yy.ravel(),
                        z=zz.ravel(),
                        value=fields["cancer"].ravel(),
                        isomin=float(fields["cancer"].min()),
                        isomax=float(fields["cancer"].max()),
                        opacity=0.12,
                        surface_count=14,
                        caps={"x_show": False, "y_show": False, "z_show": False},
                    )
                ]
            )
            volume.update_layout(
                height=650,
                scene={
                    "xaxis_title": "x (dimensionless)",
                    "yaxis_title": "y (dimensionless)",
                    "zaxis_title": "z (dimensionless)",
                    "aspectmode": "cube",
                },
                margin={"l": 0, "r": 0, "t": 20, "b": 0},
            )
            st.plotly_chart(volume, use_container_width=True)
            st.caption("This is a true 3D numerical grid, but its biological parameters remain illustrative.")
            slice_index = grid_size // 2
            st.subheader("Central tissue slice for the other state fields")
        else:
            st.subheader("Computed cancer-cell field — 2D spatial map")
            figure, axis = plt.subplots(figsize=(6, 4))
            image = axis.imshow(fields["cancer"], origin="lower")
            axis.set_xlabel("Grid x")
            axis.set_ylabel("Grid y")
            figure.colorbar(image, ax=axis, label="Dimensionless cancer field")
            st.pyplot(figure, clear_figure=True)
            plt.close(figure)

            st.subheader("3D surface rendering of the same 2D field")
            st.caption("Surface height represents the computed 2D field; this is not a 3D tissue solver.")
            surface = go.Figure(
                data=[go.Surface(z=fields["cancer"], colorscale="Viridis")]
            )
            surface.update_layout(
                height=520,
                scene={
                    "xaxis_title": "Grid x",
                    "yaxis_title": "Grid y",
                    "zaxis_title": "Dimensionless cancer field",
                },
                margin={"l": 0, "r": 0, "t": 20, "b": 0},
            )
            st.plotly_chart(surface, use_container_width=True)
            slice_index = None

        st.subheader("Other computed spatial fields")
        columns = st.columns(3)
        display_fields = [
            ("Immune-cell density (I)", "immune"),
            ("Oxygen field (O)", "oxygen"),
            ("Fibroblast field (F)", "fibroblast"),
            ("Drug field (D)", "drug"),
            ("Vasculature proxy (V)", "vasculature"),
            ("Immune-evasion proxy (E)", "evasion"),
        ]
        for index, (title, key) in enumerate(display_fields):
            field = fields[key]
            image_field = field if dimensions == 2 else field[:, :, slice_index]
            with columns[index % 3]:
                figure, axis = plt.subplots(figsize=(4, 3))
                image = axis.imshow(image_field, origin="lower")
                axis.set_title(title)
                axis.set_xticks([])
                axis.set_yticks([])
                figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
                st.pyplot(figure, clear_figure=True)
                plt.close(figure)

        st.subheader("Mean state at final simulation time")
        st.json({name: round(float(field.mean()), 5) for name, field in fields.items()})
    except ValueError as error:
        st.error(str(error))
else:
    st.info("Choose a 2D slice or 3D tissue volume, then run the simulation or connect the ESP32.")
