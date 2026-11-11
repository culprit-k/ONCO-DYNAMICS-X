"""Streamlit research console for the 2D and 3D dimensionless model."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import streamlit as st

from onco_dynamics_x.cellular_simulation.agent_based import (
    CANCER,
    FIBROBLAST,
    IMMUNE,
    initialize_agents,
)
from onco_dynamics_x.cellular_simulation.coupling import simulate_coupled
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
    coupled_mode = st.checkbox(
        "Enable bidirectional PDE–agent coupling",
        value=False,
        help="Agents sample PDE fields and feed phenotype density proxies back to the grid.",
    )
    if coupled_mode:
        coupling_strength = st.slider("PDE–agent coupling strength", 0.0, 1.0, 0.15, 0.05)
        agents_per_unit = st.number_input(
            "Agents per grid node for density=1", min_value=0.1, value=4.0, step=0.5
        )
        initial_cancer_agents = st.number_input("Initial cancer agents", 0, 1000, 80, 10)
        initial_immune_agents = st.number_input("Initial immune agents", 0, 1000, 50, 10)
        initial_fibroblast_agents = st.number_input("Initial fibroblast agents", 0, 1000, 35, 5)
        deposition_method = st.selectbox(
            "Agent-to-grid deposition",
            ["nearest_grid", "cloud_in_cell"],
            help=(
                "Cloud-in-cell distributes each agent across adjacent nodes with "
                "multilinear weights; both methods are numerical options, not calibrated biology."
            ),
        )
    else:
        coupling_strength = 0.0
        agents_per_unit = 4.0
        initial_cancer_agents = 0
        initial_immune_agents = 0
        initial_fibroblast_agents = 0
        deposition_method = "nearest_grid"
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
        final_population = None
        if coupled_mode:
            initial_population = initialize_agents(
                cancer_count=int(initial_cancer_agents),
                immune_count=int(initial_immune_agents),
                fibroblast_count=int(initial_fibroblast_agents),
                seed=42,
                spatial_dimensions=dimensions,
            )
            _, snapshots, population_snapshots = simulate_coupled(
                initial_state(grid_size, spatial_dimensions=dimensions),
                initial_population,
                parameters,
                dx=dx,
                dt=float(dt),
                steps=steps,
                coupling_strength=float(coupling_strength),
                agents_per_unit=float(agents_per_unit),
                deposition_method=deposition_method,
                infusion_schedule=schedule,
                save_every=steps,
                seed=42,
                maximum_agents=5000,
            )
            final_population = population_snapshots[-1]
        else:
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
            if final_population is not None and final_population.z is not None:
                for cell_type, label in (
                    (CANCER, "Cancer agents"),
                    (IMMUNE, "Immune agents"),
                    (FIBROBLAST, "Fibroblast agents"),
                ):
                    mask = final_population.cell_type == cell_type
                    volume.add_trace(
                        go.Scatter3d(
                            x=final_population.x[mask],
                            y=final_population.y[mask],
                            z=final_population.z[mask],
                            mode="markers",
                            name=label,
                            marker={"size": 3},
                        )
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
            if final_population is not None:
                columns_xy = (final_population.x + 1.0) * 0.5 * (grid_size - 1)
                rows_xy = (final_population.y + 1.0) * 0.5 * (grid_size - 1)
                axis.scatter(columns_xy, rows_xy, s=10, label="Cell agents")
                axis.legend(fontsize=7)
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

        if final_population is not None:
            st.subheader("Discrete cell-agent population")
            st.json(
                {
                    "total_agents": len(final_population.x),
                    "cancer_agents": final_population.count(CANCER),
                    "immune_agents": final_population.count(IMMUNE),
                    "fibroblast_agents": final_population.count(FIBROBLAST),
                    "coupling_strength": float(coupling_strength),
                    "agents_per_unit": float(agents_per_unit),
                    "biological_calibration": "not established",
                }
            )
        st.subheader("Mean state at final simulation time")
        st.json({name: round(float(field.mean()), 5) for name, field in fields.items()})
    except ValueError as error:
        st.error(str(error))
else:
    st.info("Choose a 2D slice or 3D tissue volume, then run the simulation or connect the ESP32.")
