"""Optional Streamlit interface for the dimensionless numerical demonstration."""
from __future__ import annotations

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import streamlit as st

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import simulate
from onco_dynamics_x.state import initial_state

st.set_page_config(page_title="ONCO-DYNAMICS X", layout="wide")
st.title("ONCO-DYNAMICS X | Research simulation console")
st.caption("Educational, dimensionless model only — not a diagnostic, prognostic, or treatment tool.")

with st.sidebar:
    st.header("Simulation controls")
    grid_size = st.slider("Grid size", 21, 81, 41, step=10)
    steps = st.slider("Time steps", 20, 500, 200, step=20)
    dt = st.number_input(
        "Time increment (dimensionless)", value=0.0003, min_value=0.00001, format="%.5f"
    )
    infusion = st.slider("Illustrative drug-field input", 0.0, 0.5, 0.0, 0.01)
    run = st.button("Run simulation", type="primary")

if run:
    params = ModelParameters()
    dx = 2.0 / (grid_size - 1)

    def schedule(_time: float) -> float:
        return float(infusion)

    try:
        _, snapshots = simulate(
            initial_state(grid_size),
            params,
            dx=dx,
            dt=float(dt),
            steps=steps,
            infusion_schedule=schedule,
            save_every=steps,
        )
        final = snapshots[-1]
        fields = final.as_dict()

        st.subheader("Computed tumour field — 2D spatial map")
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
            with columns[index % 3]:
                figure, axis = plt.subplots(figsize=(4, 3))
                image = axis.imshow(fields[key], origin="lower")
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
    st.info("Adjust the controls and press **Run simulation** to calculate a new state.")
