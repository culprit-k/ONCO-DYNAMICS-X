# Physical prototype build guide

## Demonstrator layout

The physical exhibit has two independent components:
1. **Breast phantom:** a removable soft-tissue teaching model with embedded removable lesion inserts and an open rear access.
2. **Research console:** a low-voltage ESP32 control panel with two rotary knobs and Run/Reset buttons, connected by USB to the laptop running ONCO-DYNAMICS X.

The phantom is a static educational prop. It is not connected to the simulation as a sensor and does not represent real patient anatomy, tumour size, stage or imaging. Simulation parameters are controlled by the ESP32 and dashboard, not measured from the phantom.

## Bill of materials

### Console electronics
- ESP32 development board (USB-powered).
- Two 10 kΩ linear potentiometers.
- Two normally-open momentary pushbuttons.
- Breadboard or small prototyping PCB.
- Jumper wires, USB data cable and insulated enclosure.

### Physical model
- 3D-printed phantom shell and stand ring from `hardware/cad/breast_phantom.scad`.
- Removable tumour inserts printed separately.
- Soft silicone casting material or another suitable non-biological teaching fill, mixed and cured strictly according to its manufacturer's safety data and instructions.
- Optional coloured silicone pigments and labels.
- Gloves, eye protection, mixing cups and disposable work surface for casting.

Do not use human tissue, clinical specimens, or biologically active materials. Confirm printer/material compatibility and curing ventilation requirements before fabrication.

## CAD export
1. Open `hardware/cad/breast_phantom.scad` in OpenSCAD.
2. Set `part = "shell"`, render (F6), then export STL.
3. Repeat with `small_insert`, `medium_insert`, `large_insert` and `stand_ring`.
4. Inspect wall thickness, opening and fit in the slicer before printing. The model is parametric and may need scaling for your printer and casting process.
5. Print the shell as a mould/form or use it as a support shell. Test a small material coupon first; silicone adhesion/release depends on the chosen formulation and surface finish.

## Console firmware and dashboard
Follow [the ESP32 wiring map](../hardware/wiring.md). Upload `hardware/esp32/ONCOConsole.ino`, connect USB, install the hardware extra, and run:

```bash
pip install -e ".[dashboard,hardware]"
streamlit run dashboard/app.py
```

Select the local serial port in the sidebar and use “Read controls from ESP32”. A hosted dashboard cannot access USB hardware connected to a different local computer; use the local dashboard for a live exhibit.

## Demonstration sequence
1. Identify the simplified breast shell and removable lesion inserts.
2. Show the eight-field mathematical model and explain which fields are proxies.
3. Rotate the growth and immune-scenario knobs; read the normalized values in the dashboard.
4. Run a 2D slice simulation and then a small 3D voxel simulation.
5. Compare the computed cancer, immune, oxygen and stromal fields.
6. Explain that all displayed values are dimensionless and synthetic, and that the phantom is not a diagnostic model.

## Electrical safety and scope
Use USB and 3.3 V logic only. Do not connect the circuit to a person or patient. Do not use medical electrodes, sensors or invasive hardware. Use the model as a classroom demonstration, not as a medical device.
