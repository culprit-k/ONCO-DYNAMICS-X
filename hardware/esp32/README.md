# ESP32 console firmware scaffold

Begin by printing button and knob readings over USB serial. Do not connect to a patient device or clinical equipment.

Suggested inputs:
- Encoder A: tumour-growth scenario control.
- Encoder B: immune-activity scenario control.
- Button 1: reset.
- Button 2: run/pause request.

Before adding host control, document message types, allowed parameter ranges, invalid-message behavior and power-off behavior.
