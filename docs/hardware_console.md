# Physical demonstration: breast phantom and research console

## Demonstration concept
A breast anatomy phantom or cutaway model provides context. A console with rotary encoders and buttons changes model parameters. A computer or Raspberry Pi runs the solver. A main display renders computed tumour, immune, oxygen, stromal, vascular and simulated drug fields; a secondary display shows equations, parameter values, traces and experiment logs.

The phantom is a teaching prop. It contains no real specimen and does not sense, screen or diagnose a person.

## Suggested electronics
ESP32 development board; two rotary encoders or potentiometers; run/pause and reset buttons; status LEDs or small OLED; USB serial connection to host; optional Raspberry Pi for the console.

## Communication example
```json
{"event":"parameter","name":"immune_activity","value":0.65}
{"event":"command","name":"reset"}
```
The host must validate allowed parameter names and ranges before changing simulation configuration. Values must be visibly labelled as illustrative model controls.

## Safety
No clinical monitors, hospital networks, infusion pumps or patient devices. No drug dispensing, stimulation or patient contact. Use low-voltage USB-powered electronics. Keep serial control local. Hardware input must never be presented as a patient treatment selector.

## Stages
1. Bench-test controls using serial output only.
2. Implement host parsing and strict validation.
3. Connect to a simulation parameter panel.
4. Add event logs and reset-to-default.
5. Integrate phantom and enclosure after software behavior is stable.
