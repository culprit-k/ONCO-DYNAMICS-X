# ESP32 console wiring

All wiring is low-voltage and powered from the ESP32 development board. Do not connect electrodes or sensors to a person. The console controls only the educational simulation.

## Pin map

| Component | ESP32 pin | Connection |
|---|---|---|
| Tumour-scenario potentiometer (10 kΩ) | GPIO 34 (ADC1) | Wiper to GPIO 34; outer legs to 3V3 and GND |
| Immune-scenario potentiometer (10 kΩ) | GPIO 35 (ADC1) | Wiper to GPIO 35; outer legs to 3V3 and GND |
| Run pushbutton | GPIO 25 | One side to GPIO 25; other side to GND |
| Reset pushbutton | GPIO 26 | One side to GPIO 26; other side to GND |
| Board power | USB | Connect ESP32 to the computer running the dashboard |

Buttons use the firmware's `INPUT_PULLUP` configuration. Do not connect a potentiometer end to 5 V: ESP32 GPIO inputs are 3.3 V logic.

## Assembly
1. Mount the two potentiometers and two momentary pushbuttons in the console enclosure.
2. Connect the wipers and button contacts using short insulated jumpers.
3. Connect the ESP32 by USB to the computer.
4. Open `hardware/esp32/ONCOConsole.ino` in Arduino IDE, select the correct ESP32 board and serial port, then upload.
5. Close Arduino Serial Monitor before launching Streamlit because only one process can normally own the serial port.
6. Start the dashboard locally and enter the serial device name (for example, `COM5`, `/dev/ttyUSB0`, or `/dev/cu.usbserial-...`).
7. Use “Read controls from ESP32”; confirm the received values update the two scenario sliders. Press Run to calculate the selected 2D or 3D simulation.

## Safety and limitations
- Use only USB/3.3 V low-voltage electronics.
- Insulate exposed conductors and prevent shorts; do not power motors or other high-current loads from GPIO.
- The knobs are arbitrary normalized simulation parameters, not clinical measurements or treatment controls.
- Reset restores the dashboard's default scenario controls; it does not reset or modify a patient/device.
- The ESP32 does not collect data from a person, patient or biological specimen.
