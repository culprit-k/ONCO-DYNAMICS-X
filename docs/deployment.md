# Dashboard deployment

## Local run
1. Install the dashboard extra: `pip install -e ".[dashboard,hardware]"`.
2. Launch: `streamlit run dashboard/app.py`.
3. The browser interface can run the dimensionless model without an ESP32.

## Streamlit Community Cloud
1. Sign in to Streamlit Community Cloud and connect the GitHub account with access to this repository.
2. Create an app from the repository `culprit-k/ONCO-DYNAMICS-X`, branch `main`, entry point `dashboard/app.py`.
3. Configure the dependency file as `requirements-dashboard.txt` if the deployment UI offers a dependency-file selector.
4. Deploy and inspect the app logs. The deployment is not created by this repository commit; the repository owner must authorize and configure it.

## Hardware note
A cloud-hosted dashboard cannot access a USB serial port attached to a separate laptop or local ESP32. For hardware demonstration, run the dashboard locally on the computer connected to the ESP32 and install the hardware extra. Do not expose an unauthenticated serial-to-network bridge.

## Troubleshooting
- If imports fail, confirm the dashboard dependency file includes the editable local package and all dashboard dependencies.
- If simulation rejects a timestep, lower dt or reduce the grid size; the error reports the explicit diffusion stability limit.
- If serial controls fail, check the device port, USB cable, board power, 115200 baud, and that no other process has opened the port.
