// ONCO-DYNAMICS X educational console — ESP32 development board.
// USB serial only. These controls alter simulation parameters, not medical inputs.

#include <Arduino.h>

constexpr uint8_t TUMOUR_KNOB_PIN = 34;   // ADC1 input
constexpr uint8_t IMMUNE_KNOB_PIN = 35;   // ADC1 input
constexpr uint8_t RUN_BUTTON_PIN = 25;    // button to GND, INPUT_PULLUP
constexpr uint8_t RESET_BUTTON_PIN = 26;  // button to GND, INPUT_PULLUP
constexpr unsigned long SEND_INTERVAL_MS = 100;

bool previousRunPressed = false;
bool previousResetPressed = false;
bool pendingRunPulse = false;
bool pendingResetPulse = false;
unsigned long lastSendMs = 0;

float readNormalizedKnob(uint8_t pin) {
  return constrain(static_cast<float>(analogRead(pin)) / 4095.0f, 0.0f, 1.0f);
}

void setup() {
  Serial.begin(115200);
  analogReadResolution(12);
  pinMode(RUN_BUTTON_PIN, INPUT_PULLUP);
  pinMode(RESET_BUTTON_PIN, INPUT_PULLUP);
}

void loop() {
  const bool runPressed = (digitalRead(RUN_BUTTON_PIN) == LOW);
  const bool resetPressed = (digitalRead(RESET_BUTTON_PIN) == LOW);
  pendingRunPulse = pendingRunPulse || (runPressed && !previousRunPressed);
  pendingResetPulse = pendingResetPulse || (resetPressed && !previousResetPressed);
  previousRunPressed = runPressed;
  previousResetPressed = resetPressed;

  const unsigned long now = millis();
  if (now - lastSendMs >= SEND_INTERVAL_MS) {
    lastSendMs = now;
    const float tumourGrowth = readNormalizedKnob(TUMOUR_KNOB_PIN);
    const float immuneActivity = readNormalizedKnob(IMMUNE_KNOB_PIN);

    Serial.printf(
      "{\"event\":\"controls\",\"tumour_growth\":%.4f,"
      "\"immune_activity\":%.4f,\"run\":%s,\"reset\":%s}\n",
      tumourGrowth,
      immuneActivity,
      pendingRunPulse ? "true" : "false",
      pendingResetPulse ? "true" : "false"
    );
    pendingRunPulse = false;
    pendingResetPulse = false;
  }
  delay(5);
}
