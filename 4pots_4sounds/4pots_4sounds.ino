#include "heartbeats.h"  // Contains newborn, normal, flatline and their *_len variables

const int potPins[4] = {33, 32, 35, 34};  // Three potentiometer pins
const int dacPin = 25;

// Sounds pointers and lengths (4 sounds)
const uint8_t* sounds[4] = {
  newborn,
  normal,
  congestive,
  flatline
};

const size_t lengths[4] = {
  newborn_len,
  normal_len,
  congestive_len,
  flatline_len
};

// Positions for playback of each sound
size_t pos[4] = {0, 0, 0, 0};

void setup() {
  analogReadResolution(12);  // ADC range: 0–4095
  Serial.begin(115200);
}

void loop() {
  float mixedSample = 0.0;

  for (int i = 0; i < 4; i++) {
    int potValue = analogRead(potPins[i]);

    // Threshold to decide if sound i is active (e.g. potValue > 2048)
    bool active = (potValue > 2048);

    if (active) {
      // Loop position reset
      if (pos[i] >= lengths[i]) pos[i] = 0;

      // Read raw sample from PROGMEM
      uint8_t rawSample = pgm_read_byte_near(sounds[i] + pos[i]);
      pos[i]++;

      // Convert unsigned 8-bit sample (0-255) to float centered [-1,1]
      float centered = ((int)rawSample - 128) / 128.0;

      // Add to mixed sample
      mixedSample += centered;
    }
  }

  // Prevent clipping by limiting mixedSample to [-1,1]
  mixedSample = constrain(mixedSample, -1.0, 1.0);

  // Convert back to unsigned 8-bit (0-255) for DAC output
  int outputValue = (int)(mixedSample * 127.0 + 128.0);
  outputValue = constrain(outputValue, 0, 255);

  dacWrite(dacPin, outputValue);

  delayMicroseconds(125);  // ~8kHz playback rate
}
