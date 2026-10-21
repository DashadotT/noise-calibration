/*
  ============================================================
  BANHA_INMP441_CALIBRATION_v1.ino
  Standalone INMP441 calibration / characterization firmware

  PURPOSE
  -------
  This is NOT the final BANHA firmware.

  It is a temporary measurement tool to establish whether the
  INMP441 + ESP32 digital signal path is behaving correctly
  before we change BANHA's calibration constant.

  SIGNAL PATH
  -----------
  INMP441
      -> I2S 48 kHz / 32-bit frame
      -> 24-bit microphone sample (>> 8)
      -> INMP441 frequency-response equalization
      -> A-weighting
      -> RMS over a 5-second measurement window
      -> datasheet-based SPL conversion

  HARDWARE
  --------
  INMP441:
    WS  = GPIO 16
    SCK = GPIO 17
    SD  = GPIO 5

  IMPORTANT
  ---------
  The INMP441 datasheet specifies:
    sensitivity = -26 dBFS typical at 94 dB SPL, 1 kHz
    valid I2S data = 24-bit two's complement

  The datasheet's -26 dBFS peak reference corresponds to:
    420426 peak digital counts

  For RMS conversion of a sine-wave reference:
    RMS reference = 420426 / sqrt(2)

  Equivalently:
    SPL = 94 + 3.0103 + 20*log10(RMS / 420426)

  A-weighting is normalized to 0 dB at 1 kHz.

  IMPORTANT LIMITATION
  --------------------
  This is NOT traceable acoustic calibration. The INMP441 has
  individual sensitivity tolerance. Without a 94 dB / 1 kHz
  acoustic calibrator, the datasheet calculation is our
  theoretical starting point. We will later validate it against
  an external reference at practical sound levels.

  Based on:
    Ivan Kostoski - ESP32-I2S-SLM
    https://hackaday.io/project/166867-esp32-i2s-slm

  Ivan's project uses 48 kHz because the IIR filters are designed
  for that sample rate. His INMP441 equalizer and A-weighting
  coefficients are reproduced below for this standalone test.

  ============================================================
*/

#include <Arduino.h>
#include <driver/i2s.h>
#include <math.h>

// ============================================================
// BLOCK RESULT TYPE
// Must be declared before processBlock() so Arduino's automatic
// function-prototype generation can resolve the return type.
// ============================================================

struct BlockResult {
  double rawSum = 0.0;
  double rawSumSquares = 0.0;
  double eqSumSquares = 0.0;
  double weightedSumSquares = 0.0;
};

// Explicit prototype prevents the Arduino preprocessor from
// generating an earlier prototype with an unknown return type.
BlockResult processBlock(float* samples);

// ============================================================
// I2S / MEASUREMENT CONFIGURATION
// ============================================================

#define I2S_PORT I2S_NUM_0

#define I2S_WS  16
#define I2S_SCK 17
#define I2S_SD   5

#define SAMPLE_RATE       48000
#define SAMPLE_BITS       32
#define MIC_BITS           24

#define BLOCK_MS           125
#define BLOCK_SAMPLES      (SAMPLE_RATE * BLOCK_MS / 1000)   // 6000

#define SETTLE_SECONDS     2
#define MEASURE_SECONDS    5
#define MEASURE_SAMPLES    (SAMPLE_RATE * MEASURE_SECONDS)   // 240000
#define TOTAL_BLOCKS       (MEASURE_SAMPLES / BLOCK_SAMPLES)  // 40

#define DMA_BUF_COUNT      8
#define DMA_BUF_LEN        512

// ============================================================
// DATASHEET-BASED MICROPHONE REFERENCE
// ============================================================

// INMP441 typical sensitivity at 1 kHz, 94 dB SPL.
static constexpr double MIC_SENSITIVITY_DBFS = -26.0;
static constexpr double MIC_REFERENCE_DB_SPL  = 94.0;

// 24-bit signed full-scale peak = 2^23 - 1.
// At -26 dBFS the datasheet gives approximately 420426.
static constexpr double MIC_REF_PEAK =
    ((double)((1UL << (MIC_BITS - 1)) - 1UL)) *
    pow(10.0, MIC_SENSITIVITY_DBFS / 20.0);

// RMS-to-peak correction for a sine-wave reference.
static constexpr double RMS_TO_PEAK_DB = 3.01029995664;

// Optional empirical correction.
// KEEP THIS AT ZERO FOR v1.
// We will determine this only after the signal chain is verified
// and external-reference tests have been completed.
static constexpr double CALIBRATION_OFFSET_DB = 0.0;

// ============================================================
// SIMPLE SECOND-ORDER-SECTION IIR FILTER
// ============================================================

struct SOSCoefficients {
  float b1;
  float b2;
  float a1;
  float a2;
};

struct SOSState {
  float w0 = 0.0f;
  float w1 = 0.0f;
};

class SOSFilter {
public:
  const SOSCoefficients* sections;
  size_t count;
  float gain;
  SOSState* state;

  SOSFilter(const SOSCoefficients* s, size_t n, float g)
    : sections(s), count(n), gain(g), state(nullptr) {}

  bool begin() {
    state = new SOSState[count];
    return state != nullptr;
  }

  void reset() {
    if (!state) return;
    for (size_t i = 0; i < count; i++) {
      state[i].w0 = 0.0f;
      state[i].w1 = 0.0f;
    }
  }

  float process(float x) {
    float y = x;

    for (size_t i = 0; i < count; i++) {
      const SOSCoefficients& c = sections[i];
      SOSState& s = state[i];

      // Same sign convention as Ivan Kostoski's SOS filter:
      // w0 = x + a1*w0(previous) + a2*w1(previous)
      float w = y + c.a1 * s.w0 + c.a2 * s.w1;

      // b0 is 1.0. b1/b2 are the stored SOS coefficients.
      y = w + c.b1 * s.w0 + c.b2 * s.w1;

      s.w1 = s.w0;
      s.w0 = w;
    }

    return y * gain;
  }

  ~SOSFilter() {
    delete[] state;
  }
};

// ============================================================
// INMP441 EQUALIZER
// Fs = 48 kHz
//
// From Ivan Kostoski's ESP32-I2S-SLM source.
// ============================================================

static const SOSCoefficients INMP441_SOS[] = {
  {
    -1.986920458344451f,
    +0.986963226946616f,
    +1.995178510504166f,
    -0.995184322194091f
  }
};

static SOSFilter inmp441EQ(
  INMP441_SOS,
  1,
  1.00197834654696f
);

// ============================================================
// A-WEIGHTING
// Fs = 48 kHz
//
// From Ivan Kostoski's ESP32-I2S-SLM source.
// Implemented as three SOS / biquad sections.
// ============================================================

static const SOSCoefficients A_WEIGHT_SOS[] = {
  {
    -2.00026996133106f,
    +1.00027056142719f,
    -1.060868438509278f,
    -0.163987445885926f
  },
  {
    +4.35912384203144f,
    +3.09120265783884f,
    +1.208419926363593f,
    -0.273166998428332f
  },
  {
    -0.70930303489759f,
    -0.29071868393580f,
    +1.982242159753048f,
    -0.982298594928989f
  }
};

static SOSFilter aWeight(
  A_WEIGHT_SOS,
  3,
  0.169994948147430f
);

// ============================================================
// RAW I2S BUFFER
// ============================================================

static int32_t i2sBuffer[DMA_BUF_LEN];

// ============================================================
// I2S INITIALIZATION
// ============================================================

void setupI2S() {

  i2s_config_t config = {};

  config.mode =
      (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_RX);

  config.sample_rate = SAMPLE_RATE;

  config.bits_per_sample =
      I2S_BITS_PER_SAMPLE_32BIT;

  config.channel_format =
      I2S_CHANNEL_FMT_ONLY_LEFT;

  config.communication_format =
      (i2s_comm_format_t)(
        I2S_COMM_FORMAT_I2S |
        I2S_COMM_FORMAT_I2S_MSB
      );

  config.intr_alloc_flags = ESP_INTR_FLAG_LEVEL1;

  config.dma_buf_count = DMA_BUF_COUNT;
  config.dma_buf_len   = DMA_BUF_LEN;

  config.use_apll = true;
  config.tx_desc_auto_clear = false;
  config.fixed_mclk = 0;

  esp_err_t err = i2s_driver_install(
    I2S_PORT,
    &config,
    0,
    nullptr
  );

  if (err != ESP_OK) {
    Serial.printf(
      "ERROR: i2s_driver_install = %d\n",
      err
    );
    while (true) delay(1000);
  }

  i2s_pin_config_t pins = {};

  pins.bck_io_num   = I2S_SCK;
  pins.ws_io_num    = I2S_WS;
  pins.data_out_num = I2S_PIN_NO_CHANGE;
  pins.data_in_num  = I2S_SD;

  err = i2s_set_pin(I2S_PORT, &pins);

  if (err != ESP_OK) {
    Serial.printf(
      "ERROR: i2s_set_pin = %d\n",
      err
    );
    while (true) delay(1000);
  }

  // Make sure DMA starts from a clean state.
  i2s_zero_dma_buffer(I2S_PORT);
}

// ============================================================
// READ ONE MICROPHONE SAMPLE
// ============================================================

bool readMicSample(int32_t& sample24) {

  size_t bytesRead = 0;

  esp_err_t err = i2s_read(
    I2S_PORT,
    i2sBuffer,
    sizeof(i2sBuffer),
    &bytesRead,
    portMAX_DELAY
  );

  if (err != ESP_OK || bytesRead < sizeof(int32_t)) {
    return false;
  }

  // This function reads a block because the I2S driver is DMA based.
  // Keep the first sample for the block-level reader below.
  sample24 = i2sBuffer[0] >> 8;
  return true;
}

// ============================================================
// BLOCK READER
//
// Reads exactly BLOCK_SAMPLES microphone samples.
// Returns false if I2S fails.
// ============================================================

bool readBlock(
  float* input,
  size_t sampleCount,
  double& rawSum,
  double& rawSumSquares
) {

  size_t collected = 0;

  while (collected < sampleCount) {

    size_t bytesRead = 0;

    esp_err_t err = i2s_read(
      I2S_PORT,
      i2sBuffer,
      sizeof(i2sBuffer),
      &bytesRead,
      portMAX_DELAY
    );

    if (err != ESP_OK) {
      return false;
    }

    size_t got = bytesRead / sizeof(int32_t);

    if (got == 0) {
      continue;
    }

    if (got > sampleCount - collected) {
      got = sampleCount - collected;
    }

    for (size_t i = 0; i < got; i++) {

      // INMP441 provides 24-bit two's-complement data.
      // ESP32 I2S receives it in a 32-bit word.
      int32_t sample24 = i2sBuffer[i] >> 8;

      float x = (float)sample24;

      input[collected++] = x;

      rawSum += (double)x;
      rawSumSquares += (double)x * (double)x;
    }
  }

  return true;
}

// ============================================================
// PROCESS ONE BLOCK
// ============================================================

BlockResult processBlock(float* samples) {

  BlockResult result;

  // Raw statistics before digital equalization.
  for (size_t i = 0; i < BLOCK_SAMPLES; i++) {
    double x = samples[i];

    result.rawSum += x;
    result.rawSumSquares += x * x;
  }

  // Equalization then A-weighting.
  for (size_t i = 0; i < BLOCK_SAMPLES; i++) {

    float eq = inmp441EQ.process(samples[i]);

    result.eqSumSquares +=
      (double)eq * (double)eq;

    float weighted = aWeight.process(eq);

    result.weightedSumSquares +=
      (double)weighted * (double)weighted;
  }

  return result;
}

// ============================================================
// PRINT HEADER
// ============================================================

void printHeader() {

  Serial.println();
  Serial.println();
  Serial.println("============================================================");
  Serial.println(" BANHA INMP441 CALIBRATION v1");
  Serial.println("============================================================");
  Serial.println("Standalone characterization firmware");
  Serial.println();
  Serial.println("I2S:");
  Serial.printf("  WS  = GPIO %d\n", I2S_WS);
  Serial.printf("  SCK = GPIO %d\n", I2S_SCK);
  Serial.printf("  SD  = GPIO %d\n", I2S_SD);
  Serial.printf("  Rate = %d Hz\n", SAMPLE_RATE);
  Serial.printf("  Word = %d bits\n", SAMPLE_BITS);
  Serial.printf("  Valid microphone bits = %d\n", MIC_BITS);
  Serial.println();
  Serial.println("Measurement:");
  Serial.printf("  Settle = %d s\n", SETTLE_SECONDS);
  Serial.printf("  Window = %d s\n", MEASURE_SECONDS);
  Serial.printf("  Samples/window = %d\n", MEASURE_SAMPLES);
  Serial.printf("  Block = %d samples (%d ms)\n",
                BLOCK_SAMPLES, BLOCK_MS);
  Serial.printf("  Blocks/window = %d\n", TOTAL_BLOCKS);
  Serial.println();
  Serial.println("INMP441 datasheet reference:");
  Serial.printf("  Sensitivity = %.1f dBFS @ %.1f dB SPL, 1 kHz\n",
                MIC_SENSITIVITY_DBFS,
                MIC_REFERENCE_DB_SPL);
  Serial.printf("  Reference peak = %.0f counts\n",
                MIC_REF_PEAK);
  Serial.printf("  RMS correction = +%.6f dB\n",
                RMS_TO_PEAK_DB);
  Serial.printf("  Empirical calibration offset = %.3f dB\n",
                CALIBRATION_OFFSET_DB);
  Serial.println();
  Serial.println("Theoretical SPL formula:");
  Serial.println("  SPL = 94 + 3.0103 + 20*log10(RMS/420426)");
  Serial.println();
  Serial.println("IMPORTANT:");
  Serial.println("  This is not traceable calibration.");
  Serial.println("  We will validate it against an external reference later.");
  Serial.println("============================================================");
}

// ============================================================
// SETTLE FILTER STATES
// ============================================================

void settleFilters() {

  Serial.println();
  Serial.printf(
    "Settling filters for %d seconds...\n",
    SETTLE_SECONDS
  );

  const uint32_t settleSamples =
    SAMPLE_RATE * SETTLE_SECONDS;

  uint32_t processed = 0;

  while (processed < settleSamples) {

    size_t want = BLOCK_SAMPLES;

    if (settleSamples - processed < want) {
      want = settleSamples - processed;
    }

    // For simplicity, the settle period uses a local buffer.
    static float settleBuffer[BLOCK_SAMPLES];

    size_t collected = 0;

    while (collected < want) {

      size_t bytesRead = 0;

      if (i2s_read(
            I2S_PORT,
            i2sBuffer,
            sizeof(i2sBuffer),
            &bytesRead,
            portMAX_DELAY
          ) != ESP_OK) {
        continue;
      }

      size_t got =
        bytesRead / sizeof(int32_t);

      if (got > want - collected) {
        got = want - collected;
      }

      for (size_t i = 0; i < got; i++) {
        settleBuffer[collected++] =
          (float)(i2sBuffer[i] >> 8);
      }
    }

    for (size_t i = 0; i < want; i++) {
      float eq = inmp441EQ.process(settleBuffer[i]);
      (void)aWeight.process(eq);
    }

    processed += want;
  }

  Serial.println("Filter settling complete.");
}

// ============================================================
// ONE COMPLETE 5-SECOND MEASUREMENT
// ============================================================

void runMeasurement() {

  static float block[BLOCK_SAMPLES];

  double totalRawSum = 0.0;
  double totalRawSumSquares = 0.0;

  double totalEqSumSquares = 0.0;
  double totalWeightedSumSquares = 0.0;

  Serial.println();
  Serial.println("------------------------------------------------------------");
  Serial.println("Starting 5-second measurement...");
  Serial.println("Keep the microphone in a stable acoustic environment.");
  Serial.println("------------------------------------------------------------");

  uint32_t measurementStart = millis();

  for (int blockIndex = 0;
       blockIndex < TOTAL_BLOCKS;
       blockIndex++) {

    double rawSum = 0.0;
    double rawSumSquares = 0.0;

    // Read exactly 125 ms of audio.
    size_t collected = 0;

    while (collected < BLOCK_SAMPLES) {

      size_t bytesRead = 0;

      if (i2s_read(
            I2S_PORT,
            i2sBuffer,
            sizeof(i2sBuffer),
            &bytesRead,
            portMAX_DELAY
          ) != ESP_OK) {
        continue;
      }

      size_t got =
        bytesRead / sizeof(int32_t);

      if (got > BLOCK_SAMPLES - collected) {
        got = BLOCK_SAMPLES - collected;
      }

      for (size_t i = 0; i < got; i++) {

        int32_t sample24 =
          i2sBuffer[i] >> 8;

        float x = (float)sample24;

        block[collected++] = x;

        rawSum += (double)x;
        rawSumSquares +=
          (double)x * (double)x;
      }
    }

    // Filter this block.
    double eqSumSquares = 0.0;
    double weightedSumSquares = 0.0;

    for (size_t i = 0; i < BLOCK_SAMPLES; i++) {

      float eq =
        inmp441EQ.process(block[i]);

      eqSumSquares +=
        (double)eq * (double)eq;

      float weighted =
        aWeight.process(eq);

      weightedSumSquares +=
        (double)weighted * (double)weighted;
    }

    totalRawSum += rawSum;
    totalRawSumSquares += rawSumSquares;

    totalEqSumSquares += eqSumSquares;
    totalWeightedSumSquares += weightedSumSquares;

    // Show progress every 125 ms.
    double blockWeightedRMS =
      sqrt(weightedSumSquares / BLOCK_SAMPLES);

    double blockDBA =
      MIC_REFERENCE_DB_SPL +
      RMS_TO_PEAK_DB +
      20.0 * log10(
        max(blockWeightedRMS, 1e-12) /
        MIC_REF_PEAK
      ) +
      CALIBRATION_OFFSET_DB;

    Serial.printf(
      "Block %02d/%02d | A-weighted RMS: %10.2f | "
      "theoretical: %7.2f dBA\n",
      blockIndex + 1,
      TOTAL_BLOCKS,
      blockWeightedRMS,
      blockDBA
    );
  }

  uint32_t elapsed =
    millis() - measurementStart;

  const double N =
    (double)MEASURE_SAMPLES;

  double rawMean =
    totalRawSum / N;

  // RMS around zero.
  double rawRMS =
    sqrt(totalRawSumSquares / N);

  // RMS after removing measured DC mean.
  double rawACRMS =
    sqrt(
      max(
        0.0,
        (totalRawSumSquares / N) -
        (rawMean * rawMean)
      )
    );

  double eqRMS =
    sqrt(totalEqSumSquares / N);

  double weightedRMS =
    sqrt(totalWeightedSumSquares / N);

  double rawDBFS =
    20.0 * log10(
      max(rawACRMS, 1e-12) /
      (double)((1UL << 23) - 1UL)
    );

  double theoreticalDBA =
    MIC_REFERENCE_DB_SPL +
    RMS_TO_PEAK_DB +
    20.0 * log10(
      max(weightedRMS, 1e-12) /
      MIC_REF_PEAK
    ) +
    CALIBRATION_OFFSET_DB;

  Serial.println();
  Serial.println("============================================================");
  Serial.println("5-SECOND RESULT");
  Serial.println("============================================================");

  Serial.printf(
    "Elapsed time:                 %.3f s\n",
    elapsed / 1000.0
  );

  Serial.printf(
    "Samples processed:            %d\n",
    MEASURE_SAMPLES
  );

  Serial.println();

  Serial.printf(
    "Raw mean (DC estimate):       %.3f\n",
    rawMean
  );

  Serial.printf(
    "Raw RMS (zero reference):     %.3f\n",
    rawRMS
  );

  Serial.printf(
    "Raw AC RMS (DC removed):      %.3f\n",
    rawACRMS
  );

  Serial.printf(
    "Raw AC level:                 %.3f dBFS RMS\n",
    rawDBFS
  );

  Serial.println();

  Serial.printf(
    "Equalized RMS:                %.3f\n",
    eqRMS
  );

  Serial.printf(
    "A-weighted RMS:               %.3f\n",
    weightedRMS
  );

  Serial.printf(
    "Reference peak:               %.3f\n",
    MIC_REF_PEAK
  );

  Serial.println();

  Serial.printf(
    "THEORETICAL LAeq,5s:          %.2f dBA\n",
    theoreticalDBA
  );

  Serial.printf(
    "Calibration offset:           %.2f dB\n",
    CALIBRATION_OFFSET_DB
  );

  Serial.println();

  Serial.println("Interpretation:");
  Serial.println("  This value is the datasheet-based theoretical result.");
  Serial.println("  Do NOT adjust the offset yet.");
  Serial.println("  First verify raw signal behavior and repeatability.");
  Serial.println("============================================================");
}

// ============================================================
// SETUP
// ============================================================

void setup() {

  Serial.begin(115200);
  delay(1000);

  printHeader();

  Serial.println();
  Serial.println("Initializing filters...");

  if (!inmp441EQ.begin()) {
    Serial.println("ERROR: INMP441 EQ allocation failed.");
    while (true) delay(1000);
  }

  if (!aWeight.begin()) {
    Serial.println("ERROR: A-weighting allocation failed.");
    while (true) delay(1000);
  }

  inmp441EQ.reset();
  aWeight.reset();

  Serial.println("Filters initialized.");

  Serial.println();
  Serial.println("Initializing I2S...");

  setupI2S();

  Serial.println("I2S initialized.");
  Serial.println();
  Serial.println("Waiting for microphone startup...");

  // INMP441 requires startup time after power-up.
  delay(500);

  Serial.println("Ready.");
}

// ============================================================
// LOOP
// ============================================================

void loop() {

  // Reset filter states before each independent measurement.
  inmp441EQ.reset();
  aWeight.reset();

  // Discard startup/filter transient.
  settleFilters();

  runMeasurement();

  Serial.println();
  Serial.println("Next measurement will begin in 3 seconds...");
  Serial.println();

  delay(3000);
}
