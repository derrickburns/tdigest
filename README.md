# tdigest

This is an implementation of Ted Dunning's [Merging T-Digest](https://github.com/tdunning/t-digest/).

This implementation batches all inserts, including merges.  This is an improvement over the original implementation where merges require sorting of all points on each merge.

This implementation does not support storing the incoming data with each centroid, (since obviously that is for testing).

## Usage

The T-Digest is a data structure for accurate online computation of quantiles and cumulative distribution functions (CDFs) on streaming data.

### Basic Example

```cpp
#include "tdigest2/TDigest.h"

using namespace tdigest;

int main() {
    // Create a T-Digest with compression factor (default: 1000)
    // Higher compression = more accuracy but more memory
    TDigest digest(1000);

    // Add data points
    digest.add(1.0);
    digest.add(2.0);
    digest.add(3.0);

    // Or add with weights
    digest.add(4.0, 2.5);  // value, weight

    // Compute quantiles (0.0 to 1.0)
    double median = digest.quantile(0.5);      // 50th percentile
    double p95 = digest.quantile(0.95);        // 95th percentile
    double p99 = digest.quantile(0.99);        // 99th percentile

    // Compute cumulative distribution function
    double cdf_at_2 = digest.cdf(2.0);         // P(X <= 2.0)

    return 0;
}
```

### Merging Multiple T-Digests

```cpp
TDigest digest1(1000);
TDigest digest2(1000);

// Add data to both digests
digest1.add(1.0);
digest2.add(2.0);

// Merge digest1 into digest2
digest2.merge(&digest1);

// Now digest2 contains data from both
double combined_median = digest2.quantile(0.5);
```

### Batch Processing

```cpp
TDigest digest(1000);

// Add many values
for (int i = 0; i < 1000000; i++) {
    digest.add(some_value);
}

// Optional: compress to ensure processing
digest.compress();

// Compute quantiles
double q25 = digest.quantile(0.25);
double q50 = digest.quantile(0.50);
double q75 = digest.quantile(0.75);
```

## Building

This project uses Bazel for building and testing:

```bash
bazel test //...
```

## API Reference

### Key Methods

- `TDigest(compression)` - Constructor with compression factor (default: 1000)
- `add(value)` - Add a single data point
- `add(value, weight)` - Add a weighted data point
- `quantile(q)` - Get the q-th quantile (0.0 to 1.0)
- `cdf(x)` - Get the cumulative probability P(X <= x)
- `merge(other)` - Merge another T-Digest into this one
- `compress()` - Force compression of unprocessed data

