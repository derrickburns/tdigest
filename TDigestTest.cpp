#include "gtest/gtest.h"
#include "glog/logging.h"
#include "tdigest2/TDigest.h"

namespace stesting {

class TDigestTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  TDigestTest() {
    // You can do set-up work for each test here.
  }

  virtual ~TDigestTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  static void SetUpTestCase() {
    static bool initialized = false;
    if (!initialized) {
      FLAGS_logtostderr = true;
      google::InstallFailureSignalHandler();
      google::InitGoogleLogging("testing::TDigestTest");
      initialized = true;
    }
  }

  // Objects declared here can be used by all tests in the test case for Foo.
};


TEST_F(TDigestTest, TDigestMergeTest2) {
  tdigest::TDigest tdigest1(1000);
  tdigest::TDigest tdigest2(1000);
  for (int i = 0; i <= 10 * 500 * 1000; i++) {
    tdigest1.add(std::rand() % 1001);
  }
  for (int i = 0; i <= 10 * 500 * 1000; i++) {
    tdigest2.add(std::rand() % 1001);
  }
  bool success;

  tdigest1.merge(tdigest2);

  EXPECT_NEAR(tdigest1.quantile(0.5, &success), 500.0, 2.0);
  EXPECT_NEAR(tdigest1.quantile(0.95, &success), 950.0, 1.0);
  EXPECT_NEAR(tdigest1.quantile(0.99, &success), 990.0, 0.5);
  EXPECT_NEAR(tdigest1.quantile(0.995, &success), 995.0, 0.75);
  EXPECT_NEAR(tdigest1.quantile(0.999, &success), 999.0, 0.60);
}

TEST_F(TDigestTest, TDigestCompressTest) {
  tdigest::TDigest tdigest1(1000);
  for (int i = 0; i <= 10 * 500 * 1000; i++) {
    tdigest1.add(std::rand() % 1001);
  }
  tdigest1.compress();

  bool success;
  EXPECT_NEAR(tdigest1.quantile(0.5, &success), 500.0, 2.0);
  EXPECT_NEAR(tdigest1.quantile(0.95, &success), 950.0, 1.0);
  EXPECT_NEAR(tdigest1.quantile(0.99, &success), 990.0, 0.5);
  EXPECT_NEAR(tdigest1.quantile(0.995, &success), 995.0, 0.75);
  EXPECT_NEAR(tdigest1.quantile(0.999, &success), 999.0, 0.60);
}

}  // namespace stesting

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
